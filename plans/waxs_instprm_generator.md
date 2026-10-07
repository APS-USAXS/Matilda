# PLAN: GSAS-II `.instprm` generator from a WAXS LaB6 standard

Status: **design only, nothing implemented** (written 2026-10-02)
Owner: Jan Ilavsky

## 1. Goal

Let a user on the **WAXS tab** of `matilda-gui` pick a LaB6 (or CeO2) standard
measurement, fit its 1D pattern, and write a GSAS-II instrument-parameter file
(`.instprm`) describing the peak profile of the USAXS-station WAXS detector.
Later: expose the same function as a script / CLI (no Qt).

Background: collaborators periodically ask for an `.instprm` for the WAXS
device. GSAS-II does not derive it from hardware; it is **fitted from a
standard** and saved (GSAS-II: Powder Data > Instrument Parameters >
Operations > Save Profile). We reproduce that fit ourselves, so GSAS-II code
is NOT required, only its file format and its width formulas.

## 2. Non-goals (first version)

- No Rietveld/whole-pattern refinement, no lattice-parameter refinement
  beyond a sanity check.
- No 2D-image calibration (distance/tilt/beam centre) — that is done by the
  existing Nika→pyFAI path in `convertSWAXS.py`.
- No TOF / energy-dispersive profile types. Constant-wavelength X-ray only
  (`Type:PXC`).
- No `.imctrl` export (possible later; GSAS-II image-control format).

## 3. What the file is

Plain text, key:value, first line a fixed comment. Expected content
(**verify against a real GSAS-II-written file and against `G2fil.py` /
`G2IO.py` readers before coding; do not trust this listing from memory**):

```
#GSAS-II instrument parameter file; do not add/delete items!
Type:PXC
Bank:1.0
Lam:<wavelength A>
Polariz.:<0.99 synchrotron>
Azimuth:0.0
Zero:<2theta zero, deg>
U:  V:  W:      Gaussian widths
X:  Y:  Z:      Lorentzian widths
SH/L:           low-angle asymmetry
```

Open items about the format:
- Exact key set/order/number formatting GSAS-II writes (also `I(L2)/I(L1)`,
  `Lam1`/`Lam2` for doublet sources; synchrotron uses single `Lam`).
- Task: save one from GSAS-II using a real WAXS LaB6 pattern and keep it in
  `TestData/` as the reference for round-trip tests.

## 4. Profile model (GSAS-II convention — VERIFY)

Thompson–Cox–Hastings pseudo-Voigt in 2θ:

- Gaussian variance  σ²(centideg²) = U·tan²θ + V·tanθ + W
- Lorentzian width   γ(centideg)   = X/cosθ + Y·tanθ + Z   (Z term: check)
- θ = half of 2θ. Units are **centidegrees** in the file.
- Convert σ² → Gaussian FWHM = 2·sqrt(2 ln2)·sqrt(σ²); check whether γ is
  FWHM or HWHM in GSAS-II. Authoritative source: GSAS-II `GSASIIpwd.py`
  (`getgamFW`, `getFWHM`, `getWidthsCW`). Read it in a clone of
  https://github.com/AdvancedPhotonSource/GSAS-II before writing the
  converter. Do NOT import GSAS-II; only mirror the formulas.
- Wavelength in file = the one in the Matilda HDF5
  (`instrument/monochromator/wavelength`, Å).

## 5. Algorithm

Input: 1D I(Q) of a standard + wavelength + standard identity.

1. **Get 1D data.** Use Matilda's reduced (uncalibrated/Normalized is fine;
   absolute scale irrelevant) I(Q) from `reduceADData(...)`. Do not use
   calibrated/blank-subtracted data unless the user supplies a blank; raw
   normalized is the default. Q in Å⁻¹.
2. **Identify the standard and predict peak positions.** The standard is
   chosen **from the sample name** (case-insensitive match on `LaB6` /
   `CeO2`, with a manual override combo in the GUI when the name is
   ambiguous or matches neither). LaB6 is the default; CeO2 is used by some
   users at specific energies. Lattice parameter is a preset **plus an
   editable field** — **take the exact certified value for the specific SRM
   lot from the NIST certificate (LaB6 SRM 660x, CeO2 SRM 674x); never
   hard-code a guess from memory.** Q_hkl = 2π·sqrt(h²+k²+l²)/a, with the
   allowed (hkl) set per structure (LaB6 primitive cubic: all hkl; CeO2
   fluorite fcc: h,k,l all even or all odd). Keep only peaks inside the
   measured Q range.
   Sanity check (not refinement): compare fitted Q0 to prediction; report
   residual and implied 2θ zero / wavelength error.
3. **Convert Q → 2θ** with the file wavelength: 2θ = 2·asin(Qλ/4π).
4. **Fit peaks.** Matilda stays independent of pyirena (decision, see §8):
   write a small self-contained fitter in `matilda/instprm.py`, copying and
   trimming what is needed from `pyirena/core/waxs_peakfit.py` (Voigt =
   Gaussian `FWHM` + Lorentzian `FWHM_L`, `voigt_fwhm`, polynomial/linear
   background). Keep only what this feature uses; add a header comment saying
   the code derives from pyirena (same author) so the two can be
   cross-checked, with no import link. Fit in Q or 2θ; keep per-peak Gaussian and Lorentzian FWHM as
   outputs. Fit all peaks in one pass with a shared background; seed
   positions from step 2 (positions are known, so constrain windows tightly).
   Reject peaks that overlap, sit on the beam stop/detector gap, or fail the
   fit (flag in report, don't silently drop).
5. **Fit width functions vs 2θ.** Weighted least squares (scipy
   `curve_fit`/`least_squares`) of the per-peak σ² and γ to the U,V,W /
   X,Y,Z forms. U/V/W and X/Y/Z are strongly correlated on a limited 2θ
   range → use bounds (σ²≥0, γ≥0), let the user **fix** parameters, and
   show parameter uncertainties.

   **Number of usable peaks depends on energy** (Jan, 2026-10-02): at 28 keV
   the WAXS range is wide enough for plenty of peaks; at 12 keV there are at
   most ~3. So the number of free parameters must be chosen **automatically
   from the number of fitted peaks**, never fixed in advance. Proposed tiers
   (to be tuned on real data; "DoF" = peaks minus free parameters per width
   function):

   | usable peaks | Gaussian σ² | Lorentzian γ |
   |---|---|---|
   | ≥ 8 | U, V, W free | X, Y free (Z fixed 0 unless justified) |
   | 4–7 | W (+U) free, V=0 | X free, Y free if DoF allows |
   | 3 | W free, U=V=0 | X free, Y=0 |
   | 1–2 | refuse: report widths of the peaks found, write no `.instprm` unless user forces (then W, X from the mean) |

   The chosen tier, the fixed parameters and a "low-confidence" warning must
   appear in the report and in the GUI. The user can override the free/fixed
   choice. Exact thresholds are an open question (see §8) and should be
   decided by fitting real 12/21/28 keV LaB6 patterns.
6. **Zero offset.** Fit `Zero` from Q0(measured) − Q0(predicted), reported
   in degrees 2θ. Cross-check that it is not just compensating for a
   wrong distance/wavelength in the Nika geometry.
7. **SH/L asymmetry.** pyirena has no asymmetric profile. v1: write a fixed
   default (e.g. 0.002 — **check GSAS-II's default**) and say so in the
   report; v2: add an asymmetric (Finger–Cox–Jephcoat) fit if the data
   demand it. Flag this in the output so users know it was not fitted.
8. **Write file** + a sidecar report (JSON or `.txt`): data file, standard,
   a, wavelength, peaks used with fitted Q0/FWHM, residuals, parameter
   values ± σ, software version, date.

## 6. Where it goes in the code

Follow the layering rule: core math has no Qt.

| Piece | Location | Notes |
|---|---|---|
| Fit + width model + file writer | new `matilda/instprm.py` (core, numpy/scipy only) | pure functions, no GUI, no HDF5 side effects |
| Standard presets | same module (dict of name → a, hkl list) | user-editable lattice parameter |
| Peak fitting | inside `matilda/instprm.py` (copied/trimmed from pyirena `waxs_peakfit.py`) | **no pyirena dependency**; scipy/numpy only, already Matilda deps |
| GUI | `matilda/gui/data_reduction/parameter_tabs.py`, `_WAXSTab` | new "Instrument profile (GSAS-II)" group box |
| Worker | `matilda/gui/data_reduction/reduction_worker.py` pattern | run fit off the GUI thread |
| Script/CLI | later, thin wrapper calling `instprm.py` | e.g. `matilda-instprm file.hdf --standard LaB6 -o waxs.instprm` |
| Tests | `tests/test_instprm.py` | see §7 |

Existing hooks to reuse:
- `_WAXSTab` already has blank selection, thickness, calibration sections and
  a placeholder label ("Additional parameters … future version") where the
  new group box can go.
- File tree (`file_tree.py`) has a context menu that already special-cases
  WAXS ("Set as WAXS Blank"). Add **"Use as LaB6 standard"** here, or a
  dedicated button in the new group box that uses the currently selected
  WAXS file.
- `process2Ddata` / `reduceADData` provide Q, Intensity, Error; wavelength
  via `instrument_dict["monochromator"]["wavelength"]`. Mask and integrator
  caching already handled in `convertSWAXS.py`.

## 7. GUI sketch (WAXS tab)

Group box **"Instrument profile (GSAS-II .instprm)"**:
- Standard: combo, **auto-set from the sample name** (LaB6 / CeO2 / custom a),
  with the detected choice shown and overridable
- Lattice parameter (Å): editable, defaults from preset
- Q range / peak list: auto from standard; table with include checkboxes
- Fix/free checkboxes for U, V, W, X, Y, Z; SH/L value
- Button **"Fit profile"** → plot fitted peaks over data (reuse
  `sas_plot.py`) and a table of fitted parameters ± σ
- Button **"Save .instprm…"** + report path
- Warnings area (too few peaks, correlated parameters, big zero offset)

Keep the panel thin: the GUI only collects parameters and displays results;
all numerical work in `matilda/instprm.py`.

## 8. Decisions

Decided (Jan, 2026-10-02):
1. **No interdependency between Matilda and pyirena.** Copy/trim or write
   new code inside Matilda; no pyirena import, no new dependency.
2. **Standards:** LaB6 normally; CeO2 for some users at specific energies.
   Distinguish by **sample name**. Still needed: the certified lattice
   parameter of the SRM lot(s) actually used (fill in the preset when known).
3. **Input:** Matilda-reduced 1D data for now. Other inputs (ASCII, NXcanSAS
   from elsewhere) can be added later without changing the core function —
   keep `fit_instrument_profile(q, I, err, wavelength, standard, ...)`
   array-based so only the loader differs.
4. **Fixed/free parameters:** unknown a priori. Strategy is the automatic
   peak-count tiers in §5; thresholds to be tuned on real data.
5. **Target format:** TCH-compatible GSAS-II `.instprm` is sufficient, since
   requesters will use GSAS-II.

Still open:
- Certified lattice parameters (LaB6 SRM 660 letter, CeO2 SRM 674 letter).
- Tier thresholds in §5 (needs real patterns at ≥ 12 keV and 28 keV, ideally
  21 keV too).
- How to behave at ≤ 2 peaks (refuse vs force). Current proposal: refuse by
  default, allow force with a loud warning.
- Whether to save/reuse the `.instprm` per energy+distance automatically
  (e.g. a small registry), since the profile is only valid for one geometry.

## 9. Validation (scientific correctness gate)

- Fit a real WAXS LaB6 pattern in GSAS-II (peak-fit mode, refine U,V,W,X,Y,
  SH/L). Compare with our result: U,V,W,X,Y within uncertainty, and — more
  meaningful given correlations — compare the **FWHM vs 2θ curve** from both.
- Round trip: write our `.instprm`, load it in GSAS-II, overlay simulated
  LaB6 pattern on the data; residual should be at noise level.
- Unit tests: synthetic pattern from known U,V,W,X,Y → recovered within
  tolerance; writer output parses with the GSAS-II reader format; Q↔2θ
  conversions; behaviour with 1–2 usable peaks (must refuse/warn).
- No change to existing reduction output: this feature must not alter any
  numbers produced by `convertSWAXS.py`.

## 10. Implementation order

1. Get a reference `.instprm` + LaB6 WAXS patterns into `TestData/` —
   **at least one at 28 keV (many peaks) and one at 12 keV (≈3 peaks)**, plus
   a CeO2 pattern if available, so the peak-count tiers can be tested.
2. Read GSAS-II `GSASIIpwd.py` width functions; pin down §4 formulas.
3. `matilda/instprm.py`: standards, peak prediction, width model, writer.
4. Fitter (per §8 decision) + tests on synthetic data.
5. Validate against the GSAS-II fit of the real pattern.
6. WAXS-tab group box + worker thread.
7. Script/CLI wrapper; docs in `docs/matilda-gui.md`; CHANGELOG entry.

## 11. Risks

- Strong U/V/W and X/Y/Z correlation on a limited 2θ range → unstable
  parameters. Mitigate with fixed defaults, bounds, reported uncertainties.
- Peak count scales with energy: ~3 peaks at 12 keV vs plenty at 28 keV.
  Low-energy profiles will be less well determined; the tier logic and
  warnings in §5 exist for this. Check actual Q coverage per energy before
  fixing the tier thresholds.
- Profile of an azimuthally integrated area detector may depend on detector
  distance/energy: the file is valid only for that geometry; record it in
  the report and filename.
- GSAS-II format details from memory may be off: verify (§3, §4).
