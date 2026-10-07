# Desmearing methods — developer reference

Technical reference for the selectable slit-desmearing methods added in
`matilda/desmearing_methods.py` (branch `feature/desmearing-methods`). Read this
first when you need to change, fix, or extend the desmearing options.

Companion documents: `docs/desmearing-paper.md` (the theory + citations for the
instrument paper) and `docs/new_desmearing.md` (the derivation, numerical
validation and Igor prototype of the truncated-Abel method). This file is the
"how the code works and where to touch it" reference.

---

## 1. What was added and why

Matilda historically desmeared with a single algorithm, the **Lake** iterative
method (`matilda/desmearing.py`, `desmearData`). It works well on strong data
but amplifies noise and, on weak / near-zero / over-subtracted intensities,
produces non-physical spikes and negative values. A batch test over 80 random
datasets found this on ~39% of files (see the paper doc for numbers).

Two families of alternatives were added, both in `desmearing_methods.py` and
both selectable from the GUI:

- **Gaussian-process (GP, "Huang-style") desmearing** — a regularised Bayesian
  fit. Stays smooth and positive on weak data and returns a posterior
  uncertainty. Added first; still available.
- **Truncated-Abel inversion** — an analytical, non-iterative inversion of the
  slit integral. This is **now the default** (see §6). It is linear in the data,
  so its uncertainties are well defined; it needs no iteration and no global
  solve; and on synthetic tests it is several times more accurate than Lake at
  the same noise level. It was developed and validated in Igor first
  (`IN3_DesmearDataAbel`) — `docs/new_desmearing.md` has the derivation and the
  error tables this implementation reproduces.

### References

- G.-R. Huang, Y. Wang, C. Do, W.-R. Chen *et al.*, "Analytical desmearing of
  Bonse–Hart USANS data via truncated Abel inversion", *J. Appl. Cryst.*
  (ORNL / NTHU). Source of Eq. 1/7/9–15 and of the χ²≈1 smoothing criterion.
- J. A. Lake, "An iterative method of slit-correcting small angle X-ray data",
  *Acta Cryst.* **23** (1967) 191–194. The historical `desmearData` method.
- J. Ilavsky & P. R. Jemian, "Irena: tool suite for modeling and analysis of
  small-angle scattering", *J. Appl. Cryst.* **42** (2009) 347–353. The Igor
  implementation whose `SlitLength` convention and `IN3_ExtendData` behaviour
  this code follows.

---

## 2. File map / touch points

| File | Change |
|------|--------|
| `matilda/desmearing_methods.py` | **New.** Slit forward operator, GP core, Abel core, `desmear_dispatch`, GUI string→key maps. numpy-only. |
| `matilda/desmearing.py` | Unchanged. `desmearData` (Lake) is called by the dispatcher for `method='lake'`; `extendData` and `calculateErrors` are reused by the Abel path. |
| `matilda/convertFlyscan.py` | `processFlyscan` gains `desmear_method`, `gp_length_scale`, `gp_kernel`, `abel_auto_smooth`, `abel_smooth_w`, `abel_num_mc` kwargs (defaults `'abel'`, `0.5`, `'matern32'`, `True`, `0.05`, `20`); the `desmearData(...)` call became `desmear_dispatch(...)`. |
| `matilda/convertUSAXS.py` | Same change for `processStepscan`. |
| `matilda/gui/data_reduction/reduction_worker.py` | Forwards the six new params from the GUI `params` dict to both converters. |
| `matilda/gui/data_reduction/parameter_tabs.py` | USAXS tab (`_USAXSTab`, shared by Flyscan+StepScan): Method / Abel smoothing / Abel MC / GP kernel / GP length-scale widgets + `get_params()` keys; imports `GUI_METHOD_MAP`, `GUI_KERNEL_MAP`. |

Data flow (unchanged structure, one extra hop):

```
GUI (_USAXSTab.get_params) → params dict
   → reduction_worker (forwards desmear_method / abel_* / gp_*)
      → processFlyscan / processStepscan
         → desmear_dispatch(method=…)          ← single entry point
            → _abel_core (default)  OR  desmearData (Lake)  OR  _gp_core (GP)
```

---

## 3. Public API — `desmear_dispatch`

```python
desmear_dispatch(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=None, *,
                 method="abel", length_scale_decades=0.5, kernel="matern32",
                 sigma_log=4.0, extrap_method="PowerLaw w flat",
                 extrap_qstart=None, max_iter=20,
                 abel_auto_smooth=True, abel_smooth_w=0.05, abel_num_mc=20)
    -> (DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ)
```

Same call signature and 4-tuple return as `desmearData`, so it is a drop-in
replacement at the call site. Returns `(None, None, None, None)` on empty input
or missing slit length (matches `desmearData`).

`method` keys:

| key | meaning |
|-----|---------|
| `abel` | **Default.** Truncated-Abel inversion with the exact finite-slit correction (§4a). `DSM_Error` from Monte Carlo. |
| `lake` | Calls `desmearData` unchanged. **Bit-for-bit identical** to the old path (regression-verified). |
| `lake_smooth` | Lake, then geometric-mean smoothing (window = `length_scale_decades`). Best for very noisy flat-plateau data where the GP would ring; no band. |
| `gp` | GP desmearing, kernel from `kernel` arg. Resolution-aware prior on by default. |
| `gp_rbf` | GP desmearing forced to the RBF kernel (shortcut). |
| `gp_lake_mean` | GP anchored on a smoothed Lake curve as prior mean (robust low-q edge). |

For GP methods `DSM_Error` is the **posterior 1σ** (in intensity units); `DSM_dQ`
is `SMR_dQ` carried through (interpolated if the point count changed).

---

## 4a. Algorithm internals (`_abel_core`) — the default

The Irena/Matilda slit smearing is, with `L` ≡ `slitLength` (the half-length,
same normalisation as `IN3_SmearData` and `smearIntensityArray`):

```
Ismr(q) = (1/L) ∫_0^L I( √(q² + y²) ) dy
```

Substituting `s = √(q²+y²)` turns this into an Abel-type integral, which can be
inverted in closed form. The paper's Eq. 15 does so but treats the upper limit
`Qm(q) = √(q²+L²)` as independent of `q`; that is fine for USANS (`Q ≪ L`) and
**wrong by up to ~5% near `Q ≈ L` for USAXS**. The exact finite-slit relation
carries one extra term and simplifies to what the code implements:

```
I(Q) = −(2L/π) · T1(Q)  +  T2(Q)

T1(Q) = ∫_Q^{√(Q²+L²)}  Ismr′(x) / √(x² − Q²) dx
T2(Q) = (2/π) ∫_0^{π/2} I( √(Q² + 2L²/(1+sin u)) ) du
```

`T2` needs the *desmeared* `I` only at `s ≥ √(Q²+L²) > Q`, so the whole curve
follows from **one sweep from the highest Q downwards** — no iteration, no
global solve. The derivation, the sanity checks (constant `I` → `T1=0, T2=C`;
`Q ≫ L` → the first-order result) and the accuracy tables are in
`docs/new_desmearing.md` §4.

Code pieces, in call order:

- **`_abel_find_smooth_width(q, y, err)`** — the inversion differentiates the
  data and the `1/√(x²−Q²)` kernel weights the neighbourhood of `Q` heavily, so
  smoothing is mandatory: 1% noise becomes ~7% scatter without it. This is the
  paper's Step 4 — bisect in `ln w` over `[0.002, 0.5]` for the **largest** width
  whose reduced χ² of `(smoothed − measured)/SMR_Error` is still ≤ 1. χ² rises
  monotonically with `w`, so bisection is safe.
  The statistic is the **robust** reduced χ², `median(r²)/median(χ²₁)` with
  `median(χ²₁) = 0.45494`, *not* the mean. This matters: real USAXS curves carry
  a few outliers (bad points, or a feature too sharp for the local-linear fit),
  and under a mean-based χ² two or three of them reach χ² = 1 on their own. On
  `TestData/NXcanSAS_SMR.h5` that collapsed the width from 0.16 to 0.016 — no
  smoothing at all — and the desmeared curve came back as noisy as Lake's.
  Regression test: `test_abel_auto_width_survives_outliers`.
  The criterion still trusts the *scale* of `SMR_Error`: inflated errors
  over-smooth, underestimated errors under-smooth. The chosen width and χ² are
  logged at INFO level; uncheck "Auto" in the GUI and set the width by hand when
  they look wrong.
- **`_abel_smooth_local_linear(q, y, w, use_log)`** — Gaussian-weighted
  local-**linear** regression of `ln I` (or `I`, if any point is ≤ 0, which
  happens at high Q after blank subtraction) against `ln q`. Width `w` is in
  `ln q` units, so `w = 0.1` is roughly ±10% in Q. Local-linear rather than the
  paper's plain Gaussian kernel because a plain average is biased at the ends of
  the range and on steep power laws — exactly where USAXS data live; it roughly
  halves the lowest-Q error. The window is truncated at `4w`, so cost is
  `O(N·window)`, not `O(N²)`.
- **`_abel_extend(...)`** — the inversion needs `Ismr` up to `√(Qmax²+2L²)`,
  about 1% past `Qmax`. It reuses Lake's `extendData` so the user-chosen
  extrapolation function and "Extrap Q start" mean the same thing for both
  methods (`extendData` reaches `√(Qmax²+(1.5L)²)`, which is enough), then
  resamples that very sparse, linearly-spaced extension onto `n_ext=24`
  log-spaced points so the `T2` quadrature near `Qmax` has something to
  interpolate. A failed extension falls back to a flat tail with a warning.
- **`_abel_invert(qe, Is, sl, n_orig, nu=64)`** — the sweep.
  - *Start values* in the extension (`Q > Qmax`, i.e. `Q ≳ 10L`): the paper's
    first-order `I(Q) = Ismr(Q) − L²/(6Q)·dIsmr/dQ`, good to ~0.01% there.
  - *T1* is integrated **analytically per segment** with `Ismr` taken piecewise
    linear, using `∫ dx/√(x²−Q²) = acosh(x_hi/Q) − acosh(x_lo/Q)`. This handles
    the `1/√` singularity at `x = Q` exactly, with no special-casing.
  - *T2* is a 64-node midpoint rule in `u`, with `I` interpolated linearly in
    `q` between already-computed points (log-log interpolation made no
    measurable difference and misbehaves on negative intensities).
  - *Implicit first bracket.* At high Q the range `[Qm, Qmm]` is only
    `≈ L²/(2Q)` wide — narrower than the point spacing — so some `T2` nodes land
    between `qe[i]` and `qe[i+1]`, where `I(qe[i])` is still unknown. Clamping
    them to `Id[i+1]` causes a systematic drift of several percent. Instead the
    dependence is kept implicit: `T2 = α·Id[i] + β`, solved as
    `Id[i] = (−(2L/π)T1 + β)/(1−α)`. `α < 1` always, because every node has
    `s > Q`.
- **Uncertainties** — the inversion is *linear* in `Ismr`, unlike Lake, so error
  propagation is well defined. `abel_num_mc` (default 20) realizations add
  `N(0, SMR_Error)` to the **measured** points only (the high-Q extension is held
  fixed), re-smooth with the same width, re-invert; `DSM_Error` is the standard
  deviation over realizations. This includes the smoothing, which no analytic
  shortcut does. `abel_num_mc < 2` falls back to the Lake-style
  `desmearing.calculateErrors` for comparison.

Numerics: `O(N²)` in the worst case (the `T1` window spans most of the data at
low Q) but vectorised per point — ~0.25 s for 500 points including 20 MC draws.
Un-rebinned fly scans with several thousand points should be rebinned first.

Validation against `docs/new_desmearing.md` §1, on the 3-level synthetic model
with `L = 0.03` and 400 log-spaced points (this implementation / the Igor
prototype): noise-free 0.18% / 0.3% rms; 1% noise with auto width 1.1% / 0.7%;
3% noise 2.9% / 1.7%. In all cases the lowest 1–3 points are the weakest (the
smoothing window is one-sided there and `T1` integrates over the most points) —
on USAXS those points sit next to the beam anyway.

---

## 4. Algorithm internals (`_gp_core`)

The GP model is multiplicative: the desmeared intensity is
`x(q) = μ0(q) · exp(s(q))`, with `s` a smooth log-correction under a GP prior.
The forward (smeared) model is linear, `y = M x`, so the MAP is found by a
damped Gauss-Newton iteration. Full theory is in `desmearing-paper.md`; the code
pieces:

- **`_SlitOperator.build(q, L, n_extra, n_slit)`** — builds the dense slit
  matrix `M` (`Nm × Nw`). Working grid = measured `q` plus `n_extra=40`
  log-spaced points up to `√(qmax²+L²)` (replaces Lake's ad-hoc high-q
  extrapolation). Slit integral by trapezoid over `n_slit=200` nodes on `[0,L]`
  with linear interpolation in `log q`.
- **Error model** (the important robustness bit):
  `σ = max(Idev, err_floor_frac·median(Idev), rel_err_floor·|y|)`
  with `err_floor_frac=0.01` (absolute floor) and **`rel_err_floor=0.01`
  (1% relative floor — the instrument's true error floor).** The relative floor
  both matches physics and regularizes files that carry absurd ~1e-7 reported
  errors (which otherwise blow up the ill-posed solve).
- **Base trend `μ0`** — `_smooth_positive_trend`: log-space (geometric-mean)
  moving average of the data after clipping to the noise floor. Log-space
  smoothing preserves multi-decade power laws; the floor lets near-zero /
  negative points map to a small positive value. Window ≈ `0.5·ls·points/decade`.
  For `gp_lake_mean`, `μ0` instead comes from `_smooth_loglog(lake)`.
- **Kernel** `_kernel(x, ell, σ, kind)` over `x = ln q`, `ell = ls·ln10`:
  RBF `σ² exp(−d²/2)` or Matérn-3/2 `σ²(1+√3 d)exp(−√3 d)`, `d=|Δx|/ell`.
- **Low-q guard** (`low_q_guard=True`, default): after the fit, over the
  under-constrained low-q block (small column weight in `M`, i.e. slit ≥ q at
  those points), clamp the desmeared intensity so it does not fall below the
  local **smeared** plateau level. Fixes the downward low-q ramp that the fit
  produces on weak flat-plateau data; leaves rising low-q power laws untouched
  (there the desmeared exceeds the smeared, so the floor never binds). This is
  the "weak/flat vs strong/rising" auto-discriminator.
- **Resolution-aware prior** (`resolution_aware=True`, default): a non-stationary
  length scale `ell(q)=max(ell0, α·W(q))` with slit resolution
  `W(q)=½·ln(1+(L/q)²)`, applied by warping the coordinate
  `u=∫dx/ell(q)` and using a unit-length kernel in `u`. This smooths the low-q
  null space on flat plateaus with slit≥q_min (kills Gibbs ringing) without
  touching data-constrained features or high-q resolution. Verified not to
  reduce a real q≈1e-3 peak (D15) or change T6 agreement.
- **Gauss-Newton loop** — Jacobian `J = M·diag(x)`; solve
  `(JᵀΣ⁻¹J + K⁻¹)δ = JᵀΣ⁻¹(y−f) − K⁻¹s`; backtracking line search (halving,
  ≤30 tries) so it cannot diverge; hard clip `|s| ≤ ln(1e8)`; stop on step
  `< tol=1e-4` or `max_iter` (GP uses `max(max_iter,30)`).
- **Uncertainty** — posterior covariance `(JᵀΣ⁻¹J + K⁻¹)⁻¹`; 1σ in `s` maps to a
  multiplicative band `I·exp(±2σ_s)` (`n_sigma_band=2`).

Numerics: dense `O(Nw³)` solve per iteration, `Nw ≈ Npts + 40`. Fine for
rebinned USAXS (~500 pts, well under a second). Step scans that keep many points
would be slow — see §7.

---

## 5. GUI wiring

`_USAXSTab` (in `parameter_tabs.py`) builds these widgets in the Desmearing
form, after "Extrap Q start":

- `_desmear_method` (QComboBox) — items from `GUI_METHOD_MAP.keys()`, in order:
  `"Truncated Abel"` (the default, because it is first), `"Lake"`,
  `"Lake (smoothed)"`, `"Huang GP"`, `"Huang GP + Lake mean"`.
- `_abel_auto_smooth` (QCheckBox, checked) + `_abel_smooth_w` (QDoubleSpinBox,
  0.002–0.5 in `ln Q` units, default 0.05) — on one row, "Abel smoothing".
- `_abel_num_mc` (QSpinBox) — 0–200, default 20.
- `_gp_kernel` (QComboBox) — `"Matérn-3/2"`, `"RBF"`.
- `_gp_length_scale` (QDoubleSpinBox) — 0.1–2.0 decades, default 0.5.

`_on_desmear_method_changed` enables the GP widgets only when the method name
contains "GP" and the Abel widgets only when it contains "Abel"; the manual
width is additionally greyed out while "Auto" is checked. `get_params()` maps
the display strings to keys via `GUI_METHOD_MAP` (→ `desmear_method`) and
`GUI_KERNEL_MAP` (→ `gp_kernel`) and adds `gp_length_scale`,
`abel_auto_smooth`, `abel_smooth_w`, `abel_num_mc`. Those keys are forwarded
verbatim by `reduction_worker`.

To add a method to the dropdown: add an entry to `GUI_METHOD_MAP` in
`desmearing_methods.py` (display string → `(method_key, kernel_key)`) and handle
`method_key` in `desmear_dispatch`. No GUI code change needed.

---

## 6. Invariants to preserve

- **The default is truncated Abel.** `desmear_method='abel'` is the default of
  `desmear_dispatch`, of both converters and of every `params.get(...)` in the
  worker — which means **the production daemon uses it too**, since `matilda.py`
  does not pass the argument. This was a deliberate promotion (the method is
  under active tuning); keep the four defaults in step if it changes again.
- **`method='lake'` must stay bit-for-bit identical to `desmearData`.** The
  dispatcher just forwards to it. Keep it that way (there is/should be a
  regression test asserting equality).
- **4-tuple contract.** Every branch returns `(Q, I, Error, dQ)`; all callers
  unpack exactly four values.

---

## 7. Known limitations / TODO

- **No provenance yet.** The desmeared NXcanSAS entry does not record which
  method/params produced it. Re-processing overwrites it (SMR is retained) but
  the method is not labelled. TODO: write `desmear_method` plus the method's
  parameters (`abel_smooth_w` and the fitted χ², `abel_num_mc`, or
  `gp_length_scale`/`gp_kernel`) as attributes on the desmeared `sasdata` group
  in `hdf5code.saveNXcanSAS`. Downstream readers must default the field to
  `"lake"` so old files keep working. The Igor side records the same thing in
  the `DSM_Int` wave note.
- **Abel: sharp features.** Data whose features (form-factor minima, structure
  peaks) are as narrow as the slit length cannot be desmeared reliably by *any*
  method, and noise is amplified at the minima. Fit SMR data with a smeared
  model instead — the same advice Irena already gives.
- **Abel: the auto width trusts the scale of `SMR_Error`.** See §4a. Watch the
  logged width and χ²; on the test file it lands around 0.16 in `ln Q`. A width
  near the `0.002` floor means the criterion has failed and the result will be
  as noisy as Lake.
- **Abel: `O(N²)`.** Rebin fly scans before desmearing (the converters already
  do this above 800 points). `T1` could be vectorised into a banded matrix of
  `acosh` weights if this ever matters.
- **Performance on large step scans.** `O(Nw³)`; if `Npts > ~1500` the GP is
  slow/memory-heavy. TODO: subsample the working grid for GP, or cap.
- **Forward-operator vs Matilda smearing.** `_SlitOperator` is internally
  self-consistent but differs from `desmearing.smearIntensityArray` by ~20% at
  the high-q tail (different tail handling). Fine for the GP fit (self-consistent)
  but reconcile before trusting cross-method χ² numerically.
- **RBF erases sharp features.** RBF (C^∞ prior) smooths real correlation peaks;
  Matérn-3/2 is the safe default. Keep RBF for featureless samples only.
- **`sigma_log`** is a secondary amplitude knob; it mostly affects the
  under-constrained edges and band width, not the well-determined range. Not
  exposed in the GUI on purpose.
- **Auto-selection** of kernel/length scale by GP marginal likelihood (evidence)
  is a natural future addition; not implemented.

---

## 8. How to test / verify

- Lake-equivalence (must pass):
  ```python
  a = desmearData(Q, I, E, dQ, slitLength=L, ExtrapMethod="PowerLaw w flat",
                  ExtrapQstart=None, MaxNumIter=20)
  b = desmear_dispatch(Q, I, E, dQ, L, method="lake",
                       extrap_method="PowerLaw w flat", max_iter=20)
  assert all(np.array_equal(x, y) for x, y in zip(a, b))
  ```
- GP smoke: `desmear_dispatch(..., method="gp")` returns positive, finite
  intensity with no giant point-to-point log jumps, on a weak/negative dataset.
- Abel closure (both in `tests/test_desmearing.py`): smear a known model by
  direct quadrature, desmear it noise-free with smoothing off, and require
  `rms(I/I_true − 1) < 1%`; then add 2% noise and require the Abel rms to be
  less than half the Lake rms on the same input.
- Real data, side by side: reduce the same sample with each method and compare.
  Useful samples are glassy carbon (smooth calibration standard), a
  Porod-dominated powder, a sample with a weak broad peak, one with a sharp
  form-factor minimum (expect every method to struggle) and a noisy weak
  scatterer. The real acceptance criterion is downstream: Unified / size
  distribution parameters from the desmeared curve should agree with a
  slit-smeared model fitted to the SMR data.
- Sandbox: the full method development, batch runner (`batch_compare.py`),
  negative-injection test, and figures live in the separate **`Better desmearing`**
  project; `desmear_lab.huang_gp` there is the reference implementation this
  module was vendored from.
