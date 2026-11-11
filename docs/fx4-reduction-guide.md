# FX4 Data Reduction — Implementation Guide

**Purpose.** This document describes exactly how Matilda recognises and reduces
data from the FX4 counting chain, for all four USAXS/SAXS/WAXS techniques, so
that the Igor Pro reduction can be reimplemented to follow the same path.

**Audience.** Someone writing the Igor Pro implementation. It assumes
familiarity with the existing Igor USAXS reduction (Indra/Nika) and with the
old counting chain, and spells out only what changed and what the new
arithmetic is.

**Status of the Python side.** The reduction described here runs and has been
validated: an FX4 fly scan of glassy carbon SRM 3600 against its blank gives
transmission 0.936 and an absolute plateau of 29.75 cm²/cm³ over
0.01 < q < 0.1 Å⁻¹, against a certified value near 30. Two constants are
deliberately **not** yet calibrated and are called out in §7.

Source of truth in the Python code:
`matilda/fx4support.py`, `matilda/supportFunctions.py`,
`matilda/convertUSAXS.py`, `matilda/convertSWAXS.py`.
Test data: `TestData/FX4Set/` (one blank + one glassy carbon per technique).

---

## 1. What changed, in one paragraph

On 2026-09-26 the 12-ID-E counting chain changed from
*Femto amplifier + V/F converter + Struck 3820 scaler* to **FX4 electrometers**.
All four data formats changed with it. The single consequence that drives
every formula below:

> **Detector values are now gain-independent mean PICOAMPS, not counts.**
> There is no amplifier gain to divide by, no per-range gain table, and
> no dwell-time divisor.

| aspect | **scaler** (until 2026-09-26) | **FX4** (from 2026-09-26) |
| --- | --- | --- |
| hardware | Femto + V/F + Struck 3820 | FX4 electrometers |
| detector value | counts accumulated over a dwell | mean current, pA |
| gain | per-range table, divided out | none — reading is gain-independent |
| dwell | divided out (1e6 Hz fly MCA / 1e7 Hz Joerger) | none — value is already a mean |
| `I(q)` | `(counts − dark·t)/(f·gain) / (I0/I0gain)` | `(upd_current − dark) / I0_current` |

The reduction *downstream* of the normalised `I` vs `AR` signal — peak finding,
Q construction, blank subtraction, K-factor calibration, desmearing — is
**unchanged**. Only the front end (file reading and normalisation) and the
uncertainty estimate differ. This is the key point for the Igor port: you are
replacing the gain-correction block, not the pipeline.

---

## 2. Detecting an FX4 file

### 2.1 The rule

Every FX4 format writes a `counting_chain` marker whose value is the string
`"FX4"`. **Absent means the old scaler chain.** That is the whole test.

The marker lives in a different place in each format. Check all of them, first
hit wins:

| # | format | location | kind |
|---|---|---|---|
| 1 | fly scan (documented) | `/entry/flyScan/counting_chain` | dataset |
| 2 | fly scan (as actually written) | `/entry/program_name` attribute `counting_chain` | **attribute** |
| 3 | SAXS / WAXS frame | `/entry/counting_chain` | dataset |
| 4 | step scan | `/entry/instrument/bluesky/metadata/counting_chain` | dataset |

Note #2: `saveFlyData.xml` v2.0 writes the fly-scan marker as an *attribute* of
`/entry/program_name`, alongside `config_version`, not as the dataset the
format spec describes. Both spellings are accepted so either works. In the
bundled test file `TestData/FX4Set/flyscan/GC_SRM3600_0131.h5`:

```
/entry/program_name = "saveFlyData.py"
/entry/program_name@config_version  = "2.0"
/entry/program_name@counting_chain  = "FX4"
```

Python reference: `fx4support.detect_counting_chain()`.

### 2.2 Two traps

**Do not branch on field names.** `UPD` and `I0` in the step-scan file kept
their names across the conversion and changed only their units and meaning.
A name-based test silently mis-reduces every step scan.

**Do not branch on `config_version`.** The step-scan file's
`program_name@config_version` is the NeXus writer's schema version and did not
move with this conversion.

### 2.3 Igor sketch

```igor
Function/S FX4_DetectCountingChain(String fileID)
    // Returns "FX4" or "scaler". Absent marker => "scaler".
    // fileID is an open HDF5 file reference (HDF5OpenFile).
    String chain = ""
    // 1) fly scan, documented location
    chain = FX4_TryReadString(fileID, "/entry/flyScan/counting_chain", "")
    if (strlen(chain) == 0)
        // 2) fly scan, as actually written: ATTRIBUTE of /entry/program_name
        chain = FX4_TryReadAttrString(fileID, "/entry/program_name", "counting_chain", "")
    endif
    if (strlen(chain) == 0)
        // 3) SAXS / WAXS frame
        chain = FX4_TryReadString(fileID, "/entry/counting_chain", "")
    endif
    if (strlen(chain) == 0)
        // 4) step scan
        chain = FX4_TryReadString(fileID, "/entry/instrument/bluesky/metadata/counting_chain", "")
    endif
    if (CmpStr(TrimString(chain), "FX4") == 0)
        return "FX4"
    endif
    return "scaler"
End
```

`FX4_TryReadString` / `FX4_TryReadAttrString` must return the default when the
path or attribute is missing rather than erroring — every one of these four
probes will miss on three of the four formats.

---

## 3. USAXS fly scan

Python: `supportFunctions._importFlyscanFX4()`, `_calculatePD_FlyFX4()`,
`_calculatePDErrorFlyFX4()`.
Test file: `TestData/FX4Set/flyscan/GC_SRM3600_0131.h5` (4000 points).

### 3.1 Arrays read

All under `/entry/flyScan/`. Lengths are the raw dataset lengths; see §3.2 for
trimming.

| dataset | shape | units | role |
|---|---|---|---|
| `upd_current` | (4000,) | pA | **detector signal.** Mean current over each PSO interval |
| `I0_current` | (4000,) | pA | **monitor.** Mean current, same gate |
| `AR_PulsePositions` | (8000,) | degrees | analyser angle, zero-padded, one extra leading point |
| `upd_sigma` | (4000,) | pA | spread of samples *within* each interval (not error of mean) |
| `I0_sigma` | (4000,) | pA | same, for I0 |
| `upd_total` | (4000,) | — | sum of raw samples in interval; gives N via `total/mean` |
| `I0_total` | (4000,) | — | same, for I0 |
| `channel_time` | (4000,) | s | per-point measurement time, derived by the writer |
| `sample_time` | (1,) | s | FX4 raw sample period (1e-4 s in the test file) |
| `values_per_read` | (1,) | — | raw readings averaged per FX4 sample (10) |
| `upd_lurange` | (1,) | — | autoranger range **once, after the scan** (4 in test file) |
| `ring_overflows` | (1,) | — | must be 0; non-zero ⇒ dropped samples |
| `ring_overflows_I0` | (1,) | — | same, for I0 |
| `upd_bkg0..4` | (1,) each | pA | dark current per range |
| `upd_bkgErr0..4` | (1,) each | pA | its measured error |

Not present on this chain, and not needed: `mca1`, `mca2`, `mca3`,
`mca_clock_frequency`, `changes_DDPCA300_ampGain`,
`changes_DDPCA300_ampReqGain`, `changes_DDPCA300_mcsChan`,
`DDPCA300_gain0..4`.

Geometry and transmission metadata stay in `/entry/metadata/`:
`detector_distance` (1073 mm), `UPDsize` (4.95 mm), `AR_center`, `ARenc_0`,
`DCM_energy`, `trans_pin_counts`, `trans_I0_counts`, `trans_pin_time`.
Wavelength: `/entry/instrument/monochromator/wavelength` (0.5904016 Å).
Thickness: `/entry/sample/thickness` (mm).

### 3.2 Trimming and alignment — unchanged from the scaler chain

`AR_PulsePositions` is 8k-padded with zeros and carries **one extra leading
point**. Handle exactly as before:

1. Drop the first element.
2. Strip trailing zeros.
3. `n = min(len(upd_current), len(I0_current), len(AR))`.
4. Truncate every array to `n`.

### 3.3 Metadata precedence — this one bites

The FX4 `saveFlyData` configuration writes the counting-chain fields inside
`/entry/flyScan` (sourced from `usxFX4:FX4:seq01`). `/entry/metadata` may still
carry **stale Femto-era copies** of `upd_bkg*` and `upd_bkgErr*` (from
`usxLAX:pd01:seq02`). So precedence differs by field:

| field group | search order |
|---|---|
| `upd_bkg0..4`, `upd_bkgErr0..4`, `upd_amp_change_mask_time0..4`, `upd_lurange`, `sample_time`, `values_per_read`, `ring_overflows*`, `counting_chain` | **`/entry/flyScan` first**, then `/entry/metadata` |
| everything else (geometry, transmission, `AR_center`, …) | **`/entry/metadata` first**, then `/entry/flyScan` |

Getting this backwards subtracts the wrong dark current. Some deployments write
no `/entry/metadata` group at all, which is why both locations are searched in
each direction.

Also set, by construction on this chain: `I0AmpGain = 1.0`, `I0Gain = 1.0`.

### 3.4 Per-point measurement time

Used only for smoothing (§3.8) — **nothing in the FX4 intensity divides by
it**. Read `channel_time` directly; it is already in seconds. The writer
derives it as `N · sample_time` with `N = I0_total / I0_current`.

Fallbacks, in order:
1. `channel_time` if present and long enough.
2. `N · sample_time` with `N = |I0_total / I0_current|`.
3. Uniform 1 ms, with a warning.

Contrast with the scaler chain, where `TimeInSec = mca1 / 1e6`.

### 3.5 Dark current and the range index

The fly-scan file records the autoranger range **only once, after the scan**
(`upd_lurange`). There is no per-point range. Consequences:

* The dark current for that single range is subtracted at **every** point:
  `updBkg[i] = upd_bkg{r}` for all `i`, where `r = upd_lurange`.
* Likewise `updBkgErr[i] = upd_bkgErr{r}`.
* The range-index array handed to the smoother is constant (`= r` everywhere).
* **There is no amplifier dead-time masking.** The scaler chain masked points
  after each gain change using `changes_DDPCA300_mcsChan` and
  `upd_amp_change_mask_time{n}`; the FX4 chain has no equivalent and no
  points are masked or set to NaN.

In the test file `upd_lurange = 4`, so `upd_bkg4 = 8.528 pA` and
`upd_bkgErr4 = 0.0562 pA` are used throughout.

### 3.6 Normalised intensity

```
I[i] = (upd_current[i] − updBkg[i]) / I0_current[i]
```

No gain term. No dwell term. Both electrometers share the same gate, so the
interval cancels in the ratio.

For comparison, the scaler chain computed:
`I = ((UPD − t·bkg)/(f·gain)) / (Monitor/I0Gain)` with `f = 1e6`.

**No negative-intensity fix here.** Weak FX4 points scatter around zero after
dark subtraction, so a lift applied now is set by the single most negative
noise point (0.3–0.7 pA on the 2026-10 data, as large as the net high-Q
signal), and sample and blank get different random offsets. The FX4 chain
applies the same Igor `IN3_FixNegativeIntensities` **after smoothing** (§3.8),
where the minimum follows the real background. The scaler chain still lifts
here, unchanged.

### 3.7 Uncertainty

Measured noise model; the data behind it are in
`docs/fx4-flyscan-uncertainties.md`. Every term depends on **signal level
only** — never on Q, never on photon energy.

```
N_upd = |upd_total / upd_current|      // reads per interval, clamped to >= 1
N_i0  = |I0_total  / I0_current|
vpr   = values_per_read                // 10; 1 if missing (conservative)
t     = channel_time                   // s
q_eff = 6.5e-4                         // pC, FX4_SHOT_CHARGE_PC

net   = upd_current − updBkg
s_upd = sqrt( upd_sigma^2/(N_upd·vpr) + q_eff·max(net,0)/t + updBkgErr^2 )
s_i0  = sqrt( I0_sigma^2 /(N_i0 ·vpr) + q_eff·I0_current/t )

R     = net / I0_current
Error = sqrt( s_upd^2 + R^2·s_i0^2 ) / I0_current
```

* **Electronic floor**, `sigma/sqrt(N·vpr)`: the recorded per-interval spread
  behaves like that of N·values_per_read independent raw samples. Dividing by
  `sqrt(N)` alone overstated weak-signal errors by ~√10.
* **Shot-like term**, `q_eff·net/t`: invisible to `upd_sigma`, dominant from
  ~1 pA net upward. `q_eff` fitted to 5–9e-4 pC at 12, 21, 24 and 28 keV alike
  (E/3.65 eV would change 2.3× over that range), so it is one constant. Diode
  efficiency does not enter; the term is written in measured current.
* The ratio is propagated in the form above, not as `|R|·sqrt(rel_num² +
  rel_den²)`: the relative form gives a point sitting exactly on the dark
  current zero error.

Where a sigma array is absent, fall back to a flat relative error of
**1 %** (`FX4_RELATIVE_CURRENT_ERROR`) of the current. Real files carry a few
NaNs in the sigma arrays; replace those points with the same 1 % fallback —
otherwise they propagate to zero error, which reads as an infinitely
well-known point.

This is the error of **one** interval. Smoothing (§3.8) and rebinning (§6.4)
reduce it for the averaging they do.

### 3.8 Smoothing

Same windows as the scaler chain; four FX4 differences (`fx4_noise_model=True`):

* The time base is **1.0**, because `channel_time` is already in seconds
  (the scaler chain passed mca1 counts with a base of 1e6 Hz).
* The range index is derived from the signal, not `upd_lurange` (§3.10).
* **Average in linear intensity**, not ln(I). Weak FX4 points include
  negatives that ln(I) cannot hold (they were NaN-replaced from neighbours,
  biasing upward), and averaging ln(I) of noisy data is biased low: 11–12 % at
  high Q on a strong ceramic, 1–2 % on glassy carbon. On strong, steep data the
  two differ by < 0.07 %.
* **Propagate the averaging into the error.** For the trapezoid boxcar of n
  points, `Error /= sqrt(n_eff)` with `n_eff = (n−1)^2/(n−1.5)` (the end points
  carry half weight). For the linear-fit branch at range boundaries,
  `Error /= sqrt(n_used)` instead of the scaler chain's fixed `/3`. Keep
  `n_eff` per point as **SmoothWindow** (1 where not smoothed): rebinning needs
  it.
* **Then** apply the negative-intensity fix of §3.6 (Igor
  `IN3_FixNegativeIntensities`, minimum over the last 80 % of points,
  `1.1·|Vmin| + 0.3e-10·max`). It is one flat constant, absorbed by the flat
  background of any model; done after smoothing it follows the background
  level instead of a noise spike.

Per-range minimum smoothing times (seconds, indexed by range 0–4), unchanged:
`[0.02, 0.02, 0.03, 0.1, 0.4]`. At t ≈ 0.019 s that is a 41–47 point window
below 30 pA and 7–15 points above.

### 3.10 Signal-derived smoothing range

The FX4 reading is gain-independent and `upd_lurange` is written **once,
after the scan**, so there is no per-point range: every point would take
`rwave_smooth_times[4]`, the slowest. On the scaler chain the amplifier range
*was* a proxy for signal level, and those smoothing times were tuned against
it, so derive a pseudo-range from the raw `upd_current` and reuse the
constants unchanged (`FX4_SIGNAL_DERIVED_SMOOTHING_RANGE`, on):

| `upd_current` | range index | window |
|---|---|---|
| > 1e6 | 0 | 0.02 s |
| > 12000 | 1 | 0.02 s |
| > 800 | 2 | 0.03 s |
| > 30 | 3 | 0.10 s |
| otherwise | 4 | 0.40 s |

Thresholds are descending on the **raw** `upd_current`, before dark-current
subtraction. A separate "start smoothing only below 1e-5 of the peak" guard is
redundant once this is in.

The two other switches of the first FX4 release are **retired** (2026-11):
"error of the mean on rebin" and "effective samples N/values_per_read". The
first double-counted the smoothing once that is propagated; the second had the
direction wrong (the floor is `N·values_per_read`, §3.7). Their apparent
cancellation on the commissioning pair was a coincidence of its 20–100 pA
signal level; on weak samples the SMR errors came out 3–7× too large.

### 3.11 Igor implementation of the noise model

**Index base differs between the two codes.** Igor's `PD_range` is 1–5 and
indexes `rwave_smooth_times` 1-based; Python's `UPD_gainsIndx` is 0–4. Same
five times, `[0.02, 0.02, 0.03, 0.1, 0.4]`, same order. Everything below is
written in Igor's 1-based convention, and applies **only when the file is FX4**;
the scaler chain stays as it is, apart from §6.3.

**Constants.**

```igor
Constant FX4_ShotChargePC      = 6.5e-4   // pC, energy independent (12-28 keV)
Constant FX4_RelErrorFloor     = 1e-3     // added once to SMR_Error after rebin
```

**1. Per-interval error** (replaces the `upd_sigma/sqrt(N)` step):

```igor
// updWvSigma, updTotal, USAXS_PD = upd_current (raw), Bkg = dark, BkgErr
// I0 = I0_current, I0Sigma, I0Total, MeasTime = channel_time, vpr from file
Duplicate/FREE USAXS_PD, Nupd, NetUPD, SigUPD2, SigI02
Nupd   = max(abs(updTotal/USAXS_PD), 1)
NetUPD = USAXS_PD - Bkg
SigUPD2 = updWvSigma^2/(Nupd*vpr) + FX4_ShotChargePC*max(NetUPD,0)/MeasTime + BkgErr^2
SigI02  = I0Sigma^2/(max(abs(I0Total/I0),1)*vpr) + FX4_ShotChargePC*I0/MeasTime
PD_Error = sqrt(SigUPD2 + (NetUPD/I0)^2*SigI02) / I0
```

**2. Pseudo-range** (switch 3, unchanged):

```igor
Static Function FX4_SignalDerivedRange(WAVE PD_range, WAVE USAXS_PD)
    // USAXS_PD is the RAW upd_current, BEFORE dark-current subtraction.
    Redimension/S PD_range
    PD_range = 5
    PD_range = (USAXS_PD[p] >    30) ? 4 : PD_range[p]
    PD_range = (USAXS_PD[p] >   800) ? 3 : PD_range[p]
    PD_range = (USAXS_PD[p] > 12000) ? 2 : PD_range[p]
    PD_range = (USAXS_PD[p] >   1e6) ? 1 : PD_range[p]
End
```

**3. Smoothing** (`IN3_SmoothRData` FX4 branch): average `R_Int` directly
rather than `log(R_Int)`, divide `R_Error` by `sqrt(n_eff)` where the boxcar
is used (`n_eff = (n-1)^2/(n-1.5)`) and by `sqrt(n_used)` in the linear-fit
branch, store `n_eff` in a wave `SmoothWindow` (1 elsewhere), and call
`IN3_FixNegativeIntensities` **after** the loop instead of before it.

**4. Rebin** (the log-binning routine above Q = 0.0002), with `SmoothWindow`
rebinned alongside as a sum:

```igor
    Rebinned_Sdev[i] = sqrt(Rebinned_Sdev[i] / cntPoints)          // RMS, as now
    if (isFX4)
        Variable meanW = max(Rebinned_W[i] / cntPoints, 1)
        Rebinned_Sdev[i] /= sqrt(max(cntPoints / meanW, 1))        // independent points only
    endif
```

**5. After rebin:** `SMR_Error = sqrt(SMR_Error^2 + (FX4_RelErrorFloor*SMR_Int)^2)`.

**6. SMR_Error and SampleTransmission** — see §6.3; this one affects both chains.

**What to look at when comparing.** Steps 1, 4 and 5 change **only error bars**.
Step 3 changes intensities: by < 0.1 % where the signal is strong, and by a flat
offset (the lift) plus the removal of the ln(I) bias where it is weak. Use the
halves test of `tests/test_fx4_noise.py` (even vs odd raw points, peak pinned to
the full scan's) to check the error bars; target true/quoted ≈ 1.

### 3.9 Ring-overflow check

If `ring_overflows` or `ring_overflows_I0` is non-zero, the FX4 driver's ring
buffer overflowed: the oldest samples in an averaging interval were discarded,
so **every mean in the file is biased toward the end of its interval**. The
data are still readable; the means are not trustworthy. Log a warning — do not
refuse the file. Both are 0 in the bundled test data.

---

## 4. USAXS step scan (uascan)

Python: `convertUSAXS.importStepScan()`, `_CorrectUPDStepFX4()`,
`_calculatePDErrorStepFX4()`, `createUPDGainsAndBkgErrArrays()`.
Test file: `TestData/FX4Set/usaxs/GC_SRM3600_0044.h5` (300 points).

**This is the format most likely to be mis-reduced**, because the array names
did not change — only their units did.

### 4.1 Arrays read

| dataset | shape | FX4 meaning | scaler meaning |
|---|---|---|---|
| `/entry/data/a_stage_r` | (300,) | analyser angle, deg | same |
| `/entry/data/UPD` | (300,) | **mean pA** | counts over `seconds` ticks |
| `/entry/data/I0` | (300,) | **mean pA** | counts |
| `/entry/data/fx4_autorange_lurange` | (300,) | **per-point range 0–4** | *absent* |
| `/entry/data/seconds` | — | **absent** | dwell, 1e7 Hz Joerger ticks |
| `/entry/data/upd_autorange_controls_gain` | — | *present but stale — do not read* | Femto gain |
| `/entry/data/I0_autorange_controls_gain` | — | *stale — do not read* | I0 gain |

Legacy name fallbacks still apply: `I0_USAXS` for `I0`, `PD_USAXS` for `UPD`.

Unlike the fly scan, the step scan **does** record the range at every point
(`fx4_autorange_lurange`), so the dark current is looked up per point.

### 4.2 Dark currents

From the baseline stream, `/entry/instrument/bluesky/streams/baseline/`,
taking element `[0]` of each `value` array:

```
fx4_autorange_ranges_range{i}_background/value        -> Bkg_map[i]       (pA)
fx4_autorange_ranges_range{i}_background_error/value  -> Bkg_err_map[i]   (pA)
```

for `i = 0..4`, defaulting to 0.0 when absent. In the test file only range 4 is
non-zero: background 8.3498, error 0.08216.

> The Femto `upd_autorange_controls_ranges_gain{i}_background` records still
> exist in the baseline but are **stale** on this chain. Reading them is the
> single most likely porting error.

Note the HDF5 `units` attribute on these datasets says `c/s`, left over from
the scaler era. The values are picoamps. Ignore the attribute.

### 4.3 Per-point count time

The FX4 reading is a mean current, so the writer records no `seconds` column.
The count time is reconstructed from the plan arguments **for reporting only** —
nothing in the FX4 reduction divides by it. If you do not need to report it,
you may skip this entirely.

```
count_time  <- regex "count_time:\s*([0-9.eE+-]+)" applied to
               /entry/instrument/bluesky/metadata/plan_args
dynamic     <- /entry/instrument/bluesky/metadata/useDynamicTime == "True"
intervals   <- /entry/instrument/bluesky/metadata/intervals
```

With `useDynamicTime`, `uascan` scales the dwell by thirds across the scan:
`count_time/3` over the first third (`fraction < 0.33`), `count_time` over the
middle, `2·count_time` over the last (`fraction >= 0.66`), where
`fraction = i / intervals`. Without it, `count_time` everywhere. If
`count_time` cannot be parsed, report 1 s per point with a warning.

Test file: `count_time = 0.3`, `useDynamicTime = True`, `intervals = 300`.

### 4.4 Normalised intensity

```
Bckg[i] = Bkg_map[ fx4_autorange_lurange[i] ]      // 0.0 where range unknown

I[i] = (UPD[i] − Bckg[i]) / I0[i]
```

No gain, no dwell. **Two further differences from the scaler path:**

* **No one-point gain shift.** The scaler code duplicated `AmpGain[0]` and
  dropped the last element, to undo the Femto gain readback lagging its own
  range change. The FX4 range is read in the same event document as the
  current, so there is nothing to undo. Do not shift.
* **No background × time.** The scaler path multiplied the background by
  `TimePerPoint/1e7`. The FX4 dark current is already a current in pA.

For reference, the scaler path was:
`I = ((UPD − bkg·t/1e7) · I0gain) / (AmpGain · I0)`.

### 4.5 Gains and background-error arrays

```
UPD_gains[i]  = 1.0                                   // gain-independent chain
UPD_bkgErr[i] = Bkg_err_map[ fx4_autorange_lurange[i] ]    // pA, no dwell factor
```

### 4.6 Uncertainty

The uascan file records **no per-point sigma** for the electrometer currents —
the counting statistics the scaler chain relied on simply are not available. So:

```
sigma_upd = sqrt( (0.01 · UPD)^2 + UPD_bkgErr^2 )
sigma_i0  = 0.01 · |I0|
Error     = ratio_error(UPD − Bckg, sigma_upd, I0, sigma_i0)
```

using the same `ratio_error` as §3.7. The 0.01 is
`FX4_RELATIVE_CURRENT_ERROR`, a placeholder (§7). This sets error bars only;
intensities are untouched.

### 4.7 Transmission inputs

Unchanged in location from the scaler chain — baseline stream, element `[0]`:

| baseline dataset | variable | test value |
|---|---|---|
| `terms_USAXS_transmission_diode_counts/value` | `trans_pin_counts` | 16578155.149 |
| `terms_USAXS_transmission_diode_gain/value` | `trans_pin_gain` | 1.0 |
| `terms_USAXS_transmission_I0_counts/value` | `trans_I0_counts` | 67729.942 |
| `terms_USAXS_transmission_I0_gain/value` | `trans_I0_gain` | 1.0 |
| `terms_USAXS_transmission_count_time/value` | `trans_I0_time` | 2.0 |

The FX4 step-scan writer still emits both gains, as exactly 1.0. Note the
diode stream maps to `trans_pin_*` and the I0 stream to `trans_I0_*` — an
earlier Python version had these swapped, which computed the **reciprocal**
transmission. Get this right.

### 4.8 Other metadata

| quantity | path | test value |
|---|---|---|
| `detector_distance` | `/entry/instrument/bluesky/metadata/SDD_mm` | 1073 mm |
| (SAD) | `/entry/instrument/bluesky/metadata/SAD_mm` | 228 mm |
| `UPDsize` | baseline `terms_USAXS_diode_upd_size/value[0]` | 4.95 mm |
| thickness | `/entry/instrument/bluesky/metadata/sample_thickness_mm` | 1.0 mm |
| wavelength | `/entry/instrument/monochromator/` | |
| timestamp | `/entry/start_time` | |

### 4.9 Known quirk in the commissioning step-scan data

In `TestData/FX4Set/usaxs/`, `I0` is **constant across all 300 points** —
`fx42` was not re-read per point in that scan. The reduction is unaffected (it
is a division by a constant) but the monitor normalisation is doing nothing.
Also, TRD shares one `Range` with UPD on `usxFX4`, so the transmission terms
give T = 1.35 for this sample/blank pair, which is unphysical. **Do not use
this particular pair to validate transmission**; use the fly-scan pair, which
gives a correct 0.936.

---

## 5. SAXS and WAXS area-detector frames

Python: `convertSWAXS.importADData()`, `_fx4Monitor()`, `calibrateAD2DData()`,
`_correction_factor()`.
Test files: `TestData/FX4Set/saxs/`, `TestData/FX4Set/waxs/`.

Marker: `/entry/counting_chain` = `"FX4"` (a dataset at the entry level).

The 2D image, the geometry, the mask, the pyFAI/Nika azimuthal integration and
the solid-angle correction are **all unchanged**. Only the monitor and
transmission terms change, plus one caveat on absolute scale.

### 5.1 The calibration formula — unchanged in form

```
solidAngle = pixel_size^2 / detector_distance^2
preFactor  = I_scaling / I0s / thickness_cm / solidAngle

calib2D = preFactor · ( sample2D / transmission  −  (I0s/I0b) · blank2D )
```

with

```
transmission = ((sTRdiode/sTRdiodeGain) / (sTRI0/sTRI0gain))
             / ((bTRdiode/bTRdiodeGain) / (bTRI0/bTRI0gain))

I0s = sampleI0 / sampleI0gain
I0b = blankI0  / blankI0gain
```

`thickness_cm = thickness_mm · 0.1`. All metadata is under `/entry/Metadata/`.

### 5.2 Monitor — both SAXS and WAXS

```
sampleI0     = I0_cts_gated
sampleI0gain = 1.0          // gain-independent
```

`I0_cts_gated` is the sum of FX4 samples over the software-gated exposure:

```
I0_cts_gated = I0_current [pA] · Exp_time_gated [s] / FX4_SampleTime [s]
```

`/entry/control/integral` is documented as a hardlink to it and may be read as
a fallback, **but it is absent from the commissioning files** — only
`/entry/control/mode` exists there. Prefer `I0_cts_gated` and treat the
control-group path as optional.

> `I0_cts` and `I0_gain` are still written on FX4 frames but are **stale**
> scaler/Femto records. Do not read them when the chain is FX4. On a
> commissioning frame `I0_cts/I0_gain = 1.3e-2` while `I0_cts_gated = 1.1e8`.

### 5.3 Transmission — SAXS

The SAXS transmission still comes from the pin-diode measurement made before
the frame. Field names are unchanged; on FX4 both values are pA and both
recorded gains are exactly 1:

```
sTRdiode = Pin_TrPD        sTRdiodeGain = Pin_TrPDgain   (= 1.0)
sTRI0    = Pin_TrI0        sTRI0gain    = Pin_TrI0gain   (= 1.0)
```

and the same four from the blank. Test values: `Pin_TrPD = 2473771.35`,
`Pin_TrI0 = 17174.90`, both gains 1.0 → transmission 0.946.

### 5.4 Transmission — WAXS

```
sTRdiode = TR_current      sTRdiodeGain = 1.0
sTRI0    = I0_current      sTRI0gain    = 1.0
```

and the same two from the blank.

> **Use `TR_current`, not `TR_cts_gated`.** `TR_cts_gated` is the leftover
> `Total_RBV` of whatever `fx4` last acquired — in practice the 0.05 s
> autoscale read — not the exposure. `TR_current` and `I0_current` are both
> means in pA over the same gate, so the exposure time cancels in the ratio.
> This is a real trap: `TR_cts_gated` is present and looks plausible.

Test values: `TR_current = 11651372.26`, `I0_current = 54912.71` →
transmission 0.558.

Contrast with the scaler WAXS path, which used `TR_cts` / `TR_gain` against
`control/integral` / `I0_gain`.

### 5.5 Field summary

| role | scaler SAXS | scaler WAXS | **FX4 SAXS** | **FX4 WAXS** |
|---|---|---|---|---|
| monitor | `I0_cts` | `control/integral` | `I0_cts_gated` | `I0_cts_gated` |
| monitor gain | `I0_gain` | `I0_gain` | 1.0 | 1.0 |
| trans. diode | `Pin_TrPD` | `TR_cts` | `Pin_TrPD` | `TR_current` |
| trans. diode gain | `Pin_TrPDgain` | `TR_gain` | `Pin_TrPDgain` (1.0) | 1.0 |
| trans. I0 | `Pin_TrI0` | monitor | `Pin_TrI0` | `I0_current` |
| trans. I0 gain | `Pin_TrI0gain` | monitor gain | `Pin_TrI0gain` (1.0) | 1.0 |

SAXS vs WAXS is distinguished by the presence of `pin_ccd_tilt_x` in
`/entry/Metadata` (present ⇒ SAXS).

### 5.6 Chain mismatch between sample and blank

A sample counted on one chain and a blank on the other are **not comparable**:
one is in counts, the other in picoamps. Detect the chain for both files and
raise an error rather than producing a silently wrong result.

### 5.7 Ring overflows

`/entry/Metadata/FX4_RingOverflows` — same meaning and same warn-don't-refuse
treatment as §3.9.

### 5.8 Absolute intensity — resolved at the beamline (2026-11)

**Status:** `I_scaling` has been recalibrated against SRM 3600 for the FX4
monitor and stored in every instrument energy profile, so FX4 frames written
since then carry the correct constant and both codes simply use the recorded
`I_scaling`. No warning is emitted any more. The commissioning frames in
`TestData/FX4Set/` predate the recalibration and are **not** on absolute
scale. The background below is kept for reference.

`I_scaling` (`/entry/Metadata/I_scaling`) is the SAXS/WAXS absolute-intensity
constant. Its value was determined against the **old** monitor, `I0_cts/I0_gain`
— a V/F count divided by a Femto gain. The FX4 chain normalises by
`I0_cts_gated` instead, a sum of picoamp samples, which is a completely
different unit. The recorded `I_scaling` therefore under-scales FX4 frames by
roughly **10 orders of magnitude**.

The curve **shape is unaffected** — diode/I0 is a ratio and the error is common
to every pixel — so only the absolute level is wrong. Until a new constant is
determined from glassy carbon SRM 3600, reduce FX4 frames with the recorded
`I_scaling`. `convertSWAXS.FX4_I_SCALING` remains as an optional override
hook (e.g. to re-reduce commissioning data); it is `None` in normal use.

---

## 6. Downstream pipeline — unchanged, listed for completeness

Everything from here on is chain-independent and should already exist in the
Igor code. It is listed so the port can confirm the FX4 front end hands over
the same quantities.

### 6.1 Beam-centre correction and Q

Fit a modified Gaussian to the top of the rocking curve
(`a·exp(−(|x−x0|/(2σ))^p)`), using points above `max/2.3`; compute FWHM from
the half-maximum crossings of the fitted curve evaluated over points above
`max/5`. Then

```
Q = −(4π · sin(radians(AR − x0)/2)) / wavelength
```

Outputs: `Q`, `Center` (x0), `Maximum` (amplitude), `FWHM`, `wavelength`.

### 6.2 Transmission terms for USAXS

```
MeasuredTransmission = ((SaPin/SaPinGain)/(SaI0/SaI0Gain))
                     / ((BlPin/BlPinGain)/(BlI0/BlI0Gain))
```

On FX4, all four gains are **1.0 by construction**. In the fly-scan format
`saveFlyData.xml` v2.0 stops writing `trans_pin_gain` / `trans_I0_gain`
altogether — their absence is **normal, not an error**. (The step-scan writer
still emits them, as 1.0.) If a gain *is* present and is not 1.0 on an FX4
scan, warn and use 1.0 anyway.

A missing *count* is a real problem — the transmission cannot be measured —
so fall back to T = 1.0 with a warning.

```
PeakToPeakTransmission = SamplePeakMax / BlankPeakMax
MSAXSCorrection        = MeasuredTransmission / PeakToPeakTransmission
```

### 6.3 Subtraction, Qmin selection, calibration

Log-interpolate the blank onto the sample Q grid, subtract, propagate errors in
quadrature, and form `IntRatio = |I_sample / I_blank_interp|`. Lift negative
intensities (`1.1·|min|` over the second half, plus `3e-11·max`). Choose the
start index as `max(indexSample, indexBlank, indexRatio)` where the first two
come from the sample and blank FWHM converted to Q, and `indexRatio` is the
first point where `IntRatio / QCorrection > 1.05` with

```
QCorrection = min(1 + (|QminBlank/Q|)^3, 2)
```

Then

```
slitLength  = 0.5 · ((4π)/wavelength) · sin(UPDsize/(2·SDD))
OmegaFactor = (UPDsize/SDD) · radians(FWHM_blank)
Kfactor     = BlankPeakMax · OmegaFactor · thickness_cm

SMR_Int   = SMR_Int   / (Kfactor · MSAXSCorrection)
SMR_Error = SMR_Error / (Kfactor · MSAXSCorrection)
```

**No `· PeakToPeakTransmission` on SMR_Error (both chains, since 2026-11).**
Igor's 2/2014 line `SMR_Error *= SampleTransmission` undoes an extra
`/SampleTransmission` earlier in Igor's own formula
(`SMR_Error = sqrt(R_error^2/T^2 + BL^2)/Kfactor`, `IN3_SupportFnct.ipf`).
Matilda divides by PeakToPeakTransmission only once (in
`normalizeByTransmission`), so carrying the multiplication over shrank the
*relative* error by T_pp: ~10 % on ordinary samples, 18× on a strongly
multiple-scattering one (T_pp = 0.055, MSAXS = 7). Igor itself also divides
SMR_Error by `Kfactor` only, not `Kfactor·MSAXSCorrection`, and has a `/3` from
12/2013; revisit both when porting.

Note the K-factor is derived from the blank's **own** peak, not from a stored
constant — which is why USAXS absolute calibration carried across the FX4
conversion unchanged, while SAXS/WAXS (§5.8) did not.

### 6.4 Rebinning and desmearing

Fly scan: rebin to 500 points when more than 800 survive (log binning above
Q = 0.0002, linear below), then desmearing with `slitLength` and
`PowerLaw w flat` extrapolation from Q = 0.15. The default method is now the
truncated-Abel inversion (`docs/desmearing-methods.md` §4a); Lake/Strobl with
20 iterations is still selectable.
Step scan: **no rebinning** (if many points were collected, they are wanted),
then the same desmearing.

On FX4 a bin of n points keeps the RMS of its inputs divided by
`sqrt(max(n / mean SmoothWindow, 1))` — the number of *independent* points it
holds; smoothed neighbours share their raw data. Then
`SMR_Error = sqrt(SMR_Error² + (1e-3·SMR_Int)²)`: ~0.1 % of I survives all
averaging on strong signal. The scaler chain keeps the plain RMS.
This matters downstream: the Abel desmearing picks its smoothing width from a
chi^2 against `SMR_Error`, so an inflated error makes it over-smooth.

---

## 7. Deliberately uncalibrated constants

One number is still a placeholder (`FX4_RELATIVE_CURRENT_ERROR`); `FX4_I_SCALING`
is resolved by recalibrating `I_scaling` at the beamline and kept only as an
override hook.

| constant | Python location | current | what it affects |
|---|---|---|---|
| `FX4_I_SCALING` | `convertSWAXS.py` | `None` ⇒ use recorded `I_scaling` (recalibrated, correct since 2026-11) | override hook only (§5.8) |
| `FX4_RELATIVE_CURRENT_ERROR` | `fx4support.py` | `0.01` | **Error bars only**, never intensities (§3.7, §4.6) |

The fly-scan noise constants (`FX4_SHOT_CHARGE_PC`, `FX4_RELATIVE_ERROR_FLOOR`)
are **measured**, not placeholders: see §3.7 and
`docs/fx4-flyscan-uncertainties.md`.

---

## 8. Validation targets

Reference numbers produced by the current Python code on `TestData/FX4Set/`
(glassy carbon SRM 3600 `_0131`/`_0044` against blank `_0130`/`_0043`). An Igor
implementation following this document should reproduce them.

**Fly scan** — the primary validation case:

| quantity | value |
|---|---|
| MeasuredTransmission | 0.93616 |
| PeakToPeakTransmission | 0.93737 |
| MSAXSCorrection | 0.99870 |
| Blank peak maximum | 1607.994 |
| Blank FWHM | 4.7668e-4 deg |
| OmegaFactor | 3.83808e-8 |
| Kfactor | 6.17160e-6 |
| slitLength | 0.0245475 Å⁻¹ |
| **desmeared plateau, 0.01 < q < 0.1** | **29.71 cm²/cm³** (certified ≈ 30) |

A repeat pair (`_R_0132`/`_R_0133`, not bundled) agrees to ~1 % over the same
range. The plateau rose from 29.52 to 29.71 with the 2026-11 noise changes: the
blank is no longer lifted before smoothing, a flat +0.17 cm⁻¹ at every Q.

**Error bars:** halves test on this pair, true/quoted per Q band
1e-3 / 3e-3 / 1e-2 / 3e-2 / 0.1 → 1.05 / 0.82 / 0.72 / 0.96 / 0.77.

**Step scan:** Kfactor 1.42487e-5, slitLength 0.0245475. Transmission comes out
1.35 — unphysical, and a property of the commissioning data, not the code
(§4.9).

**SAXS:** transmission 0.94607. Absolute level provisional (§5.8).
**WAXS:** transmission 0.55833. Absolute level provisional (§5.8).

---

## 9. Porting checklist

1. Chain detection: all four marker locations, including the fly-scan
   **attribute** spelling; absent ⇒ scaler (§2).
2. Fly scan: read `upd_current`/`I0_current`; trim `AR_PulsePositions`
   (drop leading, strip trailing zeros); `/entry/flyScan` wins for
   counting-chain metadata (§3.3).
3. Fly scan: single `upd_lurange` for dark current everywhere; **no dead-time
   masking** (§3.5).
4. Step scan: `fx4_autorange_lurange` per point; read dark currents from
   `fx4_autorange_ranges_*`, **not** `upd_autorange_controls_ranges_*` (§4.2).
5. Step scan: **no** one-point gain shift, **no** background × time (§4.4).
6. Both USAXS: `I = (UPD − dark)/I0`, no gain, no dwell.
7. SAXS/WAXS: monitor is `I0_cts_gated`, gain 1 (§5.2).
8. WAXS: transmission uses `TR_current`, **not** `TR_cts_gated` (§5.4).
9. All: ignore `I0_cts`, `I0_gain`, `TR_cts`, `TR_gain`, Femto gains and
   `DDPCA300_*` on FX4 files — present but stale.
10. Refuse sample/blank pairs from different chains (§5.6).
11. Warn on non-zero ring overflows; do not refuse (§3.9, §5.7).
12. Carry both uncalibrated constants as overridable globals (§7).
