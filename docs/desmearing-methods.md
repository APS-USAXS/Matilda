# Desmearing methods — developer reference

Technical reference for the selectable slit-desmearing methods added in
`matilda/desmearing_methods.py` (branch `feature/desmearing-methods`). Read this
first when you need to change, fix, or extend the desmearing options.

Companion document: `docs/desmearing-paper.md` (the theory + citations for the
instrument paper). This file is the "how the code works and where to touch it"
reference.

---

## 1. What was added and why

Matilda historically desmeared with a single algorithm, the **Lake** iterative
method (`matilda/desmearing.py`, `desmearData`). It works well on strong data
but amplifies noise and, on weak / near-zero / over-subtracted intensities,
produces non-physical spikes and negative values. A batch test over 80 random
datasets found this on ~39% of files (see the paper doc for numbers).

The new module adds **Gaussian-process (GP, "Huang-style") desmearing** as
selectable alternatives, chosen from the GUI. They stay smooth and positive on
weak data and return a posterior uncertainty. **Lake remains the default in all
script calls**, so batch/automated reduction is unchanged; only the GUI exposes
the new options today.

---

## 2. File map / touch points

| File | Change |
|------|--------|
| `matilda/desmearing_methods.py` | **New.** Slit forward operator, GP core, `desmear_dispatch`, GUI string→key maps. numpy-only. |
| `matilda/desmearing.py` | Unchanged. `desmearData` (Lake) is called by the dispatcher for `method='lake'`. |
| `matilda/convertFlyscan.py` | `processFlyscan` gains `desmear_method`, `gp_length_scale`, `gp_kernel` kwargs (default `'lake'`, `0.5`, `'matern32'`); the `desmearData(...)` call became `desmear_dispatch(...)`. |
| `matilda/convertUSAXS.py` | Same change for `processStepscan`. |
| `matilda/gui/data_reduction/reduction_worker.py` | Forwards the three new params from the GUI `params` dict to both converters. |
| `matilda/gui/data_reduction/parameter_tabs.py` | USAXS tab (`_USAXSTab`, shared by Flyscan+StepScan): Method / GP kernel / GP length-scale widgets + `get_params()` keys; imports `GUI_METHOD_MAP`, `GUI_KERNEL_MAP`. |

Data flow (unchanged structure, one extra hop):

```
GUI (_USAXSTab.get_params) → params dict
   → reduction_worker (forwards desmear_method/gp_length_scale/gp_kernel)
      → processFlyscan / processStepscan
         → desmear_dispatch(method=…)          ← NEW single entry point
            → desmearData (Lake)   OR   _gp_core (GP)
```

---

## 3. Public API — `desmear_dispatch`

```python
desmear_dispatch(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=None, *,
                 method="lake", length_scale_decades=0.5, kernel="matern32",
                 sigma_log=4.0, extrap_method="PowerLaw w flat",
                 extrap_qstart=None, max_iter=20)
    -> (DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ)
```

Same call signature and 4-tuple return as `desmearData`, so it is a drop-in
replacement at the call site. Returns `(None, None, None, None)` on empty input
or missing slit length (matches `desmearData`).

`method` keys:

| key | meaning |
|-----|---------|
| `lake` | Calls `desmearData` unchanged. **Bit-for-bit identical** to the old path (regression-verified). |
| `lake_smooth` | Lake, then geometric-mean smoothing (window = `length_scale_decades`). Best for very noisy flat-plateau data where the GP would ring; no band. |
| `gp` | GP desmearing, kernel from `kernel` arg. Resolution-aware prior on by default. |
| `gp_rbf` | GP desmearing forced to the RBF kernel (shortcut). |
| `gp_lake_mean` | GP anchored on a smoothed Lake curve as prior mean (robust low-q edge). |

For GP methods `DSM_Error` is the **posterior 1σ** (in intensity units); `DSM_dQ`
is `SMR_dQ` carried through (interpolated if the point count changed).

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

`_USAXSTab` (in `parameter_tabs.py`) builds three widgets in the Desmearing
form, after "Extrap Q start":

- `_desmear_method` (QComboBox) — items from `GUI_METHOD_MAP.keys()`:
  `"Lake"`, `"Huang GP"`, `"Huang GP + Lake mean"`.
- `_gp_kernel` (QComboBox) — `"Matérn-3/2"`, `"RBF"`.
- `_gp_length_scale` (QDoubleSpinBox) — 0.1–2.0 decades, default 0.5.

`_on_desmear_method_changed` enables the two GP widgets only when the method name
contains "GP". `get_params()` maps the display strings to keys via
`GUI_METHOD_MAP` (→ `desmear_method`) and `GUI_KERNEL_MAP` (→ `gp_kernel`) and
adds `gp_length_scale`. Those keys are forwarded verbatim by `reduction_worker`.

To add a method to the dropdown: add an entry to `GUI_METHOD_MAP` in
`desmearing_methods.py` (display string → `(method_key, kernel_key)`) and handle
`method_key` in `desmear_dispatch`. No GUI code change needed.

---

## 6. Invariants to preserve

- **Script default is Lake.** `desmear_method='lake'` is the default of both
  converters and every `params.get(..., "lake")` in the worker. Do not change
  these until the default is deliberately promoted.
- **`method='lake'` must stay bit-for-bit identical to `desmearData`.** The
  dispatcher just forwards to it. Keep it that way (there is/should be a
  regression test asserting equality).
- **4-tuple contract.** Every branch returns `(Q, I, Error, dQ)`; all callers
  unpack exactly four values.

---

## 7. Known limitations / TODO

- **No provenance yet.** The desmeared NXcanSAS entry does not record which
  method/params produced it. Re-processing overwrites it (SMR is retained) but
  the method is not labelled. TODO: write `desmear_method`, `gp_length_scale`,
  `gp_kernel` as attributes on the desmeared `sasdata` group in
  `hdf5code.saveNXcanSAS`.
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
- Sandbox: the full method development, batch runner (`batch_compare.py`),
  negative-injection test, and figures live in the separate **`Better desmearing`**
  project; `desmear_lab.huang_gp` there is the reference implementation this
  module was vendored from.
