# Matilda Changelog

All notable changes to this project are documented here.
Branch `fix/phase1-phase3-robustness` (2026-07-08 →).

---

## [Unreleased] — branch `fix/phase1-phase3-robustness`

### Fixed — 2026-07-09
- **ASCII export crash** (`ascii_exporter.py`): `_read_group_data` used Python `or`
  between two `_read_arr()` calls to pick the dQ array (`Qdev` or `dQw`). When
  `_read_arr` returns a multi-element numpy array, Python's `or` operator calls
  `bool(array)`, raising *"The truth value of an array with more than one element
  is ambiguous."* Fixed by replacing the `or` with an explicit `None` check.
  Also added `_hdf5_str()` helper to safely convert HDF5 attributes (which can
  return `np.ndarray`, `np.bytes_`, `bytes`, or `str` depending on file vintage
  and h5py version) to plain Python strings before comparison, preventing the
  same class of error in `_find_nxcansas_groups`. Thickness attribute is now
  normalised to `float` or `""` at read time.
  Worker now logs full tracebacks via `logging.exception` instead of bare
  `logging.warning` so future failures are easier to diagnose.

---

### 2026-07-08 — Phase 4: test suite, ruff lint, CI

- **pytest suite** (`tests/`, 50 tests, 48 pass / 2 skipped pending pyFAI):
  - `test_support_functions.py` — blank matching, subtraction, rebinning,
    smoothing with range indices, filename parsing.
  - `test_desmearing.py` — return signatures, synthetic desmear run, all four
    extrapolation methods.
  - `test_hdf5code.py` — None handling, tuple bug regression, NXcanSAS
    round-trip, malformed files.
  - `test_matilda_helpers.py` — filename regex, FIFO tracker, partner matching
    (skips without pyFAI).
  - `test_technique_detector.py` — `/entry/Metadata` path, folder fallbacks.
  - `test_smoke_reduction.py` — full flyscan and step-scan end-to-end on
    `TestData/TestSet` copies; asserts finite calibrated I(Q) and physically
    sensible transmission (0 < T < 1).
- **ruff** added to `pyproject.toml`; ~50 unused imports removed (incl. Qt
  imports in both PySide6/PyQt6 branches), duplicate `beamCenterCorrection`
  import, 12 unused variables. `rebinData` import moved to its real home in
  `supportFunctions.py`.
- `tifffile` removed from dependencies.
- **GitHub Actions CI** (`.github/workflows/ci.yml`): ruff + pytest on Python
  3.11 and 3.12 for pushes and PRs.

---

### 2026-07-08 — Phase 2: science fixes (validated against Igor)

> ⚗️ These fixes change numerical results. Validated against Igor
> Pro reference data by JIL on 2026-07-08.

- **1.1 Duplicate `"UPD_gains"` key** (`supportFunctions.py`): `calculatePD_Fly`
  now returns both `UPD_gainsIndx` (range index 0–4) and `UPD_gains` (gain
  values). `smooth_r_data` receives the index array. Previously the duplicate key
  silently dropped the index array and every point was smoothed with the 0.4 s
  time constant (range 4). Index mapping corrected from 1-based to 0-based.
- **1.2 Wrong array for requested-gain tail** (`supportFunctions.py`): tail of
  `AmpGainReq_array` now filled from `AmpReqGain` (was `AmpGain`).
- **1.3 Step-scan transmission pin ↔ I0 swap** (`convertUSAXS.py`):
  `importStepScan` now correctly maps `trans_pin_*` ← diode stream and
  `trans_I0_*` ← I0 stream. Previously `MeasuredTransmission` was the
  reciprocal of the intended value.
- **1.4 Blank cache never invalidated for step scans** (`convertUSAXS.py`):
  `getBlankStepscan` now receives the caller's `recalculateAllData` flag
  (was hardcoded `False`).

**Validation results (2026-07-08, JIL):**
- Transmissions agree with Igor (validates 1.3).
- Reduced calibrated data agree with Igor end-to-end.
- Error estimates agree; error bars match measurement noise.
- Blank R intensity is one decade higher in Igor than Matilda — confirmed to be
  a prefactor-convention difference (different V-to-F clock constants for MCA vs
  Joerger scaler geometry; see below). This difference cancels in calibrated
  data because Kfactor is anchored to the blank peak. No action required.

**Frequency constants confirmed correct (2026-07-08, JIL):**
Flyscans use the MCA analyser with a dedicated 1 MHz clock; step scans use the
Joerger scaler internal 10 MHz clock. Both constants are correct for their
geometry and must not be unified. Code comments updated at all sites
(`supportFunctions.py`, `convertUSAXS.py`). Diagnostic script
`TestData/check_clock_frequency.py` can re-verify from data.

---

### 2026-07-08 — Phase 1: crash fixes

- **2.1** `desmearData` failure path now returns 4 values (was 5, callers unpack 4).
- **2.2** `oneDesmearIteration` no longer returns bare `1`; `ExtensionFailed` flag
  now truthful (True only when the flat fallback is also impossible).
- **2.3** `find_crossing_index` returns `int` (was `float`, broke `range()` call).
- **2.4** `reduceFlyscanToQR` guards group deletion with existence check.
- **2.5** `processUSAXSFolder` skips missing technique folders with a warning
  instead of raising `IndexError`; helper `_list_sorted_files` added.
- **2.6** `forceFirstBlank` with empty blank list now logs one clear error and
  returns early instead of silently failing per-scan.
- **2.8** `save_dict_to_hdf5` skips `None` values, uses `require_group`,
  overwrites existing datasets.
- **2.9** Multiple edge-case guards:
  - `subtract_data`: `Y2_min == 0` guard (avoids `log(0)`).
  - `calibrateAndSubtractFlyscan`: `nanmax`/`nanmin` (NaN-safe).
  - `rebin_QRSdata`: `< 2`-point guard; threshold changed to `>=`.
  - `desmearData`: guard for `endme == 0` / NaN.
  - `readGenericNXcanSAS`: defensive rewrite — returns `None` on malformed files.
  - `saveNXcanSAS` / `readMyNXcanSAS`: `None`-attribute guards.
- **1.5** `readMyNXcanSAS` trailing-comma `(None,)` tuple bug fixed.
- **1.6** `extract_number_from_filename` regex now matches `.h5`, `.hdf`, `.hdf5`,
  `.nxs` extensions.

---

### 2026-07-08 — Phase 3: robustness

- **2.7** Duplicate Tiled filter keys removed; `title` (FindScanDataByName) and
  `hdf5_path` (FindLastBlankScan) now matched client-side; module docstring
  documents the one-condition-per-filter-type limitation.
- **N+1 HTTP** fixed: `exit_status:stop.exit_status` added to `select_metadata`;
  `convert_results` falls back to per-uid request only if absent.
- **2.10** `technique_detector` now checks `/entry/Metadata` (with Bluesky path
  as fallback).
- **FIFO eviction**: `_remember_file()` helper + ordered dicts replace
  `set.pop()` (which removed an arbitrary element). Single bound
  `_MAX_TRACKED_FILES = 100`.
- **Logging**: setup moved into `_setup_logging()` called from `main()` —
  importing `matilda.matilda` no longer creates directories or reconfigures the
  root logger.
- **`imagePath`** overridable via `MATILDA_IMAGE_PATH` env var (`"none"` disables
  image saving).
- **Docstring drift** fixed: 15 s → 5 s polling; `extrap_qstart` 0.15
  everywhere; desmear `MaxNumIter` doc 50 → 20.
- **`_find_matching_partner`** disambiguates multiple matches by full stem and
  warns on ambiguity.
- **`Bkg_map`/gain matching** now uses `np.isclose` (rtol 1e-3) with warnings
  on unknown gains (both `createUPDGainsAndBkgErrArrays` and
  `CorrectUPDGainsStep`).
- **Path handling**: `os.path.join` replaces `path + "/" + filename` in
  `importFlyscan`, `importStepScan`, `ImportAndReduceAD`, `reduceFlyscanToQR`.
- **`thickness_override`** file mutation documented with NOTE comments in all
  three converters.
- Stray `print(Totaltime)` → `logging.debug`; stale TODOs in `readfromtiled`
  removed.
