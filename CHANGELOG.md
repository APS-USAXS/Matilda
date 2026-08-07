# Matilda Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased] — 0.3.0.dev0

### Added
- **Selectable desmearing methods** (`desmearing_methods.py`): `desmear_dispatch()`
  offers `lake` (default, unchanged behaviour), `gp` (Huang Gaussian-process,
  Matérn/RBF, resolution-aware prior, posterior 1σ returned as `DSM_Error`) and
  `gp_lake_mean`; the converters take `desmear_method` / `gp_length_scale` /
  `gp_kernel`, the reduction GUI exposes them. See `docs/desmearing-methods.md`.

### Changed
- **Single Qt import point, `matilda/gui/_qt.py`** (house standard, matches
  pyirena and MailToVault). Nine `try: from PySide6 … except ImportError: from
  PyQt6 …` blocks across the GUI package — plus two inline ones in the middle
  of functions — collapse into one shim that also normalises `Signal`
  (`pyqtSignal` under PyQt6). GUI modules now do `from .._qt import QWidget, Qt,
  Signal`; changing binding is a one-file edit. Behaviour is unchanged: PySide6
  first, PyQt6 fallback, same names. `tests/test_gui_qt_shim.py` fails on a new
  direct binding import, and on any Qt or `matilda.gui` import from the daemon
  path — by scanning source, so it works on a machine with no Qt installed.
- **Build backend is `setuptools>=77` with a static `version`** in
  `pyproject.toml`, replacing `hatchling` + `hatch-vcs`. This matches the house
  standard used by pyirena/MailToVault and removes a real failure mode: the
  version string is stamped into every NXcanSAS file by `hdf5code.py`, and a
  git-derived version reported garbage from a shallow CI checkout or a source
  tarball with no tags. Bump `version` in `pyproject.toml`, then tag `v<version>`.
- Wheel contents now come from `[tool.setuptools.packages.find]`
  (`include = ["matilda*"]`); `LICENSE.txt` ships via PEP 639 `license-files`.
- `IMPROVEMENT_PLAN.md` renamed to `PLAN.md` — the plan-file name used across
  the other USAXS Python repos.
- The `>=3.11` floor (house baseline is `>=3.10`) is now annotated in
  `pyproject.toml` with the reason: server and CI run 3.11/3.12 only.
- `build/`, `dist/` and `*.egg-info/` are gitignored and excluded from ruff —
  setuptools leaves them in the tree where hatchling did not.

### Fixed
- Dead `from .desmearing import desmearData` imports in `convertFlyscan.py` and
  `convertUSAXS.py`, left over when the dispatcher replaced the direct call.
  These were failing `ruff check` (F401) in CI. Module docstrings updated to
  name `desmear_dispatch()`.

---

## [0.2.0] — 2026-07-10

### Added
- **Live tune plots** (`convertTune.py`): new module with `getTuneResults()` that
  searches Tiled/DataBroker for the last N tune scans and derives axes from
  metadata (no per-plan hardcoding). `FindLastTuneScans()` and
  `tiled_get_primary_data()` added to `readfromtiled.py`; `plotTuneResults()`
  added to `plotData.py`. The service loop now publishes `tune_ar.jpg`,
  `tune_mr.jpg`, and `tune_a2rp.jpg` alongside SAXS/WAXS/USAXS images.
  Failure-tolerant: offline or dead-server installations never crash the loop.
- **pytest suite** (54 tests, all pass):
  - `test_support_functions.py` — blank matching, subtraction, rebinning,
    smoothing with range indices, filename parsing.
  - `test_desmearing.py` — return signatures, synthetic desmear run, all four
    extrapolation methods.
  - `test_hdf5code.py` — None handling, tuple bug regression, NXcanSAS
    round-trip, malformed files.
  - `test_matilda_helpers.py` — filename regex, FIFO tracker, partner matching.
  - `test_technique_detector.py` — `/entry/Metadata` path, folder fallbacks.
  - `test_smoke_reduction.py` — full flyscan and step-scan end-to-end;
    asserts finite calibrated I(Q) and physically sensible transmission (0 < T < 1).
  - `test_readfromtiled.py` — URI construction and one-condition-per-filter-type
    constraint.
- **GitHub Actions CI** (`.github/workflows/ci.yml`): ruff + pytest on Python
  3.11 and 3.12 for all pushes and pull requests.
- **ruff** linting added to `pyproject.toml`; ~50 unused imports cleaned;
  `rebinData` import moved to its correct home in `supportFunctions.py`.
- `MATILDA_IMAGE_PATH` environment variable: overrides default image-output
  directory; set to `"none"` to disable image saving entirely.
- `_setup_logging()` helper: logging setup moved out of module-level code so
  `import matilda.matilda` no longer creates directories or reconfigures the
  root logger.
- `_list_sorted_files()` helper in `matilda.py` for sorted folder listings.
- Diagnostic script `TestData/check_clock_frequency.py` to verify V-to-F clock
  constants from real data.

### Fixed

**Service / data-access**
- Retry interval between live-page attempts shortened.
- pyIrena merging script path was hardcoded, forcing the wrong directory for all
  data; now resolved correctly at runtime.

**Science fixes (Phase 2, validated against Igor Pro reference data by JIL 2026-07-08)**
- Duplicate `"UPD_gains"` key in `calculatePD_Fly` silently dropped the range-index
  array; `smooth_r_data` smoothed every point with the 4 s time constant (range 4).
  Return dict now uses `UPD_gainsIndx` (index 0–4) and `UPD_gains` (values) as
  distinct keys; index mapping corrected 1-based → 0-based.
- Tail of `AmpGainReq_array` filled from `AmpGain` instead of `AmpReqGain`.
- Step-scan transmission: `trans_pin_*` ← diode and `trans_I0_*` ← I0 were
  swapped; `MeasuredTransmission` was the reciprocal of the intended value.
- Blank cache for step scans was never invalidated: `getBlankStepscan` now
  receives the caller's `recalculateAllData` flag (was hardcoded `False`).

**Crash fixes (Phase 1)**
- `desmearData` failure path returned 5 values; callers unpack 4 — fixed.
- `oneDesmearIteration` returned bare `1`; `ExtensionFailed` flag now truthful.
- `find_crossing_index` returned `float`, breaking downstream `range()` call.
- `reduceFlyscanToQR` raised `KeyError` on group deletion without existence check.
- `processUSAXSFolder` raised `IndexError` on missing technique folders; now skips
  with a warning and continues.
- `forceFirstBlank` with empty blank list silently failed; now logs one clear error
  and returns early.
- `save_dict_to_hdf5` crashed on `None` values and duplicate datasets.
- `subtract_data`: guard against `Y2_min == 0` (avoids `log(0)`).
- `calibrateAndSubtractFlyscan`: `nanmax`/`nanmin` replacing `max`/`min` (NaN-safe).
- `rebin_QRSdata`: added `< 2`-point guard; rebinning threshold changed to `>=`.
- `desmearData`: guard for `endme == 0` / NaN.
- `readGenericNXcanSAS`: defensive rewrite — returns `None` on malformed files.
- `saveNXcanSAS` / `readMyNXcanSAS`: `None`-attribute guards.
- `readMyNXcanSAS` trailing-comma `(None,)` tuple bug fixed.
- `extract_number_from_filename` regex now matches `.h5`, `.hdf`, `.hdf5`, `.nxs`.

**ASCII export (Phase 4)**
- `_read_group_data` used Python `or` between two `_read_arr()` calls to pick the
  dQ array; multi-element numpy arrays raised *"The truth value of an array is
  ambiguous"*. Fixed with explicit `None` check.
- `_hdf5_str()` helper safely converts HDF5 attributes (`np.ndarray`, `np.bytes_`,
  `bytes`, or `str`) to plain Python strings before comparison in
  `_find_nxcansas_groups`.
- Thickness attribute normalised to `float` or `""` at read time.
- Worker now logs full tracebacks via `logging.exception`.

**Robustness (Phase 3)**
- Duplicate Tiled filter keys removed; `title` and `hdf5_path` now matched
  client-side; module docstring documents the one-condition-per-filter-type limit.
- N+1 HTTP: `exit_status` added to `select_metadata`; per-uid fallback only when
  absent from metadata.
- `technique_detector` now checks `/entry/Metadata` first (Bluesky path fallback).
- FIFO eviction: `_remember_file()` helper + ordered dicts replace `set.pop()`
  (which removed an arbitrary element). Bound: `_MAX_TRACKED_FILES = 100`.
- `_find_matching_partner` disambiguates multiple matches by full stem and warns.
- `Bkg_map`/gain matching uses `np.isclose` (rtol 1e-3) with warnings on unknown
  gains in both `createUPDGainsAndBkgErrArrays` and `CorrectUPDGainsStep`.
- Path handling: `os.path.join` replaces manual string concatenation in four
  converters (`importFlyscan`, `importStepScan`, `ImportAndReduceAD`,
  `reduceFlyscanToQR`).
- Stray `print(Totaltime)` converted to `logging.debug`.

**GUI fixes (Phase 6)**
- `EXTRAP_METHODS` in `parameter_tabs.py` listed method names unknown to
  `desmearing.extendData`; the extension region silently recycled garbage.
  GUI list now matches real method names; unknown methods fall back to flat with
  a warning.
- STOP MOTORS was gated by the instrument-busy check; emergency stop now always
  sends `allstop`.
- Beamline-survey saved positions were overwritten by stale table contents on
  the next export/save; the table now reloads on dialog close.
- `main_window.closeEvent` cancels and joins running worker threads (destroying
  a running `QThread` was crashing the app on exit).
- `_set_slits` guards against `caget` timeout returning `None` before `caput`.

### Changed / Refactored (Phase 5, behavior-preserving)
- `read_group_to_dict` / `filter_nested_dict`: canonical copies in `hdf5code`,
  re-exported from `supportFunctions` for backwards compatibility.
- `clearAndCheckCachedReduction()` + `writeThicknessOverride()` in `hdf5code`
  replace ~50 duplicated lines across `processFlyscan`, `processStepscan`, and
  `process2Ddata`.
- `empty_calibrated_data()` factory replaces 4 pasted all-`None` dicts.
- `_build_search_uri()` in `readfromtiled` replaces 8 near-identical URI blocks;
  sort standardised to `-time`; 6 new unit tests cover URI construction.
- `convertSWAXS`: `_geometry_from_dicts()` + `_build_mask()` unify diverged
  geometry/mask copies; `ImportAndReduceAD` now gets detector-aware masks and
  `StartTime` metadata for the SAXS year branch.
- Large commented-out legacy blocks removed (`results_to_dataset`,
  `reduceStepScanToQR`, `find_NXcanSAS_entries`, tilt-test scaffolding).
- `tifffile` removed from runtime dependencies (only referenced in commented-out
  debug code).

### Documentation
- Clock-frequency invariants documented at all sites: flyscans use MCA 1 MHz
  clock; step scans use Joerger scaler 10 MHz clock. Both values confirmed
  correct against measurement data; must not be unified.
- Docstring drift corrected: 15 s → 5 s polling interval; `extrap_qstart` 0.15
  everywhere; desmear `MaxNumIter` 50 → 20.
- `thickness_override` file-mutation side-effect documented with NOTE comments
  in all three converters.

---

## [0.1.1] — 2026-05-11

Minor bug-fix release.

- Step-scan plotting bug fixes (two rounds).
- Blank plotting bug fix.
- Thickness override: honour user choice in HDF5 output.
- Flyscan: fix occasional wrong number of points.
- Old detector masks added.
- Thickness leakage in batch reduction fixed.
- Old step-scan structure compatibility fix.
- Blank search with name filter fixed.
- Input-variable stepping fixed.
- Graph-scaling behaviour fix.
- Installation instructions updated (clone-first flow for end users).

---

## [0.1.0] — 2026-04 (initial release)

First packaged release of Matilda. Core live-processing service for
USAXS/SAXS/WAXS at APS beamline 9-ID: flyscan and step-scan reduction,
Tiled/DataBroker integration, pyFAI-based 2-D integration, NXcanSAS output,
ASCII export, GUI tools (sample-plate setup, data-reduction monitor).
