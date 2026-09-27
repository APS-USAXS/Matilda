# Matilda Changelog

All notable changes to this project are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

---

## [Unreleased]

### Added — FX4 counting chain

The 12-ID-E counting chain changed on 2026-09-26 from Femto amplifier + V/F
converter + Struck scaler to **FX4 electrometers**. All four data formats
changed with it. Matilda now reads both chains; which one produced a file is
taken from the `counting_chain` marker each format carries, never guessed
from field names (`UPD` and `I0` in the step-scan file kept their names and
changed their units). Reference: `bits_usaxs/docs/FX4_data_formats.md`.

- **`matilda/fx4support.py`** (new): `detect_counting_chain()` / `is_fx4()`
  and the FX4 arithmetic helpers (`range_indexed_array`, `mean_of_samples`,
  `ratio_error`, `warn_ring_overflows`). Chain detection covers all four
  marker locations, including the fly scan's — which `saveFlyData.xml` writes
  as an *attribute* of `/entry/program_name`, not as the dataset the format
  document describes. Both spellings are accepted.
- **Fly scan** (`supportFunctions.importFlyscan`, `calculatePD_Fly`,
  `calculatePDErrorFly`): FX4 files are read from `upd_current` / `I0_current`
  (mean pA over each PSO interval). `I(q) = (upd_current − dark) / I0_current`
  — no gain term and no dwell divisor. Uncertainties come from the recorded
  `upd_sigma` / `I0_sigma` divided by √N (N = total/mean), combined with the
  dark-current error and propagated through the ratio. `ring_overflows` and
  `ring_overflows_I0` are checked and a non-zero value is reported.
- **Step scan** (`convertUSAXS.importStepScan`, `CorrectUPDGainsStep`,
  `createUPDGainsAndBkgErrArrays`, `calculatePDErrorStep`): FX4 `UPD` / `I0`
  are picoamps; gains are 1, dark current is looked up per point from
  `/entry/data/fx4_autorange_lurange` against the FX4 sequence program's
  `fx4_autorange_ranges_range*_background`. There is no `seconds` column on
  this chain, so the per-point count time is reconstructed from `plan_args`
  and `useDynamicTime` (reporting only — nothing divides by it).
- **SAXS / WAXS** (`convertSWAXS.calibrateAD2DData`): the monitor is
  `I0_cts_gated` (equivalently `/entry/control/integral`) with gain 1, not the
  stale `I0_cts` / `I0_gain` that FX4 frames still carry. WAXS transmission
  uses `TR_current` against `I0_current`; `TR_cts_gated` is deliberately
  ignored — only I0 is gated, so that field holds the leftover total of the
  0.05 s autoscale read rather than the exposure.
- **Transmission** (`supportFunctions.transmissionTerms`): `saveFlyData.xml`
  v2.0 stops recording `trans_pin_gain` / `trans_I0_gain` because the diode
  and I0 readings are gain-independent. Their absence is now the expected case
  on the FX4 chain and both are taken as 1; a recorded value that is not 1 is
  reported and ignored. Only a missing *count* falls back to T = 1.0 with a
  warning, and the scaler chain still divides by the recorded gains.
- Reducing an FX4 frame against a scaler-chain blank (or vice versa) now
  raises instead of silently mixing counts with picoamps.
- **`TestData/FX4Set/`**: real commissioning data (glassy carbon SRM 3600 and
  blank) for all four formats, with the known quirks of those particular files
  documented. The synthetic fly scan in `tests/test_fx4support.py` is kept
  alongside it for the cases the real files cannot exercise — a known
  transmission to recover, and a file with no `/entry/metadata` group.
- 34 new tests. The legacy scaler-chain reduction is bit-for-bit unchanged.

### Validation

`GC_SRM3600_0131` / `Blank_0130` (fly scan) reduces to transmission 0.936 and
an absolute plateau of 29.8 cm²/cm³ over 0.01 < q < 0.1, against a certified
SRM 3600 value near 30; a repeat pair agrees to ~1% there. USAXS absolute
calibration therefore carried across the conversion unchanged, which is what
the self-referencing K-factor predicts.

### Known limitations of the FX4 path

- **SAXS/WAXS absolute intensity is not calibrated.** `I_scaling` was
  determined against the old `I0_cts/I0_gain` monitor; the FX4 monitor is a
  different quantity, so FX4 frames come out roughly ten orders of magnitude
  low. The curve *shape* is correct (diode/I0 is a ratio). Set
  `convertSWAXS.FX4_I_SCALING` from a standard reference material. A warning
  is logged on every FX4 frame until that is done. USAXS is unaffected — its
  K-factor is derived from the blank's own peak.
- **FX4 step-scan uncertainties are modelled, not measured.** The uascan file
  records no per-point sigma for the currents, so
  `fx4support.FX4_RELATIVE_CURRENT_ERROR` (1%) stands in. Error bars only.
- **FX4 fly-scan uncertainties are optimistic by 3-4x.** `upd_sigma`/√N came
  out consistently below the observed point-to-point scatter of I(q) on the
  first real fly scans. √N assumes independent samples; the FX4's are filtered
  averages of `ValuesPerRead` raw readings, so the effective N is smaller.
  Calibrate against repeat scans before relying on the bars.
- **The FX4 fly scan records the amplifier range once, after the scan.** There
  is no per-point range, so one dark current is subtracted throughout and the
  amplifier dead-time masking the scaler chain did from
  `changes_DDPCA300_mcsChan` has no equivalent. In practice the range-4
  smoothing window this implies only applies above q = 3e-4, so the rocking
  curve itself is untouched.
- The first deployed `saveFlyData.xml` v2.0 wrote no `/entry/metadata` group
  (an over-long XML comment swallowed the tag); fixed at the beamline on
  2026-09-27. Matilda still falls back to reading those fields from
  `/entry/flyScan`, prefers the FX4 `upd_bkg*` over any Femto-era copies in
  `/entry/metadata`, and raises a specific error if the geometry fields are
  missing entirely.

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
