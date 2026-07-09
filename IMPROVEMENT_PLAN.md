# Matilda — Code Review Findings & Improvement Plan

Date: 2026-07-08. Review scope: all core modules (`matilda/*.py`), GUI worker + technique detector, `pyproject.toml`, tests. Large GUI files (`sample_plate_setup.py`, `parameter_tabs.py`, `graph_panel.py`, `main_window.py`) received only a light pass and deserve a dedicated follow-up review.

Items marked ⚗️ affect scientific results and should be verified against Igor/known-good data before fixing.

---

## Progress log

**2026-07-08 — Phases 1 and 3 DONE** on branch `fix/phase1-phase3-robustness` (branched from `feature/tune-plots`).

Phase 1 (crash fixes) — all complete:
- ✅ 2.1 `desmearData` failure path now returns 4 values
- ✅ 2.2 `oneDesmearIteration` no longer returns bare `1`; `ExtensionFailed` flag made truthful (True only if even the flat fallback is impossible)
- ✅ 2.3 `find_crossing_index` returns `int`
- ✅ 2.4 `reduceFlyscanToQR` guards group deletion
- ✅ 2.5 `processUSAXSFolder` skips missing technique folders with a warning (helper `_list_sorted_files` added)
- ✅ 2.6 `forceFirstBlank` + empty blank list → single clear error, early return
- ✅ 2.8 `save_dict_to_hdf5` skips `None` values, uses `require_group`, overwrites existing datasets
- ✅ 2.9 `subtract_data` Y2_min==0 guard; `nanmax`/`nanmin` in `calibrateAndSubtractFlyscan`; `rebin_QRSdata` <2-point guard + `>=` threshold; `desmearData` endme==0/NaN guard; `readGenericNXcanSAS` defensive rewrite (returns None on malformed files); `saveNXcanSAS`/`readMyNXcanSAS` None-attribute guards
- ✅ 1.5 `(None,)` tuple bug fixed in `readMyNXcanSAS`
- ✅ 1.6 `extract_number_from_filename` regex now matches .h5/.hdf/.hdf5/.nxs

Phase 3 (robustness) — all complete:
- ✅ 2.7 Duplicate Tiled filter keys removed; title (FindScanDataByName) and hdf5_path (FindLastBlankScan) now matched client-side; module docstring documents the one-condition-per-filter-type limitation
- ✅ N+1 HTTP fixed: `exit_status:stop.exit_status` added to `select_metadata`; `convert_results` falls back to per-uid request only if absent
- ✅ 2.10 `technique_detector` checks `/entry/Metadata` (and Bluesky path as fallback)
- ✅ FIFO eviction: `_remember_file()` helper + ordered dicts replace `set.pop()`; single bound `_MAX_TRACKED_FILES = 100`
- ✅ Logging setup moved into `_setup_logging()` called from `main()` — importing `matilda.matilda` no longer creates directories or reconfigures the root logger
- ✅ `imagePath` overridable via `MATILDA_IMAGE_PATH` env var ("none" disables image saving)
- ✅ Docstring drift fixed: 15 s → 5 s polling (matilda.py, __init__.py, README); `extrap_qstart` 0.15 everywhere (incl. GUI worker defaults); desmear MaxNumIter doc 50 → 20
- ✅ `_find_matching_partner` disambiguates multiple matches by full stem and warns
- ✅ `Bkg_map`/gain matching tolerant (`np.isclose`, rtol 1e-3) with warnings on unknown gains (both `createUPDGainsAndBkgErrArrays` and `CorrectUPDGainsStep`)
- ✅ `os.path.join` replaces `path+"/"+filename` in importFlyscan, importStepScan, ImportAndReduceAD, reduceFlyscanToQR
- ✅ `thickness_override` file mutation documented with NOTE comments in all three converters; `filter_nested_dict` depth limitation documented (both copies)
- ✅ Stray `print(Totaltime)` → `logging.debug`; stale TODOs in readfromtiled removed

Verification: `py_compile` clean on all modules; functional sanity tests passed for `desmearData` (early return + full synthetic run), `find_crossing_index`, `subtract_data` (Y2_min==0), `rebin_QRSdata` (guard + normal path), `extract_number_from_filename`, `_remember_file` FIFO, `save_dict_to_hdf5` (None + overwrite), `readMyNXcanSAS` (SMR keys are None, not tuples), `saveNXcanSAS` (None attrs), `readGenericNXcanSAS` (malformed file); full package imports cleanly with no import-time side effects.

**2026-07-08 — Phase 2 (science fixes ⚗️) DONE**, same branch, except two open items below. ⚗️ **These fixes change numerical results and MUST be validated against Igor/known-good data before deploying:**

- ✅ 1.1 Duplicate `"UPD_gains"` key fixed: `calculatePD_Fly` now returns both `UPD_gainsIndx` (range index 0-4) and `UPD_gains` (gain values); `smooth_r_data` receives the index array (both call sites updated). Additionally the index mapping in `smooth_r_data` was corrected from 1-based to 0-based (ranges are 0-4 in this code base: `DDPCA300_gain0..4`). Behavioral impact: previously every point was smoothed with `times[4]=0.4 s` because gain values never equaled 1-4; now each range gets its intended smoothing time. **Validate flyscan smoothing vs Igor.**
- ✅ 1.2 `AmpGainReq_array` now filled from `AmpReqGain` (was `AmpGain`) — affects the gain-change/deadtime mask tail.
- ✅ 1.3 Step-scan transmission pin↔I0 swap fixed in `importStepScan`: diode stream → `trans_pin_*`, I0 stream → `trans_I0_*`. Previously `MeasuredTransmission` was the reciprocal of the intended value. **Validate step-scan transmission & absolute calibration vs Igor.**
- ✅ 1.4 `getBlankStepscan` now receives the caller's `recalculateAllData` (was hardcoded False) — forced reprocessing invalidates cached blank data.

Open Phase 2 items (deliberately NOT changed):
- ⏳ **Frequency 1e6 vs 1e7**: scaler is 1e7, but the clock signal source must be physically traced at the instrument before unifying (flyscan code uses 1e6). TODO comments added at all three 1e7 sites in convertUSAXS.py. Do not change blindly.
- ⏳ **`Error = SigmaRwave / 5`** factor in both error calculators — author-intentional approximation, left as is pending review.

Verification: compile clean; `smooth_r_data` exercised with a proper 0-4 index array (incl. NaN masked points) through both the no-smoothing and averaging/fit branches; source-level assertions confirm all four fixes are in place. Full pipeline validation against real HDF5 data still required (⚗️).

**NEXT: validate Phase 2 against beamline data, resolve frequency question on-site.** Then Phase 4 (tests + lint), Phase 5 (refactor), Phase 6 (GUI deep review).

---

## Priority 1 — Bugs that produce wrong results ⚗️

### 1.1 Duplicate `"UPD_gains"` key silently drops gain indices (Flyscan smoothing broken)
`supportFunctions.py` lines 560–564, `calculatePD_Fly()` return dict:
```python
result = {"Intensity":PD_Intensity,
          "Error":PD_error,
          "UPD_gains":GainsIndx,   # <- silently overwritten by next line
          "UPD_gains":Gains,
          "UPD_bkgErr":updBkgErr}
```
`smooth_r_data()` compares `UPD_gains[i] == 1 … 4` — it expects the gain **index** (`GainsIndx`), but receives the actual gain **values** (1e4–1e12). Every point therefore falls into the `else` branch and gets smoothed with the longest time constant (0.4 s). Fix: return both under distinct keys (`UPD_gainsIndx`, `UPD_gains`) and pass the index array to `smooth_r_data`. Verify smoothed output against Igor.

### 1.2 Wrong array used to build requested-gain array
`supportFunctions.py` line 456:
```python
AmpGainReq_array = np.full(num_elements, AmpGain[len(AmpReqGain)-1])
```
Should be `AmpReqGain[...]`. As written, the tail of the requested-gain array is filled from `AmpGain`, which can corrupt the `AmpGain == AmpReqGain` mask used for deadtime/gain masking.

### 1.3 Step-scan transmission metadata appears swapped (pin ↔ I0) ⚗️
`convertUSAXS.py` lines 406–420, `importStepScan()`: values read from `terms_USAXS_transmission_I0_*` are stored as `trans_pin_*`, and values from `terms_USAXS_transmission_diode_*` are stored as `trans_I0_*`. Because both counts and gains are swapped consistently, `MeasuredTransmission` in `calibrateAndSubtractFlyscan()` computes the **inverse** of the intended value, which propagates into `MSAXSCorrection` and absolute calibration of step scans. Verify against a known sample; if the swap is intentional (naming quirk in the Bluesky stream), rename the local variables and add a comment.

### 1.4 Blank cache never invalidated for step scans
`convertUSAXS.py` line 172: `getBlankStepscan(blankPath, blankFilename, recalculateAllData=False)` — hardcoded `False` instead of passing the caller's `recalculateAllData`. Forced reprocessing of a step scan still reuses stale cached blank data. (Flyscan path passes it correctly.)

### 1.5 `readMyNXcanSAS` stores `(None,)` tuples instead of `None`
`hdf5code.py` lines 654–658 — trailing commas:
```python
Sample["CalibratedData"]["SMR_Qvec"] = None,   # tuple (None,) !
```
Downstream `if x is not None:` checks then pass and code attempts to plot/save a `(None,)` tuple. Remove the trailing commas.

### 1.6 USAXS files not sorted by scan number in `processUSAXSFolder`
`matilda.py` line 698: `extract_number_from_filename` regex is `_(\d+)\.hdf` — never matches `.h5` flyscan files, so every key is 0 and the USAXS list keeps `os.listdir` order. Blank matching by "nearest preceding number" still works (it parses numbers separately) but processing order and any order-dependent logic is wrong. Fix regex to `_(\d+)\.(h5|hdf5?)$` or reuse `_parse_filename_info` from supportFunctions.

---

## Priority 2 — Bugs that crash or hang under real conditions

### 2.1 `desmearData` failure path returns 5 values, callers unpack 4
`desmearing.py` line 477: `return None, None, None, None, None`. Both `convertFlyscan.py` (line 208) and `convertUSAXS.py` (line 184) unpack 4 → `ValueError` whenever desmearing is skipped (empty input or missing slit length).

### 2.2 `oneDesmearIteration` returns bare `1` on extension failure
`desmearing.py` line 396 — caller unpacks 4 values → `TypeError`. Currently unreachable because `extendData` always resets `ExtensionFailed = False` after applying the flat fallback (line 171), which itself makes the flag meaningless. Clean up both: remove the flag or make it truthful, and never return a scalar.

### 2.3 `find_crossing_index` returns a float
`supportFunctions.py` line 829: returns `0.1*len(array)` when the target is never crossed; used as `range()` start in `smooth_r_data` → `TypeError`. Return `int(0.1*len(array))` (and fix the docstring that claims it returns `None`).

### 2.4 `reduceFlyscanToQR` deletes group without existence check
`convertFlyscan.py` lines 266–269: with `recalculateAllData=True`, `del hdf_file[location]` raises `KeyError` if the group doesn't exist yet. Guard with `if location in hdf_file`.

### 2.5 `processUSAXSFolder` crashes if a technique folder is missing
`matilda.py` lines 524–528: `[f for f in folders if f.endswith('_usaxs')][0]` → `IndexError` when any of the three folders is absent (docstring says "silently returns" — it doesn't). Handle each technique independently and skip missing ones.

### 2.6 `forceFirstBlank` with an empty blank list
`matilda.py` (`processFlyscans` / `processStepscans` / `processADscans`): `ListOfBlanks[0]` raises `IndexError` per scan; caught-and-logged, so *every* scan silently fails. Validate up front and log one clear message.

### 2.7 Duplicate Tiled filter keys — filters silently overwrite each other
`readfromtiled.py`: `FindScanDataByName` uses `filter[eq]` twice (plan_name and title); `FindLastBlankScan` with `path` uses `filter[regex]` twice (title and hdf5_path). Tiled keeps only one condition per key, so one filter is ignored (already flagged in the module TODO). Use distinct keys (`filter[eq][condition]…&filter[eq2]…` per Tiled syntax) or the documented multi-condition form, and add a test against the live server.

### 2.8 `save_dict_to_hdf5` cannot handle `None` values or pre-existing groups
`hdf5code.py` line 822+: `h5file[path + key] = item` raises on `None` (present in many result dicts, e.g. the all-None `CalibratedData`) and `create_group` raises if the group exists. Skip `None`s (or write empty datasets with an attribute) and use `require_group`.

### 2.9 Edge cases in numeric helpers
- `subtract_data` (`supportFunctions.py` ~856): offset applied only when `Y2_min < 1e-30`; `Y2_min == 0` exactly gives `log(0) = -inf`. Also `np.max(SMR_Int)` / `np.min` in `calibrateAndSubtractFlyscan` are not NaN-safe → use `nanmax/nanmin`.
- `rebin_QRSdata` (~1027): `Wx_greater[1]` → `IndexError` with <2 points above Q=0.0002; points exactly at 0.0002 are dropped (`<` and `>` both exclusive).
- `desmearData` (~499): `difff = 1 - (oldendme / endme)` — degenerate if `endme` is 0/NaN; guard.
- `readGenericNXcanSAS` (`hdf5code.py`): `Int_attributes`, `Q`, `Error`, `dQ` may be referenced undefined if any expected dataset is missing → `NameError` instead of a clear error message.
- `saveNXcanSAS`: `ds.attrs['blankname'] = blankname` fails with `TypeError` if `blankname`/`thickness` is `None`; guard like `Kfactor`.

### 2.10 GUI technique detection reads the wrong metadata path
`gui/data_reduction/technique_detector.py` `_detect_area_detector()` looks for `pin_ccd_*`/`waxs_ccd_*` under `/entry/instrument/bluesky/metadata`, but SAXS/WAXS files store these under `/entry/Metadata` (see `convertSWAXS.py`). Content-based detection never fires; only the folder-name fallback works. Fix the path (check both).

---

## Priority 3 — Robustness & correctness hygiene

- **Bounded-set eviction is random**: `set.pop()` in `matilda.py` (`_runPynikaCalibration`, `_runPyirenaAnalysis`, `_runMergeData`) removes an *arbitrary* element — possibly the one just added — so files can be re-calibrated/re-analyzed. Use an ordered structure (`collections.OrderedDict` / `deque`) for FIFO eviction. Also `calibrated_set` is bounded by `NumberOfImagesInGraphs` (10) while the others use `_MAX_PYIRENA_ANALYZED` (100) — unify.
- **Import-time side effects in `matilda/matilda.py`**: `logging.basicConfig(...)` + `os.makedirs(log_dir)` run on import, and `imagePath` is hardcoded to a beamline path at module level. Move logging setup and `imagePath` resolution into `main()` (env-var override like `MATILDA_LOG_DIR` already exists — add `MATILDA_IMAGE_PATH`).
- **Docstring vs. behavior drift**: main loop sleeps 5 s but docstrings/README say 15 s; `extrap_qstart` documented as 0.1 but defaults to 0.15 (and the GUI worker defaults to 0.1 — three different values for the same knob); `desmearData` docstring says max 50 iterations, code uses 20. Pick one truth per knob.
- **`_extract_sample_key` / `_find_matching_partner`** (`matilda.py`): key is `(first_token, last_token)` so `SampleA_10min_0044` and `SampleA_20min_0044` collide; `matches[0]` of the glob is arbitrary. Consider matching on the full stem minus scan number, and warn on multiple matches.
- **N+1 HTTP requests**: `convert_results()` calls `successful_run(uid)` (one extra request) per run, every 5 s cycle. Include `exit_status:stop.exit_status` in `select_metadata` instead and drop the per-uid round trips.
- **`Frequency` inconsistency (1e6 vs 1e7)** ⚗️: `calculatePDErrorFly` uses 1e6, `calculatePDErrorStep` and `createUPDGainsAndBkgErrArrays` use 1e7, `CorrectUPDGainsStep` has its own TODO about it. Confirm the actual MCA clock per scan type and centralize as named constants.
- **`Error = SigmaRwave / 5`** ⚗️: magic factor in both error calculators, flagged "close enough for now" — document or fix.
- **Stray debug output**: `print(f"{Totaltime}")` in `calculatePD_Fly` (supportFunctions.py line 512) runs on every scan — remove or `logging.debug`.
- **Bkg_map float-equality matching** (`convertUSAXS.py`): `Bkg_map_float_keys.get(gain, 0)` and the `gain == 1e4 …` ladder rely on exact float equality with EPICS-sourced values; unknown gains silently get background 0 and `UPD_gains 0` (→ divide-by-zero later). Match with tolerance and log unknown gains.
- **`thickness_override` permanently rewrites raw data files** (all three converters): original is kept as `thickness_original`, but the raw file is still mutated even when the user only experiments with values. Consider keeping overrides in the NXcanSAS output only, or documenting this prominently.
- **Path handling**: `path+"/"+filename` in `importFlyscan`, `importStepScan`, `ImportAndReduceAD`, `reduceFlyscanToQR` → `os.path.join` everywhere (matters on Windows, which the test scripts clearly target).
- **`filter_nested_dict` drops nested containers** whose keys aren't in `keys_to_keep` even when they contain wanted keys deeper (`monochromator` works only because it's explicitly listed) — document or make it depth-aware.

---

## Priority 4 — Refactoring / de-duplication (no behavior change)

1. **Deduplicate helpers**: `read_group_to_dict` and `filter_nested_dict` exist identically in both `hdf5code.py` and `supportFunctions.py`. Keep one (hdf5code) and import it.
2. **Merge the twin pipelines**: `processFlyscan` vs `processStepscan` and `getBlankFlyscan` vs `getBlankStepscan` are ~90 % identical, including a 15-line "all-None CalibratedData" dict pasted four times. Extract shared cache-check/delete logic and an `EMPTY_CALIBRATED_DATA` factory.
3. **Factor Tiled URI building** (`readfromtiled.py`): six near-identical 20-line URI blocks → one builder function taking `plan_name, title_regex, path_regex, NumScans, lastNdays`. This also makes fixing 2.7 a one-place change.
4. **Unify mask/geometry code in `convertSWAXS.py`**: `ImportAndReduceAD` and `reduceADData` duplicate geometry setup and masks — and they have *diverged* (`ImportAndReduceAD` lacks the old-detector/year branching). Extract `_build_mask(my2DData, usingWAXS, metadata)` and `_geometry_from_dicts(...)`.
5. **Remove dead code / unused imports**: unused `plan_name` variable in `process2Ddata`; unused imports (`socket`, `tifffile`, `FindLastBlankScan` in convertSWAXS; `curve_fit`, `interp1d`, `subtract_data` in convertUSAXS; `pprint`, `copy` in several files); large commented-out blocks (`results_to_dataset`, `reduceStepScanToQR`, old test scaffolding) — either delete or move to CodeFragments/.
6. **Module layout**: follow through on the planned `matilda/reduction/`, `matilda/io/` split described in `matilda/__init__.py` and `reduction/__init__.py` once the above is stable.

---

## Priority 5 — Testing & tooling

- **No automated tests exist** (`tests/manual_test.py` has hardcoded `C:/Users/ilavsky/...` paths). The repo ships `TestData/` — add pytest smoke tests that: reduce one flyscan + blank, one step scan, one SAXS and one WAXS file end-to-end on copies of the test files (converters mutate their inputs!), and assert array shapes, positive Q, finite intensities, and round-trip via `readMyNXcanSAS`.
- Unit tests for the pure functions where the P1/P2 bugs live: `_parse_filename_info`, `extract_number_from_filename`, `find_crossing_index`, `subtract_data`, `rebin_QRSdata`, `desmearData` early-return, `_extract_sample_key`.
- Add `ruff` (or flake8) — it would have caught the duplicate dict key (1.1), unused imports, and f-string issues mechanically. Small `pyproject.toml` addition + CI via GitHub Actions.
- Consider a mocked-`requests` test for `readfromtiled` filter construction (2.7).

---

## Suggested execution order

| Phase | Items | Risk |
|---|---|---|
| 1. Crash fixes | 2.1–2.6, 2.8, 2.9, 1.5, 1.6 | Low — mechanical, testable |
| 2. Science fixes (verify vs Igor) ⚗️ | 1.1, 1.2, 1.3, 1.4, Frequency + `/5` questions | Needs beamline validation data |
| 3. Robustness | P3 items, 2.7, 2.10 | Low–medium |
| 4. Tests + lint | P5 | None — do alongside phases 1–2 |
| 5. Refactor | P4 | Medium — after tests exist |
| 6. GUI deep review | `sample_plate_setup.py`, `parameter_tabs.py`, `graph_panel.py`, `main_window.py`, `file_tree.py`, `sas_plot.py`, `ascii_exporter.py`, `viewer_2d.py` | Separate pass |
