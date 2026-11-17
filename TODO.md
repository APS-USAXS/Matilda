# Matilda TODO

Open work only. Finished items go to `CHANGELOG.md`; design notes for larger
pieces live in `PLAN.md` / `plans/`. Igor-side FX4 work is tracked separately.

## FX4

- **Step-scan error bars.** FX4 step scans still use the flat 1 %
  `FX4_RELATIVE_CURRENT_ERROR` (fx4support). It now matters, because the Abel
  auto-smoothing picks its width from χ² against these errors. Revisit after
  the upcoming FX4 support changes: apply the fly-scan shot term
  (`FX4_SHOT_CHARGE_PC`) with the count time from `_fx4StepCountTime`, and get
  a per-point sigma recorded by the Bluesky step-scan plan for the
  electronic floor.
- **Smoothing-range flicker (fix belongs in Bluesky/FX4, not here).** On
  `GC_SRM3600_0131` the signal-derived pseudo-range changes 49 times along the
  scan (5 with an 11-point median filter), because noisy points cross the
  30 pA / 800 pA thresholds back and forth; neighbouring points then get
  0.1 s vs 0.4 s smoothing and alternate between the trapezoid and the
  linear-fit branch. Intended fix: change the FX4 range overlap so ranges
  do not switch at those levels (tracked in the Bluesky code). Revisit
  `FX4_SMOOTHING_RANGE_THRESHOLDS` in Matilda afterwards. Side cleanup:
  `np.polyfit` can replace `curve_fit` for the straight-line fit.

## Backend

1. **Daemon subprocesses off the main loop.** pynika / pyirena / merge run
   synchronously with a 300 s timeout each and can stall live reduction.
   (Per-block exception isolation and bounded retries are done.)
   Develop on a branch and test on the beamline machine before merging.
2. **One pipeline for fly and step.** `processFlyscan` / `processStepscan`,
   the two `getBlank*` functions and the two scaler error functions are near
   copies; the counting chain adds a third axis.
3. **Package split** into `reduction/` + `io/` (placeholder
   `matilda/reduction/__init__.py`).
4. **Live-server check** of the July Tiled query changes (`readfromtiled.py`).
5. **Transmission uncertainty.** `sastransmission_spectrum/Tdev` is written as
   NaN ("not estimated"). Estimate it per chain: Poisson on the scaler
   counts, and the FX4 current sigma for picoamps. Then store it in `Tdev`
   and in `sasprocess`.

## GUI

1. **SAXS azimuth controls are dead** (`parameter_tabs.py` returns
   `az_min`/`az_max`, nothing consumes them): wire through `process2Ddata` and
   re-validate against Nika, or remove.
2. **Persist reduction parameters** per folder (session file stores only the
   folder and splitter sizes).
3. **Warnings vs errors** in the log panel and final status; list failed files.
4. **Confirm batch overrides**: "Process All" with Override ticked reduces
   every file with the same thickness/transmission (recorded in each reduced
   entry; raw files untouched) — ask first, naming the count.
5. **Blank overlay and chain**: re-enable the blank curve overlay by reading
   `entry/Blank_data` directly; show sample and blank counting chain.
6. **GUI smoke test in CI**: one offscreen test that builds the main window
   and runs the worker on a bundled file.
7. **`.mac` re-import** in the sample-plate editor.

## Features to consider

1. Per-folder reduction recipe honoured by the daemon (`matilda_config.json`
   next to `pyirena_config.json` / `merge_config.json`): thickness, blank,
   desmearing method, Qmin.
2. SAXS/WAXS absolute-calibration tool: reduce an SRM 3600 frame, compare with
   the certified curve, report and store the scaling.
3. Data-quality flags in the output and web plots: ring-buffer overflows,
   transmission fallback to 1.0, strong MSAXS correction, auto Qmin not found.
4. GSAS-II `.instprm` generator (designed in `plans/waxs_instprm_generator.md`).
5. Command-line batch reduction, `matilda-reduce <folder> --recalculate`, on
   top of `processUSAXSFolder` (which now detects fly vs step scans).
6. Per-experiment summary report: curves, transmissions, thicknesses, blanks
   and flags per sample folder.

## Conventions to remember (not tasks)

- A file is a blank if its name contains "blank" (case-insensitive); that
  string is not allowed in sample names.
- HDF5 files grow on each re-reduction (deleted groups are not reclaimed);
  accepted, re-reduction is rare.
- FX4 fly scans record one dark current (`upd_lurange`, the highest-gain
  range); the other ranges' dark currents are 0. Tested and sufficient — a
  per-point dark current is not needed.
