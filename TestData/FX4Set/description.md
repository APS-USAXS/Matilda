# FX4 commissioning data — 2026-09-26

First data from the FX4 electrometer counting chain, taken during the
hardware commissioning session on 2026-09-26 (glassy carbon SRM 3600 and its
blank).  Source:
`/share1/USAXS_data/2026-09/09_26_Ilavsky/Commissioning`.

| folder | format | `counting_chain` marker |
|---|---|---|
| `usaxs/` | step scan (uascan), NXWriterUascan | `/entry/instrument/bluesky/metadata/counting_chain` |
| `saxs/`  | Pilatus frame | `/entry/counting_chain` |
| `waxs/`  | Eiger frame | `/entry/counting_chain` |

Detector values are gain-independent picoamps.  See
`bits_usaxs/docs/FX4_data_formats.md` and `matilda/fx4support.py`.

**There is no FX4 fly-scan file here** — stage 8 of the commissioning runbook
had not run when these were taken.  The fly-scan FX4 branch is covered by a
synthetic file built in `tests/test_fx4support.py`; replace that with a real
file when one exists.

## Known quirks of these particular files

* `I0` in the step scan is constant across all 300 points — `fx42` was not
  re-read per point in this scan.  The reduction is unaffected (it is a
  division by a constant) but the monitor normalisation is doing nothing.
* The transmission terms give T = 1.35 for the sample against this blank,
  which is unphysical.  TRD shares one `Range` with UPD on `usxFX4`, so the
  transmission measurements in this session are not trustworthy.
* SAXS/WAXS absolute intensity is ~10 orders of magnitude low because
  `I_scaling` was calibrated against the old `I0_cts/I0_gain` monitor.  See
  `convertSWAXS.FX4_I_SCALING`.
