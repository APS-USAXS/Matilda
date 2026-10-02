# FX4 commissioning data — 2026-09-26 / 27

First data from the FX4 electrometer counting chain, taken during the hardware
commissioning session (glassy carbon SRM 3600 and its blank).  All four
formats are represented.

| folder | format | source | `counting_chain` marker |
| --- | --- | --- | --- |
| `flyscan/` | fly scan, saveFlyData v2.0 | `09_26_Ilavsky/FlyScanning/FlyScanning_usaxs` (2026-09-27) | `/entry/program_name@counting_chain` |
| `usaxs/` | step scan (uascan), NXWriterUascan | `09_26_Ilavsky/Commissioning` | `/entry/instrument/bluesky/metadata/counting_chain` |
| `saxs/` | Pilatus frame | `09_26_Ilavsky/Commissioning` | `/entry/counting_chain` |
| `waxs/` | Eiger frame | `09_26_Ilavsky/Commissioning` | `/entry/counting_chain` |

Detector values are gain-independent picoamps.  See
`bits_usaxs/docs/FX4_data_formats.md` and `matilda/fx4support.py`.

The fly-scan files have had the groups a previous Matilda run wrote
(`QRS_data`, `Blank_data`, the NXcanSAS subentries) stripped, so they are the
raw writer output.

## The fly scan reduces correctly

`GC_SRM3600_0131` against `Blank_0130` gives transmission 0.936 and an
absolute plateau of 29.8 cm²/cm³ over 0.01 < q < 0.1, against a certified
SRM 3600 value near 30.  A repeat pair (`_R_0132` / `_R_0133`, not bundled)
agrees to ~1% over the same range.  USAXS absolute calibration therefore
carried across the conversion unchanged, as expected — the K-factor is
derived from the blank's own peak, not from a stored constant.

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
