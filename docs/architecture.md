# Matilda — Architecture Overview

## Big picture

```
Tiled server (HTTP/REST)
        │
        │  FindLastScanData()
        │  FindLastBlankScan()
        ▼
matilda.py  ─── polling loop (5 s sleep between cycles)
        │
        ├── convertFlyscan.py   (USAXS flyscan → I(Q))
        ├── convertUSAXS.py     (USAXS step-scan → I(Q))
        ├── convertSWAXS.py     (SAXS/WAXS 2-D → I(Q) via pyFAI)
        │
        ├── pynika (subprocess, conda run)
        │       triggered by AgBehenateLaB6 filenames
        │
        ├── pyirena (subprocess, conda run)
        │       fit_pyirena  — auto-analysis per technique folder
        │       merge_data   — USAXS+SAXS auto-merge
        │
        ├── plotData.py         (JPEG summary plots → web directory)
        └── hdf5code.py         (NXcanSAS read/write)
```

---

## Module reference

### Orchestrator

**`matilda/matilda.py`**
Main polling loop (5 s sleep between cycles).  Queries Tiled for new scans, compares to the
previous list, dispatches reduction and plotting.

Entry point: `main()` — called by the `matilda` console script and by
`python -m matilda.matilda`.

Key module-level constants.  The first seven are read from the site config
(`matilda/config.py`, `~/.matilda/config.json`, see
[configuration.md](configuration.md)); the built-in defaults are the
12-ID-E values:

| Constant | Config key | Purpose |
|---|---|---|
| `imagePath` | `image_path` (env `MATILDA_IMAGE_PATH` wins) | Directory for JPEG web plots |
| `NumberOfDaysToLookBack` | `days_to_look_back` | Tiled query window for scans |
| `NumberOfDaysToLookBackBlanks` | `days_to_look_back_blanks` | Tiled query window for blanks |
| `NumberOfImagesInGraphs` | `images_in_graphs` | Max datasets per plot |
| `CONDA_EXECUTABLE` | `conda_executable` | Path to conda for pynika/pyirena subprocesses |
| `PYNIKA_CONDA_ENV_PATH` | `pynika_conda_env_path` | pynika environment path |
| `PYIRENA_CONDA_ENV_PATH` | `pyirena_conda_env_path` | pyirena environment path |

| Constant | Purpose |
|---|---|
| `_PYIRENA_CONFIG_FILENAME` | Per-folder config that triggers pyirena analysis (`pyirena_config.json`) |
| `_MERGE_CONFIG_FILENAME` | Per-experiment config that triggers USAXS+SAXS merging (`merge_config.json`) |
| `_MERGED_FOLDER_NAME` | Output subfolder for merged data (`data_usaxs_merged`) |

---

### Counting chains

12-ID-E has had two counting chains. Every reduction path branches on which
one wrote the file; nothing else about the pipeline differs.

| aspect | **scaler** (until 2026-09-26) | **FX4** (from 2026-09-26) |
| --- | --- | --- |
| hardware | Femto amplifier + V/F converter + Struck 3820 | FX4 electrometers |
| detector value | counts accumulated over a dwell | mean current, picoamps |
| gain | per-range table, divided out | none — reading is gain-independent |
| dwell | divided out (1e6 Hz fly MCA / 1e7 Hz Joerger) | none — the value is already a mean |
| `I(q)` | `(counts − dark·t)/(f·gain) / (I0/I0gain)` | `(upd_current − dark) / I0_current` |

**`matilda/fx4support.py`**
Chain detection and the FX4-specific arithmetic.

Key functions: `detect_counting_chain`, `is_fx4`, `range_indexed_array`,
`mean_of_samples`, `ratio_error`, `warn_ring_overflows`

Every format declares `counting_chain`; absent means the old scaler chain.
That marker is the only correct test — `UPD` and `I0` in the step-scan file
kept their names across the conversion and changed only their units, and the
step-scan `program_name@config_version` is the NeXus writer's schema version,
unrelated to the chain. The marker lives in a different place in each format,
and `detect_counting_chain` checks all of them:

| format | location |
|---|---|
| fly scan | `/entry/program_name@counting_chain` (also accepts `/entry/flyScan/counting_chain`) |
| SAXS / WAXS frame | `/entry/counting_chain` |
| step scan | `/entry/instrument/bluesky/metadata/counting_chain` |

Reference: `bits_usaxs/docs/FX4_data_formats.md`. One knob is still a
placeholder: `fx4support.FX4_RELATIVE_CURRENT_ERROR` (step-scan error bars;
see TODO.md). SAXS/WAXS absolute intensity uses the recorded `I_scaling`,
recalibrated at the beamline for FX4; `convertSWAXS.FX4_I_SCALING` is only an
override hook.

For the full per-technique description of the FX4 reduction — dataset paths,
arrays, formulas, transmission terms and validation numbers, written so the
Igor Pro reduction can follow the same path — see
[`fx4-reduction-guide.md`](fx4-reduction-guide.md).

---

### Data reduction

**`matilda/convertFlyscan.py`**
USAXS flyscan: HDF5 NXsas → normalised I(Q) → desmeared calibrated I(Q).

Pipeline:
`importFlyscan` → `calculatePD_Fly` → `beamCenterCorrection` →
`getBlankFlyscan` → `normalizeByTransmission` →
`calibrateAndSubtractFlyscan` → `desmear_dispatch` → `saveNXcanSAS`

**`matilda/convertUSAXS.py`**
USAXS step-scan: identical pipeline to flyscan.
Also exports `rebinData()` used by convertFlyscan.

**`matilda/convertSWAXS.py`**
SAXS/WAXS area-detector: HDF5 → pyFAI azimuthal integration → I(Q).
Converts Nika geometry (beam-centre, tilts) to pyFAI PONI via Fit2D helper.

Pipeline:
`importADData` → `calibrateAD2DData` → `reduceADData` → `saveNXcanSAS`

**`matilda/desmearing.py`**
Lake/Strobl iterative slit-smearing correction, ported from Igor Pro.
Entry point: `desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength, …)`

**`matilda/desmearing_methods.py`**
Selectable desmearing methods behind one dispatcher with the same signature and
4-tuple return as `desmearData`. Default is the truncated-Abel inversion;
Lake and the Gaussian-process variants remain selectable.
Entry point: `desmear_dispatch(…, method='abel'|'lake'|'gp'|…)`.
See `docs/desmearing-methods.md`.

---

### Support

**`matilda/supportFunctions.py`**
Shared numerical helpers used by all three reduction modules.

Key functions: `importFlyscan`, `calculatePD_Fly`, `beamCenterCorrection`,
`smooth_r_data`, `getBlankFlyscan`, `normalizeByTransmission`,
`calibrateAndSubtractFlyscan`, `rebinData`, `findProperBlankScan`,
`read_group_to_dict`, `filter_nested_dict`, `subtract_data`

**`matilda/supportNikaFunctions.py`**
Nika → Fit2D → pyFAI PONI geometry conversion.

Key function: `convert_Nika_to_Fit2D(*, SSD, pix_size, BCX, BCY, HorTilt, VertTilt, wavelength)`

**`matilda/technique_detector.py`**
Flyscan / StepScan / SAXS / WAXS detection from file content (h5py only, no
Qt).  Used by the GUI and by `processUSAXSFolder`, which sends each scan in a
`_usaxs` folder to the fly-scan or step-scan reducer accordingly.

---

### I/O

**`matilda/hdf5code.py`**
HDF5/NXcanSAS read-write helpers.

Key functions: `saveNXcanSAS`, `readMyNXcanSAS`, `readGenericNXcanSAS`,
`save_dict_to_hdf5`, `load_dict_from_hdf5`, `find_matching_groups`

The `Matilda_version` attribute in saved files is read from
`importlib.metadata.version('Matilda')` at write time.

**Cached results.** Each group `saveNXcanSAS` writes carries a
`matilda_reduction_version` attribute (`hdf5code.REDUCTION_VERSION`).
`clearAndCheckCachedReduction` (shared by all three converters) treats a
missing or different value as a cache miss, so results written by an older
reduction are recomputed rather than served.  Bump `REDUCTION_VERSION`
whenever a change alters reduced numbers.  Cached USAXS blanks
(`entry/blankData`, fly and step) are versioned the same way by
`supportFunctions.BLANK_CACHE_VERSION` and record their counting chain; a
sample is never calibrated against a blank from the other chain.

**Entry name.** The NXcanSAS subentry is `/entry/<sample name>` (and
`<sample name>_SMR`).  If that path is an existing raw-data group (a sample
named `sample`, `metadata`, `data`…), a reserved name, or the sample name is
empty, the entry is written under a safe name instead (`<name>_NXcanSAS`, or
the file stem) and a warning is logged.  Raw data are never overwritten.

**Reduction provenance.** Both subentries (`<name>` and `<name>_SMR`) carry
an NXcanSAS `sasprocess` group (`NX_class` NXprocess, `canSAS_class`
SASprocess) holding `name` (`Matilda`), `date`, `description`,
`Matilda_version` and `reduction_version`.  Its `sasnote` group (NXnote)
holds one dataset per reduction parameter.  The values come from
`supportFunctions.reduction_provenance`, called by each converter after
calibration.  None values are skipped, flags are stored as 0/1, and the
`units` attributes come from `hdf5code.PROVENANCE_UNITS`.

| Group of keys | Datasets |
|---|---|
| What was reduced | `technique` (Flyscan/StepScan/SAXS/WAXS), `counting_chain` (FX4/scaler), `blank_file`, `units` |
| Transmission | `transmission_used`, `transmission_measured` (before any override), `transmission_override`, `transmission_overridden`, `PeakToPeakTransmission`, `MSAXSCorrection` (USAXS) |
| Thickness | `thickness_used`, `thickness_file`, `thickness_override`, `thickness_overridden` (true only when the override is what was used) |
| μ mode | `mu_mode`, `mu`, `calculated_mu`, `per_gram`, `density` |
| Qmin (USAXS) | `qmin_used`, `calculated_qmin` (auto), `qmin_override`, `qmin_overridden`, `minQMinFindRatio` |
| Desmearing (USAXS) | `desmeared`, `desmear_method`, `desmear_extrap_method`, `desmear_extrap_qstart`; Abel: `abel_auto_smooth`, `abel_smooth_w` (the width actually used, ln q), `abel_chi2`, `abel_num_mc`; Lake: `desmear_max_iter`; GP: `gp_length_scale`, `gp_kernel`; smoothed Lake: `lake_smooth_decades` |
| Binning | `num_points` |

`readMyNXcanSAS` reads the group back into `Sample["provenance"]`, and copies
`MSAXSCorrection`, `calculated_qmin` and `calculated_mu` into
`CalibratedData`, `PeakToPeakTransmission` into `reducedData`, and the
beam-centre fit (`FWHM`, `Center`, `Chi-Square`, `wavelength`) from
`sasinstrument` into `reducedData`.  As a result, a cached Sample has the same
keys as a fresh one.  The exception is the per-point arrays of
`reducedData`, which are not stored.  Files written before this change have
no `sasprocess` and read as before.  The record is additive only: never
rename a key, because pyirena, DataReporter and Igor Irena read these files.

The NXcanSAS `run` field is the Bluesky uid (`/entry/entry_identifier`) when
the raw file has one; otherwise it is the file name without its extension.
`sastransmission_spectrum/Tdev` is NaN with the attribute
`note = "not estimated by Matilda"`, because the transmission uncertainty is
not computed yet.

**Overrides never touch raw data.** A thickness, transmission or Qmin
override is used in the reduction and recorded only in the reduced entry:
`sasdata/I@thickness`, `sassample/thickness` and `sasprocess`.  Before
v0.3.0 a thickness override rewrote `/entry/sample/thickness` (or the
step-scan `sample_thickness_mm`) in the raw file and kept the original once
as `<path>_original`.  Files edited that way keep the overridden value.
Matilda no longer writes or reads `_original`.

**`matilda/readfromtiled.py`**
Tiled server REST API queries.

Key functions: `FindLastScanData`, `FindLastBlankScan`, `FindScanDataByName`

Default server: `http://usaxscontrol.xray.aps.anl.gov:8000` / catalog: `usaxs_MongoDB`.
Server, port, catalog and timeout come from the site config (`tiled_*` keys).

**`matilda/config.py`**
Site configuration: `~/.matilda/config.json` (or `$MATILDA_CONFIG`), merged
over the built-in 12-ID-E defaults (`DEFAULTS`).  Importing it only reads the
file.  `ensure_config_file()` writes the defaults when the file is missing;
`main()` in the daemon calls it.  Bad keys, bad types or bad JSON are logged
and fall back to the defaults.  Standard library only.  User guide:
[configuration.md](configuration.md).

---

### Plotting

**`matilda/plotData.py`**
Headless matplotlib export to JPEG files in `imagePath`.

Key functions: `plotUSAXSResults`, `plotSWAXSResults`

> matplotlib is intentional here for headless/server output only.
> All interactive GUI plotting uses pyqtgraph (see `matilda/gui/`).

---

### GUI subpackage

**`matilda/gui/_qt.py`**
Single Qt import point for everything under `matilda/gui/`. PySide6, with a
PyQt6 fallback for environments that only have that binding; `Signal` is
normalised across the two. GUI modules do `from .._qt import QWidget, Qt,
Signal` — a direct `from PySide6… import` in a panel is a test failure
(`tests/test_gui_qt_shim.py`), as is any Qt import on the daemon path.

**`matilda/gui/sample_plate_setup.py`**
Sample Plate Setup GUI — replaces Igor Pro "Setup Sample Plates".

Launched via:
```bash
matilda-sample-plates        # installed console script
python -m matilda.gui.sample_plate_setup
```

Features: plate geometry templates, sample table, clickable plate canvas,
`.mac` export, EPICS motor survey, HDF5 save/load.
Full user guide: [docs/sample-plate-setup.md](sample-plate-setup.md).

**`matilda/gui/data_reduction/`**
Interactive data reduction GUI (USAXS/SAXS/WAXS).

Launched via:
```bash
matilda-gui                  # installed console script
```

Features: file tree with auto blank selection, per-technique parameter tabs
(thickness, μ-based calibration, per-gram normalisation, Qmin override),
pyqtgraph dual-axis log-log plot, batch processing, ASCII export.
Full user guide: [docs/matilda-gui.md](matilda-gui.md).

---

## Data dictionary schema

All `processXxx()` functions return a list of result dicts:

```python
Sample = {
    "RawData": {
        "filename": str,
        "samplename": str,      # from HDF5 sample/name
        "metadata": dict,
        "instrument": dict,
        "sample": dict,
        "control": dict,        # SAXS/WAXS only
        "data": np.ndarray,     # raw 2-D image (SAXS/WAXS only)
    },
    "reducedData": {
        "Q": np.ndarray,        # 1/Å
        "Intensity": np.ndarray,
        "Error": np.ndarray,
        "dQ": np.ndarray,       # USAXS only
    },
    "CalibratedData": {
        "Q": np.ndarray,        # 1/Å
        "Intensity": np.ndarray,  # cm²/cm³ (absolute scale)
        "Error": np.ndarray,
        "dQ": np.ndarray,
        "units": str,           # '[cm2/cm3]'
        "blankname": str,
        "thickness": float,     # mm, the value actually used (override / file / from μ)
        "MeasuredTransmission": float,  # USAXS: transmission used
        "MSAXSCorrection": float,       # USAXS
        "calculated_qmin": float,       # USAXS: auto Qmin, 1/Å
        "calculated_mu": float,         # 1/cm
        "Kfactor": float,       # USAXS only
        "OmegaFactor": float,   # USAXS only
        "SMR_Int": np.ndarray,  # USAXS: slit-smeared data
        "SMR_Qvec": np.ndarray,
        "SMR_Error": np.ndarray,
        "SMR_dQ": np.ndarray,
        "slitLength": float,
    },
    "BlankData": {              # only when blank was subtracted
        "Q": np.ndarray,
        "Intensity": np.ndarray,
        "Error": np.ndarray,
    },
    "calib2DData": {            # SAXS/WAXS only
        "data": np.ndarray,
        "blankname": str,
        "transmission": float,
        "thickness": float,     # mm, actually used
        "calculated_mu": float,
    },
    "provenance": dict,         # flat scalars written to <entry>/sasprocess/sasnote
}
```

---

## File layout

```
Matilda/
├── matilda/                    # Python package
│   ├── __init__.py             # Architecture overview, bug register
│   ├── matilda.py              # Main polling loop + entry point
│   ├── convertFlyscan.py       # USAXS flyscan reduction
│   ├── convertUSAXS.py         # USAXS step-scan reduction
│   ├── convertSWAXS.py         # SAXS/WAXS 2-D reduction
│   ├── desmearing.py           # Slit-smearing correction (Lake)
│   ├── desmearing_methods.py   # Method dispatcher (Abel default, Lake, GP)
│   ├── fx4support.py           # Counting-chain detection, FX4 helpers
│   ├── supportFunctions.py     # Shared numerical helpers
│   ├── supportNikaFunctions.py # Geometry conversion (Nika → pyFAI)
│   ├── technique_detector.py   # Fly/step/SAXS/WAXS detection (GUI + batch)
│   ├── hdf5code.py             # NXcanSAS HDF5 I/O
│   ├── plotData.py             # Headless JPEG plot export
│   └── gui/
│       ├── __init__.py
│       ├── _qt.py                  # Single Qt import point (PySide6)
│       ├── sample_plate_setup.py   # Sample Plate Setup GUI
│       └── data_reduction/         # matilda-gui data reduction application
├── tests/
│   ├── __init__.py
│   └── manual_test.py
├── docs/
│   ├── installation.md         # Conda install for end users
│   ├── matilda-gui.md          # Data reduction GUI user guide
│   ├── sample-plate-setup.md   # Sample plate setup GUI user guide
│   ├── service.md              # Beamline systemd service (staff only)
│   ├── operations.md           # pynika / pyirena scripting integrations
│   └── architecture.md         # This file
├── matilda_server.service      # systemd user unit file
├── serv_matilda.sh             # Service launch script
├── environment.yml             # conda environment spec
├── pyproject.toml              # Package metadata and entry points
└── README.md
```
