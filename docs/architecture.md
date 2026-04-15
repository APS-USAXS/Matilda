# Matilda — Architecture Overview

## Big picture

```
Tiled server (HTTP/REST)
        │
        │  FindLastScanData()
        │  FindLastBlankScan()
        ▼
matilda.py  ─── 15-second polling loop
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
Main 15-second polling loop.  Queries Tiled for new scans, compares to the
previous list, dispatches reduction and plotting.

Entry point: `main()` — called by the `matilda` console script and by
`python -m matilda.matilda`.

Key module-level constants (configure near top of file):

| Constant | Purpose |
|---|---|
| `imagePath` | Directory for JPEG web plots |
| `NumberOfDaysToLookBack` | Tiled query window for scans |
| `NumberOfDaysToLookBackBlanks` | Tiled query window for blanks |
| `NumberOfImagesInGraphs` | Max datasets per plot |
| `CONDA_EXECUTABLE` | Path to conda for pynika/pyirena subprocesses |
| `PYNIKA_CONDA_ENV_PATH` | pynika environment path |
| `PYIRENA_CONDA_ENV_PATH` | pyirena environment path |
| `_PYIRENA_CONFIG_FILENAME` | Per-folder config that triggers pyirena analysis (`pyirena_config.json`) |
| `_MERGE_CONFIG_FILENAME` | Per-experiment config that triggers USAXS+SAXS merging (`merge_config.json`) |
| `_MERGED_FOLDER_NAME` | Output subfolder for merged data (`data_usaxs_merged`) |

---

### Data reduction

**`matilda/convertFlyscan.py`**
USAXS flyscan: HDF5 NXsas → normalised I(Q) → desmeared calibrated I(Q).

Pipeline:
`importFlyscan` → `calculatePD_Fly` → `beamCenterCorrection` →
`getBlankFlyscan` → `normalizeByTransmission` →
`calibrateAndSubtractFlyscan` → `desmearData` → `saveNXcanSAS`

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

---

### I/O

**`matilda/hdf5code.py`**
HDF5/NXcanSAS read-write helpers.

Key functions: `saveNXcanSAS`, `readMyNXcanSAS`, `readGenericNXcanSAS`,
`save_dict_to_hdf5`, `load_dict_from_hdf5`, `find_matching_groups`

The `Matilda_version` attribute in saved files is read from
`importlib.metadata.version('Matilda')` at write time.

**`matilda/readfromtiled.py`**
Tiled server REST API queries.

Key functions: `FindLastScanData`, `FindLastBlankScan`, `FindScanDataByName`

Target server: `http://usaxscontrol.xray.aps.anl.gov:8000` / catalog: `usaxs_MongoDB`

---

### Plotting

**`matilda/plotData.py`**
Headless matplotlib export to JPEG files in `imagePath`.

Key functions: `plotUSAXSResults`, `plotSWAXSResults`

> matplotlib is intentional here for headless/server output only.
> All interactive GUI plotting uses pyqtgraph (see `matilda/gui/`).

---

### GUI subpackage

**`matilda/gui/sample_plate_setup.py`**
PyQt6/pyqtgraph GUI tool replacing the Igor Pro "Setup Sample Plates" panel.

Launched via:
```bash
matilda-sample-plates        # installed console script
python -m matilda.gui.sample_plate_setup
```

Features: define sample positions, visualise on plate image, export
Bluesky `.mac` command files, save/load state in HDF5.
See the module docstring for full feature list.

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
        "thickness": float,     # mm
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
    },
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
│   ├── desmearing.py           # Slit-smearing correction
│   ├── supportFunctions.py     # Shared numerical helpers
│   ├── supportNikaFunctions.py # Geometry conversion (Nika → pyFAI)
│   ├── hdf5code.py             # NXcanSAS HDF5 I/O
│   ├── plotData.py             # Headless JPEG plot export
│   └── gui/
│       ├── __init__.py
│       └── sample_plate_setup.py   # Sample Plate Setup GUI
├── tests/
│   ├── __init__.py
│   └── manual_test.py
├── docs/
│   ├── installation.md         # Environment setup
│   ├── operations.md           # Service and CLI operations
│   └── architecture.md         # This file
├── matilda_server.service      # systemd user unit file
├── serv_matilda.sh             # Service launch script
├── environment.yml             # conda environment spec
├── pyproject.toml              # Package metadata and entry points
└── README.md
```
