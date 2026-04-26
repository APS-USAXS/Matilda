# Matilda

Data-reduction tools for USAXS / SAXS / WAXS measurements collected at
APS beamline 12-ID-E (Argonne National Laboratory).

Matilda provides two desktop GUI tools for users who collect data at the
beamline (or work with data from the beamline) and a headless daemon for
beamline staff. End users normally only need the GUIs.

| Tool | Description |
|---|---|
| **matilda-gui** | Interactive USAXS/SAXS/WAXS data reduction. Open an HDF5 file, tune parameters, see calibrated I(Q) immediately, save back to NXcanSAS. |
| **matilda-sample-plates** | Sample-plate position editor. Place samples on a plate image, fill in names/thicknesses, export a Bluesky `.mac` command file for the beamline. |
| **matilda** (daemon) | Beamline service that polls the Tiled server every 15 s and reduces new data automatically. Runs only on `usaxscontrol.xray.aps.anl.gov`. |

---

## Install

```bash
conda create -n matilda python=3.12
conda activate matilda
pip install "matilda[gui] @ git+https://github.com/jilavsky/Matilda.git"
```

Full instructions, including platform notes and developer install, are in
[docs/installation.md](docs/installation.md).

---

## Launch a GUI

```bash
conda activate matilda

matilda-gui              # data reduction GUI
matilda-sample-plates    # sample plate setup GUI
```

---

## Documentation

| Document | Audience | Description |
|---|---|---|
| [docs/installation.md](docs/installation.md) | All users | Conda install, updates, platform notes |
| [docs/matilda-gui.md](docs/matilda-gui.md) | End users | Data reduction GUI user guide |
| [docs/sample-plate-setup.md](docs/sample-plate-setup.md) | End users | Sample plate setup GUI user guide |
| [docs/service.md](docs/service.md) | Beamline staff | Daemon, systemd service, log files |
| [docs/operations.md](docs/operations.md) | Beamline staff | Auto-calibration (pynika), auto-analysis and merging (pyirena) |
| [docs/architecture.md](docs/architecture.md) | Developers | Module map, data dictionary, file layout |

---

## Data reduction GUI — at a glance

`matilda-gui` is the main tool for interactive reduction. From an HDF5
file you can:

- Detect the technique (Flyscan, StepScan, SAXS, WAXS) automatically.
- Pick a blank manually or let the auto blank-finder choose the
  nearest-preceding blank scan.
- Override sample thickness, transmission, or Qmin per file.
- Switch the calibration model: standard (mm thickness),
  μ-based (`t = −ln(T) / μ`), or per-gram (output in cm²/g for powders).
- Tune desmearing parameters (Lake/Strobl, extrapolation method, iteration count).
- Adjust SAXS/WAXS azimuthal-integration bin count and angular range.
- See the slit-smeared and desmeared (or 2D-integrated) curves on a
  log-log plot with a D-spacing top axis.
- Re-reduce a folder full of files in batch (Process Selected / Process All).
- Export ASCII (4-column, NXcanSAS-aware) for plotting in other tools.

Reduced data are saved back into the original HDF5 file as NXcanSAS
groups, ready for downstream analysis.

Full reference: [docs/matilda-gui.md](docs/matilda-gui.md).

---

## Sample plate setup GUI — at a glance

`matilda-sample-plates` replaces the Igor Pro "Setup Sample Plates" tool.
It lets you:

- Choose a plate geometry (9×9 acrylic, NMR acrylic, old-style aluminium,
  NMR tubes, or generic grid) and click on hole positions in the plate
  image to drop samples into the table.
- Enter sample names, X/Y positions, thicknesses, and per-sample
  USAXS/SAXS/WAXS toggles.
- Estimate total run time from sample count, scan times, and stage move speed.
- Export a Bluesky `.mac` command file in any of four scan orderings
  (USAXS-SAXS-WAXS, USAXS-WAXS-SAXS, SAXS-WAXS-USAXS, USWAXS).
- Save and reload the full state (sample sets, plate, hooks, timing
  constants) as a portable HDF5 file.
- Optionally drive the sample stage live via EPICS (Beamline Survey),
  available when pyepics is installed and the beamline is reachable.

Full reference: [docs/sample-plate-setup.md](docs/sample-plate-setup.md).

---

## What the data reduction does

Calibrated absolute intensity from each technique:

- **USAXS Flyscan / StepScan**: HDF5 → normalisation → blank
  subtraction → absolute calibration (Kfactor / transmission) →
  Lake/Strobl iterative slit-smearing correction → I(Q) on a
  desmeared grid.
- **SAXS / WAXS area-detector**: HDF5 → 2D calibration (transmission,
  thickness, solid angle) → pyFAI azimuthal integration → I(Q).
- **NXcanSAS** output written back into the original HDF5 file with
  intensity-unit attributes (`[cm2/cm3]` or `[cm2/g]`) that depend on
  the calibration mode used.

---

## Beamline service (USAXS staff only)

When run as a service on `usaxscontrol.xray.aps.anl.gov`, Matilda polls
the Tiled server every 15 seconds, reduces new data automatically, runs
pynika to recalibrate detector geometry from `AgBehenateLaB6` standards,
and (when configured per folder) calls pyirena to fit models and merge
USAXS+SAXS curves.

```bash
# Service control on the beamline server
systemctl --user status matilda_server
systemctl --user restart matilda_server
journalctl --user -q -f -u matilda_server
```

See [docs/service.md](docs/service.md) for installation and operation of
the daemon, and [docs/operations.md](docs/operations.md) for the
pynika and pyirena integrations.

---

## Repository layout

```
matilda/                Python package
├── matilda.py          Daemon entry point and polling loop
├── convertFlyscan.py   USAXS flyscan reduction
├── convertUSAXS.py     USAXS step-scan reduction
├── convertSWAXS.py     SAXS / WAXS 2D reduction (pyFAI)
├── desmearing.py       Lake/Strobl slit desmearing
├── supportFunctions.py Shared numerical helpers
├── hdf5code.py         NXcanSAS read/write
├── plotData.py         Headless JPEG plot export
└── gui/
    ├── data_reduction/ matilda-gui application
    └── sample_plate_setup.py   matilda-sample-plates application

docs/                   User and developer documentation
tests/                  pytest suite and manual integration tests
matilda_server.service  systemd unit file (beamline service)
serv_matilda.sh         Service launch script
environment.yml         conda environment specification
pyproject.toml          Package metadata and console-script entry points
```

---

## Issues and feedback

Please file bug reports and feature requests on the
[GitHub issue tracker](https://github.com/jilavsky/Matilda/issues).
