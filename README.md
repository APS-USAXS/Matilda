# Matilda

Live data-processing package for the APS USAXS/SAXS/WAXS instrument at
sector 12-ID-E, Argonne National Laboratory.

Matilda runs as a systemd user service on `usaxscontrol.xray.aps.anl.gov`.
Every 15 seconds it polls the Tiled server for new scan data, reduces any new
files to calibrated 1-D I(Q) curves, optionally runs pynika for detector-geometry
recalibration, and saves summary JPEG plots to a web-visible directory.

## Documentation

| Document | Description |
|---|---|
| [docs/installation.md](docs/installation.md) | Environment setup and pip install steps |
| [docs/operations.md](docs/operations.md) | Running as a service and from the command line |
| [docs/architecture.md](docs/architecture.md) | Module overview and data dictionary |

## Quick start

```bash
# On the beamline server — the service manages itself:
systemctl --user status matilda_server
systemctl --user restart matilda_server
journalctl --user -q -f -u matilda_server

# Manual one-shot run from the repo root:
conda activate matilda
python -m matilda.matilda
```

## Features

- Automated 15-second polling loop running as a systemd user service
- USAXS flyscan and step-scan: HDF5 → normalised and desmeared I(Q)
- SAXS/WAXS area-detector: HDF5 → pyFAI azimuthal integration → I(Q)
- Automatic blank (background) selection and subtraction
- Absolute intensity calibration (Kfactor / transmission)
- Iterative slit-smearing correction (Lake/Strobl) for USAXS data
- pynika auto-calibration triggered by AgBehenateLaB6 standard scans
- JPEG summary plots saved to web-visible directory for live monitoring
- NXcanSAS-compliant HDF5 output

## GUI tools

```bash
# Sample Plate Setup (requires pip install matilda[gui])
matilda-sample-plates
```

## Repository layout

```
matilda/            Python package (reduction, I/O, plotting, GUI)
tests/              pytest suite and manual integration tests
docs/               Operations and architecture documentation
matilda_server.service   systemd unit file
serv_matilda.sh     Service launch script
environment.yml     conda environment specification
pyproject.toml      Package metadata and console-script entry points
```
