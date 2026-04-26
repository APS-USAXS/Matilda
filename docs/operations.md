# Matilda — Scripting Integrations Guide

> **Audience:** Beamline staff who want to enable automatic detector
> calibration, model fitting, and USAXS+SAXS merging in the daemon.
>
> For end-user GUI guides see [matilda-gui.md](matilda-gui.md) and
> [sample-plate-setup.md](sample-plate-setup.md). For the service itself
> see [service.md](service.md).

## Overview

Matilda runs as a 15-second polling loop on `usaxscontrol.xray.aps.anl.gov`.
Every cycle it queries the Tiled server for new scan data, reduces any new files
to calibrated I(Q) curves, and saves JPEG summary plots to a web-visible directory.

Three optional scripting integrations extend the daemon without any code changes —
they are activated by configuration files placed in the data folders:

| Integration | Trigger | What it does |
|---|---|---|
| **pynika** auto-calibration | Filename contains `AgBehenateLaB6` | Fits detector geometry from the calibrant scan and pushes PONI parameters to EPICS PVs before reducing the file |
| **pyirena** auto-analysis | `pyirena_config.json` present in a technique folder | Runs model fitting on each newly reduced non-blank file |
| **pyirena** USAXS+SAXS merge | `merge_config.json` present in the parent data folder | Merges matching USAXS and SAXS I(Q) curves into a single dataset |

---

> For starting / stopping the service, log files, and configuration constants,
> see [service.md](service.md).

---

## pynika auto-calibration

When Matilda detects a scan whose filename contains `AgBehenateLaB6`
(case-insensitive) it automatically runs pynika to recalibrate the SAXS or
WAXS detector geometry and push results to EPICS PVs before reducing the file.

Calibration runs at most once per file per service session.
Separate tracking sets are kept for SAXS and WAXS.

```
Trigger: filename contains "AgBehenateLaB6"
Command: conda run --no-capture-output -p <PYNIKA_CONDA_ENV_PATH>
         pynika --file <path> --instrument SAXS|WAXS --auto-fit --save-to-pvs
Timeout: 300 seconds
```

---

## pyirena auto-analysis

After each data-reduction step (Flyscan, Step scan, SAXS, WAXS), Matilda checks
whether a `pyirena_config.json` file exists in the same folder as the reduced
data.  If it does, `pyirena.batch.fit_pyirena()` is called on every non-blank
file that has not already been analyzed this session.

This allows automatic model fitting (e.g., size distributions, form factors)
immediately after reduction, with no user intervention.

### Enabling auto-analysis

Place a `pyirena_config.json` file in each technique folder that should be
analyzed:

```
.../05_02_UserName/data/
    data_usaxs/
        pyirena_config.json    ← enables analysis for USAXS files
    data_saxs/
        pyirena_config.json    ← enables analysis for SAXS files
    data_waxs/
        pyirena_config.json    ← enables analysis for WAXS files
```

The config file format is defined by pyirena.  See the
[pyirena documentation](https://github.com/jilavsky/pyirena) for details.

### Subprocess call

```
Command: conda run --no-capture-output -p <PYIRENA_CONDA_ENV_PATH>
         python -c "from pyirena.batch import fit_pyirena;
                     fit_pyirena('<data_file>', '<config_file>',
                                 save_to_nexus=True, with_uncertainty=False, n_mc_runs=10)"
Timeout: 300 seconds
```

Blank files (filename containing "blank") are always skipped.  Each file is
analyzed at most once per service session, tracked by bounded in-memory sets
(max 100 entries per technique).

---

## Automatic USAXS+SAXS merging

When both a USAXS and a matching SAXS dataset exist for the same sample,
Matilda can automatically merge them into a single combined I(Q) curve using
`pyirena.batch.merge_data()`.

### Prerequisites

1. Both the USAXS and SAXS files must be reduced (present in their respective
   folders).
2. A `merge_config.json` file must exist in the **parent data folder** (one
   level above `data_usaxs` and `data_saxs`).

```
.../05_02_UserName/data/
    merge_config.json          ← required to enable merging
    data_usaxs/                ← reduced USAXS files (*.h5)
    data_saxs/                 ← reduced SAXS files (*.hdf)
    data_usaxs_merged/         ← output (created automatically)
```

### File matching logic

Two files are considered to represent the same sample when:

- The **prefix** before the first `_` is identical.
- The **scan number** between the last `_` and the file extension is identical.

The middle portion (timestamps, temperatures, etc.) is ignored because it may
differ between USAXS and SAXS measurements of the same sample.

| USAXS file | SAXS file | Match? |
|---|---|---|
| `Sample1_55C_10min_1234.h5` | `Sample1_58C_11min_1234.hdf` | Yes (`Sample1` + `1234`) |
| `Sample1_55C_10min_1234.h5` | `Sample2_55C_10min_1234.hdf` | No (prefix differs) |
| `Sample1_55C_10min_1234.h5` | `Sample1_55C_10min_1235.hdf` | No (scan number differs) |

### When merging is triggered

Merging is attempted after pyirena analysis in every processing cycle for
Flyscans, Step scans, and SAXS data.  It fires regardless of which technique
arrives first — if USAXS data arrive before SAXS, the merge is attempted when
SAXS data appear (and vice versa).

### Merge order

USAXS data are always **file 1** (lower Q, absolute intensity scale) and SAXS
data are always **file 2** (higher Q).  This is required by `merge_data()`.

### Output

Merged files are written to the `data_usaxs_merged` subfolder inside the parent
data directory.  The output filename matches the USAXS input filename.

### Post-merge pyirena analysis

After a successful merge, Matilda checks whether a `pyirena_config.json` file
exists in the `data_usaxs_merged` folder.  If it does, `fit_pyirena()` is called
on the merged file, enabling automatic model fitting on the combined dataset.

### Subprocess call

```
Command: conda run --no-capture-output -p <PYIRENA_CONDA_ENV_PATH>
         python -c "from pyirena.batch import merge_data;
                     merge_data('<usaxs_file>', '<saxs_file>',
                                config_file='<merge_config.json>',
                                output_folder='<data_usaxs_merged/>',
                                save_to_nexus=True, verbose=True)"
Timeout: 300 seconds
```

### Tracking and memory

Each USAXS+SAXS pair is merged at most once per service session, tracked by a
bounded in-memory set (max 100 entries).  Failed merges are also recorded — they
are not retried.  When the set exceeds 100 entries, the oldest entry is evicted;
re-merging an evicted pair is harmless (the output file is overwritten).

### merge_config.json format

The config file controls the overlap region and optimization parameters.
See the [pyirena merge documentation](https://github.com/jilavsky/pyirena/blob/main/docs/data_merge_gui.md)
for the full specification.  A minimal example:

```json
{
  "version": "1.0",
  "q_overlap_min": 0.08,
  "q_overlap_max": 0.25,
  "fit_scale": true,
  "scale_dataset": 2,
  "fit_qshift": false
}
```

---

## serv_matilda.sh

The service launch script performs three actions in order:

1. Sources the conda shell integration.
2. Activates the `matilda` conda environment.
3. Sets `MATILDA_LOG_DIR` and launches the polling loop via `python -m matilda.matilda`.

```bash
#!/bin/bash
source /APSshare/miniconda/x86_64/etc/profile.d/conda.sh
conda activate matilda
export MATILDA_LOG_DIR=/share1/log/matilda
cd /home/beams/USAXS/Apps/Matilda
exec python -m matilda.matilda
```

`exec` replaces the shell process with the Python process so that systemd
tracks the correct PID and `SIGTERM` reaches Python directly on `systemctl stop`.

---

## matilda_server.service

Key settings in the systemd unit file:

```ini
[Unit]
ConditionHost=usaxscontrol.xray.aps.anl.gov   # only starts on this host

[Service]
Restart=on-failure       # auto-restart if the process crashes
RestartSec=120           # wait 2 minutes before restarting
StartLimitBurst=5        # max 5 restarts within StartLimitInterval
StartLimitInterval=600   # … within 10 minutes, then give up
TimeoutStopSec=30        # kill after 30 s if it doesn't stop cleanly
```

The unit file lives at:
```
~/.config/systemd/user/matilda_server.service
```
The source copy in the repo is `matilda_server.service` at the repo root.

---

## Troubleshooting

### Service will not start

```bash
# Check for startup errors
systemctl --user status matilda_server
journalctl --user -xe -u matilda_server
```

Common causes:
- `matilda` conda environment not found — run `pip install -e .` inside the env
- `MATILDA_LOG_DIR` directory does not exist and cannot be created — check permissions
- Tiled server unreachable — logged as an error, service keeps running

### Reducing data manually (one-off re-run)

```bash
conda activate matilda
cd /home/beams/USAXS/Apps/Matilda
python - <<'EOF'
from matilda.matilda import processUSAXSFolder
processUSAXSFolder("/path/to/date/folder")
EOF
```

### Checking the version written into HDF5 files

```bash
python -c "from importlib.metadata import version; print(version('Matilda'))"
```
