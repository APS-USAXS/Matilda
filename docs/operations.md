# Matilda — Operations Guide

## Overview

Matilda runs as a 15-second polling loop on `usaxscontrol.xray.aps.anl.gov`.
Every cycle it queries the Tiled server for new scan data, reduces any new files
to calibrated I(Q) curves, and saves JPEG summary plots to a web-visible directory.

---

## Running from the command line

All of the following must be run from the repo root with the `matilda` conda
environment active.

```bash
cd /home/beams/USAXS/Apps/Matilda
conda activate matilda

# Option 1 — installed console script (preferred after pip install -e .)
matilda

# Option 2 — module invocation (equivalent, no install needed beyond pip install -e .)
python -m matilda.matilda
```

> **Note:** The old direct-script launch (`python matilda/matilda.py`) no longer
> works because module imports are now relative.  Use one of the two forms above.

Stop with **Ctrl-C** (KeyboardInterrupt is caught and logged cleanly).

---

## Running as a systemd user service

The service is managed by systemd as a **per-user** service — no `sudo` required
for day-to-day operations.

### Install the service unit (one-time setup)

```bash
# Copy the unit file to the systemd user config directory
cp /home/beams/USAXS/Apps/Matilda/matilda_server.service \
   ~/.config/systemd/user/matilda_server.service

# Tell systemd to read the new file
systemctl --user daemon-reload

# Enable the service to start automatically at login / after reboot
systemctl --user enable matilda_server
```

### Allow the service to survive logout (one-time, requires admin)

By default a user service stops when the user logs out.
To keep Matilda running continuously:

```bash
# An administrator must run this once for the USAXS user
sudo loginctl enable-linger USAXS

# Verify linger is active
loginctl show-user $USER --property=Linger
```

### Daily operations

| Action | Command |
|---|---|
| Start | `systemctl --user start matilda_server` |
| Stop | `systemctl --user stop matilda_server` |
| Restart | `systemctl --user restart matilda_server` |
| Status | `systemctl --user status matilda_server` |
| Follow live logs | `journalctl --user -q -f -u matilda_server` |
| Reload unit file after editing | `systemctl --user daemon-reload` |
| Enable auto-start | `systemctl --user enable matilda_server` |
| Disable auto-start | `systemctl --user disable matilda_server` |

---

## Log files

Matilda writes a rotating log file.  The location is controlled by the
`MATILDA_LOG_DIR` environment variable.

| Environment | Log path |
|---|---|
| Beamline service | `/share1/log/matilda/matilda.log` |
| Dev machine (default) | `~/.local/share/matilda/log/matilda.log` |
| Custom | Set `export MATILDA_LOG_DIR=/your/path` before launching |

The rotating log keeps **4 files × 1 MB = 4 MB maximum** on disk.
Older entries are automatically discarded when the active log reaches 1 MB.

```bash
# Tail the live log on the beamline server
tail -f /share1/log/matilda/matilda.log

# Or via journalctl (also captures stdout/stderr from the service)
journalctl --user -q -f -u matilda_server
```

To override the log directory without editing any file, set the environment
variable before starting Matilda:

```bash
export MATILDA_LOG_DIR=/my/custom/log
matilda
```

For the service, the variable is set in [serv_matilda.sh](../serv_matilda.sh).

---

## Configuration

Runtime parameters are constants near the top of `matilda/matilda.py`.
Edit them and restart the service for changes to take effect.

| Variable | Default | Description |
|---|---|---|
| `imagePath` | `/home/joule/WEBUSAXS/www_live/` | Directory where JPEG summary plots are saved. Set to `None` to disable plotting. |
| `NumberOfDaysToLookBack` | `1` | How far back (days) to search for new scans. |
| `NumberOfDaysToLookBackBlanks` | `5` | How far back (days) to search for blank scans. |
| `NumberOfImagesInGraphs` | `10` | Maximum number of data sets shown in each summary plot. |
| `CONDA_EXECUTABLE` | `/APSshare/miniconda/x86_64/bin/conda` | Full path to conda — used when invoking pynika calibration. |
| `PYNIKA_CONDA_ENV_PATH` | `/home/beams/USAXS/.conda/envs/pynika` | pynika conda environment used for auto-calibration. |

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
