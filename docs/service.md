# Matilda — Beamline Service Guide

This document covers installing and operating the Matilda headless daemon
on the beamline control server `usaxscontrol.xray.aps.anl.gov`.

End users who only need the desktop GUI tools should see
[installation.md](installation.md) instead.

---

## What the daemon does

The `matilda` daemon runs as a systemd user service. Every 15 seconds it:

1. Queries the Tiled server for new scan data.
2. Reduces any new files to calibrated 1-D I(Q) curves.
3. Optionally runs pynika to recalibrate detector geometry
   (triggered by `AgBehenateLaB6` filenames).
4. Optionally runs pyirena for automatic model fitting and USAXS+SAXS
   data merging (configured per-folder with JSON files).
5. Saves JPEG summary plots to a web-visible directory for live monitoring.

For the pynika and pyirena scripting details see [operations.md](operations.md).

---

## Install (server)

```bash
# Conda paths on the beamline server
# conda:           /APSshare/miniconda/x86_64/bin/conda
# Environment:     /home/beams/USAXS/.conda/envs/matilda

# Build the environment from scratch
cd /home/beams/USAXS/Apps/Matilda
/APSshare/miniconda/x86_64/bin/conda env create -f environment.yml

# Activate and install Matilda (headless — no GUI extras needed)
conda activate matilda
pip install -e .

# Update an existing environment
conda env update -f environment.yml --prune
pip install -e .
```

---

## Install the systemd service (one-time)

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

# Verify
loginctl show-user $USER --property=Linger
```

---

## Daily operations

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

## Running from the command line (one-off)

All of the following must be run from the repo root with the `matilda`
conda environment active.

```bash
cd /home/beams/USAXS/Apps/Matilda
conda activate matilda

# Option 1 — installed console script (preferred)
matilda

# Option 2 — module invocation (no install beyond pip install -e . needed)
python -m matilda.matilda
```

> **Note:** The old direct-script launch (`python matilda/matilda.py`) no
> longer works because imports are now relative. Use one of the two forms above.

Stop with **Ctrl-C** (KeyboardInterrupt is caught and logged cleanly).

---

## Log files

Matilda writes a rotating log file. The location is controlled by the
`MATILDA_LOG_DIR` environment variable.

| Environment | Log path |
|---|---|
| Beamline service | `/share1/log/matilda/matilda.log` |
| Dev machine (default) | `~/.local/share/matilda/log/matilda.log` |
| Custom | `export MATILDA_LOG_DIR=/your/path` before launching |

The rotating log keeps **4 files × 1 MB = 4 MB maximum** on disk.

```bash
# Tail the live log on the beamline server
tail -f /share1/log/matilda/matilda.log

# Or via journalctl (captures stdout/stderr from the service)
journalctl --user -q -f -u matilda_server
```

---

## Configuration

Runtime parameters are constants near the top of `matilda/matilda.py`.
Edit them and restart the service for changes to take effect.

| Variable | Default | Description |
|---|---|---|
| `imagePath` | `/home/joule/WEBUSAXS/www_live/` | JPEG summary plot directory. Set to `None` to disable. |
| `NumberOfDaysToLookBack` | `1` | How far back (days) to search for new scans. |
| `NumberOfDaysToLookBackBlanks` | `5` | How far back (days) to search for blank scans. |
| `NumberOfImagesInGraphs` | `10` | Maximum datasets per summary plot. |
| `CONDA_EXECUTABLE` | `/APSshare/miniconda/x86_64/bin/conda` | Full path to conda — used when invoking pynika and pyirena. |
| `PYNIKA_CONDA_ENV_PATH` | `/home/beams/USAXS/.conda/envs/pynika` | pynika conda environment for auto-calibration. |
| `PYIRENA_CONDA_ENV_PATH` | `/home/beams/USAXS/.conda/envs/pyirena` | pyirena conda environment for auto-analysis and merging. |

---

## serv_matilda.sh

The service launch script:

1. Sources the conda shell integration.
2. Activates the `matilda` conda environment.
3. Sets `MATILDA_LOG_DIR` and launches the polling loop.

```bash
#!/bin/bash
source /APSshare/miniconda/x86_64/etc/profile.d/conda.sh
conda activate matilda
export MATILDA_LOG_DIR=/share1/log/matilda
cd /home/beams/USAXS/Apps/Matilda
exec python -m matilda.matilda
```

`exec` replaces the shell so systemd tracks the correct PID and
`SIGTERM` reaches Python directly on `systemctl stop`.

---

## matilda_server.service — key settings

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

The unit file lives at `~/.config/systemd/user/matilda_server.service`.
The source copy in the repo is `matilda_server.service` at the repo root.

---

## Troubleshooting

### Service will not start

```bash
systemctl --user status matilda_server
journalctl --user -xe -u matilda_server
```

Common causes:
- `matilda` conda environment not found — run `pip install -e .` inside the env.
- `MATILDA_LOG_DIR` directory does not exist and cannot be created — check permissions.
- Tiled server unreachable — logged as an error, service keeps running.

### Re-reducing data manually

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

---

## See also

- [operations.md](operations.md) — pynika auto-calibration, pyirena auto-analysis,
  automatic USAXS+SAXS merging
- [installation.md](installation.md) — end-user GUI install
- [architecture.md](architecture.md) — module map and data dictionary
