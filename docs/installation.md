# Matilda — Installation

Matilda ships as a Python package with two GUI tools:

- **matilda-gui** — interactive USAXS/SAXS/WAXS data reduction
- **matilda-sample-plates** — sample-plate position editor and Bluesky command-file generator

This guide installs both for end users using conda.

---

## Requirements

- **Python:** 3.11 or 3.12
- **conda:** Miniconda or Anaconda
  ([download Miniconda](https://docs.conda.io/en/latest/miniconda.html))
- **Operating system:** Linux, macOS, or Windows
- **Internet access** for downloading dependencies

You do **not** need beamline access, EPICS, or pyirena/pynika to use the
GUI tools on your own data.

---

## Quick install (recommended)

```bash
# 1. Create a fresh conda environment for matilda
conda create -n matilda python=3.12
conda activate matilda

# 2. Install matilda with all GUI dependencies, straight from GitHub
pip install "matilda[gui] @ git+https://github.com/jilavsky/Matilda.git"
```

That installs the package and the two console scripts:

```bash
matilda-gui              # data reduction GUI
matilda-sample-plates    # sample plate setup GUI
```

> **zsh users (macOS):** quote the bracket expression so the shell does not
> try to expand it: `pip install "matilda[gui] @ git+https://github.com/jilavsky/Matilda.git"`

---

## Install from a local clone (developers and beamline staff)

```bash
git clone https://github.com/jilavsky/Matilda.git
cd Matilda

# Build the conda environment from the pinned spec
conda env create -f environment.yml
conda activate matilda

# Install matilda in editable mode — pick the right extras:
pip install -e .[gui]    # GUI tools (data reduction + sample plates)
pip install -e .         # headless / server only (no GUI)
pip install -e .[dev]    # adds pytest, pytest-cov
pip install -e .[all]    # everything
```

`environment.yml` pins the conda dependencies (numpy, scipy, h5py, pyFAI,
matplotlib, requests, tifffile). The `pip install -e .[gui]` step adds
PySide6, pyqtgraph, and pyepics on top.

> **zsh users (macOS):** quote the bracket expression: `pip install -e '.[gui]'`

---

## Updating an existing install

```bash
conda activate matilda

# From a clone:
git pull
conda env update -f environment.yml --prune
pip install -e .[gui]

# From the GitHub install (no clone):
pip install --upgrade --force-reinstall \
    "matilda[gui] @ git+https://github.com/jilavsky/Matilda.git"
```

---

## Verifying the installation

```bash
conda activate matilda

# Print the installed version
python -c "from importlib.metadata import version; print(version('Matilda'))"

# Quick import check (no network)
python -c "import matilda; print('OK')"

# Launch the GUIs (close the window when done)
matilda-gui
matilda-sample-plates
```

If the GUIs open, the install is complete.

---

## Console scripts reference

After `pip install`, the following commands are available in the active
conda environment:

| Command                  | Requires    | Description                                           |
|--------------------------|-------------|-------------------------------------------------------|
| `matilda-gui`            | `[gui]`     | Interactive data reduction GUI                        |
| `matilda-sample-plates`  | `[gui]`     | Sample plate setup GUI (Bluesky `.mac` generator)     |
| `matilda`                | core        | Headless polling daemon (beamline service only)       |

---

## Optional: scripting integrations

The GUI tools work standalone. If you also want **automatic data analysis
and merging**, install pyirena and pynika in their own conda environments
and configure the per-folder `pyirena_config.json` / `merge_config.json`
files. See [operations.md](operations.md) for details.

| Tool      | Purpose                                                      | Repo                                          |
|-----------|--------------------------------------------------------------|-----------------------------------------------|
| pyirena   | Model fitting and USAXS+SAXS data merging                    | https://github.com/jilavsky/pyirena           |
| pynika    | Detector geometry calibration from AgBehenate / LaB6 standards | https://github.com/jilavsky/pynika           |

These are only used by the **headless daemon** today; the GUIs do not call
them directly.

---

## Beamline server install (USAXS staff only)

For installing the headless daemon on `usaxscontrol.xray.aps.anl.gov`,
see [service.md](service.md).

---

## Platform notes

### Linux (RHEL 8 / Rocky 8 / CentOS Stream 8)

PySide6 6.8+ requires GLIBC 2.32; RHEL 8 ships GLIBC 2.28. The pin in
`pyproject.toml` (`PySide6>=6.4,<6.8`) handles this automatically.

### macOS

Do **not** install both PySide6 and PyQt6 in the same environment — they
conflict ("cocoa platform plugin not found"). The `[gui]` extra installs
PySide6 only.

### Windows

The default conda installation works; no extra steps required. Install
Miniconda for Windows from
[docs.conda.io](https://docs.conda.io/en/latest/miniconda.html) if you do
not already have conda.
