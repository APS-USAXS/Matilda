# Matilda — Installation

## Requirements

- Python 3.11 or 3.12
- conda (Miniconda or Anaconda)

---

## Beamline server (usaxscontrol.xray.aps.anl.gov)

The server runs headless — no GUI stack needed.

```bash
# Conda install location
/APSshare/miniconda/x86_64/bin/conda

# Environment location
/home/beams/USAXS/.conda/envs/matilda

# Build the environment from scratch
conda env create -f environment.yml

# Activate and install Matilda in editable mode (headless only)
conda activate matilda
cd /home/beams/USAXS/Apps/Matilda
pip install -e .

# Update an existing environment
conda env update -f environment.yml --prune
pip install -e .
```

---

## Local development (Windows)

```bash
# Conda base: C:\ProgramData\miniconda3
# Environments: C:\Users\ilavsky\.conda\envs\

# Build fresh
conda env create -f environment.yml
conda activate matilda

# Install in editable mode — choose one:
pip install -e .          # core (no GUI)
pip install -e .[gui]     # + PyQt6 and pyqtgraph
pip install -e .[dev]     # + pytest
pip install -e .[all]     # everything

# Note: zsh users must quote the bracket expression:
pip install -e '.[all]'
```

---

## Installed console scripts

After any of the `pip install` variants above the following commands are
available in the active environment:

| Command | Requires | Description |
|---|---|---|
| `matilda` | core | Start the 15-second polling service loop |
| `matilda-sample-plates` | `[gui]` | Launch the Sample Plate Setup GUI tool |

---

## Verifying the installation

```bash
# Check the installed version
python -c "from importlib.metadata import version; print(version('Matilda'))"

# Quick import check (no network needed)
python -c "import matilda; print('OK')"
```
