# Matilda — Installation

Matilda ships as a Python package with two GUI tools:

- **matilda-gui** — interactive USAXS/SAXS/WAXS data reduction
- **matilda-sample-plates** — sample-plate position editor and Bluesky command-file generator

Follow the steps below from the top. Each step builds on the previous one.

---

## Step 1 — Install conda

Matilda uses conda to manage its Python environment and dependencies. If you
already have Miniconda or Anaconda installed, skip to Step 2.

**Download Miniconda** (the small, recommended installer):
https://docs.conda.io/en/latest/miniconda.html

Choose the installer for your operating system (Windows, macOS, or Linux)
and run it. Accept the defaults — there is no need to add conda to PATH
system-wide; the installer creates a "Miniconda" shortcut / terminal that
has conda ready to use.

> **Windows users:** after installation, open the
> **Anaconda Prompt (Miniconda)** from the Start menu for all commands below.
>
> **macOS / Linux users:** open a regular terminal. After installation run
> `conda init` once to enable the `conda activate` command in your shell,
> then restart the terminal.

---

## Step 2 — Install Git

Git is used to download the Matilda source code. If `git --version` works in
your terminal, skip to Step 3.

- **Windows:** download from https://git-scm.com/download/win and install
  with default options.
- **macOS:** run `xcode-select --install` in a terminal (installs the Xcode
  command-line tools which include Git), or install via Homebrew: `brew install git`.
- **Linux:** `sudo apt install git` (Debian/Ubuntu) or
  `sudo dnf install git` (RHEL/Rocky/Fedora).

---

## Step 3 — Clone the repository

Open a terminal (or Anaconda Prompt on Windows) and choose a folder where
you want to keep the Matilda source. Then run:

```bash
git clone https://github.com/APS-USAXS/Matilda.git
cd Matilda
```

This creates a `Matilda/` folder containing the source code. All remaining
steps must be run from inside that folder.

---

## Step 4 — Create the conda environment

The repository includes an `environment.yml` file that lists all required
packages. Create the environment with:

```bash
conda env create -f environment.yml
```

This takes a few minutes on first run (it downloads numpy, scipy, h5py,
pyFAI, and other packages). You only need to do this once.

Activate the new environment:

```bash
conda activate matilda
```

Your terminal prompt should now show `(matilda)` at the start.

---

## Step 5 — Install Matilda

With the environment active and the terminal still inside the `Matilda/`
folder, run:

```bash
pip install -e .[gui]
```

The `-e` flag installs Matilda in **editable mode** — the source folder
you cloned is the live package, so any future `git pull` updates take
effect immediately without reinstalling.

The `[gui]` part adds the GUI-specific packages (PySide6 and pyqtgraph).

> **macOS / Linux (zsh shell):** quote the bracket so the shell does not
> try to interpret it: `pip install -e '.[gui]'`

---

## Step 6 — Verify and launch

```bash
# Confirm the install
python -c "from importlib.metadata import version; print(version('Matilda'))"

# Launch the data reduction GUI
matilda-gui

# Launch the sample plate setup GUI
matilda-sample-plates
```

If both windows open, the installation is complete.

---

## Updating to a newer version

When a new version of Matilda is released, update your local copy:

```bash
conda activate matilda
cd /path/to/Matilda       # the folder you cloned in Step 3

git pull                  # download the latest changes
conda env update -f environment.yml --prune   # update packages if needed
pip install -e .[gui]     # reinstall (picks up any new entry points)
```

> **macOS / Linux (zsh):** `pip install -e '.[gui]'`

---

## Console scripts reference

After Step 5, the following commands are available whenever the `matilda`
environment is active:

| Command | Description |
|---|---|
| `matilda-gui` | Interactive data reduction GUI |
| `matilda-sample-plates` | Sample plate position editor and `.mac` file generator |
| `matilda` | Headless polling daemon (beamline service only — not needed by end users) |

---

## Optional: automatic data analysis integrations

The GUI tools work completely standalone. If you also want **automatic model
fitting and USAXS+SAXS data merging** (daemon mode only), install pyirena
and pynika in their own conda environments. See [operations.md](operations.md)
for details.

| Tool | Purpose | Repository |
|---|---|---|
| pyirena | Model fitting and USAXS+SAXS merging | https://github.com/jilavsky/pyirena |
| pynika | Detector geometry calibration from calibrant scans | https://github.com/jilavsky/pynika |

---

## Beamline server install (USAXS staff only)

For the headless daemon on `usaxscontrol.xray.aps.anl.gov`, see
[service.md](service.md).

---

## Platform notes

### Windows

Use the **Anaconda Prompt (Miniconda)** from the Start menu for all commands.
PowerShell and regular Command Prompt work too once conda is initialised, but
the Anaconda Prompt is the easiest starting point.

### macOS

Do **not** install both PySide6 and PyQt6 in the same environment — they
conflict at runtime ("cocoa platform plugin not found"). The `[gui]` extra
installs PySide6 only, so this is handled automatically.

### Linux (RHEL 8 / Rocky 8 / CentOS Stream 8)

PySide6 6.8+ requires GLIBC 2.32; RHEL 8 ships GLIBC 2.28. The version
pin in `pyproject.toml` (`PySide6>=6.4,<6.8`) handles this automatically —
no manual action needed.
