# Matilda GUI — User Guide

> Other guides: [sample-plate-setup.md](sample-plate-setup.md) |
> [installation.md](installation.md) | [operations.md](operations.md) |
> [service.md](service.md)

## Overview

The Matilda GUI (`matilda-gui`) is a desktop application for interactive
reduction of USAXS, SAXS, and WAXS data collected at APS beamline 12-ID-E.
It provides the same data reduction pipeline as the Matilda daemon but with
user-controlled parameters, visual feedback, and the ability to re-reduce
data with different settings.

---

## Launching

```bash
conda activate matilda
matilda-gui
```

Or from Python:

```python
from matilda.gui.data_reduction import main
main()
```

---

## Window layout

The main window is divided into three panels with a resizable splitter:

```
┌────────────────────────────────────────────────────────────────────┐
│  📁 Select Folder…     /path/to/data/folder                       │
├──────────────┬──────────────────┬──────────────────────────────────┤
│              │                  │        D (Å)  (top axis)         │
│  File Tree   │  Parameter Tabs  │                                  │
│              │                  │    ┌──────────────────────┐      │
│  (left)      │  Flyscan |       │    │                      │      │
│              │  StepScan |      │    │    Log-log plot       │      │
│              │  SAXS |          │    │                      │      │
│              │  WAXS            │    │    I (cm⁻¹)    I (arb)│     │
│              │                  │    └──────────────────────┘      │
│              │                  │        Q (Å⁻¹) (bottom axis)    │
├──────────────┴──────────────────┴──────────────────────────────────┤
│ ▶▶ Process All          Status: Idle          ▌▌▌▌▌  ✕ Cancel     │
│ Errors will appear here…                                          │
└────────────────────────────────────────────────────────────────────┘
```

---

## Quick start

1. Click **Select Folder** and choose a directory containing HDF5 data files.
2. The file tree populates with all `.h5`, `.hdf`, `.hdf5`, and `.nxs` files.
3. Click a file — the GUI auto-detects the technique (Flyscan, StepScan,
   SAXS, or WAXS) and switches to the matching parameter tab.
4. Adjust parameters if needed (or leave defaults).
5. **Double-click** a file to process it immediately, or select multiple files
   and click **Process Selected** on the parameter tab.
6. Results appear in the graph. Reduced data are saved to the HDF5 file
   automatically as NXcanSAS groups.

---

## File tree (left panel)

### Opening a folder

Click the **Select Folder** button at the top of the window. The tree
shows the folder hierarchy with HDF5 files. Subfolders load lazily when
expanded.

### Sorting and filtering

- **Sort dropdown**: Sort files by Name, Date modified, Scan number, or
  Temperature (extracted from filename patterns like `_25C_` or `_1234.h5`).
- **Filter text box**: Enter a substring or regex to show only matching
  filenames. Click the clear button to reset.

### Selecting files

- **Click** a file to select it and auto-detect the technique.
- **Ctrl+Click** or **Shift+Click** to select multiple files.
- The status bar shows the number of selected files and detected technique.

### Blank assignment

Blanks can be assigned in two ways:

**Automatic (default):** When the blank mode dropdown is set to
"auto (nearest preceding)", the system scans all loaded files for blank
files (those with "blank" in the filename) and selects the one with the
largest scan number that is still smaller than the sample's scan number.
This is the same logic the Matilda daemon uses.

**Manual:** Right-click a file in the tree to assign it as a blank:

- *Set as Flyscan Blank* / *Set as StepScan Blank* — for USAXS files
- *Set as SAXS Blank* — for SAXS files
- *Set as WAXS Blank* — for WAXS files

When a blank is manually assigned, the blank mode dropdown automatically
switches to "manual". Assigned blanks are highlighted in light blue in
the tree and listed in the summary at the bottom.

USAXS blanks are cross-compatible: a Flyscan blank can be used for
StepScan data and vice versa.

---

## Parameter tabs (center panel)

Each technique has its own tab. The GUI switches to the matching tab
when a file is selected. All parameters have sensible defaults that match
the Matilda daemon's behavior.

### Common parameters (all techniques)

#### Blank mode

- **auto (nearest preceding)**: Automatically finds the best blank from
  the loaded files (default).
- **manual**: Uses only the explicitly assigned blank.
- **Browse**: Opens a file dialog to select a blank from outside the tree.

#### Sample thickness

The HDF5 thickness value is displayed in grey when a file is selected.

- **Override checkbox**: Check to enter a custom thickness value (in mm)
  instead of the HDF5 value.

#### Calibration mode

Two checkboxes control alternative calibration methods:

- **Use μ for thickness**: Calculate thickness from the measured
  transmission and a user-provided linear absorption coefficient μ
  (in 1/cm). The formula is `t = −ln(T) / μ`. After processing, the
  calculated thickness is displayed next to the μ input.

  Use the Scattering Contrast Calculator in pyirena to obtain μ.

- **Normalize per gram**: Divide calibrated intensity by the solid-frame
  density (in g/cm³) to get intensity in [cm²/g] instead of [cm²/cm³].
  This is used for powder samples where geometric thickness is not
  meaningful. Requires "Use μ" to be checked.

  When active, the graph Y-axis label changes to "Intensity (cm²/g)"
  and the HDF5 file stores the correct unit attribute.

> **Note:** The μ-based method requires that the measured transmission is
> between 0 and 1. If transmission is outside this range (can happen with
> certain samples), the system falls back to the standard thickness and
> displays a warning in the error log.

**Calculated μ (after processing):** Whenever data is processed (in any
mode), the GUI computes μ = −ln(T) / thickness from the actual
transmission and thickness used, and loads the result into the μ
spinbox for reference. The spinbox stays disabled (read-only) until
"Use μ for thickness" is checked. This is informational — useful when
you want to know the absorption coefficient of the current sample.

#### Transmission override

- **Override transmission** (off by default): Replace the measured
  diode transmission with a manual value (0–1, dimensionless). Rarely
  needed — useful when the diode-measured transmission is bad or
  outside the valid range. The label next to the spinbox shows the
  transmission value actually used after processing (`used: 0.8523`).
  An info message is logged to the error panel each time the override
  is applied. The override persists across files until you uncheck it.

### Flyscan (USAXS) parameters

| Parameter | Default | Description |
|---|---|---|
| Min Q ratio | 1.05 | Threshold for Q-minimum selection after blank subtraction. Higher values trim more low-Q data. |
| Output points | 500 | Number of points after rebinning. |
| Max iterations | 20 | Maximum Lake/Strobl desmearing iterations. |
| Extrap method | PowerLaw w flat | High-Q extrapolation method for desmearing. Options: PowerLaw w flat, PowerLaw, Flat, Linear. |
| Extrap Q start | 0.1 Å⁻¹ | Q value above which the extrapolation is applied. |

#### Qmin override (USAXS only)

- **Override Qmin** (off by default): Manually truncate SMR data
  below a user-specified Q value (in 1/Å). This bypasses the
  auto-calculated Qmin (which uses Min Q ratio, instrument FWHM,
  and other experimental factors) and just removes points below the
  threshold. It is an emergency override for cases where blank
  subtraction goes bad at low Q and you want to salvage data at
  higher Q.

  After each processing run, the auto-calculated Qmin from the current
  dataset is displayed next to the spinbox in grey (`calc'd: 1.234e-04`).
  When the override is OFF, the spinbox value is automatically seeded
  with the calculated Qmin so you can see what the auto value is and
  enable the override to tweak it. When the override is ON, your
  entered value is preserved across processing runs (the spinbox is
  not overwritten). The override is automatically unchecked on every
  file selection because Qmin varies a lot between samples.

### StepScan (USAXS) parameters

Same as Flyscan (including the Qmin override) except no Output Points
control — step scans retain all measured points.

### SAXS parameters

| Parameter | Default | Description |
|---|---|---|
| Max number of points | unchecked | When checked, uses the detector dimension as bin count (same as WAXS). |
| Output Q points | 200 | Number of azimuthal integration bins. Disabled when "Max" is checked. |
| Azimuth range | −30° to +30° | Angular range for azimuthal integration. |

### WAXS parameters

WAXS always uses the maximum number of output points (detector dimension).
Additional parameters (pixel mask, normalization factors) will be added in
a future version.

---

## Processing data

### Single file

**Double-click** a file in the tree. The file is processed immediately with
the current parameter settings.

### Multiple files

1. Select files in the tree (Ctrl+Click or Shift+Click).
2. Click **Process Selected** on the active parameter tab.

### All files

Click **Process All Files** at the bottom of the window to process every
visible file in the tree.

### During processing

- A progress bar shows the current file count.
- The status label updates with each file.
- Click **Cancel** to stop after the current file completes.
- Warnings and errors appear in the error log at the bottom.

### Results

- Reduced data are **saved immediately** to the original HDF5 file as
  NXcanSAS groups. There is no separate "Save" step.
- If the file is processed again with `recalculateAllData=True` (always
  true in the GUI), old NXcanSAS groups are deleted and recreated.

---

## Graph (right panel)

The graph shows a dual-axis log-log plot of the reduction results.

### Axes

| Axis | Shows | Scale |
|---|---|---|
| Bottom (X) | Q in Å⁻¹ | Log (toggleable) |
| Top (X) | D-spacing in Å (d = 2π/Q) | Linked to Q axis |
| Left (Y) | Calibrated intensity in cm⁻¹ or cm²/g | Log (toggleable) |
| Right (Y) | Raw intensity in arbitrary units | Log |

### Curves

| Curve | Color | Style | Axis |
|---|---|---|---|
| Sample raw (R_data / Normalized) | Grey | Dashed | Right |
| Blank raw | Teal | Dotted | Right |
| Slit-smeared (SMR) | Blue | Solid | Left |
| Desmeared (DSM) / Calibrated | **Red** | Solid | Left |

The final result (Desmeared for USAXS, Calibrated for SAXS/WAXS) is
always shown in red.

### Right axis scaling

- USAXS (Flyscan/StepScan): up to 10 decades displayed.
- SAXS/WAXS: up to 5 decades displayed.
- Full data range is shown from top to bottom.

### Toolbar

| Button | Description |
|---|---|
| X log | Toggle X-axis log scale (default: on) |
| Y log | Toggle Y-axis log scale (default: on) |
| Error bars | Show/hide error bars on calibrated curves (default: off). Error bars are colored to match their curve. |
| Clear | Remove all curves |
| View 2D Image | Open 2D detector image viewer (SAXS/WAXS only) |
| Export CSV | Save displayed curves to a CSV file |
| Export PNG | Save the plot as a PNG image |

### 2D image viewer

For SAXS and WAXS files, clicking **View 2D Image** opens a separate
window showing the raw 2D detector image. Features:

- Log intensity toggle (default: on)
- Draggable histogram for adjusting color levels
- Reset levels button

---

## Data flow

The reduction pipeline in the GUI follows the same steps as the daemon:

### USAXS (Flyscan / StepScan)

```
HDF5 file
  ↓ import raw data
Raw detector arrays (R_data)          → plotted on right axis
  ↓ beam center correction, smoothing
  ↓ blank subtraction + calibration
  ↓   thickness: from HDF5, override, or μ-based
  ↓   blank: auto-selected or manually assigned
Slit-smeared I(Q) [SMR]              → plotted on left axis (blue)
  ↓ rebinning (Flyscan: N points; StepScan: all points kept)
  ↓ Lake/Strobl desmearing
Desmeared I(Q) [DSM]                 → plotted on left axis (red)
  ↓ optional: divide by density for [cm²/g]
  ↓ save NXcanSAS to HDF5
```

### SAXS / WAXS

```
HDF5 file
  ↓ import raw 2D image
  ↓ blank subtraction + 2D calibration
  ↓   thickness: from HDF5, override, or μ-based
  ↓   blank: auto-selected or manually assigned
Calibrated 2D image
  ↓ azimuthal integration (pyFAI)
  ↓   SAXS: N bins or max, azimuth range
  ↓   WAXS: always max bins
Calibrated I(Q)                       → plotted on left axis (red)
  ↓ optional: divide by density for [cm²/g]
  ↓ save NXcanSAS to HDF5
```

### Intensity units

| Mode | Unit | HDF5 attribute |
|---|---|---|
| Standard (thickness) | cm²/cm³ (= cm⁻¹) | `[cm2/cm3]` |
| μ-based | cm²/cm³ (= cm⁻¹) | `[cm2/cm3]` |
| Per gram (powder) | cm²/g | `[cm2/g]` |

---

## Error log

The error log at the bottom of the window shows:

- **Processing errors**: exceptions raised during file reduction, prefixed
  with the filename.
- **Warnings**: messages from the converter code (e.g., transmission out of
  range for μ calculation, insufficient data points for desmearing).

The log accumulates messages across processing runs and holds up to 200
lines. It is cleared at the start of each new processing batch.

---

## Keyboard and mouse reference

| Action | Effect |
|---|---|
| Double-click file in tree | Process that file immediately |
| Right-click file in tree | Context menu for blank assignment |
| Click file in tree | Select file, auto-detect technique, show HDF5 thickness |
| Ctrl+Click / Shift+Click | Multi-select files |
| Scroll wheel on graph | Zoom in/out |
| Click+drag on graph | Pan |
| Right-click on graph | pyqtgraph context menu (export, transform, etc.) |

---

## Session persistence

The GUI remembers the last-used folder and panel splitter sizes between
sessions. This state is stored in `~/.matilda_gui_session.json`.

---

## Tips

- **Re-reducing data**: The GUI always forces full recalculation
  (`recalculateAllData=True`). Old cached NXcanSAS groups in the HDF5
  file are deleted and recreated with the new parameters.

- **Checking the blank**: When using auto blank mode, the blank label
  on the parameter tab updates to show "auto: filename" after processing.
  Verify this is the expected blank file.

- **pyFAI integrator caching**: For SAXS/WAXS files in the same folder
  (same detector geometry), the pyFAI integrator is cached and reused,
  significantly speeding up batch processing.

- **Large datasets**: Use the filter box to narrow down the file list
  before processing. Sort by scan number or date to find specific files.

- **Comparing results**: Process a file, note the graph, then change
  parameters and process again. The graph clears and redraws with each
  new result.

---

## See also

- [sample-plate-setup.md](sample-plate-setup.md) — prepare sample positions and export Bluesky `.mac` command files
- [installation.md](installation.md) — conda install instructions
- [operations.md](operations.md) — pynika auto-calibration, pyirena auto-analysis and merging (daemon only)
- [service.md](service.md) — beamline systemd service (staff only)
