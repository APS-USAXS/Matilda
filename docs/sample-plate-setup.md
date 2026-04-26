# Sample Plate Setup — User Guide

`matilda-sample-plates` is a desktop tool for laying out USAXS/SAXS/WAXS
sample plates and exporting a Bluesky `.mac` command file that the
beamline runs to collect data. It replaces the Igor Pro
"Setup Sample Plates" tool (`IN3_SamplePlate.ipf`).

The tool runs anywhere — you don't need beamline access — except for the
optional **Beamline Survey** dialog, which requires EPICS connectivity.

---

## Launching

```bash
conda activate matilda
matilda-sample-plates
```

Or from Python:

```python
from matilda.gui.sample_plate_setup import run_sample_plate_setup
run_sample_plate_setup()
```

---

## Window layout

The main window is organised into three tabs:

| Tab | Purpose |
|---|---|
| **Sample Table** | Sample list + clickable plate image |
| **Option Controls** | Scan times, timing constants, hook function (offset re-runs) |
| **Export Controls** | Export order, command file name, save mode (current set or list of sets) |

A toolbar above the tabs holds the set name, global USAXS/SAXS/WAXS
checkboxes, runtime estimate, and save/load-settings buttons.

A menu bar provides **File → New / Open HDF5 / Save HDF5 /
Save HDF5 As / Import 12IDB command file / Exit** and **Help → About**.

---

## Quick start

1. Pick a **Template** (e.g., `9x9 Acrylic/magnetic plate`) and click
   **Create New Set**. The plate image appears on the right.
2. Select a row in the table (left panel).
3. Click a hole on the plate image — the row's SX / SY are filled in.
4. Type a **Sample Name** and **Thickness** for the row.
5. Repeat for the rest of your samples.
6. Click **Export usaxs.mac** to write the command file.
7. Save your work with **File → Save HDF5 As…** so you can reload later.

---

## Sample Table tab

### Template and plate selection

The **Template** dropdown lists the available plate geometries:

| Template | Holes | Layout |
|---|---|---|
| 9×9 Acrylic/magnetic plate | 81 | 9 columns × 9 rows, 20 mm spacing |
| NMR Acrylic plate | 100 | 20 columns × 5 rows |
| Old Style Al Plate | 60 | 4 columns × 15 holes per column |
| NMR Tubes holder | 20 | One row of 20 NMR tubes |
| Generic Grid | user-defined | Set rows, columns, spacing in dialog |
| AgBehenateLaB6 | 1 | Single calibration position |

Set the desired **#** (sample count) and click **Create New Set** to
populate the table.

**Add Lines** appends additional empty rows without overwriting existing data.

### Sample table columns

| Column | Description |
|---|---|
| Sample Name | Free text. Sanitised to `[A-Za-z0-9_]+` on entry (`%` becomes `pct`, other unsupported characters become `_`). |
| X [mm] | Sample stage X position. |
| Y [mm] | Sample stage Y position. |
| Thick [mm] | Sample thickness. Defaults to **1.0 mm** if blank or zero. |
| USAXS / SAXS / WAXS | Per-row checkboxes (overridden by the global "all" toggles in the toolbar). |
| MetaData | Free text — written into the command file as a comment. |

Empty `Sample Name` rows are **skipped** during export, so you can leave
gaps in the table.

### Plate canvas (right panel)

The plate image renders the chosen template's hole positions. The currently
selected row is highlighted on the plate.

- **Click** a hole → fills the SX/SY of the selected row.
- **Mouse-over** → status bar shows the position in mm.

### Right-click menu (sample table)

Right-click a row to access:

- Insert Row Above / Below
- Delete Row, Duplicate Row
- Move Row Up / Down
- Cut / Copy / Paste
- **Increment SX / SY (stepped)** — auto-fill positions in a sequence
- **Add to SX / SY** — apply an offset to a range of rows
- **Write same SX / SY / Thickness** — broadcast a value to selected rows
- **Write same name** — broadcast a name to selected rows
- **Set as Blank** — name the row "Blank"
- **Set as AgBehenateLaB6** — name the row "AgBehenateLaB6" (triggers
  pynika auto-calibration when the daemon picks it up later)
- **Clear Row**

### Buttons under the table

| Button | What it does |
|---|---|
| Save Set | Save the current sample set as a named entry in the in-memory list of saved sets. |
| Load Set | Load a previously saved set into the table. |
| Preview Commands | Show the generated `.mac` content in a dialog before exporting. |
| **Export usaxs.mac** | Write the command file to disk (default name `usaxs.mac`). |
| Beamline Survey | Open the EPICS motor control dialog (greyed out if pyepics is unavailable). |

---

## Option Controls tab

### Scan Times

| Field | Default | Used for |
|---|---|---|
| USAXS scan time | 90 s | Per-sample USAXS exposure time |
| SAXS scan time | 1 s | Per-sample SAXS exposure time |
| WAXS scan time | 3 s | Per-sample WAXS exposure time |
| Default sample thickness | 1.0 mm | Used when a row's thickness is blank / zero |

These values feed the **Est. time** label in the toolbar.

### Timing Constants

Used only by the runtime estimator (the values written to the `.mac`
file are unaffected):

| Field | Default | Description |
|---|---|---|
| USAXS / SAXS / WAXS overhead | 15 / 5 / 3 s | Per-scan overhead beyond the exposure time |
| Sample move speed | 30 mm/s | Aerotech stage average speed |
| Geometry switch time | 15 s | USAXS↔SAXS↔WAXS swap |
| Retune interval | 600 s | Time between SAXS/WAXS retunes |
| USAXS retune every N scans | 3 | Frequency of USAXS retune |
| USAXS retune time | 20 s | Per retune |
| SAXS/WAXS retune time | 0 s | Disabled on 12-ID-E |

### Export Hook Function

Repeats the export at offset positions to build edge / centre maps.

The four default rows (Right `_R`, Top `_T`, Left `_L`, Bottom `_B`)
each shift by ±1 mm in SX or SY and append a suffix to the sample name.

**Enable hook function in exports** turns the feature on; individual
rows can be enabled or disabled separately. Each enabled row generates
an extra block in the command file with the offset applied to every
sample's SX / SY and the suffix appended to every sample name.

---

## Export Controls tab

### Export order

Choose how the scan blocks are interleaved in the `.mac` file:

| Order | Effect |
|---|---|
| USAXS-SAXS-WAXS | Default — three blocks, one per technique |
| USAXS-WAXS-SAXS | USAXS, then WAXS, then SAXS |
| SAXS-WAXS-USAXS | SAXS, WAXS, then USAXS — adds `preUSAXStune` calls before the USAXS block when needed |
| USWAXS | Combined block — for each sample run USAXS / SAXS / WAXS together |

### Command file name

The default is `usaxs.mac`. Change it here or use **Export with custom
name…** to set it interactively.

### Export mode

- **Export current position set only** — only the table you can see now.
- **Export list of saved sets** — drag saved sets from the left list to
  the right "Export order" list to chain multiple sets in one file.

### Export buttons

| Button | What it does |
|---|---|
| **Export as usaxs.mac** | Write the command file to the last-used export folder. |
| Export with custom name… | Choose folder and filename. |
| Append to existing file… | Append the new commands to an existing `.mac` file. |

### Generated command-file format

```
        CURRENT_EXPERIMENT_NAME "MySamples"
        # USAXS measurements
      USAXSscan      45.070      98.300      1.000      "Sample1"
      USAXSscan      65.070      98.300      1.000      "Sample2"

        # SAXS measurements
      saxsExp        45.070      98.300      1.000      "Sample1"
      saxsExp        65.070      98.300      1.000      "Sample2"

        # WAXS measurements
      waxsExp        45.070      98.300      1.000      "Sample1"
      waxsExp        65.070      98.300      1.000      "Sample2"
        # END of batch of measurements
```

Header comments document the syntax. Run the file at the beamline with
`USAXS> CollectData usaxs.mac`.

### Blank-ratio warning

If you have more than one sample but fewer than 1 blank per 15 samples
in any technique, a warning appears in the status bar before export.
Use **Right-click → Set as Blank** to mark blank rows.

---

## Saving and loading

The tool persists its state to a single HDF5 file. Saved files include:

- All saved sample sets (named).
- The currently loaded set.
- Export order, command file name, export mode.
- Timing constants.
- Hook function rows and enable state.

| Action | Menu / button |
|---|---|
| Save the full tool state | **File → Save HDF5 As…** then **File → Save HDF5** to overwrite |
| Load a saved file | **File → Open HDF5…** |
| Save / load just the configuration | **Save Settings…** / **Load Settings…** in the toolbar |
| Reset everything to defaults | **Reset Tool** (red button in toolbar) |

The toolbar `Settings: …` label shows the active settings file path and
a `*` when there are unsaved changes.

---

## Beamline Survey (EPICS motor control)

Click **Beamline Survey** on the Sample Table tab to open the live
sample-stage control dialog. This requires `pyepics` and a connection
to the beamline EPICS gateway; otherwise the buttons are greyed out
and the dialog reports "EPICS not available".

Features:

- Live readback of the sample stage SX / SY motors (2 Hz polling).
- Drive to the SX / SY stored in the selected table row.
- Save the current motor positions back into the selected table row
  (handy when surveying a real plate by eye).
- **SX +/– / SY +/–** step buttons with a configurable step size.
- **Go to 0, 0** to send the stage home.
- **STOP MOTORS** (red) — emergency stop using the EPICS allstop PV.

PVs used: `usxAERO:m8` (SX), `usxAERO:m9` (SY), `usxLAX:allstop`,
plus aperture-slit PVs for the survey display.

---

## Importing legacy 12IDB command files

**File → Import 12IDB command file…** parses an existing `.mac` file
written by the legacy 12-ID-B beamline tools and populates the table
with its samples (name, SX, SY, thickness, technique flags). Use this
to migrate previous experiments to Matilda.

---

## Tips

- Click on the plate image to fill positions — typing them by hand
  invites typos.
- **Sample names** are sanitised on entry (digits, letters, underscore
  only). `%` becomes `pct`. If you paste an invalid name the tool
  silently cleans it.
- Save your work as an HDF5 file (`File → Save HDF5 As…`) before
  exporting — that file is the portable record of the experiment.
- Export the `.mac` file straight to the beamline's spec macros folder
  if you have access; otherwise scp it across.
- Use the **Hook Function** for systematic edge measurements — set
  ±1 mm in SX/SY for `_L`/`_R`/`_T`/`_B` and the tool generates a
  full edge map automatically.
- The **Beamline Survey** dialog is the safe way to visit each sample
  position in turn before running a long batch — verify nothing is
  going to crash.

---

## Where to next

- For interactive reduction of the data the beamline collects, see
  [matilda-gui.md](matilda-gui.md).
- For the headless reduction service that runs at the beamline, see
  [service.md](service.md).
- For the auto-analysis (pyirena) and auto-calibration (pynika) hooks,
  see [operations.md](operations.md).
