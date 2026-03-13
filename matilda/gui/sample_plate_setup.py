"""
matilda/gui/sample_plate_setup.py
==================================
GUI tool for setting up sample plate positions for USAXS/SAXS/WAXS measurements.

Replaces the Igor Pro "Setup Sample Plates" tool (IN3_SamplePlate.ipf).

The tool allows users to:
- Define sample positions (name, SX, SY, thickness, which detectors to use)
- Visualise positions on a plate image (generated from plate geometry)
- Click on the plate image to assign positions to table rows
- Save/load the sample set to/from HDF5 files for portability
- Export an ASCII command file (.mac) for the Bluesky CollectData command
- Optionally connect to EPICS PVs to drive the sample stage (Beamline Survey)

Usage::

    from matilda.gui.sample_plate_setup import run_sample_plate_setup
    run_sample_plate_setup()

Or as a standalone script::

    python -m matilda.gui.sample_plate_setup
"""

from __future__ import annotations

import os
import re
import math
import datetime
from dataclasses import dataclass, field
from typing import Optional

import numpy as np
import h5py

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QTabWidget,
    QTableWidget, QTableWidgetItem, QHeaderView,
    QVBoxLayout, QHBoxLayout, QGridLayout, QFormLayout,
    QLabel, QPushButton, QLineEdit, QDoubleSpinBox, QSpinBox,
    QCheckBox, QComboBox, QFileDialog, QMessageBox, QDialog,
    QDialogButtonBox, QGroupBox, QSplitter, QStatusBar,
    QAbstractItemView, QMenu, QListWidget, QListWidgetItem,
    QTextEdit, QSizePolicy, QScrollArea,
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal, QThread, QObject
from PyQt6.QtGui import QAction, QFont, QColor, QIcon

import pyqtgraph as pg

# ---------------------------------------------------------------------------
# Constants — scan times and overhead (seconds)
# ---------------------------------------------------------------------------

USAXS_SCAN_TIME_DEFAULT = 90   # seconds per USAXS flyscan
SAXS_SCAN_TIME_DEFAULT = 1     # seconds per SAXS exposure
WAXS_SCAN_TIME_DEFAULT = 3     # seconds per WAXS exposure

USAXS_OVERHEAD = 25            # seconds overhead per USAXS scan
SAXS_OVERHEAD = 10             # seconds overhead per SAXS scan
WAXS_OVERHEAD = 3              # seconds overhead per WAXS scan
SAMPLE_MOVE_SPEED = 8          # mm/second average sample stage speed
GEOMETRY_SWITCH_TIME = 20      # seconds to switch between USAXS/SAXS/WAXS geometry
USAXS_RETUNE_INTERVAL = 600    # seconds between SAXS/WAXS retunes
USAXS_RETUNE_EVERY_N = 3       # retune every N USAXS flyscans
USAXS_RETUNE_TIME = 40         # seconds per USAXS retune
SWAXS_RETUNE_TIME = 14         # seconds per SAXS/WAXS retune

DEFAULT_THICKNESS = 1.0        # mm
DEFAULT_CMD_FILENAME = "usaxs.mac"

BLANK_RATIO_WARN = 1 / 15      # warn if fewer than 1 blank per 15 samples

EXPORT_ORDERS = [
    "USAXS-SAXS-WAXS",
    "USAXS-WAXS-SAXS",
    "SAXS-WAXS-USAXS",
    "USWAXS",
]

# EPICS PV names (usxAERO / usxLAX prefixes)
PV_SX_RBV = "usxAERO:m8.RBV"
PV_SX_VAL = "usxAERO:m8.VAL"
PV_SX_SSET = "usxAERO:m8.SSET"
PV_SX_SUSE = "usxAERO:m8.SUSE"
PV_SX_VELO = "usxAERO:m8.VELO"
PV_SY_RBV = "usxAERO:m9.RBV"
PV_SY_VAL = "usxAERO:m9.VAL"
PV_SY_SSET = "usxAERO:m9.SSET"
PV_SY_SUSE = "usxAERO:m9.SUSE"
PV_SY_VELO = "usxAERO:m9.VELO"
PV_DATA_COLLECTING = "usxLAX:dataColInProgress"
PV_ALL_STOP = "usxLAX:allstop"
# Slits
PV_USAXS_HSLIT = "usxLAX:USAXS_hslit_ap"
PV_USAXS_VSLIT = "usxLAX:USAXS_vslit_ap"
PV_USAXS_HGSLIT = "usxLAX:USAXS_hgslit_ap"
PV_USAXS_VGSLIT = "usxLAX:USAXS_vgslit_ap"
PV_SAXS_HSLIT = "usxLAX:SAXS_hslit_ap"
PV_SAXS_VSLIT = "usxLAX:SAXS_vslit_ap"
PV_SAXS_HGSLIT = "usxLAX:SAXS_hgslit_ap"
PV_SAXS_VGSLIT = "usxLAX:SAXS_vgslit_ap"
PV_GSLIT1H = "usxLAX:GSlit1H:size"
PV_GSLIT1V = "usxLAX:GSlit1V:size"
PV_C1M8 = "usxLAX:m58:c1:m8.VAL"
PV_C1M7 = "usxLAX:m58:c1:m7.VAL"

# Try importing pyepics; feature-flag it
try:
    import epics
    EPICS_AVAILABLE = True
except ImportError:
    EPICS_AVAILABLE = False


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------

@dataclass
class SampleRow:
    """One row in the sample table."""
    name: str = ""
    sx: float = 0.0
    sy: float = 0.0
    thickness: float = DEFAULT_THICKNESS
    usaxs: bool = True
    saxs: bool = True
    waxs: bool = True
    metadata: str = ""


@dataclass
class SampleSet:
    """A named collection of sample rows plus global options."""
    name: str = "MySamples"
    rows: list = field(default_factory=list)
    usaxs_all: bool = True
    saxs_all: bool = True
    waxs_all: bool = True

    def __post_init__(self):
        if not self.rows:
            self.rows = [SampleRow() for _ in range(20)]

    def copy(self) -> "SampleSet":
        import copy
        return copy.deepcopy(self)


# ---------------------------------------------------------------------------
# Plate geometry definitions
# ---------------------------------------------------------------------------

def _make_9x9_acrylic_centers():
    """81 positions in a 9x9 grid, 20 mm spacing."""
    centers = []
    for i in range(81):
        col = i // 9 + 1       # 1..9
        row = i - 9 * (i // 9) # 0..8
        centers.append((20 * col, 20 + 20 * row))
    return centers


def _make_nmr_acrylic_centers():
    """100 positions in a 20x5 grid (NMR acrylic plate)."""
    centers = []
    for i in range(100):
        col = i - 20 * (i // 20)           # x: 0..19 -> 10..190 step 10
        row = i // 20                       # y: 0..4  -> 24.181 step 40
        centers.append((10 + 10 * col, 24.181 + 40 * row))
    return centers


def _make_oldstyle_al_centers():
    """60 positions - 4 columns, 15 per column (two sub-grids per column)."""
    centers = []
    for col in range(4):
        for j in range(8):
            centers.append((12.5 + col * 25, 12.5 + j * 25))
        for j in range(7):
            centers.append((25 + col * 25, 25 + j * 25))
    return centers


def _make_nmr_tubes_centers():
    """20 NMR tube positions in a row."""
    return [(15 + i * 10, 43) for i in range(20)]


PLATE_DEFINITIONS = {
    "9x9 Acrylic/magnetic plate": {
        "centers": _make_9x9_acrylic_centers(),
        "hole_radius": 4.0,   # mm
        "plate_w": 300.0,     # mm
        "plate_h": 200.0,
        "shape": "circle",
    },
    "NMR Acrylic plate": {
        "centers": _make_nmr_acrylic_centers(),
        "hole_radius": 2.0,
        "plate_w": 330.0,
        "plate_h": 200.0,
        "shape": "circle",
    },
    "Old Style Al Plate": {
        "centers": _make_oldstyle_al_centers(),
        "hole_radius": 4.0,
        "plate_w": 250.0,
        "plate_h": 200.0,
        "shape": "circle",
    },
    "NMR Tubes holder": {
        "centers": _make_nmr_tubes_centers(),
        "hole_radius": 4.0,
        "plate_w": 300.0,
        "plate_h": 75.0,
        "shape": "nmr_tube",
    },
    "Generic Grid": None,      # filled dynamically
    "AgBehenateLaB6": None,    # special: 1 position
}


def generate_plate_image(plate_name: str, generic_params: dict | None = None,
                          pix_per_mm: float = 2.0) -> np.ndarray | None:
    """
    Generate a greyscale NumPy array (uint8) showing the plate with hole positions.

    Returns None for plates without defined geometry.
    pix_per_mm: resolution - 2 pixels per mm gives a good display.
    """
    defn = PLATE_DEFINITIONS.get(plate_name)
    if defn is None:
        return None

    w = int(defn["plate_w"] * pix_per_mm)
    h = int(defn["plate_h"] * pix_per_mm)
    img = np.full((w, h), 128, dtype=np.uint8)  # grey background

    centers = defn["centers"]
    radius = defn["hole_radius"]

    if defn["shape"] == "circle":
        # Draw filled white circles for each hole
        xs = np.arange(w)
        ys = np.arange(h)
        XX, YY = np.meshgrid(xs, ys, indexing='ij')
        for cx, cy in centers:
            px = cx * pix_per_mm
            py = cy * pix_per_mm
            pr = radius * pix_per_mm
            mask = (XX - px) ** 2 + (YY - py) ** 2 < pr ** 2
            img[mask] = 255

    elif defn["shape"] == "nmr_tube":
        # Rectangular slot for NMR tubes
        img[int(10 * pix_per_mm):int(210 * pix_per_mm),
            int(20 * pix_per_mm):int(66 * pix_per_mm)] = 255
        for cx, cy in centers:
            px = cx * pix_per_mm
            r = 2 * pix_per_mm
            l_span = 23 * pix_per_mm
            x0 = max(0, int(px - r))
            x1 = min(w - 1, int(px + r))
            y0 = max(0, int(cy * pix_per_mm - l_span))
            y1 = min(h - 1, int(cy * pix_per_mm + l_span))
            img[x0:x1, y0:y1] = 64   # darker fill
            img[x0, y0:y1] = 255     # left wall
            img[x1, y0:y1] = 255     # right wall

    return img


# ---------------------------------------------------------------------------
# HDF5 persistence
# ---------------------------------------------------------------------------

def save_sets_to_hdf5(path: str, sets: dict[str, SampleSet]) -> None:
    """Save all named SampleSets to an HDF5 file."""
    with h5py.File(path, "w") as f:
        f.attrs["creator"] = "Matilda SamplePlateSetup"
        f.attrs["version"] = "1.0"
        f.attrs["saved"] = datetime.datetime.now().isoformat()
        grp = f.require_group("saved_sets")
        for name, ss in sets.items():
            sg = grp.require_group(name)
            sg.attrs["name"] = ss.name
            sg.attrs["usaxs_all"] = int(ss.usaxs_all)
            sg.attrs["saxs_all"] = int(ss.saxs_all)
            sg.attrs["waxs_all"] = int(ss.waxs_all)
            n = len(ss.rows)
            names_arr = np.array([r.name.encode() for r in ss.rows])
            sx_arr = np.array([r.sx for r in ss.rows])
            sy_arr = np.array([r.sy for r in ss.rows])
            th_arr = np.array([r.thickness for r in ss.rows])
            u_arr = np.array([int(r.usaxs) for r in ss.rows])
            s_arr = np.array([int(r.saxs) for r in ss.rows])
            w_arr = np.array([int(r.waxs) for r in ss.rows])
            md_arr = np.array([r.metadata.encode() for r in ss.rows])
            sg.create_dataset("name", data=names_arr)
            sg.create_dataset("sx", data=sx_arr)
            sg.create_dataset("sy", data=sy_arr)
            sg.create_dataset("thickness", data=th_arr)
            sg.create_dataset("usaxs", data=u_arr)
            sg.create_dataset("saxs", data=s_arr)
            sg.create_dataset("waxs", data=w_arr)
            sg.create_dataset("metadata", data=md_arr)


def load_sets_from_hdf5(path: str) -> dict[str, SampleSet]:
    """Load all saved SampleSets from an HDF5 file."""
    result = {}
    with h5py.File(path, "r") as f:
        grp = f.get("saved_sets")
        if grp is None:
            return result
        for name in grp:
            sg = grp[name]
            ss = SampleSet(name=sg.attrs.get("name", name))
            ss.usaxs_all = bool(sg.attrs.get("usaxs_all", 1))
            ss.saxs_all = bool(sg.attrs.get("saxs_all", 1))
            ss.waxs_all = bool(sg.attrs.get("waxs_all", 1))
            names_arr = sg["name"][:] if "name" in sg else []
            sx_arr = sg["sx"][:] if "sx" in sg else []
            sy_arr = sg["sy"][:] if "sy" in sg else []
            th_arr = sg["thickness"][:] if "thickness" in sg else []
            u_arr = sg["usaxs"][:] if "usaxs" in sg else []
            s_arr = sg["saxs"][:] if "saxs" in sg else []
            w_arr = sg["waxs"][:] if "waxs" in sg else []
            md_arr = sg["metadata"][:] if "metadata" in sg else []
            rows = []
            for i in range(len(names_arr)):
                row = SampleRow(
                    name=names_arr[i].decode() if hasattr(names_arr[i], "decode") else str(names_arr[i]),
                    sx=float(sx_arr[i]) if i < len(sx_arr) else 0.0,
                    sy=float(sy_arr[i]) if i < len(sy_arr) else 0.0,
                    thickness=float(th_arr[i]) if i < len(th_arr) else DEFAULT_THICKNESS,
                    usaxs=bool(u_arr[i]) if i < len(u_arr) else True,
                    saxs=bool(s_arr[i]) if i < len(s_arr) else True,
                    waxs=bool(w_arr[i]) if i < len(w_arr) else True,
                    metadata=md_arr[i].decode() if i < len(md_arr) and hasattr(md_arr[i], "decode") else "",
                )
                rows.append(row)
            ss.rows = rows
            result[name] = ss
    return result


# ---------------------------------------------------------------------------
# Command file generation
# ---------------------------------------------------------------------------

def check_for_blanks(sample_set: SampleSet) -> str:
    """Return a warning string if there are not enough blank measurements."""
    warnings = []
    usaxs_names = [r.name for r in sample_set.rows
                   if r.name and (sample_set.usaxs_all or r.usaxs)]
    saxs_names = [r.name for r in sample_set.rows
                  if r.name and (sample_set.saxs_all or r.saxs)]
    waxs_names = [r.name for r in sample_set.rows
                  if r.name and (sample_set.waxs_all or r.waxs)]

    blank_re = re.compile(r"(blank|empty)", re.IGNORECASE)

    def blank_ratio(names):
        if not names:
            return 1.0
        blanks = sum(1 for n in names if blank_re.search(n))
        return blanks / len(names)

    if len(usaxs_names) > 1 and blank_ratio(usaxs_names) < BLANK_RATIO_WARN:
        warnings.append("USAXS")
    if len(saxs_names) > 1 and blank_ratio(saxs_names) < BLANK_RATIO_WARN:
        warnings.append("SAXS")
    if len(waxs_names) > 1 and blank_ratio(waxs_names) < BLANK_RATIO_WARN:
        warnings.append("WAXS")

    if warnings:
        return f"Warning: not enough blanks for {', '.join(warnings)} (need >=1 per 15 samples)"
    return ""


def generate_command_file(sets_to_export: list[SampleSet],
                          export_order: str = "USAXS-SAXS-WAXS",
                          include_header: bool = True) -> str:
    """Generate the ASCII command file content for one or more SampleSets."""
    lines = []
    if include_header:
        exp_name = sets_to_export[0].name if sets_to_export else "experiment"
        lines += [
            f'        CURRENT_EXPERIMENT_NAME "{exp_name}"',
            "        # This file runs USAXS, SAXS and WAXS scans according to the syntax shown below",
            "        #",
            "        # Scan Type      sx         sy   Thickness  Sample Name",
            "        # ------------------------------------------------------",
            '        # USAXSscan    45.07       98.3     0      "Water Blank"',
            '        # saxsExp      45.07       98.3     0      "Water Blank"',
            '        # waxsExp      45.07       98.3     0      "Water Blank"',
            "        #      Use a space (not a tab) to separate arguments",
            "",
            "        # Run this file by typing in the spec window:   USAXS> CollectData usaxs.mac",
            "",
            "        # Stop using the 'Stop after this scan?' checkbox in USAXS user main intf",
            "",
            "        ############ PLACE ALL USER COMMANDS AFTER THIS LINE ############",
            "",
        ]

    for ss in sets_to_export:
        lines.append("")
        lines += _write_commands_for_set(ss, export_order)

    return "\n".join(lines)


def _write_commands_for_set(ss: SampleSet, export_order: str) -> list[str]:
    lines = []

    def usaxs_rows():
        return [r for r in ss.rows
                if r.name and r.sx != "" and r.sy != ""
                and (ss.usaxs_all or r.usaxs)]

    def saxs_rows():
        return [r for r in ss.rows
                if r.name and r.sx != "" and r.sy != ""
                and (ss.saxs_all or r.saxs)]

    def waxs_rows():
        return [r for r in ss.rows
                if r.name and r.sx != "" and r.sy != ""
                and (ss.waxs_all or r.waxs)]

    def fmt_usaxs(r: SampleRow):
        t = r.thickness if r.thickness > 0 else DEFAULT_THICKNESS
        return f'      USAXSscan      {r.sx:.3f}      {r.sy:.3f}      {t:.3f}      "{r.name}"'

    def fmt_saxs(r: SampleRow):
        t = r.thickness if r.thickness > 0 else DEFAULT_THICKNESS
        return f'      saxsExp        {r.sx:.3f}      {r.sy:.3f}      {t:.3f}      "{r.name}"'

    def fmt_waxs(r: SampleRow):
        t = r.thickness if r.thickness > 0 else DEFAULT_THICKNESS
        return f'      waxsExp        {r.sx:.3f}      {r.sy:.3f}      {t:.3f}      "{r.name}"'

    if export_order == "USAXS-SAXS-WAXS":
        lines.append("        #USAXS measurements")
        lines += [fmt_usaxs(r) for r in usaxs_rows()]
        lines.append("")
        lines.append("        #SAXS measurements")
        lines += [fmt_saxs(r) for r in saxs_rows()]
        lines.append("")
        lines.append("        #WAXS measurements")
        lines += [fmt_waxs(r) for r in waxs_rows()]
        lines.append("        #END of batch of measurements")

    elif export_order == "USAXS-WAXS-SAXS":
        lines.append("        #USAXS measurements")
        lines += [fmt_usaxs(r) for r in usaxs_rows()]
        lines.append("")
        lines.append("        #WAXS measurements")
        lines += [fmt_waxs(r) for r in waxs_rows()]
        lines.append("")
        lines.append("        #SAXS measurements")
        lines += [fmt_saxs(r) for r in saxs_rows()]
        lines.append("        #END of batch of measurements")

    elif export_order == "SAXS-WAXS-USAXS":
        lines.append("        #SAXS measurements")
        lines += [fmt_saxs(r) for r in saxs_rows()]
        lines.append("")
        lines.append("        #WAXS measurements")
        lines += [fmt_waxs(r) for r in waxs_rows()]
        # Add preUSAXStune lines if needed
        num_sw = len(saxs_rows()) + len(waxs_rows())
        ur = usaxs_rows()
        if ur:
            lines.append("")
            lines.append("        preUSAXStune")
            if num_sw > 15:
                lines.append("        preUSAXStune")
            if num_sw > 30:
                lines.append("        preUSAXStune")
        lines.append("")
        lines.append("        #USAXS measurements")
        lines += [fmt_usaxs(r) for r in ur]
        lines.append("        #END USAXS measurements")

    elif export_order == "USWAXS":
        lines.append("        #Combined USAXS - SAXS - WAXS measurements")
        lines.append("")
        for r in usaxs_rows():
            t = r.thickness if r.thickness > 0 else DEFAULT_THICKNESS
            lines.append(f'      USAXSscan      {r.sx:.3f}      {r.sy:.3f}      {t:.3f}      "{r.name}"')
            lines.append(f'      saxsExp        {r.sx:.3f}      {r.sy:.3f}      {t:.3f}      "{r.name}"')
            lines.append(f'      waxsExp        {r.sx:.3f}      {r.sy:.3f}      {t:.3f}      "{r.name}"')
            lines.append("")
        lines.append("        #END of batch of measurements")

    return lines


# ---------------------------------------------------------------------------
# Run-time estimation
# ---------------------------------------------------------------------------

def estimate_run_time(ss: SampleSet, usaxs_time: float, saxs_time: float,
                      waxs_time: float) -> tuple[int, int, int, int]:
    """
    Returns (num_usaxs, num_saxs, num_waxs, total_minutes).
    """
    usaxs_rows = [r for r in ss.rows
                  if r.name and (ss.usaxs_all or r.usaxs)]
    saxs_rows = [r for r in ss.rows
                 if r.name and (ss.saxs_all or r.saxs)]
    waxs_rows = [r for r in ss.rows
                 if r.name and (ss.waxs_all or r.waxs)]

    def total_travel(rows):
        dist = 0.0
        for i in range(1, len(rows)):
            dx = rows[i].sx - rows[i - 1].sx
            dy = rows[i].sy - rows[i - 1].sy
            dist += math.sqrt(dx * dx + dy * dy)
        return dist

    total = 0.0
    n_usaxs = len(usaxs_rows)
    n_saxs = len(saxs_rows)
    n_waxs = len(waxs_rows)

    total += n_usaxs * (usaxs_time + USAXS_OVERHEAD)
    total += total_travel(usaxs_rows) / SAMPLE_MOVE_SPEED + n_usaxs

    if n_usaxs > 0:
        total += GEOMETRY_SWITCH_TIME

    total += n_saxs * (saxs_time + SAXS_OVERHEAD)
    total += total_travel(saxs_rows) / SAMPLE_MOVE_SPEED + n_saxs
    sw_time = n_saxs * (saxs_time + SAXS_OVERHEAD) + total_travel(saxs_rows) / SAMPLE_MOVE_SPEED + n_saxs

    if n_saxs > 0:
        total += GEOMETRY_SWITCH_TIME

    total += n_waxs * (waxs_time + WAXS_OVERHEAD)
    total += total_travel(waxs_rows) / SAMPLE_MOVE_SPEED + n_waxs
    sw_time += n_waxs * (waxs_time + WAXS_OVERHEAD) + total_travel(waxs_rows) / SAMPLE_MOVE_SPEED + n_waxs

    if n_waxs > 0:
        total += GEOMETRY_SWITCH_TIME

    num_sw_tunes = sw_time / USAXS_RETUNE_INTERVAL
    num_u_tunes = n_usaxs / USAXS_RETUNE_EVERY_N
    total += num_sw_tunes * SWAXS_RETUNE_TIME + num_u_tunes * USAXS_RETUNE_TIME

    total_min = max(0, round(total / 60))
    return n_usaxs, n_saxs, n_waxs, total_min


# ---------------------------------------------------------------------------
# Sample table widget
# ---------------------------------------------------------------------------

COL_NAME = 0
COL_SX = 1
COL_SY = 2
COL_THICK = 3
COL_USAXS = 4
COL_SAXS = 5
COL_WAXS = 6
COL_META = 7
COLUMNS = ["Sample Name", "X [mm]", "Y [mm]", "Thick [mm]", "USAXS", "SAXS", "WAXS", "MetaData"]


class SampleTable(QTableWidget):
    """
    Table widget for sample positions.  Columns: Name, SX, SY, Thickness,
    USAXS-checkbox, SAXS-checkbox, WAXS-checkbox, MetaData.
    """
    dataChanged = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__(0, 8, parent)
        self.setHorizontalHeaderLabels(COLUMNS)
        hdr = self.horizontalHeader()
        hdr.setSectionResizeMode(COL_NAME, QHeaderView.ResizeMode.Stretch)
        for col in (COL_SX, COL_SY, COL_THICK):
            hdr.setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents)
        for col in (COL_USAXS, COL_SAXS, COL_WAXS):
            hdr.setSectionResizeMode(col, QHeaderView.ResizeMode.Fixed)
            self.setColumnWidth(col, 50)
        hdr.setSectionResizeMode(COL_META, QHeaderView.ResizeMode.Stretch)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        self.cellChanged.connect(self._on_cell_changed)
        self._block_signals = False

    def _on_cell_changed(self, row, col):
        if not self._block_signals:
            self.dataChanged.emit()

    def load_from_sample_set(self, ss: SampleSet):
        self._block_signals = True
        self.setRowCount(0)
        for row in ss.rows:
            self._append_row(row)
        self._block_signals = False

    def _append_row(self, row: SampleRow):
        r = self.rowCount()
        self.insertRow(r)
        self.setItem(r, COL_NAME, QTableWidgetItem(row.name))
        self.setItem(r, COL_SX, QTableWidgetItem(f"{row.sx:.3f}" if row.sx else ""))
        self.setItem(r, COL_SY, QTableWidgetItem(f"{row.sy:.3f}" if row.sy else ""))
        self.setItem(r, COL_THICK, QTableWidgetItem(f"{row.thickness:.3f}"))
        for col, val in ((COL_USAXS, row.usaxs), (COL_SAXS, row.saxs), (COL_WAXS, row.waxs)):
            chk = QCheckBox()
            chk.setChecked(val)
            chk.setStyleSheet("margin-left:12px;")
            chk.stateChanged.connect(self.dataChanged)
            self.setCellWidget(r, col, chk)
        self.setItem(r, COL_META, QTableWidgetItem(row.metadata))

    def get_sample_set(self) -> list[SampleRow]:
        rows = []
        for r in range(self.rowCount()):
            def txt(col):
                item = self.item(r, col)
                return item.text().strip() if item else ""

            def num(col, default=0.0):
                try:
                    return float(txt(col))
                except (ValueError, TypeError):
                    return default

            def chk(col):
                w = self.cellWidget(r, col)
                return w.isChecked() if w else True

            rows.append(SampleRow(
                name=txt(COL_NAME),
                sx=num(COL_SX),
                sy=num(COL_SY),
                thickness=num(COL_THICK, DEFAULT_THICKNESS),
                usaxs=chk(COL_USAXS),
                saxs=chk(COL_SAXS),
                waxs=chk(COL_WAXS),
                metadata=txt(COL_META),
            ))
        return rows

    def _current_row(self):
        rows = self.selectedItems()
        if rows:
            return rows[0].row()
        return self.rowCount() - 1

    def insert_row_above(self):
        r = self._current_row()
        self.insertRow(r)
        row = SampleRow()
        self.setItem(r, COL_NAME, QTableWidgetItem(""))
        self.setItem(r, COL_SX, QTableWidgetItem(""))
        self.setItem(r, COL_SY, QTableWidgetItem(""))
        self.setItem(r, COL_THICK, QTableWidgetItem(f"{DEFAULT_THICKNESS:.3f}"))
        for col, val in ((COL_USAXS, True), (COL_SAXS, True), (COL_WAXS, True)):
            chk = QCheckBox()
            chk.setChecked(val)
            chk.setStyleSheet("margin-left:12px;")
            chk.stateChanged.connect(self.dataChanged)
            self.setCellWidget(r, col, chk)
        self.setItem(r, COL_META, QTableWidgetItem(""))
        self.dataChanged.emit()

    def insert_row_below(self):
        r = self._current_row() + 1
        self.insertRow(r)
        self.setItem(r, COL_NAME, QTableWidgetItem(""))
        self.setItem(r, COL_SX, QTableWidgetItem(""))
        self.setItem(r, COL_SY, QTableWidgetItem(""))
        self.setItem(r, COL_THICK, QTableWidgetItem(f"{DEFAULT_THICKNESS:.3f}"))
        for col, val in ((COL_USAXS, True), (COL_SAXS, True), (COL_WAXS, True)):
            chk = QCheckBox()
            chk.setChecked(val)
            chk.setStyleSheet("margin-left:12px;")
            chk.stateChanged.connect(self.dataChanged)
            self.setCellWidget(r, col, chk)
        self.setItem(r, COL_META, QTableWidgetItem(""))
        self.dataChanged.emit()

    def delete_row(self):
        r = self._current_row()
        if self.rowCount() > 1:
            self.removeRow(r)
            self.dataChanged.emit()

    def duplicate_row(self):
        r = self._current_row()
        rows = self.get_sample_set()
        if 0 <= r < len(rows):
            self.insertRow(r + 1)
            row = rows[r]
            ri = r + 1
            self.setItem(ri, COL_NAME, QTableWidgetItem(row.name))
            self.setItem(ri, COL_SX, QTableWidgetItem(f"{row.sx:.3f}"))
            self.setItem(ri, COL_SY, QTableWidgetItem(f"{row.sy:.3f}"))
            self.setItem(ri, COL_THICK, QTableWidgetItem(f"{row.thickness:.3f}"))
            for col, val in ((COL_USAXS, row.usaxs), (COL_SAXS, row.saxs), (COL_WAXS, row.waxs)):
                chk = QCheckBox()
                chk.setChecked(val)
                chk.setStyleSheet("margin-left:12px;")
                chk.stateChanged.connect(self.dataChanged)
                self.setCellWidget(ri, col, chk)
            self.setItem(ri, COL_META, QTableWidgetItem(row.metadata))
            self.dataChanged.emit()

    def move_row_up(self):
        r = self._current_row()
        if r > 0:
            rows = self.get_sample_set()
            rows[r - 1], rows[r] = rows[r], rows[r - 1]
            self._reload_rows(rows)
            self.selectRow(r - 1)

    def move_row_down(self):
        r = self._current_row()
        rows = self.get_sample_set()
        if r < len(rows) - 1:
            rows[r], rows[r + 1] = rows[r + 1], rows[r]
            self._reload_rows(rows)
            self.selectRow(r + 1)

    def _reload_rows(self, rows: list[SampleRow]):
        self._block_signals = True
        self.setRowCount(0)
        for row in rows:
            self._append_row(row)
        self._block_signals = False
        self.dataChanged.emit()

    def increment_sx(self, delta=1.0):
        for r in set(item.row() for item in self.selectedItems()):
            item = self.item(r, COL_SX)
            try:
                val = float(item.text()) + delta
            except (ValueError, AttributeError):
                val = delta
            self.setItem(r, COL_SX, QTableWidgetItem(f"{val:.3f}"))
        self.dataChanged.emit()

    def increment_sy(self, delta=1.0):
        for r in set(item.row() for item in self.selectedItems()):
            item = self.item(r, COL_SY)
            try:
                val = float(item.text()) + delta
            except (ValueError, AttributeError):
                val = delta
            self.setItem(r, COL_SY, QTableWidgetItem(f"{val:.3f}"))
        self.dataChanged.emit()

    def add_to_sx(self):
        val, ok = _ask_float(self, "Add to SX", "Value to add to SX [mm]:", 1.0)
        if ok:
            self.increment_sx(val)

    def add_to_sy(self):
        val, ok = _ask_float(self, "Add to SY", "Value to add to SY [mm]:", 1.0)
        if ok:
            self.increment_sy(val)

    def set_as_blank(self):
        for r in set(item.row() for item in self.selectedItems()):
            item = self.item(r, COL_NAME)
            if item:
                name = item.text()
                if "Blank" not in name and "blank" not in name:
                    item.setText(name + " Blank" if name else "Blank")
        self.dataChanged.emit()

    def set_as_agbehenaatelab6(self):
        for r in set(item.row() for item in self.selectedItems()):
            self.setItem(r, COL_NAME, QTableWidgetItem("AgBehenateLaB6"))
            self.setItem(r, COL_SX, QTableWidgetItem("20.000"))
            self.setItem(r, COL_SY, QTableWidgetItem("20.000"))
            self.setItem(r, COL_THICK, QTableWidgetItem("1.000"))
            chk = self.cellWidget(r, COL_USAXS)
            if chk:
                chk.setChecked(False)
        self.dataChanged.emit()

    def clear_row(self):
        r = self._current_row()
        for col in (COL_NAME, COL_SX, COL_SY, COL_THICK, COL_META):
            self.setItem(r, col, QTableWidgetItem(""))
        for col in (COL_USAXS, COL_SAXS, COL_WAXS):
            chk = self.cellWidget(r, col)
            if chk:
                chk.setChecked(True)
        self.dataChanged.emit()

    def copy_selection(self):
        rows = sorted(set(item.row() for item in self.selectedItems()))
        lines = []
        for r in rows:
            row_data = []
            for col in range(self.columnCount()):
                chkw = self.cellWidget(r, col)
                if chkw:
                    row_data.append("1" if chkw.isChecked() else "0")
                else:
                    item = self.item(r, col)
                    row_data.append(item.text() if item else "")
            lines.append("\t".join(row_data))
        QApplication.clipboard().setText("\n".join(lines))

    def paste_selection(self):
        text = QApplication.clipboard().text()
        if not text:
            return
        r = self._current_row()
        for line in text.split("\n"):
            parts = line.split("\t")
            if r >= self.rowCount():
                self._append_row(SampleRow())
            for col, val in enumerate(parts):
                if col >= self.columnCount():
                    break
                chkw = self.cellWidget(r, col)
                if chkw:
                    chkw.setChecked(val == "1")
                else:
                    self.setItem(r, col, QTableWidgetItem(val))
            r += 1
        self.dataChanged.emit()

    def _show_context_menu(self, pos):
        r = self.rowAt(pos.y())
        if r < 0:
            return
        menu = QMenu(self)
        menu.addAction("Insert Row Above", self.insert_row_above)
        menu.addAction("Insert Row Below", self.insert_row_below)
        menu.addAction("Delete Row", self.delete_row)
        menu.addAction("Duplicate Row", self.duplicate_row)
        menu.addSeparator()
        menu.addAction("Move Row Up", self.move_row_up)
        menu.addAction("Move Row Down", self.move_row_down)
        menu.addSeparator()
        menu.addAction("Cut", lambda: (self.copy_selection(), self.clear_row()))
        menu.addAction("Copy", self.copy_selection)
        menu.addAction("Paste", self.paste_selection)
        menu.addSeparator()
        menu.addAction("Increment SX (+1 mm)", lambda: self.increment_sx(1.0))
        menu.addAction("Increment SY (+1 mm)", lambda: self.increment_sy(1.0))
        menu.addAction("Add value to SX...", self.add_to_sx)
        menu.addAction("Add value to SY...", self.add_to_sy)
        menu.addSeparator()
        menu.addAction("Set as Blank", self.set_as_blank)
        menu.addAction("Set as AgBehenateLaB6", self.set_as_agbehenaatelab6)
        menu.addAction("Clear Row", self.clear_row)
        menu.exec(self.viewport().mapToGlobal(pos))

    def set_position_for_current_row(self, sx: float, sy: float):
        r = self._current_row()
        if 0 <= r < self.rowCount():
            self.setItem(r, COL_SX, QTableWidgetItem(f"{sx:.3f}"))
            self.setItem(r, COL_SY, QTableWidgetItem(f"{sy:.3f}"))
            self.dataChanged.emit()


def _ask_float(parent, title, label, default=0.0):
    """Simple dialog to ask for a float value."""
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    layout = QVBoxLayout(dlg)
    layout.addWidget(QLabel(label))
    spin = QDoubleSpinBox()
    spin.setRange(-9999, 9999)
    spin.setValue(default)
    spin.setDecimals(3)
    layout.addWidget(spin)
    buttons = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
    buttons.accepted.connect(dlg.accept)
    buttons.rejected.connect(dlg.reject)
    layout.addWidget(buttons)
    ok = dlg.exec() == QDialog.DialogCode.Accepted
    return spin.value(), ok


# ---------------------------------------------------------------------------
# Plate image canvas
# ---------------------------------------------------------------------------

class PlateCanvas(pg.GraphicsLayoutWidget):
    """
    pyqtgraph widget showing the plate image and sample position markers.
    Click on the image to assign position to the currently selected table row.
    """
    positionClicked = pyqtSignal(float, float)  # sx, sy in mm

    def __init__(self, parent=None):
        super().__init__(parent)
        self._plot = self.addPlot(row=0, col=0)
        self._plot.setAspectLocked(True)
        self._plot.invertY(True)
        self._img_item = pg.ImageItem()
        self._plot.addItem(self._img_item)
        self._scatter = pg.ScatterPlotItem(
            size=10, pen=pg.mkPen("r", width=2), brush=pg.mkBrush(None))
        self._plot.addItem(self._scatter)
        self._text_items: list[pg.TextItem] = []
        self._plate_w = 200.0
        self._plate_h = 200.0
        self._pix_per_mm = 2.0

        self._plot.scene().sigMouseClicked.connect(self._on_click)

    def set_plate(self, plate_name: str):
        img = generate_plate_image(plate_name, pix_per_mm=self._pix_per_mm)
        if img is None:
            self._img_item.clear()
            return
        defn = PLATE_DEFINITIONS.get(plate_name)
        if defn:
            self._plate_w = defn["plate_w"]
            self._plate_h = defn["plate_h"]
        # ImageItem: pixel axis -> mm
        self._img_item.setImage(img, autoLevels=True)
        self._img_item.setRect(0, 0, self._plate_w, self._plate_h)
        self._plot.setXRange(0, self._plate_w)
        self._plot.setYRange(0, self._plate_h)

    def update_markers(self, rows: list[SampleRow], current_row: int = -1):
        for t in self._text_items:
            self._plot.removeItem(t)
        self._text_items.clear()

        spots = []
        for i, r in enumerate(rows):
            if r.name and r.sx is not None and r.sy is not None:
                color = "r" if i == current_row else "y"
                spots.append({
                    "pos": (r.sx, r.sy),
                    "brush": pg.mkBrush(color),
                    "size": 10 if i == current_row else 7,
                })
                txt = pg.TextItem(r.name, anchor=(0, 1), color=color)
                txt.setPos(r.sx, r.sy)
                self._plot.addItem(txt)
                self._text_items.append(txt)
        self._scatter.setData(spots)

    def _on_click(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        pos = event.scenePos()
        view_pos = self._plot.vb.mapSceneToView(pos)
        sx = view_pos.x()
        sy = view_pos.y()
        if 0 <= sx <= self._plate_w and 0 <= sy <= self._plate_h:
            self.positionClicked.emit(sx, sy)


# ---------------------------------------------------------------------------
# Generic grid dialog
# ---------------------------------------------------------------------------

class GenericGridDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Generic Grid Parameters")
        layout = QFormLayout(self)
        self.cols = QSpinBox(); self.cols.setRange(1, 50); self.cols.setValue(5)
        self.rows = QSpinBox(); self.rows.setRange(1, 50); self.rows.setValue(5)
        self.x_start = QDoubleSpinBox(); self.x_start.setRange(0, 500); self.x_start.setValue(10)
        self.y_start = QDoubleSpinBox(); self.y_start.setRange(0, 500); self.y_start.setValue(10)
        self.x_step = QDoubleSpinBox(); self.x_step.setRange(1, 100); self.x_step.setValue(20)
        self.y_step = QDoubleSpinBox(); self.y_step.setRange(1, 100); self.y_step.setValue(20)
        layout.addRow("Columns:", self.cols)
        layout.addRow("Rows:", self.rows)
        layout.addRow("X start [mm]:", self.x_start)
        layout.addRow("Y start [mm]:", self.y_start)
        layout.addRow("X step [mm]:", self.x_step)
        layout.addRow("Y step [mm]:", self.y_step)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addRow(buttons)

    def get_centers(self):
        cs = []
        for col in range(self.cols.value()):
            for row in range(self.rows.value()):
                cs.append((
                    self.x_start.value() + col * self.x_step.value(),
                    self.y_start.value() + row * self.y_step.value(),
                ))
        return cs

    def total(self):
        return self.cols.value() * self.rows.value()


# ---------------------------------------------------------------------------
# Multi-set export widget (drag-from-saved, drop-to-export)
# ---------------------------------------------------------------------------

class MultiSetExportWidget(QWidget):
    """Shows saved sets on the left; user drags them to the export list on the right."""

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QHBoxLayout(self)

        lgroup = QGroupBox("Saved sets (drag to export list)")
        ll = QVBoxLayout(lgroup)
        self.source_list = QListWidget()
        self.source_list.setDragEnabled(True)
        self.source_list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        ll.addWidget(self.source_list)
        layout.addWidget(lgroup)

        rgroup = QGroupBox("Export order (drag to reorder, right-click to remove)")
        rl = QVBoxLayout(rgroup)
        self.export_list = QListWidget()
        self.export_list.setAcceptDrops(True)
        self.export_list.setDragEnabled(True)
        self.export_list.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)
        self.export_list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self.export_list.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.export_list.customContextMenuRequested.connect(self._remove_item)
        rl.addWidget(self.export_list)
        layout.addWidget(rgroup)

    def _remove_item(self, pos):
        item = self.export_list.itemAt(pos)
        if item:
            row = self.export_list.row(item)
            self.export_list.takeItem(row)

    def set_saved_sets(self, names: list[str]):
        self.source_list.clear()
        for name in names:
            self.source_list.addItem(QListWidgetItem(name))

    def get_export_order(self) -> list[str]:
        return [self.export_list.item(i).text()
                for i in range(self.export_list.count())]

    def add_selected_to_export(self):
        for item in self.source_list.selectedItems():
            if item.text() not in self.get_export_order():
                self.export_list.addItem(QListWidgetItem(item.text()))


# ---------------------------------------------------------------------------
# Beamline Survey Dialog (EPICS-dependent)
# ---------------------------------------------------------------------------

class BeamlineSurveyDialog(QDialog):
    """
    Modal-less dialog for driving sample stage motors at the beamline.
    Requires pyepics (EPICS_AVAILABLE must be True).
    """
    def __init__(self, sample_set: SampleSet, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Beamline Sample Plate Survey")
        self.setWindowFlags(self.windowFlags() & ~Qt.WindowType.WindowContextHelpButtonHint)
        self._ss = sample_set
        self._epics_timer = QTimer(self)
        self._epics_timer.setInterval(500)  # 2 Hz
        self._epics_timer.timeout.connect(self._poll_epics)

        self._build_ui()

        if EPICS_AVAILABLE:
            self._epics_timer.start()
        else:
            self._status_lbl.setText("EPICS not available - motor control disabled")
            for btn in self._motor_buttons:
                btn.setEnabled(False)

    def _build_ui(self):
        layout = QVBoxLayout(self)

        # Row selection
        top = QHBoxLayout()
        top.addWidget(QLabel("Selected row:"))
        self._row_spin = QSpinBox()
        self._row_spin.setRange(0, max(0, len(self._ss.rows) - 1))
        self._row_spin.valueChanged.connect(self._on_row_changed)
        top.addWidget(self._row_spin)
        top.addWidget(QLabel("Sample:"))
        self._sample_lbl = QLabel("")
        self._sample_lbl.setFont(QFont("Arial", 12, QFont.Weight.Bold))
        top.addWidget(self._sample_lbl)
        layout.addLayout(top)

        # Navigation
        nav = QHBoxLayout()
        btn_up = QPushButton("Row up")
        btn_up.clicked.connect(self._row_up)
        btn_dn = QPushButton("Row down")
        btn_dn.clicked.connect(self._row_down)
        nav.addWidget(btn_up)
        nav.addWidget(btn_dn)
        layout.addLayout(nav)

        # Position display
        pos_grp = QGroupBox("Motor positions")
        pos_lay = QGridLayout(pos_grp)
        pos_lay.addWidget(QLabel("RBV SX:"), 0, 0)
        self._sx_rbv = QLabel("---")
        self._sx_rbv.setFont(QFont("Courier", 14, QFont.Weight.Bold))
        pos_lay.addWidget(self._sx_rbv, 0, 1)
        pos_lay.addWidget(QLabel("RBV SY:"), 0, 2)
        self._sy_rbv = QLabel("---")
        self._sy_rbv.setFont(QFont("Courier", 14, QFont.Weight.Bold))
        pos_lay.addWidget(self._sy_rbv, 0, 3)

        pos_lay.addWidget(QLabel("Table SX:"), 1, 0)
        self._sx_tbl = QLabel("---")
        pos_lay.addWidget(self._sx_tbl, 1, 1)
        pos_lay.addWidget(QLabel("Table SY:"), 1, 2)
        self._sy_tbl = QLabel("---")
        pos_lay.addWidget(self._sy_tbl, 1, 3)
        layout.addWidget(pos_grp)

        # Target + move buttons
        tgt_grp = QGroupBox("Target position [mm]")
        tgt_lay = QGridLayout(tgt_grp)
        tgt_lay.addWidget(QLabel("SX target:"), 0, 0)
        self._sx_tar = QDoubleSpinBox()
        self._sx_tar.setRange(-300, 300)
        self._sx_tar.setDecimals(3)
        tgt_lay.addWidget(self._sx_tar, 0, 1)
        tgt_lay.addWidget(QLabel("SY target:"), 0, 2)
        self._sy_tar = QDoubleSpinBox()
        self._sy_tar.setRange(-300, 300)
        self._sy_tar.setDecimals(3)
        tgt_lay.addWidget(self._sy_tar, 0, 3)

        btn_sx_dn = QPushButton("SX -")
        btn_sx_up = QPushButton("SX +")
        btn_sy_dn = QPushButton("SY -")
        btn_sy_up = QPushButton("SY +")
        self._step_sx = QDoubleSpinBox()
        self._step_sx.setRange(0.01, 100)
        self._step_sx.setValue(1.0)
        self._step_sx.setDecimals(2)
        self._step_sy = QDoubleSpinBox()
        self._step_sy.setRange(0.01, 100)
        self._step_sy.setValue(1.0)
        self._step_sy.setDecimals(2)

        tgt_lay.addWidget(QLabel("SX step:"), 1, 0)
        tgt_lay.addWidget(self._step_sx, 1, 1)
        tgt_lay.addWidget(btn_sx_dn, 2, 0)
        tgt_lay.addWidget(btn_sx_up, 2, 1)
        tgt_lay.addWidget(QLabel("SY step:"), 1, 2)
        tgt_lay.addWidget(self._step_sy, 1, 3)
        tgt_lay.addWidget(btn_sy_dn, 2, 2)
        tgt_lay.addWidget(btn_sy_up, 2, 3)
        layout.addWidget(tgt_grp)

        # Action buttons
        act_lay = QHBoxLayout()
        btn_drive = QPushButton("Drive to table values")
        btn_drive.setStyleSheet("background-color: #88ff88")
        btn_save = QPushButton("Save current position to table")
        btn_save.setStyleSheet("background-color: #ff8888")
        btn_go00 = QPushButton("Go to 0, 0")
        btn_stop = QPushButton("STOP MOTORS")
        btn_stop.setStyleSheet("background-color: red; color: white; font-weight: bold;")
        act_lay.addWidget(btn_drive)
        act_lay.addWidget(btn_save)
        act_lay.addWidget(btn_go00)
        act_lay.addWidget(btn_stop)
        layout.addLayout(act_lay)

        # Slit presets
        slit_lay = QHBoxLayout()
        btn_slits_lg = QPushButton("Open Slits Large")
        btn_slits_us = QPushButton("USAXS Slits")
        btn_slits_sw = QPushButton("SAXS/WAXS Slits")
        slit_lay.addWidget(btn_slits_lg)
        slit_lay.addWidget(btn_slits_us)
        slit_lay.addWidget(btn_slits_sw)
        layout.addLayout(slit_lay)

        # Status
        self._status_lbl = QLabel("Ready")
        layout.addWidget(self._status_lbl)

        self._motor_buttons = [btn_sx_dn, btn_sx_up, btn_sy_dn, btn_sy_up,
                                btn_drive, btn_save, btn_go00, btn_stop,
                                btn_slits_lg, btn_slits_us, btn_slits_sw]

        btn_sx_dn.clicked.connect(lambda: self._jog("SX", -self._step_sx.value()))
        btn_sx_up.clicked.connect(lambda: self._jog("SX", +self._step_sx.value()))
        btn_sy_dn.clicked.connect(lambda: self._jog("SY", -self._step_sy.value()))
        btn_sy_up.clicked.connect(lambda: self._jog("SY", +self._step_sy.value()))
        btn_drive.clicked.connect(self._drive_to_table)
        btn_save.clicked.connect(self._save_position)
        btn_go00.clicked.connect(self._go_to_zero)
        btn_stop.clicked.connect(self._stop_motors)
        btn_slits_lg.clicked.connect(lambda: self._set_slits("large"))
        btn_slits_us.clicked.connect(lambda: self._set_slits("usaxs"))
        btn_slits_sw.clicked.connect(lambda: self._set_slits("saxswaxs"))

        self._on_row_changed(0)

    def closeEvent(self, event):
        self._epics_timer.stop()
        super().closeEvent(event)

    def _on_row_changed(self, row):
        if 0 <= row < len(self._ss.rows):
            r = self._ss.rows[row]
            self._sample_lbl.setText(r.name)
            self._sx_tbl.setText(f"{r.sx:.3f}")
            self._sy_tbl.setText(f"{r.sy:.3f}")
            self._sx_tar.setValue(r.sx)
            self._sy_tar.setValue(r.sy)

    def _row_up(self):
        v = self._row_spin.value()
        if v > 0:
            self._row_spin.setValue(v - 1)

    def _row_down(self):
        v = self._row_spin.value()
        if v < self._row_spin.maximum():
            self._row_spin.setValue(v + 1)

    def _poll_epics(self):
        if not EPICS_AVAILABLE:
            return
        try:
            sx = epics.caget(PV_SX_RBV)
            sy = epics.caget(PV_SY_RBV)
            self._sx_rbv.setText(f"{sx:.2f}" if sx is not None else "---")
            self._sy_rbv.setText(f"{sy:.2f}" if sy is not None else "---")
        except Exception:
            pass

    def _check_instrument_busy(self):
        if not EPICS_AVAILABLE:
            return False
        try:
            busy = epics.caget(PV_DATA_COLLECTING)
            if busy:
                QMessageBox.warning(self, "Instrument busy",
                                    "Instrument is collecting data - cannot move motors.")
                return True
        except Exception:
            pass
        return False

    def _move_motor(self, pv_val: str, value: float):
        if self._check_instrument_busy():
            return
        try:
            epics.caput(pv_val, value)
        except Exception as e:
            self._status_lbl.setText(f"EPICS error: {e}")

    def _jog(self, axis: str, delta: float):
        if axis == "SX":
            new_val = self._sx_tar.value() + delta
            self._sx_tar.setValue(new_val)
            self._move_motor(PV_SX_VAL, new_val)
        else:
            new_val = self._sy_tar.value() + delta
            self._sy_tar.setValue(new_val)
            self._move_motor(PV_SY_VAL, new_val)

    def _drive_to_table(self):
        row = self._row_spin.value()
        if 0 <= row < len(self._ss.rows):
            r = self._ss.rows[row]
            self._sx_tar.setValue(r.sx)
            self._sy_tar.setValue(r.sy)
            self._move_motor(PV_SX_VAL, r.sx)
            self._move_motor(PV_SY_VAL, r.sy)

    def _save_position(self):
        if not EPICS_AVAILABLE:
            return
        try:
            sx = epics.caget(PV_SX_RBV)
            sy = epics.caget(PV_SY_RBV)
            if sx is not None and sy is not None:
                row = self._row_spin.value()
                if 0 <= row < len(self._ss.rows):
                    self._ss.rows[row].sx = float(sx)
                    self._ss.rows[row].sy = float(sy)
                    self._sx_tbl.setText(f"{sx:.3f}")
                    self._sy_tbl.setText(f"{sy:.3f}")
                    self._status_lbl.setText(f"Saved SX={sx:.3f} SY={sy:.3f} to row {row}")
        except Exception as e:
            self._status_lbl.setText(f"EPICS error: {e}")

    def _go_to_zero(self):
        self._sx_tar.setValue(0.0)
        self._sy_tar.setValue(0.0)
        self._move_motor(PV_SX_VAL, 0.0)
        self._move_motor(PV_SY_VAL, 0.0)

    def _stop_motors(self):
        if self._check_instrument_busy():
            return
        try:
            epics.caput(PV_ALL_STOP, 1)
        except Exception as e:
            self._status_lbl.setText(f"EPICS error: {e}")

    def _set_slits(self, mode: str):
        if self._check_instrument_busy():
            return
        try:
            if mode == "large":
                epics.caput(PV_C1M8, 2.5)
                epics.caput(PV_GSLIT1H, 2.8)
                epics.caput(PV_C1M7, 1.2)
                epics.caput(PV_GSLIT1V, 1.4)
            elif mode == "usaxs":
                sx = epics.caget(PV_USAXS_HSLIT)
                sy = epics.caget(PV_USAXS_VSLIT)
                epics.caput(PV_C1M8, sx)
                epics.caput(PV_C1M7, sy)
            elif mode == "saxswaxs":
                sx = epics.caget(PV_SAXS_HSLIT)
                sy = epics.caget(PV_SAXS_VSLIT)
                epics.caput(PV_C1M8, sx)
                epics.caput(PV_C1M7, sy)
        except Exception as e:
            self._status_lbl.setText(f"EPICS error: {e}")


# ---------------------------------------------------------------------------
# Main window
# ---------------------------------------------------------------------------

class SamplePlateSetupWindow(QMainWindow):
    """Main window for the Setup Sample Plates GUI tool."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Matilda - Setup Sample Plates")
        self.resize(1200, 750)

        # State
        self._current_set = SampleSet()
        self._saved_sets: dict[str, SampleSet] = {}
        self._hdf5_path: str | None = None
        self._unsaved = False
        self._usaxs_time = USAXS_SCAN_TIME_DEFAULT
        self._saxs_time = SAXS_SCAN_TIME_DEFAULT
        self._waxs_time = WAXS_SCAN_TIME_DEFAULT

        self._build_ui()
        self._load_template("9x9 Acrylic/magnetic plate")
        self._update_runtime()

    # ------------------------------------------------------------------
    # UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        # Toolbar row: set name + global checkboxes + status
        toolbar = QHBoxLayout()
        toolbar.addWidget(QLabel("Set name:"))
        self._set_name_edit = QLineEdit("MySamples")
        self._set_name_edit.setMaximumWidth(200)
        self._set_name_edit.textChanged.connect(self._on_name_changed)
        toolbar.addWidget(self._set_name_edit)

        self._usaxs_all = QCheckBox("USAXS all")
        self._usaxs_all.setChecked(True)
        self._saxs_all = QCheckBox("SAXS all")
        self._saxs_all.setChecked(True)
        self._waxs_all = QCheckBox("WAXS all")
        self._waxs_all.setChecked(True)
        for cb in (self._usaxs_all, self._saxs_all, self._waxs_all):
            cb.stateChanged.connect(self._on_global_check)
            toolbar.addWidget(cb)

        toolbar.addStretch()
        self._runtime_lbl = QLabel("Est. time: - min  |  Samples: 0")
        toolbar.addWidget(self._runtime_lbl)
        main_layout.addLayout(toolbar)

        # Status message
        self._status_lbl = QLabel("Tool started - define sample positions and export command file.")
        self._status_lbl.setStyleSheet("color: darkblue; font-style: italic;")
        main_layout.addWidget(self._status_lbl)

        # Tab widget
        tabs = QTabWidget()
        main_layout.addWidget(tabs)

        # -- Tab 1: Sample Table --
        tab1 = QWidget()
        tab1_lay = QHBoxLayout(tab1)

        # Left: table + controls
        left = QVBoxLayout()

        # Template selector + create button
        tmpl_row = QHBoxLayout()
        tmpl_row.addWidget(QLabel("Template:"))
        self._template_combo = QComboBox()
        self._template_combo.addItems(list(PLATE_DEFINITIONS.keys()))
        tmpl_row.addWidget(self._template_combo)
        tmpl_row.addWidget(QLabel("#:"))
        self._n_samples = QSpinBox()
        self._n_samples.setRange(1, 500)
        self._n_samples.setValue(20)
        tmpl_row.addWidget(self._n_samples)
        btn_create = QPushButton("Create New Set")
        btn_create.clicked.connect(self._on_create_set)
        tmpl_row.addWidget(btn_create)
        btn_add_lines = QPushButton("Add Lines")
        btn_add_lines.clicked.connect(self._on_add_lines)
        tmpl_row.addWidget(btn_add_lines)
        left.addLayout(tmpl_row)

        # Table
        self._table = SampleTable()
        self._table.dataChanged.connect(self._on_table_changed)
        left.addWidget(self._table)

        # Table action buttons
        btn_row = QHBoxLayout()
        btn_save_set = QPushButton("Save Set")
        btn_save_set.clicked.connect(self._on_save_set)
        btn_load_set = QPushButton("Load Set")
        btn_load_set.clicked.connect(self._on_load_set)
        btn_preview = QPushButton("Preview Commands")
        btn_preview.clicked.connect(self._on_preview)
        btn_export = QPushButton("Export usaxs.mac")
        btn_export.setStyleSheet("background-color: #88cc88; font-weight: bold;")
        btn_export.clicked.connect(self._on_export)
        btn_survey = QPushButton("Beamline Survey")
        btn_survey.setStyleSheet("background-color: #aaccff;")
        btn_survey.clicked.connect(self._on_beamline_survey)
        if not EPICS_AVAILABLE:
            btn_survey.setToolTip("pyepics not installed - EPICS motor control unavailable")
        for btn in (btn_save_set, btn_load_set, btn_preview, btn_export, btn_survey):
            btn_row.addWidget(btn)
        left.addLayout(btn_row)

        tab1_lay.addLayout(left, 60)

        # Right: plate image
        right = QVBoxLayout()
        right.addWidget(QLabel("Plate view - click to assign SX/SY to selected row:"))
        self._canvas = PlateCanvas()
        self._canvas.positionClicked.connect(self._on_plate_clicked)
        right.addWidget(self._canvas, 1)
        tab1_lay.addLayout(right, 40)

        tabs.addTab(tab1, "Sample Table")

        # -- Tab 2: Option Controls --
        tab2 = QWidget()
        tab2_lay = QFormLayout(tab2)
        self._usaxs_time_spin = QDoubleSpinBox()
        self._usaxs_time_spin.setRange(10, 600)
        self._usaxs_time_spin.setValue(USAXS_SCAN_TIME_DEFAULT)
        self._usaxs_time_spin.setSuffix(" s")
        self._usaxs_time_spin.valueChanged.connect(self._update_runtime)
        tab2_lay.addRow("USAXS scan time:", self._usaxs_time_spin)
        self._saxs_time_spin = QDoubleSpinBox()
        self._saxs_time_spin.setRange(0.1, 3600)
        self._saxs_time_spin.setValue(SAXS_SCAN_TIME_DEFAULT)
        self._saxs_time_spin.setSuffix(" s")
        self._saxs_time_spin.valueChanged.connect(self._update_runtime)
        tab2_lay.addRow("SAXS scan time:", self._saxs_time_spin)
        self._waxs_time_spin = QDoubleSpinBox()
        self._waxs_time_spin.setRange(0.1, 3600)
        self._waxs_time_spin.setValue(WAXS_SCAN_TIME_DEFAULT)
        self._waxs_time_spin.setSuffix(" s")
        self._waxs_time_spin.valueChanged.connect(self._update_runtime)
        tab2_lay.addRow("WAXS scan time:", self._waxs_time_spin)
        self._default_thick = QDoubleSpinBox()
        self._default_thick.setRange(0, 20)
        self._default_thick.setValue(DEFAULT_THICKNESS)
        self._default_thick.setSuffix(" mm")
        tab2_lay.addRow("Default sample thickness:", self._default_thick)
        tabs.addTab(tab2, "Option Controls")

        # -- Tab 3: Export Controls --
        tab3 = QWidget()
        tab3_lay = QVBoxLayout(tab3)

        order_row = QHBoxLayout()
        order_row.addWidget(QLabel("Export order:"))
        self._order_combo = QComboBox()
        self._order_combo.addItems(EXPORT_ORDERS)
        order_row.addWidget(self._order_combo)
        order_row.addStretch()
        tab3_lay.addLayout(order_row)

        fname_row = QHBoxLayout()
        fname_row.addWidget(QLabel("Command file name:"))
        self._cmd_fname = QLineEdit(DEFAULT_CMD_FILENAME)
        fname_row.addWidget(self._cmd_fname)
        tab3_lay.addLayout(fname_row)

        # Export mode
        mode_grp = QGroupBox("Export mode")
        mode_lay = QVBoxLayout(mode_grp)
        self._export_current = QCheckBox("Export current position set only")
        self._export_current.setChecked(True)
        self._export_list = QCheckBox("Export list of saved sets (drag to order below)")
        self._export_current.stateChanged.connect(
            lambda s: self._export_list.setChecked(not bool(s)))
        self._export_list.stateChanged.connect(
            lambda s: self._export_current.setChecked(not bool(s)))
        mode_lay.addWidget(self._export_current)
        mode_lay.addWidget(self._export_list)
        tab3_lay.addWidget(mode_grp)

        # Multi-set widget
        self._multi_export = MultiSetExportWidget()
        tab3_lay.addWidget(self._multi_export)

        # Export buttons
        exp_btns = QHBoxLayout()
        btn_exp_default = QPushButton("Export as usaxs.mac")
        btn_exp_default.setStyleSheet("background-color: #88cc88; font-weight: bold;")
        btn_exp_default.clicked.connect(self._on_export)
        btn_exp_named = QPushButton("Export with custom name...")
        btn_exp_named.clicked.connect(self._on_export_named)
        btn_exp_append = QPushButton("Append to existing file...")
        btn_exp_append.clicked.connect(self._on_export_append)
        exp_btns.addWidget(btn_exp_default)
        exp_btns.addWidget(btn_exp_named)
        exp_btns.addWidget(btn_exp_append)
        tab3_lay.addLayout(exp_btns)
        tabs.addTab(tab3, "Export Controls")

        # Menubar
        menubar = self.menuBar()
        file_menu = menubar.addMenu("File")
        file_menu.addAction("New", self._on_new)
        file_menu.addAction("Open HDF5...", self._on_open_hdf5)
        file_menu.addAction("Save HDF5", self._on_save_hdf5)
        file_menu.addAction("Save HDF5 As...", self._on_save_hdf5_as)
        file_menu.addSeparator()
        file_menu.addAction("Import 12IDB command file...", self._on_import_cmdfile)
        file_menu.addSeparator()
        file_menu.addAction("Exit", self.close)

        help_menu = menubar.addMenu("Help")
        help_menu.addAction("About", self._on_about)

    # ------------------------------------------------------------------
    # Event handlers
    # ------------------------------------------------------------------

    def _on_name_changed(self, text: str):
        self._current_set.name = text
        self._unsaved = True

    def _on_global_check(self):
        self._current_set.usaxs_all = self._usaxs_all.isChecked()
        self._current_set.saxs_all = self._saxs_all.isChecked()
        self._current_set.waxs_all = self._waxs_all.isChecked()
        self._update_runtime()

    def _on_table_changed(self):
        self._unsaved = True
        self._current_set.rows = self._table.get_sample_set()
        self._update_runtime()
        self._canvas.update_markers(
            self._current_set.rows,
            self._table.currentRow(),
        )

    def _on_create_set(self):
        if self._unsaved:
            reply = QMessageBox.question(
                self, "Unsaved changes",
                "Current set has unsaved changes. Create new set anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return
        plate_name = self._template_combo.currentText()
        self._load_template(plate_name)

    def _load_template(self, plate_name: str):
        n = self._n_samples.value()
        defn = PLATE_DEFINITIONS.get(plate_name)

        if plate_name == "Generic Grid":
            dlg = GenericGridDialog(self)
            if dlg.exec() != QDialog.DialogCode.Accepted:
                return
            centers = dlg.get_centers()
            n = dlg.total()
            rows = []
            for i, (cx, cy) in enumerate(centers):
                rows.append(SampleRow(name=f"Sample{i+1:03d}", sx=cx, sy=cy))
        elif plate_name == "AgBehenateLaB6":
            rows = [SampleRow(name="AgBehenateLaB6", sx=20.0, sy=20.0, usaxs=False)]
        elif defn is not None:
            centers = defn["centers"]
            rows = []
            for i in range(min(n, len(centers))):
                cx, cy = centers[i]
                rows.append(SampleRow(name=f"Sample{i+1:03d}", sx=cx, sy=cy))
            while len(rows) < n:
                rows.append(SampleRow())
        else:
            rows = [SampleRow() for _ in range(n)]

        self._current_set = SampleSet(name=self._set_name_edit.text(), rows=rows)
        self._table.load_from_sample_set(self._current_set)
        self._canvas.set_plate(plate_name)
        self._canvas.update_markers(self._current_set.rows)
        self._status_lbl.setText(f"Created new set: {plate_name}")
        self._unsaved = True
        self._update_runtime()

    def _on_add_lines(self):
        n = self._n_samples.value()
        for _ in range(n):
            self._table.insert_row_below()
        self._status_lbl.setText(f"Added {n} new lines")

    def _on_plate_clicked(self, sx: float, sy: float):
        row = self._table.currentRow()
        self._table.set_position_for_current_row(sx, sy)
        self._status_lbl.setText(f"Set row {row}: SX={sx:.3f} SY={sy:.3f}")
        self._canvas.update_markers(self._table.get_sample_set(), row)

    def _on_save_set(self):
        """Save current set into in-memory store (and to HDF5 if file is open)."""
        name = self._set_name_edit.text().strip()
        if not name:
            QMessageBox.warning(self, "Name required", "Please enter a set name.")
            return
        self._current_set.name = name
        self._current_set.rows = self._table.get_sample_set()
        self._current_set.usaxs_all = self._usaxs_all.isChecked()
        self._current_set.saxs_all = self._saxs_all.isChecked()
        self._current_set.waxs_all = self._waxs_all.isChecked()
        self._saved_sets[name] = self._current_set.copy()
        self._multi_export.set_saved_sets(list(self._saved_sets.keys()))
        if self._hdf5_path:
            save_sets_to_hdf5(self._hdf5_path, self._saved_sets)
        self._status_lbl.setText(f"Saved set '{name}'")
        self._unsaved = False

    def _on_load_set(self):
        names = list(self._saved_sets.keys())
        if not names:
            QMessageBox.information(self, "No saved sets",
                                    "No sets saved yet. Save a set first.")
            return
        from PyQt6.QtWidgets import QInputDialog
        name, ok = QInputDialog.getItem(self, "Load set", "Select set:", names, 0, False)
        if ok and name in self._saved_sets:
            if self._unsaved:
                reply = QMessageBox.question(
                    self, "Unsaved changes",
                    "Current set has unsaved changes. Load anyway?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                )
                if reply != QMessageBox.StandardButton.Yes:
                    return
            self._current_set = self._saved_sets[name].copy()
            self._set_name_edit.setText(self._current_set.name)
            self._table.load_from_sample_set(self._current_set)
            self._canvas.update_markers(self._current_set.rows)
            self._update_runtime()
            self._status_lbl.setText(f"Loaded set '{name}'")
            self._unsaved = False

    def _collect_sets_for_export(self) -> list[SampleSet]:
        if self._export_list.isChecked():
            names = self._multi_export.get_export_order()
            if not names:
                QMessageBox.warning(self, "No sets in export list",
                                    "Drag sets from the left list to the export order list.")
                return []
            return [self._saved_sets[n] for n in names if n in self._saved_sets]
        else:
            self._current_set.rows = self._table.get_sample_set()
            return [self._current_set]

    def _on_preview(self):
        sets = self._collect_sets_for_export()
        if not sets:
            return
        warn = check_for_blanks(sets[0])
        if warn:
            self._status_lbl.setText(warn)
        content = generate_command_file(
            sets,
            export_order=self._order_combo.currentText(),
            include_header=True,
        )
        dlg = QDialog(self)
        dlg.setWindowTitle("Command File Preview")
        dlg.resize(700, 500)
        lay = QVBoxLayout(dlg)
        txt = QTextEdit()
        txt.setReadOnly(True)
        txt.setPlainText(content)
        txt.setFont(QFont("Courier", 9))
        lay.addWidget(txt)
        btn = QPushButton("Close")
        btn.clicked.connect(dlg.accept)
        lay.addWidget(btn)
        dlg.exec()

    def _on_export(self):
        sets = self._collect_sets_for_export()
        if not sets:
            return
        warn = check_for_blanks(sets[0])
        if warn:
            self._status_lbl.setText(warn)
        content = generate_command_file(
            sets,
            export_order=self._order_combo.currentText(),
            include_header=True,
        )
        fname = self._cmd_fname.text().strip() or DEFAULT_CMD_FILENAME
        desktop = os.path.join(os.path.expanduser("~"), "Desktop")
        path = os.path.join(desktop, fname)
        with open(path, "w", encoding="utf-8") as f:
            f.write(content)
        self._status_lbl.setText(f"Exported command file: {path}")

    def _on_export_named(self):
        sets = self._collect_sets_for_export()
        if not sets:
            return
        warn = check_for_blanks(sets[0])
        if warn:
            self._status_lbl.setText(warn)
        content = generate_command_file(
            sets,
            export_order=self._order_combo.currentText(),
            include_header=True,
        )
        fname = self._cmd_fname.text().strip() or DEFAULT_CMD_FILENAME
        path, _ = QFileDialog.getSaveFileName(
            self, "Save command file", fname, "Command files (*.mac);;All files (*)")
        if path:
            with open(path, "w", encoding="utf-8") as f:
                f.write(content)
            self._status_lbl.setText(f"Exported: {path}")

    def _on_export_append(self):
        sets = self._collect_sets_for_export()
        if not sets:
            return
        content_new = generate_command_file(
            sets,
            export_order=self._order_combo.currentText(),
            include_header=False,
        )
        fname = self._cmd_fname.text().strip() or DEFAULT_CMD_FILENAME
        path, _ = QFileDialog.getSaveFileName(
            self, "Append to command file", fname,
            "Command files (*.mac);;All files (*)")
        if path:
            existing = ""
            if os.path.exists(path):
                with open(path, "r", encoding="utf-8") as f:
                    existing = f.read()
            with open(path, "w", encoding="utf-8") as f:
                f.write(existing)
                f.write("\n\n###  Appended commands\n\n")
                f.write(content_new)
            self._status_lbl.setText(f"Appended to: {path}")

    def _on_beamline_survey(self):
        self._current_set.rows = self._table.get_sample_set()
        dlg = BeamlineSurveyDialog(self._current_set, self)
        dlg.show()

    # -- File menu actions --

    def _on_new(self):
        if self._unsaved:
            reply = QMessageBox.question(
                self, "Unsaved changes",
                "Create a new empty set? Unsaved changes will be lost.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                return
        self._current_set = SampleSet()
        self._table.load_from_sample_set(self._current_set)
        self._canvas.update_markers([])
        self._update_runtime()
        self._unsaved = False

    def _on_open_hdf5(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open HDF5 file", "", "HDF5 files (*.h5 *.hdf5);;All files (*)")
        if path:
            self._saved_sets = load_sets_from_hdf5(path)
            self._hdf5_path = path
            self._multi_export.set_saved_sets(list(self._saved_sets.keys()))
            if self._saved_sets:
                first = next(iter(self._saved_sets.values()))
                self._current_set = first.copy()
                self._set_name_edit.setText(self._current_set.name)
                self._table.load_from_sample_set(self._current_set)
                self._canvas.update_markers(self._current_set.rows)
                self._update_runtime()
            self._status_lbl.setText(
                f"Loaded {len(self._saved_sets)} set(s) from {os.path.basename(path)}")
            self._unsaved = False

    def _on_save_hdf5(self):
        if self._hdf5_path is None:
            self._on_save_hdf5_as()
            return
        name = self._set_name_edit.text().strip() or "MySamples"
        self._current_set.rows = self._table.get_sample_set()
        self._saved_sets[name] = self._current_set.copy()
        save_sets_to_hdf5(self._hdf5_path, self._saved_sets)
        self._status_lbl.setText(f"Saved to {self._hdf5_path}")
        self._unsaved = False

    def _on_save_hdf5_as(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Save HDF5 file", "sample_positions.h5",
            "HDF5 files (*.h5 *.hdf5);;All files (*)")
        if path:
            self._hdf5_path = path
            self._on_save_hdf5()

    def _on_import_cmdfile(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Import command file", "",
            "Command files (*.mac *.txt);;CSV files (*.csv);;All files (*)")
        if not path:
            return
        rows = []
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = [p.strip().strip('"') for p in line.split(",")]
                    if len(parts) >= 3:
                        rows.append(SampleRow(
                            name=parts[0],
                            sx=float(parts[1]) if len(parts) > 1 else 0.0,
                            sy=float(parts[2]) if len(parts) > 2 else 0.0,
                            thickness=float(parts[4]) if len(parts) > 4 else DEFAULT_THICKNESS,
                            metadata=",".join(parts[5:]) if len(parts) > 5 else "",
                        ))
        except Exception as e:
            QMessageBox.critical(self, "Import error", str(e))
            return
        if rows:
            self._current_set = SampleSet(
                name=os.path.splitext(os.path.basename(path))[0], rows=rows)
            self._set_name_edit.setText(self._current_set.name)
            self._table.load_from_sample_set(self._current_set)
            self._canvas.update_markers(self._current_set.rows)
            self._update_runtime()
            self._status_lbl.setText(f"Imported {len(rows)} rows from {os.path.basename(path)}")

    def _on_about(self):
        QMessageBox.about(
            self, "About - Setup Sample Plates",
            "Matilda - Setup Sample Plates\n\n"
            "GUI tool for defining USAXS/SAXS/WAXS sample positions,\n"
            "exporting Bluesky command files, and optionally driving\n"
            "the sample stage via EPICS (Beamline Survey).\n\n"
            "Part of the Matilda data processing package.\n"
            "APS 12-ID-E, Argonne National Laboratory.",
        )

    def _update_runtime(self):
        self._current_set.rows = self._table.get_sample_set()
        n_u, n_s, n_w, t_min = estimate_run_time(
            self._current_set,
            self._usaxs_time_spin.value(),
            self._saxs_time_spin.value(),
            self._waxs_time_spin.value(),
        )
        total = max(n_u, n_s, n_w)
        self._runtime_lbl.setText(
            f"USAXS:{n_u} SAXS:{n_s} WAXS:{n_w} | "
            f"Est. time: {t_min} min")

    def closeEvent(self, event):
        if self._unsaved:
            reply = QMessageBox.question(
                self, "Unsaved changes",
                "The current set has unsaved changes. Exit anyway?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            )
            if reply != QMessageBox.StandardButton.Yes:
                event.ignore()
                return
        event.accept()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def run_sample_plate_setup():
    """Launch the Setup Sample Plates GUI.  Creates QApplication if needed."""
    app = QApplication.instance() or QApplication([])
    win = SamplePlateSetupWindow()
    win.show()
    app.exec()


if __name__ == "__main__":
    run_sample_plate_setup()
