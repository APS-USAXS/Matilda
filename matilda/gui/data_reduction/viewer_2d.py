"""
matilda.gui.data_reduction.viewer_2d
=====================================
Non-modal dialog that shows the 2D detector image for a SAXS or WAXS file.

Data source: ``/entry/data/data`` (raw 2D pixel array) read from HDF5.

Features
--------
- pyqtgraph ImageView with built-in LUT histogram / draggable levels
- Log / linear intensity toggle (applies np.log10 to data before display)
- "Reset levels" button to auto-level the display
- Header label showing filename and array dimensions
"""

import os

import h5py
import numpy as np
import pyqtgraph as pg

from .._qt import (
    QDialog, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QCheckBox,
)

_HDF5_DATA_PATH = "/entry/data/data"


class Viewer2D(QDialog):
    """2D detector image viewer dialog (SAXS / WAXS)."""

    def __init__(
        self,
        path: str,
        filename: str,
        technique: str,
        parent=None,
    ):
        super().__init__(parent)
        self.setWindowTitle(f"2D Image — {filename}  [{technique}]")
        self.resize(720, 620)
        self.setModal(False)

        self._path     = path
        self._filename = filename
        self._data: np.ndarray | None = None

        self._build_ui()
        self._load_data()

    # ── UI ────────────────────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.setSpacing(4)

        # Info header
        self._info_label = QLabel("Loading…")
        self._info_label.setStyleSheet("font-size: 11px; color: #444;")
        layout.addWidget(self._info_label)

        # Toolbar
        tb = QHBoxLayout()

        self._chk_log = QCheckBox("Log intensity")
        self._chk_log.setChecked(True)
        self._chk_log.setToolTip(
            "Display log10(intensity); negative values shown as NaN (transparent)"
        )
        self._chk_log.toggled.connect(self._update_display)
        tb.addWidget(self._chk_log)

        self._btn_reset = QPushButton("Reset levels")
        self._btn_reset.setToolTip("Auto-scale the colour/intensity range")
        self._btn_reset.clicked.connect(self._reset_levels)
        tb.addWidget(self._btn_reset)

        tb.addStretch()
        layout.addLayout(tb)

        # ImageView (pyqtgraph)
        self._image_view = pg.ImageView()
        # Hide the ROI and menu buttons — keep only the histogram/levels bar
        self._image_view.ui.roiBtn.hide()
        self._image_view.ui.menuBtn.hide()
        layout.addWidget(self._image_view, 1)

    # ── Data loading ──────────────────────────────────────────────────────────

    def _load_data(self):
        full_path = os.path.join(self._path, self._filename)
        try:
            with h5py.File(full_path, "r") as f:
                raw = f[_HDF5_DATA_PATH][()]
            data = np.asarray(raw, dtype=float)
            # Collapse leading frame dimensions if 3-D or higher
            while data.ndim > 2:
                data = data[0]
            self._data = data
            rows, cols = data.shape
            self._info_label.setText(
                f"{self._filename}   —   {rows} × {cols} pixels"
            )
            self._update_display()
        except KeyError:
            self._info_label.setText(
                f"No data found at '{_HDF5_DATA_PATH}' in {self._filename}"
            )
        except Exception as exc:
            self._info_label.setText(f"Error loading image: {exc}")

    # ── Display update ────────────────────────────────────────────────────────

    def _update_display(self):
        if self._data is None:
            return
        data = self._data
        if self._chk_log.isChecked():
            with np.errstate(divide="ignore", invalid="ignore"):
                display = np.where(data > 0, np.log10(data), np.nan)
        else:
            display = data
        # Transpose so that row 0 is at the bottom (physical detector orientation)
        self._image_view.setImage(display.T, autoLevels=True)

    def _reset_levels(self):
        self._image_view.autoLevels()
