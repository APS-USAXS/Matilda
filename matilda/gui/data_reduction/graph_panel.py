"""
matilda.gui.data_reduction.graph_panel
=======================================
Right panel: dual-axis pyqtgraph log-log plot.

Axis layout
-----------
Left  axis (cm⁻¹): calibrated / processed curves
  USAXS: Slit-smeared (SMR), Desmeared (DSM)
  SAXS/WAXS: Calibrated

Right axis (arb. units): raw / unnormalised curve (dashed grey)
  USAXS: Raw (R_data)
  SAXS/WAXS: Normalized (raw 1D)

Both axes share the same Q (X) axis.  Log-log mode is the default.
Error bars are shown on calibrated (left-axis) curves only and can be
toggled with the toolbar button.

A "View 2D Image" button appears for SAXS/WAXS files only.
"""

import csv
import os

import numpy as np
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QPen
except ImportError:
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal
    from PyQt6.QtGui import QPen

from .sas_plot import make_sas_plot, draw_error_bars, set_robust_y_range

# White background, black foreground — set before any pg widget is created.
pg.setConfigOption("background", "w")
pg.setConfigOption("foreground", "k")

# ── Palette ──────────────────────────────────────────────────────────────────

# Right-axis raw curves: (color, Qt.PenStyle, width)
# Index 0 = sample raw, Index 1 = blank raw
_RAW_STYLES = [
    ((100, 100, 100), Qt.PenStyle.DashLine, 1.5),    # grey dashed  — sample raw
    (( 40, 180, 170), Qt.PenStyle.DotLine,  1.5),    # teal dotted  — blank raw
]


# ── Main widget ───────────────────────────────────────────────────────────────

class GraphPanel(QWidget):
    """Dual-axis log-log plot for reduction intermediate curves."""

    # Emitted when the "View 2D Image" button is clicked (SAXS/WAXS only)
    view_2d_requested = Signal(str, str, str)   # path, filename, technique

    def __init__(self, parent=None):
        super().__init__(parent)
        self._left_items: list  = []       # items on the left PlotItem
        self._error_bar_items: list = []   # subset of _left_items (error bars)
        self._right_items: list = []       # items on the right ViewBox
        self._last_result:    dict | None = None
        self._last_technique: str  | None = None
        self._last_path:      str  | None = None
        self._last_filename:  str  | None = None
        self._build_ui()

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(4, 4, 4, 4)
        layout.setSpacing(4)

        # Toolbar
        tb = QHBoxLayout()

        self._btn_xlog = QPushButton("X log")
        self._btn_xlog.setCheckable(True)
        self._btn_xlog.setChecked(True)
        self._btn_xlog.setToolTip("Toggle X-axis log scale")
        self._btn_xlog.clicked.connect(self._update_log_mode)
        tb.addWidget(self._btn_xlog)

        self._btn_ylog = QPushButton("Y log")
        self._btn_ylog.setCheckable(True)
        self._btn_ylog.setChecked(True)
        self._btn_ylog.setToolTip("Toggle Y-axis log scale")
        self._btn_ylog.clicked.connect(self._update_log_mode)
        tb.addWidget(self._btn_ylog)

        self._btn_errbar = QPushButton("Error bars")
        self._btn_errbar.setCheckable(True)
        self._btn_errbar.setChecked(False)
        self._btn_errbar.setToolTip("Show/hide error bars on calibrated curves")
        self._btn_errbar.clicked.connect(self._toggle_error_bars)
        tb.addWidget(self._btn_errbar)

        self._btn_clear = QPushButton("Clear")
        self._btn_clear.setToolTip("Remove all curves from the plot")
        self._btn_clear.clicked.connect(self.clear)
        tb.addWidget(self._btn_clear)

        tb.addStretch()

        self._btn_2d = QPushButton("View 2D Image")
        self._btn_2d.setToolTip(
            "Open 2D detector image viewer (SAXS/WAXS only)"
        )
        self._btn_2d.setVisible(False)
        self._btn_2d.clicked.connect(self._on_view_2d)
        tb.addWidget(self._btn_2d)

        self._btn_csv = QPushButton("Export CSV")
        self._btn_csv.setToolTip("Export displayed curves to CSV")
        self._btn_csv.clicked.connect(self._export_csv)
        tb.addWidget(self._btn_csv)

        self._btn_png = QPushButton("Export PNG")
        self._btn_png.setToolTip("Export plot as PNG image")
        self._btn_png.clicked.connect(self._export_png)
        tb.addWidget(self._btn_png)

        layout.addLayout(tb)

        # Graphics layout
        self._gl = pg.GraphicsLayoutWidget()
        self._plot = make_sas_plot(
            self._gl,
            row=0, col=0,
            x_label="Q  (Å⁻¹)",
            y_label="Intensity  (cm⁻¹)",
            parent_widget=self,
            jpeg_default_name="matilda_graph",
        )
        self._legend = self._plot.addLegend(offset=(10, 10))

        # ── Right Y axis (raw/arb. units) ─────────────────────────────────
        self._right_vb = pg.ViewBox()
        self._right_ax = pg.AxisItem("right")
        self._right_ax.setLogMode(True)
        self._right_ax.enableAutoSIPrefix(False)
        self._right_ax.setLabel("Intensity  (arb. units)")

        # Row 2, col 3 = right side of the plot layout grid
        self._plot.layout.addItem(self._right_ax, 2, 3)
        self._plot.scene().addItem(self._right_vb)
        self._right_ax.linkToView(self._right_vb)
        self._right_vb.setXLink(self._plot)
        self._plot.vb.sigResized.connect(self._sync_right_vb)

        layout.addWidget(self._gl, 1)

    def _sync_right_vb(self):
        """Keep right ViewBox geometry in sync with the main ViewBox."""
        self._right_vb.setGeometry(self._plot.vb.sceneBoundingRect())
        self._right_vb.linkedViewChanged(self._plot.vb, self._right_vb.XAxis)

    # ── Public API ────────────────────────────────────────────────────────────

    def clear(self):
        """Remove all curves and clear the legend."""
        # Remove each left-axis item from both the plot and the legend.
        # Note: legend.removeItem(item) searches for the item by identity;
        # it is a no-op if the item was never added to the legend (e.g. error
        # bar segments added without a name).  Do NOT use legend.setParentItem
        # (None) — pyqtgraph overrides setParentItem and calls self.anchor()
        # which requires a parent, raising "Cannot anchor; parent is not set."
        for item in self._left_items:
            self._legend.removeItem(item)
            self._plot.removeItem(item)
        self._left_items.clear()
        self._error_bar_items.clear()

        for item in self._right_items:
            self._legend.removeItem(item)
            self._right_vb.removeItem(item)
        self._right_items.clear()

        # Do NOT call setLimits(xMin=None, ...) here — passing None to
        # setLimits triggers _effectiveLimits() which does max(None, float)
        # and raises TypeError in Python 3 (pyqtgraph bug, some versions).
        # Limits are always overwritten by update_curves() when new data
        # arrives, so no explicit reset is needed.

        self._last_result    = None
        self._last_technique = None
        self._last_path      = None
        self._last_filename  = None
        self._btn_2d.setVisible(False)

    def update_curves(
        self,
        result: dict,
        technique: str,
        path: str = "",
        filename: str = "",
    ):
        """Plot reduction intermediate curves from a result dict.

        Parameters
        ----------
        result :    dict returned by processFlyscan / processStepscan / process2Ddata
        technique : 'Flyscan', 'StepScan', 'SAXS', or 'WAXS'
        path :      directory containing the file (for 2D viewer)
        filename :  filename (for 2D viewer)
        """
        self.clear()
        self._last_result    = result
        self._last_technique = technique
        self._last_path      = path
        self._last_filename  = filename

        self._btn_2d.setVisible(technique in ("SAXS", "WAXS"))

        raw_curves, calibrated = _extract_curves(result, technique)

        # Accumulate all Q arrays so we can set a shared X limit at the end.
        all_q: list[np.ndarray] = []

        # ── Right axis: raw curves (arb. units, dashed) ───────────────────
        # Includes sample raw and (if available) blank raw overlay.
        all_right_I: list[np.ndarray] = []
        for i, (q, I, _dI, label) in enumerate(raw_curves):
            mask = (q > 0) & (I > 0) & np.isfinite(q) & np.isfinite(I)
            q_, I_ = q[mask], I[mask]
            if len(q_) < 2:
                continue
            all_q.append(q_)
            all_right_I.append(I_)
            # Pre-log10 transform: right ViewBox has no setLogMode
            q_log = np.log10(q_)
            I_log = np.log10(I_)
            color, style, width = _RAW_STYLES[i % len(_RAW_STYLES)]
            pen  = pg.mkPen(color=color, width=width, style=style)
            item = pg.PlotDataItem(q_log, I_log, pen=pen, name=label)
            self._right_vb.addItem(item)
            self._right_items.append(item)
            self._legend.addItem(item, label)

        # Set right Y range: show full data top-to-bottom, limit decades.
        # USAXS data can span ~13 decades; cap display at 10.
        # SAXS/WAXS typically spans less; cap at 5.
        if all_right_I:
            combined_r = np.concatenate(all_right_I)
            valid_r = combined_r[combined_r > 0]
            if len(valid_r) >= 3:
                log_r = np.log10(valid_r)
                hi_r  = float(np.max(log_r)) + 0.3
                max_decades = 10 if technique in ("Flyscan", "StepScan") else 5
                lo_r  = max(float(np.min(log_r)) - 0.3, hi_r - max_decades)
                self._right_vb.enableAutoRange(axis=pg.ViewBox.YAxis, enable=False)
                self._right_vb.setYRange(lo_r, hi_r, padding=0)
                self._right_vb.setLimits(
                    yMin=float(np.min(log_r)) - 1,
                    yMax=float(np.max(log_r)) + 1,
                )

        # ── Left axis: calibrated / cm⁻¹ (solid colored lines) ───────────
        all_cal_I: list[np.ndarray] = []
        for i, (q, I, dI, label, color) in enumerate(calibrated):
            mask = (q > 0) & (I > 0) & np.isfinite(q) & np.isfinite(I)
            q_, I_ = q[mask], I[mask]
            if len(q_) < 2:
                continue
            all_q.append(q_)
            all_cal_I.append(I_)

            pen   = pg.mkPen(color=color, width=1.5)
            scatter = self._plot.plot(q_, I_, pen=pen, name=label)
            self._left_items.append(scatter)

            if dI is not None:
                dI_ = np.asarray(dI, dtype=float)
                if dI_.shape == q.shape:
                    dI_ = dI_[mask]
                elif dI_.shape != q_.shape:
                    dI_ = None

                if dI_ is not None:
                    valid_I = I_[I_ > 0]
                    y_global_max = None
                    if len(valid_I) >= 5:
                        y_global_max = 10.0 ** (
                            float(np.percentile(np.log10(valid_I), 99)) + 3
                        )
                    eb = draw_error_bars(
                        self._plot, q_, I_, dI_,
                        y_global_max=y_global_max,
                        color=color,
                    )
                    if eb is not None:
                        self._left_items.append(eb)
                        self._error_bar_items.append(eb)
                        eb.setVisible(self._btn_errbar.isChecked())

        # Set left Y range: percentile-based view, 1-decade hard limits.
        if all_cal_I:
            combined_c = np.concatenate(all_cal_I)
            set_robust_y_range(self._plot, combined_c)
            # Override limits to 1 decade outside actual data min/max
            valid_c = combined_c[combined_c > 0]
            if len(valid_c) >= 2:
                self._plot.getViewBox().setLimits(
                    yMin=np.log10(float(valid_c.min())) - 1,
                    yMax=np.log10(float(valid_c.max())) + 1,
                )

        # ── Shared X axis: set initial view + 1-decade hard limit ─────────
        if all_q:
            q_all   = np.concatenate(all_q)
            q_valid = q_all[q_all > 0]
            if len(q_valid) >= 2:
                xlo_lim = np.log10(float(q_valid.min())) - 1
                xhi_lim = np.log10(float(q_valid.max())) + 1
                self._plot.getViewBox().setLimits(xMin=xlo_lim, xMax=xhi_lim)
                # Also limit the right ViewBox x range
                self._right_vb.setLimits(xMin=xlo_lim, xMax=xhi_lim)
                # Set initial X view range
                self._plot.setXRange(
                    np.log10(float(q_valid.min())),
                    np.log10(float(q_valid.max())),
                    padding=0.05,
                )

    # ── Private helpers ───────────────────────────────────────────────────────

    def _toggle_error_bars(self):
        visible = self._btn_errbar.isChecked()
        for item in self._error_bar_items:
            item.setVisible(visible)

    def _update_log_mode(self):
        self._plot.setLogMode(
            x=self._btn_xlog.isChecked(),
            y=self._btn_ylog.isChecked(),
        )

    def _on_view_2d(self):
        if self._last_path and self._last_filename and self._last_technique:
            self.view_2d_requested.emit(
                self._last_path, self._last_filename, self._last_technique
            )

    def _export_png(self):
        path, _ = QFileDialog.getSaveFileName(
            self, "Export PNG", "", "PNG Files (*.png)"
        )
        if not path:
            return
        exporter = pg.exporters.ImageExporter(self._plot)
        exporter.export(path)

    def _export_csv(self):
        if self._last_result is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export CSV", "", "CSV Files (*.csv)"
        )
        if not path:
            return
        technique = self._last_technique or "Unknown"
        raw_curves, calibrated = _extract_curves(self._last_result, technique)
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            for q, I, dI, label in raw_curves:
                writer.writerow([f"# {label}"])
                writer.writerow(["Q (1/A)", "Intensity (arb.)", "Error"])
                err = dI if dI is not None else [0.0] * len(q)
                for row in zip(q, I, err):
                    writer.writerow(row)
                writer.writerow([])
            for q, I, dI, label, _color in calibrated:
                writer.writerow([f"# {label}"])
                writer.writerow(["Q (1/A)", "Intensity (1/cm)", "Error"])
                err = dI if dI is not None else [0.0] * len(q)
                for row in zip(q, I, err):
                    writer.writerow(row)
                writer.writerow([])


# ── Result-dict → curve extraction ───────────────────────────────────────────

def _safe(d: dict, key: str) -> np.ndarray | None:
    """Return float array from *d[key]*, or None if absent or unusable."""
    v = d.get(key)
    if v is None:
        return None
    try:
        arr = np.asarray(v, dtype=float).ravel()
    except (TypeError, ValueError):
        return None
    return arr if arr.size > 1 else None


# Colors for calibrated curves by role
_COLOR_INTERMEDIATE = (0,   114, 189)   # blue — SMR (intermediate)
_COLOR_FINAL        = (204,  30,  30)   # red  — DSM / Calibrated (final result)


def _extract_curves(
    result: dict,
    technique: str,
) -> tuple[list[tuple], list[tuple]]:
    """Parse a result dict into (raw_curves, calibrated_curves).

    Returns
    -------
    raw_curves : list of (Q, I, dI, label)
        Curves for the right (arb. units) axis.
        Index 0 = sample raw; index 1 = blank raw (if available).
    calibrated : list of (Q, I, dI, label, color)
        One or more curves for the left (cm⁻¹) axis.
        *color* is an (R, G, B) tuple.
    """
    raw_curves: list[tuple] = []
    calibrated: list[tuple] = []

    if technique in ("Flyscan", "StepScan"):
        rd = result.get("reducedData", {})
        q  = _safe(rd, "Q")
        I  = _safe(rd, "Intensity")
        dI = _safe(rd, "Error")
        if q is not None and I is not None:
            raw_curves.append((q, I, dI, "Raw (R_data)"))

        # Blank raw curve
        brd = result.get("BlankData")
        if brd:
            bq = _safe(brd, "Q")
            bI = _safe(brd, "Intensity")
            if bq is not None and bI is not None:
                raw_curves.append((bq, bI, None, "Blank (R_data)"))

        cd = result.get("CalibratedData", {})

        sq = _safe(cd, "SMR_Qvec")
        si = _safe(cd, "SMR_Int")
        se = _safe(cd, "SMR_Error")
        if sq is not None and si is not None:
            calibrated.append((sq, si, se, "Slit-smeared (SMR)", _COLOR_INTERMEDIATE))

        dq = _safe(cd, "Q")
        di = _safe(cd, "Intensity")
        de = _safe(cd, "Error")
        if dq is not None and di is not None:
            calibrated.append((dq, di, de, "Desmeared (DSM)", _COLOR_FINAL))

    elif technique in ("SAXS", "WAXS"):
        rd = result.get("reducedData", {})
        q  = _safe(rd, "Q")
        I  = _safe(rd, "Intensity")
        dI = _safe(rd, "Error")
        if q is not None and I is not None:
            raw_curves.append((q, I, dI, "Normalized (raw 1D)"))

        # Blank raw curve
        brd = result.get("BlankData")
        if brd:
            bq = _safe(brd, "Q")
            bI = _safe(brd, "Intensity")
            if bq is not None and bI is not None:
                raw_curves.append((bq, bI, None, "Blank (normalized)"))

        cd = result.get("CalibratedData", {})
        cq = _safe(cd, "Q")
        ci = _safe(cd, "Intensity")
        ce = _safe(cd, "Error")
        if cq is not None and ci is not None:
            calibrated.append((cq, ci, ce, "Calibrated", _COLOR_FINAL))

    return raw_curves, calibrated
