"""
matilda.gui.data_reduction.graph_panel
=======================================
Right panel: pyqtgraph log-log PlotWidget.

Displays intermediate data-reduction curves:
  USAXS (Flyscan / StepScan):  R_data (raw), SMR (slit-smeared), DSM (desmeared)
  SAXS / WAXS:                 normalized (raw 1D), subtracted, calibrated
"""

import csv

import numpy as np
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
    )
    from PySide6.QtCore import Qt
except ImportError:
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QFileDialog,
    )
    from PyQt6.QtCore import Qt

# White background, black foreground — must be set before any pg widget is created.
pg.setConfigOption("background", "w")
pg.setConfigOption("foreground", "k")


# Colorblind-friendly palette (blue, orange, green, purple, sky-blue, dark-red, yellow)
_PALETTE = [
    (0,   114, 189),
    (217,  83,  25),
    ( 32, 178,  34),
    (126,  47, 142),
    ( 77, 190, 238),
    (162,  20,  47),
    (237, 177,  32),
]

# Maps technique → ordered list of (result-dict key, display label)
_CURVE_SPECS: dict[str, list[tuple[str, str]]] = {
    "Flyscan": [
        ("R_data", "Raw (R_data)"),
        ("SMR",    "Slit-smeared (SMR)"),
        ("DSM",    "Desmeared (DSM)"),
    ],
    "StepScan": [
        ("R_data", "Raw (R_data)"),
        ("SMR",    "Slit-smeared (SMR)"),
        ("DSM",    "Desmeared (DSM)"),
    ],
    "SAXS": [
        ("normalized", "Normalized (raw 1D)"),
        ("subtracted", "Blank-subtracted"),
        ("calibrated", "Calibrated"),
    ],
    "WAXS": [
        ("normalized", "Normalized (raw 1D)"),
        ("subtracted", "Blank-subtracted"),
        ("calibrated", "Calibrated"),
    ],
}


class GraphPanel(QWidget):
    """pyqtgraph log-log plot for reduction intermediate curves."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._curves: list[pg.PlotDataItem] = []
        self._last_result: dict | None = None
        self._last_technique: str | None = None
        self._build_ui()

    # ── UI ────────────────────────────────────────────────────────────────────

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

        self._btn_clear = QPushButton("Clear")
        self._btn_clear.setToolTip("Remove all curves")
        self._btn_clear.clicked.connect(self.clear)
        tb.addWidget(self._btn_clear)

        tb.addStretch()

        self._btn_csv = QPushButton("Export CSV")
        self._btn_csv.setToolTip("Export displayed curves to CSV")
        self._btn_csv.clicked.connect(self._export_csv)
        tb.addWidget(self._btn_csv)

        self._btn_png = QPushButton("Export PNG")
        self._btn_png.setToolTip("Export plot as PNG image")
        self._btn_png.clicked.connect(self._export_png)
        tb.addWidget(self._btn_png)

        layout.addLayout(tb)

        # Plot widget
        self._plot_widget = pg.GraphicsLayoutWidget()
        self._plot = self._plot_widget.addPlot()
        self._plot.setLogMode(x=True, y=True)
        self._plot.showGrid(x=True, y=True, alpha=0.3)
        self._legend = self._plot.addLegend(offset=(10, 10))
        self._plot.setLabel("bottom", "Q", units="Å⁻¹")
        self._plot.setLabel("left", "Intensity", units="cm⁻¹")

        layout.addWidget(self._plot_widget, 1)

    # ── Public API ────────────────────────────────────────────────────────────

    def clear(self):
        """Remove all curves from the plot."""
        for c in self._curves:
            self._plot.removeItem(c)
        self._curves.clear()
        # Rebuild legend (removeItem leaves stale entries)
        self._plot.removeItem(self._legend)
        self._legend = self._plot.addLegend(offset=(10, 10))
        self._last_result = None
        self._last_technique = None

    def update_curves(self, result: dict, technique: str):
        """Plot intermediate reduction curves from a result dict.

        Parameters
        ----------
        result:    dict returned by processFlyscan / processStepscan / process2Ddata
        technique: one of 'Flyscan', 'StepScan', 'SAXS', 'WAXS'
        """
        self.clear()
        self._last_result = result
        self._last_technique = technique

        specs = _CURVE_SPECS.get(technique, [])
        curve_data = _extract_curves(result, technique)

        for i, (key, label) in enumerate(specs):
            if key not in curve_data:
                continue
            q, inten, _err = curve_data[key]
            if q is None or inten is None or len(q) == 0:
                continue
            # Strip non-positive values so log scale renders correctly.
            # pyqtgraph silently drops them, but the auto-range then fails.
            mask = (q > 0) & (inten > 0) & np.isfinite(q) & np.isfinite(inten)
            q_plot = q[mask]
            i_plot = inten[mask]
            if len(q_plot) < 2:
                continue
            color = _PALETTE[i % len(_PALETTE)]
            pen = pg.mkPen(color=color, width=1.5)
            item = self._plot.plot(q_plot, i_plot, pen=pen, name=label)
            self._curves.append(item)

    # ── Private helpers ───────────────────────────────────────────────────────

    def _update_log_mode(self):
        self._plot.setLogMode(
            x=self._btn_xlog.isChecked(),
            y=self._btn_ylog.isChecked(),
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
        curve_data = _extract_curves(self._last_result, technique)
        specs = _CURVE_SPECS.get(technique, [])
        with open(path, "w", newline="") as f:
            writer = csv.writer(f)
            for key, label in specs:
                if key not in curve_data:
                    continue
                q, inten, err = curve_data[key]
                if q is None:
                    continue
                writer.writerow([f"# {label}"])
                writer.writerow(["Q (1/A)", "Intensity (1/cm)", "Error"])
                err_vals = err if err is not None else [0.0] * len(q)
                for row in zip(q, inten, err_vals):
                    writer.writerow(row)
                writer.writerow([])


# ── Result-dict → curve extraction ───────────────────────────────────────────

def _safe(d: dict, key: str):
    """Return float array from *d[key]*, or None if absent / unusable."""
    v = d.get(key)
    if v is None:
        return None
    try:
        arr = np.asarray(v, dtype=float).ravel()
    except (TypeError, ValueError):
        return None
    return arr if arr.size > 1 else None


def _extract_curves(result: dict, technique: str) -> dict[str, tuple]:
    """Map result dict → {curve_key: (Q, Intensity, Error)} triples."""
    curves: dict[str, tuple] = {}

    if technique in ("Flyscan", "StepScan"):
        # R_data: raw, gain-corrected, beam-centre-corrected
        rd = result.get("reducedData", {})
        q = _safe(rd, "Q")
        i = _safe(rd, "Intensity")
        e = _safe(rd, "Error")
        if q is not None and i is not None:
            curves["R_data"] = (q, i, e)

        cd = result.get("CalibratedData", {})

        # SMR: slit-smeared raw (before desmearing)
        sq = _safe(cd, "SMR_Qvec")
        si = _safe(cd, "SMR_Int")
        se = _safe(cd, "SMR_Error")
        if sq is not None and si is not None:
            curves["SMR"] = (sq, si, se)

        # DSM: desmeared result
        dq = _safe(cd, "Q")
        di = _safe(cd, "Intensity")
        de = _safe(cd, "Error")
        if dq is not None and di is not None:
            curves["DSM"] = (dq, di, de)

    elif technique in ("SAXS", "WAXS"):
        # Normalised raw 1D from reducedData
        rd = result.get("reducedData", {})
        q = _safe(rd, "Q")
        i = _safe(rd, "Intensity")
        e = _safe(rd, "Error")
        if q is not None and i is not None:
            curves["normalized"] = (q, i, e)

        # Calibrated data (blank-subtracted + absolute scale)
        cd = result.get("CalibratedData", {})
        cq = _safe(cd, "Q")
        ci = _safe(cd, "Intensity")
        ce = _safe(cd, "Error")
        if cq is not None and ci is not None:
            # process2Ddata currently returns one calibrated dataset;
            # expose it as both "subtracted" and "calibrated" until the
            # converter is split into intermediate steps.
            curves["subtracted"] = (cq, ci, ce)
            curves["calibrated"] = (cq, ci, ce)

    return curves
