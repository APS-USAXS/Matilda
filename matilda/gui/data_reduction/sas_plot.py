"""
Reusable SAS I(Q) log-log plot utilities for Matilda.

Adapted from pyirena/gui/sas_plot.py (https://github.com/jilavsky/pyirena).
Modifications vs pyirena original:
  - ``_draw_error_bars`` renamed to ``draw_error_bars`` (made public for import).
  - ``plot_iq_data`` accepts optional ``brush`` and ``symbol_size`` keyword
    arguments to allow per-call style overrides.

Design notes
------------
* All public helpers work with pyqtgraph ``PlotItem`` objects.
* Data values passed to ``plot_iq_data`` / ``plot_iq_model`` must be in
  **physical (linear) units**.  ``make_sas_plot`` calls
  ``setLogMode(x=True, y=True)`` so pyqtgraph applies the log10 transform
  internally — you never pre-transform the data.
* Cursor positions and ``TextItem`` positions go into the **ViewBox
  coordinate system**, which equals log10(data) when setLogMode is active.
  Use ``get_cursor_q_range()`` to read back cursor positions in linear units.
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
import pyqtgraph as pg

try:
    from PySide6.QtWidgets import QFileDialog, QMessageBox
    from PySide6.QtCore import Qt
except ImportError:
    from PyQt6.QtWidgets import QFileDialog, QMessageBox
    from PyQt6.QtCore import Qt


# ===========================================================================
# Style constants
# ===========================================================================

class SASPlotStyle:
    """Central definition of visual style for all Matilda SAS plots."""

    DATA_BRUSH       = pg.mkBrush(44, 62, 80, 200)    # dark slate
    DATA_SIZE        = 5

    ERROR_PEN        = pg.mkPen((127, 140, 141, 180), width=1)
    ERROR_CAP_FRAC   = 0.05

    FIT_PEN          = pg.mkPen(color=(200, 30, 30), width=2)

    RESID_BRUSH      = pg.mkBrush(127, 140, 141, 200)
    RESID_SIZE       = 4

    CURSOR_A_PEN     = pg.mkPen('#e74c3c', width=2)
    CURSOR_B_PEN     = pg.mkPen('#3498db', width=2)
    CURSOR_A_COLOR   = '#e74c3c'
    CURSOR_B_COLOR   = '#3498db'

    ANNOT_COLOR      = (40, 40, 40)
    GRID_ALPHA       = 0.3


# ===========================================================================
# _SafeInfiniteLine
# ===========================================================================

class _SafeInfiniteLine(pg.InfiniteLine):
    """pg.InfiniteLine subclass that prevents PySide6/Shiboken segfaults."""

    def mouseMoveEvent(self, ev):
        try:
            super().mouseMoveEvent(ev)
        except Exception:
            ev.ignore()

    def mouseDragEvent(self, ev):
        try:
            super().mouseDragEvent(ev)
        except Exception:
            ev.ignore()


# ===========================================================================
# Plot factory
# ===========================================================================

def make_sas_plot(
    graphics_layout: pg.GraphicsLayoutWidget,
    row: int,
    col: int,
    x_label: str = 'Q  (Å⁻¹)',
    y_label: str = 'I',
    title: str | None = None,
    x_link=None,
    log_x: bool = True,
    log_y: bool = True,
    parent_widget=None,
    jpeg_default_name: str = 'matilda_graph',
) -> pg.PlotItem:
    """Create a log-log (or semi-log) SAS plot with standard style.

    Parameters
    ----------
    graphics_layout : pg.GraphicsLayoutWidget
    row, col : int
    x_label, y_label : str
    title : str or None
    x_link : PlotItem or None
    log_x, log_y : bool
    parent_widget : QWidget or None
        Parent for JPEG export dialog.  Pass ``None`` to suppress.
    jpeg_default_name : str

    Returns
    -------
    pg.PlotItem
    """
    plot = graphics_layout.addPlot(row=row, col=col)
    plot.setLogMode(x=log_x, y=log_y)
    plot.setLabel('left',   y_label)
    plot.setLabel('bottom', x_label)
    plot.showGrid(x=True, y=True, alpha=SASPlotStyle.GRID_ALPHA)
    plot.getAxis('left').enableAutoSIPrefix(False)
    plot.getAxis('bottom').enableAutoSIPrefix(False)
    if title:
        plot.setTitle(title)
    if x_link is not None:
        plot.setXLink(x_link)
    if parent_widget is not None:
        _add_jpeg_export(plot, parent_widget, jpeg_default_name)
    return plot


# ===========================================================================
# JPEG export
# ===========================================================================

def _add_jpeg_export(plot: pg.PlotItem, parent_widget, default_name: str):
    vb = plot.getViewBox()
    vb.menu.addSeparator()
    action = vb.menu.addAction("Save graph as JPEG…")
    action.triggered.connect(
        lambda checked=False, p=plot, pw=parent_widget, n=default_name:
            _save_plot_as_jpeg(p, pw, n)
    )


def _save_plot_as_jpeg(plot: pg.PlotItem, parent, default_name: str):
    from pyqtgraph.exporters import ImageExporter
    default_path = str(Path.home() / f'{default_name}.jpg')
    file_path, _ = QFileDialog.getSaveFileName(
        parent, 'Save Graph as JPEG', default_path,
        'JPEG Images (*.jpg *.jpeg);;All Files (*)',
    )
    if not file_path:
        return
    try:
        exporter = ImageExporter(plot)
        exporter.parameters()['width'] = 1600
        exporter.export(file_path)
    except Exception as exc:
        QMessageBox.warning(parent, 'Export Failed',
                            f'Could not save image:\n{exc}')


# ===========================================================================
# Data plotting
# ===========================================================================

def plot_iq_data(
    plot: pg.PlotItem,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray | None = None,
    label: str = 'Data',
    brush=None,
    symbol_size: int | None = None,
) -> tuple:
    """Add a data scatter and optional error bars to *plot*.

    Parameters
    ----------
    plot : pg.PlotItem
        Target plot with log mode active (via ``make_sas_plot``).
    q, I : array
        Q and intensity in physical (linear) units.
    dI : array or None
        Uncertainties.  Error bars are drawn when provided.
    label : str
    brush : pg.QtGui.QBrush or None
        Symbol fill brush; defaults to ``SASPlotStyle.DATA_BRUSH``.
    symbol_size : int or None
        Symbol diameter in pixels; defaults to ``SASPlotStyle.DATA_SIZE``.

    Returns
    -------
    (scatter_item, error_item)
    """
    q = np.asarray(q, dtype=float)
    I = np.asarray(I, dtype=float)
    mask = np.isfinite(q) & np.isfinite(I) & (q > 0) & (I > 0)
    q_, I_ = q[mask], I[mask]

    _brush = brush if brush is not None else SASPlotStyle.DATA_BRUSH
    _size  = symbol_size if symbol_size is not None else SASPlotStyle.DATA_SIZE

    scatter = plot.plot(
        q_, I_,
        pen=None,
        symbol='o', symbolSize=_size,
        symbolBrush=_brush, symbolPen=pg.mkPen(None),
        name=label,
    )

    error_item = None
    if dI is not None:
        dI_ = np.asarray(dI, dtype=float)
        if dI_.shape == q.shape:
            dI_ = dI_[mask]
        elif dI_.shape == q_.shape:
            pass
        else:
            dI_ = np.zeros_like(q_)
        valid_I = I_[I_ > 0]
        if len(valid_I) >= 5:
            y_global_max = 10.0 ** (float(np.percentile(np.log10(valid_I), 99)) + 3)
        else:
            y_global_max = None
        error_item = draw_error_bars(plot, q_, I_, dI_, y_global_max=y_global_max)

    set_robust_y_range(plot, I_)

    valid_q = q_[q_ > 0]
    if len(valid_q) >= 2:
        q_lo = int(np.floor(np.log10(float(valid_q.min())))) - 1
        q_hi = int(np.ceil(np.log10(float(valid_q.max())))) + 1
        plot.getViewBox().setLimits(xMin=q_lo, xMax=q_hi)
        plot.setXRange(
            np.log10(float(valid_q.min())),
            np.log10(float(valid_q.max())),
            padding=0.05,
        )

    return scatter, error_item


def set_robust_y_range(plot: pg.PlotItem, I: np.ndarray) -> None:
    """Set Y axis range from percentile-based bounds of *I*.

    Call this after adding error bars so that pyqtgraph's auto-range
    algorithm (which includes error bar segments) does not drive the axis
    to the full double-precision range.
    """
    valid = (np.asarray(I) > 0) & np.isfinite(I)
    if np.sum(valid) < 3:
        return
    log_i = np.log10(np.asarray(I)[valid])
    lo = np.percentile(log_i, 2) - 0.5
    hi = np.percentile(log_i, 99) + 0.5
    plot.setYRange(lo, hi, padding=0)
    plot.getViewBox().setLimits(yMin=lo - 3, yMax=hi + 3)


def draw_error_bars(
    plot: pg.PlotItem,
    q: np.ndarray,
    I: np.ndarray,
    dI: np.ndarray,
    y_global_max: float | None = None,
) -> pg.PlotDataItem | None:
    """Draw I(Q) error bars as NaN-separated line segments.

    Each bar: vertical line from I-dI to I+dI, top cap, bottom cap.
    Caps are ±5 % of Q in log space.  Lower bound clipped to stay positive.
    Upper bound clipped to 3 decades above local I and optionally to
    *y_global_max* (99th-percentile + 3 decades of the dataset).

    Python lists (not numpy arrays) are passed to ``plot.plot()`` to match
    the pattern that avoids edge cases in pyqtgraph's bounds computation
    with NaN-containing arrays under log mode.
    """
    cap = SASPlotStyle.ERROR_CAP_FRAC
    x_lines: list[float] = []
    y_lines: list[float] = []

    for qi, Ii, dIi in zip(q, I, dI):
        if not (np.isfinite(dIi) and dIi > 0):
            continue
        y_top = min(Ii + dIi, Ii * 1000)
        if y_global_max is not None:
            y_top = min(y_top, y_global_max)
        y_bot = max(Ii - dIi, Ii * 0.001)

        x_lines.extend([qi, qi, np.nan])
        y_lines.extend([y_bot, y_top, np.nan])

        xl, xr = qi / (1.0 + cap), qi * (1.0 + cap)
        x_lines.extend([xl, xr, np.nan])
        y_lines.extend([y_top, y_top, np.nan])
        x_lines.extend([xl, xr, np.nan])
        y_lines.extend([y_bot, y_bot, np.nan])

    if not x_lines:
        return None

    return plot.plot(
        x_lines, y_lines,
        pen=SASPlotStyle.ERROR_PEN,
        connect='finite',
    )


def plot_iq_model(
    plot: pg.PlotItem,
    q: np.ndarray,
    I_model: np.ndarray,
    label: str = 'Model',
) -> pg.PlotDataItem | None:
    """Overlay a model/fit curve on *plot*."""
    q = np.asarray(q, dtype=float)
    I_model = np.asarray(I_model, dtype=float)
    mask = np.isfinite(q) & np.isfinite(I_model) & (q > 0) & (I_model > 0)
    if mask.sum() < 2:
        return None
    return plot.plot(q[mask], I_model[mask],
                     pen=SASPlotStyle.FIT_PEN, name=label)


# ===========================================================================
# Cursor helpers
# ===========================================================================

def make_cursors(
    plot: pg.PlotItem,
    q_min: float,
    q_max: float,
) -> tuple[_SafeInfiniteLine, _SafeInfiniteLine]:
    """Create two movable vertical cursors on *plot* in log10 space."""
    log_min = np.log10(max(float(q_min), 1e-20))
    log_max = np.log10(max(float(q_max), 1e-20))
    span    = log_max - log_min

    cursor_a = _SafeInfiniteLine(
        pos=log_min + 0.1 * span, angle=90, movable=True,
        pen=SASPlotStyle.CURSOR_A_PEN,
        label='A',
        labelOpts={'position': 0.05, 'color': SASPlotStyle.CURSOR_A_COLOR},
    )
    cursor_b = _SafeInfiniteLine(
        pos=log_max - 0.1 * span, angle=90, movable=True,
        pen=SASPlotStyle.CURSOR_B_PEN,
        label='B',
        labelOpts={'position': 0.10, 'color': SASPlotStyle.CURSOR_B_COLOR},
    )
    plot.addItem(cursor_a, ignoreBounds=True)
    plot.addItem(cursor_b, ignoreBounds=True)
    return cursor_a, cursor_b


def get_cursor_q_range(
    cursor_a: _SafeInfiniteLine | None,
    cursor_b: _SafeInfiniteLine | None,
) -> tuple[float | None, float | None]:
    """Return ``(q_min, q_max)`` in linear units from cursor positions."""
    if cursor_a is None or cursor_b is None:
        return None, None
    a = 10.0 ** cursor_a.getPos()[0]
    b = 10.0 ** cursor_b.getPos()[0]
    return (min(a, b), max(a, b))


def set_cursor_q_range(
    plot: pg.PlotItem,
    cursor_a: _SafeInfiniteLine | None,
    cursor_b: _SafeInfiniteLine | None,
    q_min: float,
    q_max: float,
) -> tuple[_SafeInfiniteLine, _SafeInfiniteLine]:
    """Position cursors at *q_min* and *q_max* (linear units)."""
    if q_min is not None and q_max is not None and q_min > 0 and q_max > 0:
        if cursor_a is None or cursor_b is None:
            cursor_a, cursor_b = make_cursors(plot, q_min, q_max)
        else:
            cursor_a.setPos(np.log10(q_min))
            cursor_b.setPos(np.log10(q_max))
    return cursor_a, cursor_b


# ===========================================================================
# Annotation helpers
# ===========================================================================

def add_plot_annotation(
    plot: pg.PlotItem,
    text: str,
    corner: str = 'lower_left',
) -> pg.TextItem:
    """Place a text annotation in a corner of *plot*'s visible area."""
    vr = plot.viewRange()
    dx = vr[0][1] - vr[0][0]
    dy = vr[1][1] - vr[1][0]
    margin_x = 0.02 * dx
    margin_y = 0.03 * dy

    if corner == 'lower_left':
        x = vr[0][0] + margin_x
        y = vr[1][0] + margin_y
        anchor = (0, 1)
    elif corner == 'upper_left':
        x = vr[0][0] + margin_x
        y = vr[1][1] - margin_y
        anchor = (0, 0)
    elif corner == 'lower_right':
        x = vr[0][1] - margin_x
        y = vr[1][0] + margin_y
        anchor = (1, 1)
    else:
        x = vr[0][1] - margin_x
        y = vr[1][1] - margin_y
        anchor = (1, 0)

    item = pg.TextItem(text=text, color=SASPlotStyle.ANNOT_COLOR, anchor=anchor)
    item.setPos(x, y)
    plot.addItem(item)
    return item
