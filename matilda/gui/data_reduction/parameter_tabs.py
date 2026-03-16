"""
matilda.gui.data_reduction.parameter_tabs
==========================================
Center panel: QTabWidget with 4 technique-specific parameter tabs.

Each tab exposes reduction parameters and a "Process Selected" button.
Tabs are always visible; technique detection auto-activates the matching tab.

Thickness handling
------------------
Each tab shows the HDF5 thickness by default (read by main_window.py on
file selection and pushed via ``update_hdf5_thickness``).  When the
"Override" checkbox is checked the spinbox becomes enabled and the manually
entered value is used instead.

Adding new parameters later
---------------------------
1. Add widget to the relevant tab's ``__init__``.
2. Include the value in ``get_params()``.
3. Pass the value through ReductionWorker to the converter function.
"""

try:
    from PySide6.QtWidgets import (
        QWidget, QTabWidget, QVBoxLayout, QFormLayout, QHBoxLayout,
        QLabel, QSpinBox, QDoubleSpinBox, QComboBox, QPushButton,
        QFrame, QCheckBox,
    )
    from PySide6.QtCore import Signal
except ImportError:
    from PyQt6.QtWidgets import (
        QWidget, QTabWidget, QVBoxLayout, QFormLayout, QHBoxLayout,
        QLabel, QSpinBox, QDoubleSpinBox, QComboBox, QPushButton,
        QFrame, QCheckBox,
    )
    from PyQt6.QtCore import pyqtSignal as Signal


EXTRAP_METHODS = [
    "PowerLaw w flat",
    "PowerLaw",
    "Flat",
    "Linear",
]

# Maps technique name → tab index (must stay in sync with addTab order)
TECHNIQUE_TAB_IDX: dict[str, int] = {
    "Flyscan":  0,
    "StepScan": 1,
    "SAXS":     2,
    "WAXS":     3,
}


def _separator() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line


class ParameterTabWidget(QTabWidget):
    """4-tab parameter widget: Flyscan | StepScan | SAXS | WAXS."""

    # technique name of the tab whose "Process Selected" button was clicked
    process_selected_clicked = Signal(str)
    # technique name of the tab whose "Browse…" blank button was clicked
    blank_browse_clicked = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._tabs: dict[str, _TechniqueTab] = {}

        for technique in TECHNIQUE_TAB_IDX:
            tab = _make_tab(technique, self)
            tab.process_selected_clicked.connect(
                lambda checked=False, t=technique: self.process_selected_clicked.emit(t)
            )
            tab.blank_browse_clicked.connect(
                lambda checked=False, t=technique: self.blank_browse_clicked.emit(t)
            )
            self.addTab(tab, technique)
            self._tabs[technique] = tab

    # ── Public API ────────────────────────────────────────────────────────────

    def activate_technique(self, technique: str):
        """Switch to the tab matching *technique* (called on file selection)."""
        idx = TECHNIQUE_TAB_IDX.get(technique)
        if idx is not None:
            self.setCurrentIndex(idx)

    def set_blank(self, technique: str, label: str):
        """Update the blank display label on the matching tab."""
        tab = self._tabs.get(technique)
        if tab is not None:
            tab.set_blank_label(label)

    def update_hdf5_thickness(self, technique: str, value: float | None):
        """Push the HDF5 thickness value to the matching tab."""
        tab = self._tabs.get(technique)
        if tab is not None:
            tab.update_hdf5_thickness(value)

    def get_params(self, technique: str) -> dict:
        """Return reduction parameters for *technique* as a plain dict."""
        tab = self._tabs.get(technique)
        return tab.get_params() if tab is not None else {}

    def get_all_params(self) -> dict[str, dict]:
        """Return params for all techniques keyed by technique name."""
        return {t: tab.get_params() for t, tab in self._tabs.items()}


# ── Base class ────────────────────────────────────────────────────────────────

class _TechniqueTab(QWidget):
    process_selected_clicked = Signal()
    blank_browse_clicked = Signal()

    # Shared state for HDF5 thickness
    _hdf5_thickness: float | None = None

    def get_params(self) -> dict:           # pragma: no cover
        raise NotImplementedError

    def set_blank_label(self, label: str):  # pragma: no cover
        raise NotImplementedError

    def update_hdf5_thickness(self, value: float | None):
        """Update the HDF5 default thickness display."""
        self._hdf5_thickness = value
        if value is not None:
            self._thickness_hdf5_lbl.setText(f"HDF5: {value:.4f} mm")
        else:
            self._thickness_hdf5_lbl.setText("HDF5: (not found)")

    def _get_thickness(self) -> float:
        """Return the effective thickness value."""
        if self._thickness_override.isChecked():
            return self._thickness_spin.value()
        if self._hdf5_thickness is not None:
            return self._hdf5_thickness
        return 1.0   # safe fallback

    def _setup_thickness_section(self, form: QFormLayout):
        """Add the thickness row (HDF5 label + override checkbox + spinbox)."""
        self._hdf5_thickness = None

        thickness_outer = QWidget()
        t_layout = QVBoxLayout(thickness_outer)
        t_layout.setContentsMargins(0, 0, 0, 0)
        t_layout.setSpacing(2)

        # Row 1: HDF5 value label + override toggle
        hdf5_row = QHBoxLayout()
        self._thickness_hdf5_lbl = QLabel("HDF5: (not loaded)")
        self._thickness_hdf5_lbl.setStyleSheet("color: grey; font-size: 11px;")
        hdf5_row.addWidget(self._thickness_hdf5_lbl, 1)

        self._thickness_override = QCheckBox("Override")
        self._thickness_override.setToolTip(
            "Check to enter a custom thickness instead of the HDF5 value"
        )
        hdf5_row.addWidget(self._thickness_override)
        t_layout.addLayout(hdf5_row)

        # Row 2: spinbox (disabled until Override is checked)
        self._thickness_spin = QDoubleSpinBox()
        self._thickness_spin.setRange(0.001, 100.0)
        self._thickness_spin.setValue(1.0)
        self._thickness_spin.setDecimals(3)
        self._thickness_spin.setSuffix("  mm")
        self._thickness_spin.setEnabled(False)
        self._thickness_spin.setToolTip("Custom sample thickness in mm")
        self._thickness_override.toggled.connect(self._thickness_spin.setEnabled)
        t_layout.addWidget(self._thickness_spin)

        form.addRow("Thickness:", thickness_outer)


def _make_tab(technique: str, parent=None) -> _TechniqueTab:
    if technique in ("Flyscan", "StepScan"):
        return _USAXSTab(technique, parent)
    elif technique == "SAXS":
        return _SAXSTab(parent)
    elif technique == "WAXS":
        return _WAXSTab(parent)
    return _TechniqueTab(parent)


# ── USAXS tab (shared by Flyscan and StepScan) ────────────────────────────────

class _USAXSTab(_TechniqueTab):
    def __init__(self, technique: str, parent=None):
        super().__init__(parent)
        self._technique = technique

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        # ── Blank ─────────────────────────────────────────────────────────────
        blank_row = QHBoxLayout()
        self._blank_mode = QComboBox()
        self._blank_mode.addItems(["auto (nearest preceding)", "manual"])
        self._blank_mode.setToolTip(
            "auto: use nearest-preceding blank scan (same as matilda daemon)\n"
            "manual: use the file designated via right-click in the file tree"
        )
        blank_row.addWidget(self._blank_mode, 1)
        self._blank_browse = QPushButton("Browse…")
        self._blank_browse.setToolTip("Browse for a blank file")
        self._blank_browse.clicked.connect(self.blank_browse_clicked)
        blank_row.addWidget(self._blank_browse)
        form.addRow("Blank:", blank_row)

        self._blank_label = QLabel("(none assigned)")
        self._blank_label.setStyleSheet("color: grey; font-style: italic;")
        form.addRow("", self._blank_label)

        form.addRow(_separator())

        # ── Sample thickness (HDF5 default + override) ────────────────────────
        self._setup_thickness_section(form)

        # ── Output points ─────────────────────────────────────────────────────
        self._npts = QSpinBox()
        self._npts.setRange(0, 5000)
        self._npts.setValue(500)
        self._npts.setToolTip(
            "Number of Q points after rebinning.\n0 = no rebinning (keep all raw points)."
        )
        form.addRow("Output points:", self._npts)

        form.addRow(_separator())
        form.addRow(QLabel("<b>Desmearing</b>"))

        # ── Max iterations ────────────────────────────────────────────────────
        self._desmear_iter = QSpinBox()
        self._desmear_iter.setRange(1, 500)
        self._desmear_iter.setValue(20)
        self._desmear_iter.setToolTip("Maximum Lake/Strobl desmearing iterations")
        form.addRow("Max iterations:", self._desmear_iter)

        # ── Extrapolation method ──────────────────────────────────────────────
        self._extrap_method = QComboBox()
        self._extrap_method.addItems(EXTRAP_METHODS)
        self._extrap_method.setToolTip("High-Q extrapolation method used during desmearing")
        form.addRow("Extrap method:", self._extrap_method)

        # ── Extrap Q start ────────────────────────────────────────────────────
        self._extrap_qstart = QDoubleSpinBox()
        self._extrap_qstart.setRange(0.0001, 10.0)
        self._extrap_qstart.setValue(0.1)
        self._extrap_qstart.setDecimals(4)
        self._extrap_qstart.setSuffix("  Å⁻¹")
        self._extrap_qstart.setToolTip("Q value above which the extrapolation is applied")
        form.addRow("Extrap Q start:", self._extrap_qstart)

        layout.addLayout(form)
        layout.addStretch()

        layout.addWidget(_separator())
        self._btn_process = QPushButton(f"▶  Process Selected  ({technique})")
        self._btn_process.setToolTip(f"Process all selected files as {technique}")
        self._btn_process.clicked.connect(self.process_selected_clicked)
        layout.addWidget(self._btn_process)

    def set_blank_label(self, label: str):
        if label:
            self._blank_label.setText(label)
            self._blank_label.setStyleSheet("color: darkblue;")
        else:
            self._blank_label.setText("(none assigned)")
            self._blank_label.setStyleSheet("color: grey; font-style: italic;")

    def get_params(self) -> dict:
        return {
            "blank_mode":         self._blank_mode.currentText(),
            "thickness":          self._get_thickness(),
            "npts":               self._npts.value(),
            "desmear_iter":       self._desmear_iter.value(),
            "extrap_method":      self._extrap_method.currentText(),
            "extrap_qstart":      self._extrap_qstart.value(),
            "recalculateAllData": True,
        }


# ── SAXS tab ──────────────────────────────────────────────────────────────────

class _SAXSTab(_TechniqueTab):
    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        # ── Blank ─────────────────────────────────────────────────────────────
        blank_row = QHBoxLayout()
        self._blank_mode = QComboBox()
        self._blank_mode.addItems(["auto (nearest preceding)", "manual"])
        blank_row.addWidget(self._blank_mode, 1)
        self._blank_browse = QPushButton("Browse…")
        self._blank_browse.clicked.connect(self.blank_browse_clicked)
        blank_row.addWidget(self._blank_browse)
        form.addRow("Blank:", blank_row)

        self._blank_label = QLabel("(none assigned)")
        self._blank_label.setStyleSheet("color: grey; font-style: italic;")
        form.addRow("", self._blank_label)

        form.addRow(_separator())

        # ── Sample thickness (HDF5 default + override) ────────────────────────
        self._setup_thickness_section(form)

        # ── Output Q points ───────────────────────────────────────────────────
        self._npts = QSpinBox()
        self._npts.setRange(10, 5000)
        self._npts.setValue(200)
        self._npts.setToolTip("Number of Q points in azimuthal integration (pyFAI npt)")
        form.addRow("Output Q points:", self._npts)

        form.addRow(_separator())
        form.addRow(QLabel("<b>Azimuthal Integration</b>"))

        # ── Azimuth range ─────────────────────────────────────────────────────
        az_row = QHBoxLayout()
        self._az_min = QDoubleSpinBox()
        self._az_min.setRange(-180.0, 180.0)
        self._az_min.setValue(-30.0)
        self._az_min.setSuffix("°")
        az_row.addWidget(self._az_min)
        az_row.addWidget(QLabel("to"))
        self._az_max = QDoubleSpinBox()
        self._az_max.setRange(-180.0, 180.0)
        self._az_max.setValue(30.0)
        self._az_max.setSuffix("°")
        az_row.addWidget(self._az_max)
        form.addRow("Azimuth range:", az_row)

        layout.addLayout(form)
        layout.addStretch()

        layout.addWidget(_separator())
        self._btn_process = QPushButton("▶  Process Selected  (SAXS)")
        self._btn_process.clicked.connect(self.process_selected_clicked)
        layout.addWidget(self._btn_process)

    def set_blank_label(self, label: str):
        if label:
            self._blank_label.setText(label)
            self._blank_label.setStyleSheet("color: darkblue;")
        else:
            self._blank_label.setText("(none assigned)")
            self._blank_label.setStyleSheet("color: grey; font-style: italic;")

    def get_params(self) -> dict:
        return {
            "blank_mode":         self._blank_mode.currentText(),
            "thickness":          self._get_thickness(),
            "npts":               self._npts.value(),
            "az_min":             self._az_min.value(),
            "az_max":             self._az_max.value(),
            "recalculateAllData": True,
        }


# ── WAXS tab ──────────────────────────────────────────────────────────────────

class _WAXSTab(_TechniqueTab):
    def __init__(self, parent=None):
        super().__init__(parent)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        # ── Blank ─────────────────────────────────────────────────────────────
        blank_row = QHBoxLayout()
        self._blank_mode = QComboBox()
        self._blank_mode.addItems(["auto (nearest preceding)", "manual"])
        blank_row.addWidget(self._blank_mode, 1)
        self._blank_browse = QPushButton("Browse…")
        self._blank_browse.clicked.connect(self.blank_browse_clicked)
        blank_row.addWidget(self._blank_browse)
        form.addRow("Blank:", blank_row)

        self._blank_label = QLabel("(none assigned)")
        self._blank_label.setStyleSheet("color: grey; font-style: italic;")
        form.addRow("", self._blank_label)

        form.addRow(_separator())

        # ── Sample thickness (HDF5 default + override) ────────────────────────
        self._setup_thickness_section(form)

        # Placeholder for future parameters
        form.addRow(_separator())
        placeholder = QLabel(
            "<i>Additional parameters (pixel mask, normalisation factors)\n"
            "will be added here in a future version.</i>"
        )
        placeholder.setStyleSheet("color: grey;")
        form.addRow(placeholder)

        layout.addLayout(form)
        layout.addStretch()

        layout.addWidget(_separator())
        self._btn_process = QPushButton("▶  Process Selected  (WAXS)")
        self._btn_process.clicked.connect(self.process_selected_clicked)
        layout.addWidget(self._btn_process)

    def set_blank_label(self, label: str):
        if label:
            self._blank_label.setText(label)
            self._blank_label.setStyleSheet("color: darkblue;")
        else:
            self._blank_label.setText("(none assigned)")
            self._blank_label.setStyleSheet("color: grey; font-style: italic;")

    def get_params(self) -> dict:
        return {
            "blank_mode":         self._blank_mode.currentText(),
            "thickness":          self._get_thickness(),
            "recalculateAllData": True,
        }
