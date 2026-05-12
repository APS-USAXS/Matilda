"""
matilda.gui.data_reduction.parameter_tabs
==========================================
Center panel: QTabWidget with 4 technique-specific parameter tabs.

Each tab exposes reduction parameters and a "Process Selected" button.
Tabs are always visible; technique detection auto-activates the matching tab.

Thickness handling
------------------
Each tab shows the HDF5 thickness of the last single-clicked file (read
by main_window.py and pushed via ``update_hdf5_thickness``).  This label
is for display only.  ``get_params()`` returns the thickness override
ONLY when the "Override" checkbox is checked — otherwise it returns
None so the converter reads each file's own /entry/sample/thickness
during batch processing.  Same pattern is used for transmission and Qmin
overrides.

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

    def set_blank_mode(self, technique: str, mode: str):
        """Switch the blank-mode dropdown on the matching tab.

        *mode* should be one of the combo-box items, e.g. ``"manual"``
        or ``"auto (nearest preceding)"``.
        """
        tab = self._tabs.get(technique)
        if tab is not None and hasattr(tab, "_blank_mode"):
            idx = tab._blank_mode.findText(mode)
            if idx >= 0:
                tab._blank_mode.setCurrentIndex(idx)

    def update_hdf5_thickness(self, technique: str, value: float | None):
        """Push the HDF5 thickness value to the matching tab."""
        tab = self._tabs.get(technique)
        if tab is not None:
            tab.update_hdf5_thickness(value)

    def update_mu_thickness(self, technique: str, thickness_mm: float):
        """Show the calculated thickness on the matching tab's μ label."""
        tab = self._tabs.get(technique)
        if tab is not None and hasattr(tab, "_mu_thickness_lbl"):
            if tab._use_mu.isChecked():
                tab._mu_thickness_lbl.setText(f"t = {thickness_mm:.4f} mm")

    def update_calculated_mu(self, technique: str, mu: float | None):
        """Set the μ field on the matching tab to the calculated μ value."""
        tab = self._tabs.get(technique)
        if tab is not None and hasattr(tab, "update_calculated_mu"):
            tab.update_calculated_mu(mu)

    def update_measured_transmission(self, technique: str, t: float | None):
        """Update the 'used: T' label on the matching tab."""
        tab = self._tabs.get(technique)
        if tab is not None and hasattr(tab, "update_measured_transmission"):
            tab.update_measured_transmission(t)

    def update_calculated_qmin(self, technique: str, qmin: float | None):
        """Update the 'calc'd: Qmin' label on the matching tab (USAXS only)."""
        tab = self._tabs.get(technique)
        if tab is not None and hasattr(tab, "update_calculated_qmin"):
            tab.update_calculated_qmin(qmin)

    def reset_per_file_state(self, technique: str):
        """Reset per-file display state when a new file is selected."""
        tab = self._tabs.get(technique)
        if tab is not None and hasattr(tab, "reset_per_file_state"):
            tab.reset_per_file_state()

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

    def _get_thickness(self) -> float | None:
        """Return the thickness override to pass to the converter.

        Returns None when the user has NOT enabled override — this signals
        the converter to read each file's own /entry/sample/thickness.
        Returning the cached _hdf5_thickness here would force every file
        in a batch to use the thickness of the last single-clicked file
        (or zero, if that file had no thickness recorded), causing
        divide-by-zero downstream in calibrateAD2DData / calibrateAndSubtractFlyscan.
        The cached _hdf5_thickness is used only for the on-screen label.
        """
        if self._thickness_override.isChecked():
            return self._thickness_spin.value()
        return None

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
        self._thickness_spin.setRange(0.0, 50.0)
        self._thickness_spin.setValue(1.0)
        self._thickness_spin.setDecimals(3)
        self._thickness_spin.setSingleStep(0.02)
        self._thickness_spin.setSuffix("  mm")
        self._thickness_spin.setEnabled(False)
        self._thickness_spin.setToolTip("Custom sample thickness in mm")
        self._thickness_override.toggled.connect(self._thickness_spin.setEnabled)
        t_layout.addWidget(self._thickness_spin)

        form.addRow("Thickness:", thickness_outer)

    def _setup_calibration_section(self, form: QFormLayout):
        """Add μ-based thickness and per-gram calibration controls."""
        self._use_mu = QCheckBox("Use μ for thickness")
        self._use_mu.setToolTip(
            "Calculate thickness from measured transmission and\n"
            "linear absorption coefficient μ:  t = −ln(T) / μ"
        )
        form.addRow("", self._use_mu)

        # μ input + calculated thickness readout
        mu_row = QHBoxLayout()
        self._mu_spin = QDoubleSpinBox()
        self._mu_spin.setRange(0.001, 1000.0)
        self._mu_spin.setValue(12.5)
        self._mu_spin.setDecimals(3)
        self._mu_spin.setSuffix("  1/cm")
        self._mu_spin.setEnabled(False)
        self._mu_spin.setToolTip(
            "Linear absorption coefficient from Scattering Contrast Calculator.\n"
            "Auto-populated after processing from T and thickness for reference."
        )
        mu_row.addWidget(self._mu_spin, 1)
        self._mu_thickness_lbl = QLabel("")
        self._mu_thickness_lbl.setStyleSheet("color: grey; font-size: 11px;")
        mu_row.addWidget(self._mu_thickness_lbl)
        form.addRow("μ (1/cm):", mu_row)

        self._per_gram = QCheckBox("Normalize per gram")
        self._per_gram.setToolTip(
            "Divide intensity by solid-frame density to get [cm²/g].\n"
            "Used for powder samples where thickness is not meaningful."
        )
        self._per_gram.setEnabled(False)
        form.addRow("", self._per_gram)

        self._density_spin = QDoubleSpinBox()
        self._density_spin.setRange(0.0, 30.0)
        self._density_spin.setValue(2.2)
        self._density_spin.setDecimals(3)
        self._density_spin.setSingleStep(0.02)
        self._density_spin.setSuffix("  g/cm³")
        self._density_spin.setEnabled(False)
        self._density_spin.setToolTip("Solid-frame density of the sample material")
        form.addRow("Density:", self._density_spin)

        # ── Transmission override ─────────────────────────────────────────
        self._override_transmission = QCheckBox("Override transmission")
        self._override_transmission.setToolTip(
            "Replace the measured transmission with a manual value.\n"
            "Rarely needed — useful when the diode-measured transmission is bad."
        )
        form.addRow("", self._override_transmission)

        trans_row = QHBoxLayout()
        self._transmission_spin = QDoubleSpinBox()
        self._transmission_spin.setRange(0.0, float('inf'))
        self._transmission_spin.setValue(0.5)
        self._transmission_spin.setDecimals(4)
        self._transmission_spin.setSingleStep(0.05)
        self._transmission_spin.setEnabled(False)
        self._transmission_spin.setToolTip("Manual transmission value (0–1, dimensionless; rarely >1 due to amplifier failure)")
        trans_row.addWidget(self._transmission_spin, 1)
        self._measured_t_lbl = QLabel("")
        self._measured_t_lbl.setStyleSheet("color: grey; font-size: 11px;")
        trans_row.addWidget(self._measured_t_lbl)
        form.addRow("Transmission:", trans_row)

        form.addRow(_separator())

        # Wiring: use_mu toggles mu_spin and per_gram availability
        def _on_use_mu_toggled(checked):
            self._mu_spin.setEnabled(checked)
            self._per_gram.setEnabled(checked)
            if not checked:
                self._per_gram.setChecked(False)
                self._mu_thickness_lbl.setText("")
            # Disable thickness override when using μ
            self._thickness_override.setEnabled(not checked)
            if checked:
                self._thickness_override.setChecked(False)
                self._update_mu_thickness()

        self._use_mu.toggled.connect(_on_use_mu_toggled)
        self._per_gram.toggled.connect(self._density_spin.setEnabled)
        self._mu_spin.valueChanged.connect(lambda _: self._update_mu_thickness())
        self._override_transmission.toggled.connect(self._transmission_spin.setEnabled)

    def _update_mu_thickness(self):
        """Reset the thickness label when μ value changes (actual value shown after processing)."""
        if not self._use_mu.isChecked():
            self._mu_thickness_lbl.setText("")
            return
        self._mu_thickness_lbl.setText("t = −ln(T)/μ  (process to calculate)")

    def _update_qmin_step(self):
        """Set Qmin step to 10% of current value."""
        current = self._qmin_spin.value()
        if current > 0:
            self._qmin_spin.setSingleStep(max(current * 0.1, 1e-6))

    def update_calculated_mu(self, mu: float | None):
        """Show the μ value calculated from T and thickness during processing.

        Only seeds the μ spinbox when μ-based mode is OFF (so the user's
        entered μ in μ-based mode is not silently replaced by the value
        the converter computed from it).
        """
        if mu is None or not (0 < mu < 1e6):
            return
        if self._use_mu.isChecked():
            return  # preserve user-entered μ in μ-based mode
        try:
            self._mu_spin.blockSignals(True)
            self._mu_spin.setValue(max(self._mu_spin.minimum(),
                                        min(self._mu_spin.maximum(), float(mu))))
        finally:
            self._mu_spin.blockSignals(False)

    def update_measured_transmission(self, t: float | None):
        """Show the measured (or used) transmission next to the override field."""
        if t is None:
            self._measured_t_lbl.setText("")
            return
        self._measured_t_lbl.setText(f"used: {t:.4f}")

    def reset_per_file_state(self):
        """Reset per-file display state (called when a new file is selected)."""
        self._measured_t_lbl.setText("")

    def _get_calibration_params(self) -> dict:
        """Return calibration-mode parameters."""
        return {
            "use_mu":              self._use_mu.isChecked(),
            "mu":                  self._mu_spin.value(),
            "per_gram":            self._per_gram.isChecked(),
            "density":             self._density_spin.value(),
            "transmission_override": (
                self._transmission_spin.value()
                if self._override_transmission.isChecked() else None
            ),
        }


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
        self._setup_calibration_section(form)

        # ── Min Q ratio (MinQMinFindRatio) ───────────────────────────────────
        self._min_q_ratio = QDoubleSpinBox()
        self._min_q_ratio.setRange(0.9, 10.0)
        self._min_q_ratio.setValue(1.05)
        self._min_q_ratio.setDecimals(2)
        self._min_q_ratio.setSingleStep(0.01)
        self._min_q_ratio.setToolTip(
            "Threshold for Q-minimum selection after blank subtraction.\n"
            "The first Q point where the sample/blank intensity ratio\n"
            "exceeds this value defines the low-Q cutoff."
        )
        form.addRow("Min Q ratio:", self._min_q_ratio)

        # ── Manual Qmin override (truncates data below this Q value) ─────────
        self._override_qmin = QCheckBox("Override Qmin (truncate from below)")
        self._override_qmin.setToolTip(
            "Emergency override of the auto-calculated Qmin.\n"
            "Useful when blank subtraction goes bad at low Q and you want\n"
            "to salvage data at higher Q values."
        )
        form.addRow("", self._override_qmin)

        qmin_row = QHBoxLayout()
        self._qmin_spin = QDoubleSpinBox()
        self._qmin_spin.setRange(1e-5, 1.0)
        self._qmin_spin.setValue(1e-4)
        self._qmin_spin.setDecimals(6)
        self._qmin_spin.setSuffix("  1/Å")
        self._qmin_spin.setEnabled(False)
        self._qmin_spin.setToolTip(
            "Manual Qmin: data points with Q below this value are removed.\n"
            "Step size is 10% of the current value."
        )
        qmin_row.addWidget(self._qmin_spin, 1)
        self._calculated_qmin_lbl = QLabel("")
        self._calculated_qmin_lbl.setStyleSheet("color: grey; font-size: 11px;")
        qmin_row.addWidget(self._calculated_qmin_lbl)
        form.addRow("Qmin:", qmin_row)

        self._override_qmin.toggled.connect(self._qmin_spin.setEnabled)
        self._qmin_spin.valueChanged.connect(self._update_qmin_step)

        # ── Output points (Flyscan only) ─────────────────────────────────────
        self._npts = None
        if technique == "Flyscan":
            self._npts = QSpinBox()
            self._npts.setRange(200, 8000)
            self._npts.setValue(500)
            self._npts.setSingleStep(100)
            self._npts.setToolTip("Number of output points after rebinning")
            form.addRow("Output points:", self._npts)

        form.addRow(_separator())
        form.addRow(QLabel("<b>Desmearing</b>"))

        # ── Max iterations ────────────────────────────────────────────────────
        self._desmear_iter = QSpinBox()
        self._desmear_iter.setRange(10, 100)
        self._desmear_iter.setValue(20)
        self._desmear_iter.setSingleStep(5)
        self._desmear_iter.setToolTip("Maximum Lake/Strobl desmearing iterations")
        form.addRow("Max iterations:", self._desmear_iter)

        # ── Extrapolation method ──────────────────────────────────────────────
        self._extrap_method = QComboBox()
        self._extrap_method.addItems(EXTRAP_METHODS)
        self._extrap_method.setToolTip("High-Q extrapolation method used during desmearing")
        form.addRow("Extrap method:", self._extrap_method)

        # ── Extrap Q start ────────────────────────────────────────────────────
        self._extrap_qstart = QDoubleSpinBox()
        self._extrap_qstart.setRange(0.01, 1.0)
        self._extrap_qstart.setValue(0.1)
        self._extrap_qstart.setDecimals(4)
        self._extrap_qstart.setSingleStep(0.05)
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
        params = {
            "blank_mode":         self._blank_mode.currentText(),
            "thickness":          self._get_thickness(),
            "minQMinFindRatio":   self._min_q_ratio.value(),
            "desmear_iter":       self._desmear_iter.value(),
            "extrap_method":      self._extrap_method.currentText(),
            "extrap_qstart":      self._extrap_qstart.value(),
            "qmin_override":      (
                self._qmin_spin.value()
                if self._override_qmin.isChecked() else None
            ),
            "recalculateAllData": True,
            **self._get_calibration_params(),
        }
        if self._npts is not None:
            params["npts"] = self._npts.value()
        return params

    def update_calculated_qmin(self, qmin: float | None):
        """Always show the auto-calculated Qmin; seed the spinbox only when
        override is OFF so a user-entered value is not silently replaced."""
        if qmin is None or qmin <= 0:
            self._calculated_qmin_lbl.setText("")
            return
        self._calculated_qmin_lbl.setText(f"calc'd: {qmin:.4e}")
        if self._override_qmin.isChecked():
            return  # preserve user-entered override value
        try:
            self._qmin_spin.blockSignals(True)
            self._qmin_spin.setValue(max(self._qmin_spin.minimum(),
                                         min(self._qmin_spin.maximum(), float(qmin))))
        finally:
            self._qmin_spin.blockSignals(False)

    def reset_per_file_state(self):
        """Reset per-file display state (called when a new file is selected).
        Qmin override is unchecked so the calculated Qmin from the new file
        is what gets used by default (Qmin varies a lot between samples)."""
        super().reset_per_file_state()
        self._calculated_qmin_lbl.setText("")
        self._override_qmin.setChecked(False)


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
        self._setup_calibration_section(form)

        # ── Output Q points ───────────────────────────────────────────────────
        self._max_pts = QCheckBox("Max number of points")
        self._max_pts.setToolTip(
            "When checked, use the maximum number of output points\n"
            "(determined by the detector dimensions)."
        )
        form.addRow("", self._max_pts)

        self._npts = QSpinBox()
        self._npts.setRange(200, 8000)
        self._npts.setValue(200)
        self._npts.setSingleStep(100)
        self._npts.setToolTip("Number of output Q points for azimuthal integration")
        form.addRow("Output Q points:", self._npts)

        self._max_pts.toggled.connect(lambda checked: self._npts.setEnabled(not checked))

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
            "npts":               None if self._max_pts.isChecked() else self._npts.value(),
            "az_min":             self._az_min.value(),
            "az_max":             self._az_max.value(),
            "recalculateAllData": True,
            **self._get_calibration_params(),
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
        self._setup_calibration_section(form)

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
            **self._get_calibration_params(),
        }
