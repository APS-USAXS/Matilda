"""
matilda.gui.data_reduction.main_window
=======================================
Top-level QMainWindow for the Matilda data-reduction GUI.

Layout (horizontal QSplitter):
  Left   — FileTreeWidget  (folder browser, blank assignment)
  Centre — ParameterTabWidget  (4 tabs: Flyscan | StepScan | SAXS | WAXS)
  Right  — GraphPanel  (pyqtgraph log-log curves)

Bottom bar:
  [▶▶ Process All]  status label  progress bar  [✕ Cancel]
  Error log (QPlainTextEdit — accumulates per-file errors; always visible)
"""

import json
import os
import webbrowser

import h5py
import numpy as np

from .._qt import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QPushButton, QLabel, QProgressBar, QSizePolicy,
    QFileDialog, QMessageBox, QPlainTextEdit, Qt,
)

from .file_tree import FileTreeWidget
from .parameter_tabs import ParameterTabWidget
from .graph_panel import GraphPanel
from .reduction_worker import ReductionWorker
from .ascii_exporter import AsciiExportWorker
from .technique_detector import detect_technique

_SESSION_FILE = os.path.expanduser("~/.matilda_gui_session.json")
_WINDOW_TITLE  = "Matilda — Data Reduction"


class MatildaReductionWindow(QMainWindow):
    """Main application window."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle(_WINDOW_TITLE)
        self.resize(1400, 860)
        self._worker: ReductionWorker | None = None
        self._export_worker: AsciiExportWorker | None = None
        self._last_folder: str = os.path.expanduser("~")
        self._build_ui()
        self._restore_session()

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        # ── Toolbar ───────────────────────────────────────────────────────────
        tb = self.addToolBar("Main")
        tb.setMovable(False)

        self._btn_folder = QPushButton("📁  Select Folder…")
        self._btn_folder.setToolTip("Choose a root folder containing HDF5 data files")
        self._btn_folder.clicked.connect(self._select_folder)
        tb.addWidget(self._btn_folder)

        self._folder_label = QLabel("  (no folder selected)")
        self._folder_label.setStyleSheet("color: grey;")
        tb.addWidget(self._folder_label)

        # Spacer to push help button to the right
        spacer = QWidget()
        spacer.setMinimumWidth(0)
        spacer.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        tb.addWidget(spacer)

        btn_help = QPushButton("? Help")
        btn_help.setStyleSheet(
            "QPushButton { background-color: #e84040; color: white;"
            " font-weight: bold; padding: 4px 12px; border-radius: 3px; }"
            "QPushButton:hover { background-color: #cc3030; }"
        )
        btn_help.setToolTip("Open Matilda GUI documentation in web browser")
        btn_help.clicked.connect(
            lambda: webbrowser.open(
                "https://github.com/APS-USAXS/Matilda/blob/main/docs/matilda-gui.md"
            )
        )
        tb.addWidget(btn_help)

        # ── Central area ──────────────────────────────────────────────────────
        central = QWidget()
        self.setCentralWidget(central)
        outer = QVBoxLayout(central)
        outer.setContentsMargins(4, 4, 4, 4)
        outer.setSpacing(4)

        # Three-panel splitter
        self._splitter = QSplitter(Qt.Orientation.Horizontal)

        self._file_tree = FileTreeWidget()
        self._file_tree.setMinimumWidth(200)
        self._splitter.addWidget(self._file_tree)

        self._param_tabs = ParameterTabWidget()
        self._param_tabs.setMinimumWidth(280)
        self._splitter.addWidget(self._param_tabs)

        self._graph = GraphPanel()
        self._graph.setMinimumWidth(360)
        self._splitter.addWidget(self._graph)

        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 0)
        self._splitter.setStretchFactor(2, 1)
        self._splitter.setSizes([260, 370, 770])

        outer.addWidget(self._splitter, 1)

        # ── Bottom bar ────────────────────────────────────────────────────────
        bottom = QWidget()
        bl = QHBoxLayout(bottom)
        bl.setContentsMargins(4, 2, 4, 2)

        self._btn_process_all = QPushButton("▶▶  Process All Files")
        self._btn_process_all.setToolTip(
            "Process every visible file in the tree using its auto-detected technique"
        )
        self._btn_process_all.clicked.connect(self._on_process_all)
        bl.addWidget(self._btn_process_all)

        self._btn_export_ascii = QPushButton("📄  Export ASCII from All")
        self._btn_export_ascii.setToolTip(
            "Export processed I(Q) data from all HDF5 files to ASCII .dat files.\n"
            "Files without NXcanSAS data are skipped."
        )
        self._btn_export_ascii.clicked.connect(self._on_export_ascii)
        bl.addWidget(self._btn_export_ascii)

        bl.addStretch()

        self._status_label = QLabel("Idle")
        bl.addWidget(self._status_label)

        self._progress = QProgressBar()
        self._progress.setVisible(False)
        self._progress.setMinimumWidth(200)
        self._progress.setTextVisible(True)
        bl.addWidget(self._progress)

        self._btn_cancel = QPushButton("✕  Cancel")
        self._btn_cancel.setEnabled(False)
        self._btn_cancel.clicked.connect(self._on_cancel)
        bl.addWidget(self._btn_cancel)

        outer.addWidget(bottom)

        # ── Error log ─────────────────────────────────────────────────────────
        # Always visible; accumulates per-file errors so they are not lost
        # when the status bar text is overwritten by the next file.
        self._error_log = QPlainTextEdit()
        self._error_log.setReadOnly(True)
        self._error_log.setMaximumBlockCount(200)
        self._error_log.setFixedHeight(70)
        self._error_log.setPlaceholderText("Errors will appear here…")
        self._error_log.setStyleSheet(
            "QPlainTextEdit { font-family: 'Courier New', Courier; font-size: 11px;"
            " background: #fff8f8; color: #800; border: 1px solid #dbb; }"
        )
        outer.addWidget(self._error_log)

        # ── Signal wiring ─────────────────────────────────────────────────────
        self._file_tree.selection_changed.connect(self._on_selection_changed)
        self._file_tree.file_double_clicked.connect(self._start_reduction)
        self._file_tree.blank_assigned.connect(self._on_blank_assigned)
        self._file_tree.blank_cleared.connect(self._on_blank_cleared)
        self._param_tabs.process_selected_clicked.connect(self._on_process_selected)
        self._param_tabs.blank_browse_clicked.connect(self._on_blank_browse)
        self._graph.view_2d_requested.connect(self._open_2d_viewer)

    # ── Folder selection ──────────────────────────────────────────────────────

    def _select_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select Data Folder", self._last_folder
        )
        if folder:
            self._last_folder = folder
            self._folder_label.setText(f"  {folder}")
            self._folder_label.setStyleSheet("")
            self._file_tree.set_folder(folder)
            self._status_label.setText("Folder loaded.")

    # ── Selection ─────────────────────────────────────────────────────────────

    def _on_selection_changed(self, selected: list[tuple[str, str]]):
        if not selected:
            self._status_label.setText("Idle")
            return
        path, fname = selected[0]
        technique = detect_technique(path, fname)
        if technique != "Unknown":
            self._param_tabs.activate_technique(technique)
            thickness = self._read_hdf5_thickness(path, fname, technique)
            self._param_tabs.update_hdf5_thickness(technique, thickness)
            # Clear per-file display state (measured T, calculated Qmin)
            self._param_tabs.reset_per_file_state(technique)
        n = len(selected)
        suffix = (f"  —  detected: {technique}") if n == 1 else ""
        self._status_label.setText(f"{n} file(s) selected{suffix}")

    # ── Blank assignment ──────────────────────────────────────────────────────

    def _on_blank_assigned(self, technique: str, path: str, fname: str):
        self._param_tabs.set_blank(technique, fname)
        self._param_tabs.set_blank_mode(technique, "manual")
        self._status_label.setText(f"Blank set for {technique}: {fname}")

    def _on_blank_cleared(self, technique: str):
        self._param_tabs.set_blank(technique, "")
        self._status_label.setText(f"Blank cleared for {technique}")

    def _on_blank_auto(self, technique: str, blank_filename: str):
        """Update the blank label when the worker auto-selects a blank."""
        self._param_tabs.set_blank(technique, f"auto: {blank_filename}")

    def _on_blank_browse(self, technique: str):
        """Open a file dialog to browse for a blank outside the current tree."""
        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Select blank file for {technique}",
            self._last_folder,
            "HDF5 Files (*.h5 *.hdf *.hdf5 *.nxs);;All Files (*)",
        )
        if path:
            folder = os.path.dirname(path)
            fname  = os.path.basename(path)
            self._file_tree.set_blank_external(technique, folder, fname)
            self._param_tabs.set_blank_mode(technique, "manual")

    # ── Processing ────────────────────────────────────────────────────────────

    def _on_process_selected(self, technique: str):
        files = self._file_tree.get_selected_files()
        if not files:
            QMessageBox.information(
                self, "No selection",
                "Please select one or more files in the file tree."
            )
            return
        self._start_reduction(files)

    def _on_process_all(self):
        files = self._file_tree.get_all_files()
        if not files:
            QMessageBox.information(
                self, "No files",
                "No files are visible in the file tree.\n"
                "Select a folder first (or clear the filter)."
            )
            return
        self._start_reduction(files)

    def _start_reduction(self, file_list: list[tuple[str, str]]):
        if self._worker and self._worker.isRunning():
            QMessageBox.warning(
                self, "Busy",
                "A reduction job is already running.\n"
                "Please wait or click Cancel."
            )
            return

        self._error_log.clear()
        blanks = self._file_tree.get_blanks()
        params = self._param_tabs.get_all_params()
        # Use unfiltered file list for blank detection (filter is only for UI display)
        all_files = self._file_tree.get_all_files_unfiltered()

        self._progress.setMaximum(len(file_list))
        self._progress.setValue(0)
        self._progress.setVisible(True)
        self._btn_cancel.setEnabled(True)
        self._btn_process_all.setEnabled(False)
        self._status_label.setText(f"Processing  0 / {len(file_list)}…")
        self._graph.clear()

        self._worker = ReductionWorker(file_list, blanks, params, all_files=all_files, parent=self)
        self._worker.progress.connect(self._on_progress)
        self._worker.file_done.connect(self._on_file_done)
        self._worker.file_error.connect(self._on_file_error)
        self._worker.log_message.connect(self._on_log_message)
        self._worker.blank_auto.connect(self._on_blank_auto)
        self._worker.all_done.connect(self._on_all_done)
        self._worker.start()

    def _on_progress(self, current: int, total: int):
        self._progress.setValue(current)
        self._status_label.setText(f"Processing  {current} / {total}…")

    def _on_file_done(self, filepath: str, result: dict, technique: str):
        fname = os.path.basename(filepath)
        path  = os.path.dirname(filepath)
        self._status_label.setText(f"Done: {fname}  [{technique}]")
        self._graph.update_curves(result, technique, path=path, filename=fname)
        # Show calculated thickness in the mu label (if mu mode was used)
        cd = result.get("CalibratedData", {})
        if cd:
            thickness = cd.get("thickness")
            if thickness is not None:
                self._param_tabs.update_mu_thickness(technique, thickness)
            # Push calculated μ from T and thickness into the μ field for reference
            calc_mu = cd.get("calculated_mu")
            if calc_mu is not None:
                self._param_tabs.update_calculated_mu(technique, calc_mu)
            # Show measured/used transmission next to the override field
            t_used = cd.get("MeasuredTransmission")
            if t_used is None:
                t_used = result.get("calib2DData", {}).get("transmission")
            if t_used is not None:
                self._param_tabs.update_measured_transmission(technique, t_used)
            # USAXS only: show the auto-calculated Qmin
            calc_qmin = cd.get("calculated_qmin")
            if calc_qmin is not None:
                self._param_tabs.update_calculated_qmin(technique, calc_qmin)

    def _on_file_error(self, filepath: str, error: str):
        fname = os.path.basename(filepath)
        self._status_label.setText(f"Error processing: {fname}")
        # Append to persistent error log so errors are not lost
        self._error_log.appendPlainText(f"[{fname}]  {error}")

    def _on_log_message(self, message: str):
        """Show converter warnings/errors in the error log panel."""
        self._error_log.appendPlainText(message)

    def _on_all_done(self):
        self._progress.setVisible(False)
        self._btn_cancel.setEnabled(False)
        self._btn_process_all.setEnabled(True)
        n_errors = self._error_log.document().blockCount() - 1
        if n_errors > 0:
            self._status_label.setText(
                f"Done — {n_errors} error(s), see log below."
            )
        else:
            self._status_label.setText("Done — all files processed successfully.")
        self._save_session()

    def _on_cancel(self):
        if self._worker:
            self._worker.cancel()
        if self._export_worker:
            self._export_worker.cancel()
        self._btn_cancel.setEnabled(False)
        self._status_label.setText("Cancelling after current file…")

    # ── ASCII export ──────────────────────────────────────────────────────────

    def _on_export_ascii(self):
        all_files = self._file_tree.get_all_files()
        if not all_files:
            QMessageBox.information(
                self, "No files",
                "No files are visible in the file tree.\n"
                "Select a folder first (or clear the filter)."
            )
            return

        if self._export_worker and self._export_worker.isRunning():
            QMessageBox.warning(
                self, "Busy",
                "An export job is already running.\n"
                "Please wait or click Cancel."
            )
            return

        output_dir = QFileDialog.getExistingDirectory(
            self,
            "Select output folder for ASCII files",
            self._last_folder,
        )
        if not output_dir:
            return

        self._error_log.clear()
        self._progress.setMaximum(len(all_files))
        self._progress.setValue(0)
        self._progress.setVisible(True)
        self._btn_cancel.setEnabled(True)
        self._btn_export_ascii.setEnabled(False)
        self._btn_process_all.setEnabled(False)
        self._status_label.setText(f"Exporting  0 / {len(all_files)}…")

        self._export_worker = AsciiExportWorker(all_files, output_dir, parent=self)
        self._export_worker.progress.connect(self._on_export_progress)
        self._export_worker.file_done.connect(self._on_export_file_done)
        self._export_worker.file_skipped.connect(self._on_export_file_skipped)
        self._export_worker.file_error.connect(self._on_export_file_error)
        self._export_worker.all_done.connect(self._on_export_all_done)
        self._export_worker.start()

    def _on_export_progress(self, current: int, total: int):
        self._progress.setValue(current)
        self._status_label.setText(f"Exporting  {current} / {total}…")

    def _on_export_file_done(self, filename: str, n: int):
        self._status_label.setText(f"Exported: {filename}  ({n} dataset(s))")

    def _on_export_file_skipped(self, filename: str, reason: str):
        pass  # silently skip files with no processed data

    def _on_export_file_error(self, filename: str, error: str):
        self._error_log.appendPlainText(f"[export] {filename}:  {error}")

    def _on_export_all_done(self, n_exported: int, n_skipped: int, n_errored: int):
        self._progress.setVisible(False)
        self._btn_cancel.setEnabled(False)
        self._btn_export_ascii.setEnabled(True)
        self._btn_process_all.setEnabled(True)
        if n_errored > 0:
            self._status_label.setText(
                f"Export done — {n_exported} file(s) exported, {n_errored} error(s)."
            )
        else:
            self._status_label.setText(
                f"Export done — {n_exported} file(s) exported."
            )

    # ── HDF5 helpers ──────────────────────────────────────────────────────────

    def _read_hdf5_thickness(
        self, path: str, filename: str, technique: str
    ) -> float | None:
        """Read sample thickness from the HDF5 file. Returns mm or None."""
        if technique == "StepScan":
            hdf5_key = "/entry/instrument/bluesky/metadata/sample_thickness_mm"
        else:
            hdf5_key = "/entry/sample/thickness"
        full_path = os.path.join(path, filename)
        try:
            with h5py.File(full_path, "r") as f:
                val = f[hdf5_key][()]
            return float(np.asarray(val).ravel()[0])
        except Exception:
            return None

    def _open_2d_viewer(self, path: str, filename: str, technique: str):
        """Open a non-modal 2D detector image dialog."""
        from .viewer_2d import Viewer2D
        dlg = Viewer2D(path, filename, technique, parent=self)
        dlg.show()

    # ── Session persistence ───────────────────────────────────────────────────

    def _save_session(self):
        state = {
            "last_folder":    self._last_folder,
            "splitter_sizes": self._splitter.sizes(),
        }
        try:
            with open(_SESSION_FILE, "w", encoding="utf-8") as f:
                json.dump(state, f, indent=2)
        except Exception:
            pass  # Non-fatal

    def _restore_session(self):
        try:
            with open(_SESSION_FILE, encoding="utf-8") as f:
                state = json.load(f)
        except Exception:
            return

        folder = state.get("last_folder", "")
        if folder and os.path.isdir(folder):
            self._last_folder = folder
            self._folder_label.setText(f"  {folder}")
            self._folder_label.setStyleSheet("")
            self._file_tree.set_folder(folder)

        sizes = state.get("splitter_sizes")
        if sizes and len(sizes) == 3:
            self._splitter.setSizes(sizes)

    def closeEvent(self, event):
        # Stop background workers before the window (their parent) is
        # destroyed — destroying a running QThread crashes the application.
        for worker in (self._worker, self._export_worker):
            if worker is not None and worker.isRunning():
                worker.cancel()
                worker.wait(10_000)   # cancellation is checked between files
        self._save_session()
        super().closeEvent(event)
