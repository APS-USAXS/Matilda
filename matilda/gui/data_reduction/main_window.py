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

import h5py
import numpy as np

try:
    from PySide6.QtWidgets import (
        QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
        QToolBar, QPushButton, QLabel, QProgressBar,
        QFileDialog, QMessageBox, QPlainTextEdit,
    )
    from PySide6.QtCore import Qt
    from PySide6.QtGui import QColor
except ImportError:
    from PyQt6.QtWidgets import (
        QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
        QToolBar, QPushButton, QLabel, QProgressBar,
        QFileDialog, QMessageBox, QPlainTextEdit,
    )
    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QColor

from .file_tree import FileTreeWidget
from .parameter_tabs import ParameterTabWidget
from .graph_panel import GraphPanel
from .reduction_worker import ReductionWorker
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
        all_files = self._file_tree.get_all_files()

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
        self._btn_cancel.setEnabled(False)
        self._status_label.setText("Cancelling after current file…")

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
        self._save_session()
        super().closeEvent(event)
