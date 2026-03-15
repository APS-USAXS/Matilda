"""
matilda.gui.data_reduction.file_tree
=====================================
Left panel: folder browser and file list.

Mirrors pyirena's FileTreeWidget with additions for technique-aware
blank assignment via a right-click context menu.

Blank assignment rules:
  - SAXS blanks must be SAXS files.
  - WAXS blanks must be WAXS files.
  - Flyscan / StepScan blanks are cross-compatible (both USAXS).
  - Files whose technique is Unknown show all four options.
"""

import os
import re

try:
    from PySide6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QTreeWidget, QTreeWidgetItem,
        QPushButton, QLabel, QLineEdit, QComboBox, QMenu,
    )
    from PySide6.QtCore import Qt, Signal
    from PySide6.QtGui import QColor, QFont, QBrush
except ImportError:
    from PyQt6.QtWidgets import (
        QWidget, QVBoxLayout, QHBoxLayout, QTreeWidget, QTreeWidgetItem,
        QPushButton, QLabel, QLineEdit, QComboBox, QMenu,
    )
    from PyQt6.QtCore import Qt, pyqtSignal as Signal
    from PyQt6.QtGui import QColor, QFont, QBrush

from .technique_detector import detect_technique

# File extensions considered HDF5 data files
_HDF5_EXTS = {".h5", ".hdf", ".hdf5", ".nxs"}

# QTreeWidgetItem data roles
_ROLE_PATH     = Qt.ItemDataRole.UserRole
_ROLE_IS_FILE  = Qt.ItemDataRole.UserRole + 1
_ROLE_FILENAME = Qt.ItemDataRole.UserRole + 2

_SORT_MODES = ["Name", "Date modified", "Scan number", "Temperature"]

_BLANK_COLOR = QColor(200, 220, 255)   # light-blue tint for blank files


class FileTreeWidget(QWidget):
    """Folder browser and file selector with right-click blank assignment."""

    # Emits list of (path, filename) tuples for currently selected files
    selection_changed = Signal(list)
    # Emits (technique, path, filename) when a blank is assigned
    blank_assigned = Signal(str, str, str)
    # Emits technique name when a blank is cleared
    blank_cleared = Signal(str)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._root_folder: str | None = None
        # blanks: technique → (path, filename) | None
        self._blanks: dict[str, tuple[str, str] | None] = {
            "Flyscan":  None,
            "StepScan": None,
            "SAXS":     None,
            "WAXS":     None,
        }
        self._filter_re: re.Pattern | None = None
        self._sort_mode: str = "Name"
        self._build_ui()

    # ── Public API ────────────────────────────────────────────────────────────

    def set_folder(self, folder: str):
        """Load *folder* and populate the tree."""
        self._root_folder = folder
        self._refresh()

    def get_selected_files(self) -> list[tuple[str, str]]:
        """Return [(path, filename), …] for all selected file items."""
        result = []
        for item in self._tree.selectedItems():
            if item.data(0, _ROLE_IS_FILE):
                result.append((
                    item.data(0, _ROLE_PATH),
                    item.data(0, _ROLE_FILENAME),
                ))
        return result

    def get_all_files(self) -> list[tuple[str, str]]:
        """Return all visible file items (recursively, respecting filter)."""
        return self._collect_visible_files(self._tree.invisibleRootItem())

    def get_blanks(self) -> dict[str, tuple[str, str] | None]:
        """Return current blank assignments."""
        return dict(self._blanks)

    def set_blank_external(self, technique: str, path: str, filename: str):
        """Assign a blank from outside the tree (e.g. Browse dialog)."""
        self._set_blank(technique, path, filename)

    # ── UI construction ───────────────────────────────────────────────────────

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(2, 2, 2, 2)
        layout.setSpacing(4)

        # Sort row
        sort_row = QHBoxLayout()
        sort_row.addWidget(QLabel("Sort:"))
        self._sort_combo = QComboBox()
        self._sort_combo.addItems(_SORT_MODES)
        self._sort_combo.currentTextChanged.connect(self._on_sort_changed)
        sort_row.addWidget(self._sort_combo, 1)
        layout.addLayout(sort_row)

        # Filter row
        filter_row = QHBoxLayout()
        filter_row.addWidget(QLabel("Filter:"))
        self._filter_edit = QLineEdit()
        self._filter_edit.setPlaceholderText("regex…")
        self._filter_edit.setToolTip(
            "Regex or substring filter applied to filenames.\n"
            "Folders are shown when any child matches."
        )
        self._filter_edit.textChanged.connect(self._on_filter_changed)
        filter_row.addWidget(self._filter_edit, 1)

        self._btn_clear_filter = QPushButton("✕")
        self._btn_clear_filter.setMaximumWidth(28)
        self._btn_clear_filter.setToolTip("Clear filter")
        self._btn_clear_filter.clicked.connect(self._filter_edit.clear)
        filter_row.addWidget(self._btn_clear_filter)
        layout.addLayout(filter_row)

        # Tree widget
        self._tree = QTreeWidget()
        self._tree.setHeaderHidden(True)
        self._tree.setSelectionMode(QTreeWidget.SelectionMode.ExtendedSelection)
        self._tree.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self._tree.customContextMenuRequested.connect(self._on_context_menu)
        self._tree.itemSelectionChanged.connect(self._on_selection_changed)
        self._tree.itemExpanded.connect(self._on_item_expanded)
        layout.addWidget(self._tree, 1)

        # Blank summary label at the bottom
        self._blank_summary = QLabel("No blanks assigned")
        self._blank_summary.setWordWrap(True)
        self._blank_summary.setStyleSheet("font-size: 10px; color: #444;")
        layout.addWidget(self._blank_summary)

    # ── Tree population ───────────────────────────────────────────────────────

    def _refresh(self):
        self._tree.clear()
        if not self._root_folder or not os.path.isdir(self._root_folder):
            return
        self._populate_folder(self._tree.invisibleRootItem(), self._root_folder)
        self._tree.expandToDepth(0)
        self._apply_filter()

    def _populate_folder(self, parent_item: QTreeWidgetItem, folder: str):
        """Add subfolders and HDF5 files under *parent_item*."""
        try:
            entries = list(os.scandir(folder))
        except PermissionError:
            return

        dirs  = sorted(e.name for e in entries if e.is_dir(follow_symlinks=False))
        files = self._sort_files(
            [e.name for e in entries
             if e.is_file() and os.path.splitext(e.name)[1].lower() in _HDF5_EXTS],
            folder,
        )

        for d in dirs:
            child = QTreeWidgetItem(parent_item, [f"📁  {d}"])
            child.setData(0, _ROLE_PATH, os.path.join(folder, d))
            child.setData(0, _ROLE_IS_FILE, False)
            # Placeholder so Qt shows the expand arrow
            placeholder = QTreeWidgetItem(child, [""])
            placeholder.setData(0, _ROLE_IS_FILE, False)

        for fname in files:
            self._add_file_item(parent_item, folder, fname)

    def _add_file_item(
        self, parent: QTreeWidgetItem, folder: str, fname: str
    ) -> QTreeWidgetItem:
        item = QTreeWidgetItem(parent, [fname])
        item.setData(0, _ROLE_PATH, folder)
        item.setData(0, _ROLE_IS_FILE, True)
        item.setData(0, _ROLE_FILENAME, fname)
        self._refresh_item_appearance(item)
        return item

    def _refresh_item_appearance(self, item: QTreeWidgetItem):
        if not item.data(0, _ROLE_IS_FILE):
            return
        path  = item.data(0, _ROLE_PATH)
        fname = item.data(0, _ROLE_FILENAME)
        blank_for = self._blank_techniques_for(path, fname)

        font = item.font(0)
        if blank_for:
            font.setItalic(True)
            item.setFont(0, font)
            item.setBackground(0, QBrush(_BLANK_COLOR))
            item.setToolTip(0, f"Blank for: {', '.join(blank_for)}")
        else:
            font.setItalic(False)
            item.setFont(0, font)
            item.setBackground(0, QBrush())
            item.setToolTip(0, "")

    def _blank_techniques_for(self, path: str, fname: str) -> list[str]:
        """Return list of techniques this file is currently assigned as blank for."""
        return [
            t for t, blank in self._blanks.items()
            if blank and blank[0] == path and blank[1] == fname
        ]

    # ── Lazy subfolder loading ────────────────────────────────────────────────

    def _on_item_expanded(self, item: QTreeWidgetItem):
        if item.data(0, _ROLE_IS_FILE):
            return
        # Has exactly one empty placeholder child → load real contents
        if (item.childCount() == 1
                and not item.child(0).data(0, _ROLE_IS_FILE)
                and item.child(0).text(0) == ""):
            item.takeChild(0)
            folder = item.data(0, _ROLE_PATH)
            self._populate_folder(item, folder)
            self._apply_filter()

    # ── Sorting ───────────────────────────────────────────────────────────────

    def _sort_files(self, files: list[str], folder: str) -> list[str]:
        mode = self._sort_mode
        if mode == "Name":
            return sorted(files)
        if mode == "Date modified":
            return sorted(
                files,
                key=lambda f: os.path.getmtime(os.path.join(folder, f)),
            )
        if mode == "Scan number":
            def scan_num(f: str) -> int:
                m = re.search(r"_(\d+)\.", f)
                return int(m.group(1)) if m else 0
            return sorted(files, key=scan_num)
        if mode == "Temperature":
            def temp_key(f: str) -> float:
                m = re.search(r"_(-?\d+(?:\.\d+)?)C", f, re.IGNORECASE)
                return float(m.group(1)) if m else 0.0
            return sorted(files, key=temp_key)
        return files

    def _on_sort_changed(self, text: str):
        self._sort_mode = text
        self._refresh()

    # ── Filtering ─────────────────────────────────────────────────────────────

    def _on_filter_changed(self, text: str):
        if text:
            try:
                self._filter_re = re.compile(text, re.IGNORECASE)
            except re.error:
                self._filter_re = re.compile(re.escape(text), re.IGNORECASE)
        else:
            self._filter_re = None
        self._apply_filter()

    def _apply_filter(self):
        self._filter_item(self._tree.invisibleRootItem())

    def _filter_item(self, parent: QTreeWidgetItem) -> bool:
        """Recursively show/hide items. Returns True if any visible child exists."""
        any_visible = False
        for i in range(parent.childCount()):
            child = parent.child(i)
            if child.data(0, _ROLE_IS_FILE):
                fname = child.data(0, _ROLE_FILENAME) or child.text(0)
                visible = self._matches_filter(fname)
                child.setHidden(not visible)
                if visible:
                    any_visible = True
            else:
                child_visible = self._filter_item(child)
                child.setHidden(not child_visible)
                if child_visible:
                    any_visible = True
        return any_visible

    def _matches_filter(self, fname: str) -> bool:
        if self._filter_re is None:
            return True
        return bool(self._filter_re.search(fname))

    # ── Selection ─────────────────────────────────────────────────────────────

    def _on_selection_changed(self):
        self.selection_changed.emit(self.get_selected_files())

    # ── Context menu ──────────────────────────────────────────────────────────

    def _on_context_menu(self, pos):
        item = self._tree.itemAt(pos)
        if item is None or not item.data(0, _ROLE_IS_FILE):
            return

        path  = item.data(0, _ROLE_PATH)
        fname = item.data(0, _ROLE_FILENAME)
        technique = detect_technique(path, fname)

        menu = QMenu(self)

        # Which blank options to expose depends on technique
        if technique in ("Flyscan", "StepScan", "Unknown"):
            self._add_blank_action(menu, "Set as Flyscan Blank",  "Flyscan",  path, fname)
            self._add_blank_action(menu, "Set as StepScan Blank", "StepScan", path, fname)

        if technique in ("SAXS", "Unknown"):
            self._add_blank_action(menu, "Set as SAXS Blank", "SAXS", path, fname)

        if technique in ("WAXS", "Unknown"):
            self._add_blank_action(menu, "Set as WAXS Blank", "WAXS", path, fname)

        # Clear options for any technique this file is already blank for
        assigned = self._blank_techniques_for(path, fname)
        if assigned:
            menu.addSeparator()
            for tech in assigned:
                act = menu.addAction(f"Clear as {tech} Blank")
                act.triggered.connect(
                    lambda checked=False, t=tech: self._clear_blank(t)
                )

        if not menu.isEmpty():
            menu.exec(self._tree.viewport().mapToGlobal(pos))

    def _add_blank_action(
        self, menu: QMenu, label: str, technique: str, path: str, fname: str
    ):
        act = menu.addAction(label)
        act.triggered.connect(
            lambda checked=False, t=technique, p=path, f=fname: self._set_blank(t, p, f)
        )

    # ── Blank management ──────────────────────────────────────────────────────

    def _set_blank(self, technique: str, path: str, fname: str):
        self._blanks[technique] = (path, fname)
        self._refresh_all_item_appearances()
        self._update_blank_summary()
        self.blank_assigned.emit(technique, path, fname)

    def _clear_blank(self, technique: str):
        self._blanks[technique] = None
        self._refresh_all_item_appearances()
        self._update_blank_summary()
        self.blank_cleared.emit(technique)

    def _refresh_all_item_appearances(self):
        self._refresh_items_recursive(self._tree.invisibleRootItem())

    def _refresh_items_recursive(self, parent: QTreeWidgetItem):
        for i in range(parent.childCount()):
            child = parent.child(i)
            if child.data(0, _ROLE_IS_FILE):
                self._refresh_item_appearance(child)
            else:
                self._refresh_items_recursive(child)

    def _update_blank_summary(self):
        lines = [
            f"{t}: {blank[1]}"
            for t, blank in self._blanks.items()
            if blank is not None
        ]
        self._blank_summary.setText(
            "Blanks:\n" + "\n".join(lines) if lines else "No blanks assigned"
        )

    # ── Helpers ───────────────────────────────────────────────────────────────

    def _collect_visible_files(
        self, parent: QTreeWidgetItem
    ) -> list[tuple[str, str]]:
        result = []
        for i in range(parent.childCount()):
            child = parent.child(i)
            if child.isHidden():
                continue
            if child.data(0, _ROLE_IS_FILE):
                result.append((
                    child.data(0, _ROLE_PATH),
                    child.data(0, _ROLE_FILENAME),
                ))
            else:
                result.extend(self._collect_visible_files(child))
        return result
