"""
matilda.gui._qt — single Qt import point for the whole ``matilda.gui`` package.

PySide6 is the supported binding (it is what ``pip install matilda[gui]`` and
``environment.yml`` install); PyQt6 is accepted as a fallback so an environment
that only has PyQt6 still runs.  **Never install both into one environment** —
they fight over Qt plugin resolution and the symptom on macOS is the familiar
"could not load the Qt platform plugin cocoa".

Every GUI module imports the Qt names it needs from here instead of repeating a
``try: from PySide6 … except ImportError: from PyQt6 …`` block.  That removed
eight duplicated import blocks, and it means a binding change is a one-file
edit rather than a sweep through the package.

Two import styles are supported:

* flat names — ``from matilda.gui._qt import QWidget, Qt, Signal``
* submodules — ``from matilda.gui._qt import QtWidgets, QtCore, QtGui``

``Signal`` is normalised across bindings (it is ``pyqtSignal`` under PyQt6).

Nothing here may be imported from the daemon path: ``matilda.matilda`` and the
``convert*`` modules must stay headless (see CLAUDE.md invariant 3).
"""

try:
    from PySide6 import QtCore, QtGui, QtWidgets
    from PySide6.QtCore import Signal
    QT_BINDING = "PySide6"
except ImportError as _pyside_error:  # pragma: no cover - only without PySide6
    try:
        from PyQt6 import QtCore, QtGui, QtWidgets  # type: ignore[no-redef]
        from PyQt6.QtCore import pyqtSignal as Signal  # type: ignore[no-redef]
        QT_BINDING = "PyQt6"
    except ImportError as _pyqt_error:  # pragma: no cover
        # Do NOT report this as "not installed": ImportError also covers an
        # installed-but-unloadable binding — a shiboken6 version skew, or the
        # freetype/harfbuzz ABI mismatch that environment.yml pins against
        # ("undefined symbol: FT_Get_Colorline_Stops").  Telling the user to
        # reinstall Qt in that case sends them in a circle, so show both the
        # install hint and the underlying PySide6 error.
        del _pyqt_error  # PySide6 is the binding we report on
        raise ImportError(
            "matilda.gui could not load a Qt binding.\n"
            "  If Qt is missing:      pip install matilda[gui]\n"
            "  If Qt is installed:    the import failed for another reason — "
            "the original PySide6 error is chained below. A shiboken6 version "
            "skew or a freetype/harfbuzz ABI mismatch both look like this; see "
            "the notes in environment.yml."
        ) from _pyside_error

# --- QtWidgets ---------------------------------------------------------------
QAbstractItemView = QtWidgets.QAbstractItemView
QApplication = QtWidgets.QApplication
QCheckBox = QtWidgets.QCheckBox
QComboBox = QtWidgets.QComboBox
QDialog = QtWidgets.QDialog
QDialogButtonBox = QtWidgets.QDialogButtonBox
QDoubleSpinBox = QtWidgets.QDoubleSpinBox
QFileDialog = QtWidgets.QFileDialog
QFormLayout = QtWidgets.QFormLayout
QFrame = QtWidgets.QFrame
QGridLayout = QtWidgets.QGridLayout
QGroupBox = QtWidgets.QGroupBox
QHBoxLayout = QtWidgets.QHBoxLayout
QHeaderView = QtWidgets.QHeaderView
QInputDialog = QtWidgets.QInputDialog
QLabel = QtWidgets.QLabel
QLineEdit = QtWidgets.QLineEdit
QListWidget = QtWidgets.QListWidget
QListWidgetItem = QtWidgets.QListWidgetItem
QMainWindow = QtWidgets.QMainWindow
QMenu = QtWidgets.QMenu
QMessageBox = QtWidgets.QMessageBox
QPlainTextEdit = QtWidgets.QPlainTextEdit
QProgressBar = QtWidgets.QProgressBar
QPushButton = QtWidgets.QPushButton
QScrollArea = QtWidgets.QScrollArea
QSizePolicy = QtWidgets.QSizePolicy
QSpinBox = QtWidgets.QSpinBox
QSplitter = QtWidgets.QSplitter
QTabWidget = QtWidgets.QTabWidget
QTableWidget = QtWidgets.QTableWidget
QTableWidgetItem = QtWidgets.QTableWidgetItem
QTextEdit = QtWidgets.QTextEdit
QTreeWidget = QtWidgets.QTreeWidget
QTreeWidgetItem = QtWidgets.QTreeWidgetItem
QVBoxLayout = QtWidgets.QVBoxLayout
QWidget = QtWidgets.QWidget

# --- QtCore ------------------------------------------------------------------
Qt = QtCore.Qt
QSettings = QtCore.QSettings
QThread = QtCore.QThread
QTimer = QtCore.QTimer

# --- QtGui -------------------------------------------------------------------
QBrush = QtGui.QBrush
QColor = QtGui.QColor
QFont = QtGui.QFont

__all__ = [
    "QtWidgets", "QtCore", "QtGui", "Signal", "QT_BINDING",
    # QtWidgets
    "QAbstractItemView", "QApplication", "QCheckBox", "QComboBox", "QDialog",
    "QDialogButtonBox", "QDoubleSpinBox", "QFileDialog", "QFormLayout",
    "QFrame", "QGridLayout", "QGroupBox", "QHBoxLayout", "QHeaderView",
    "QInputDialog", "QLabel", "QLineEdit", "QListWidget", "QListWidgetItem",
    "QMainWindow", "QMenu", "QMessageBox", "QPlainTextEdit", "QProgressBar",
    "QPushButton", "QScrollArea", "QSizePolicy", "QSpinBox", "QSplitter",
    "QTabWidget",
    "QTableWidget", "QTableWidgetItem", "QTextEdit", "QTreeWidget",
    "QTreeWidgetItem", "QVBoxLayout", "QWidget",
    # QtCore
    "Qt", "QSettings", "QThread", "QTimer",
    # QtGui
    "QBrush", "QColor", "QFont",
]
