"""
matilda.gui.data_reduction
==========================
Interactive data-reduction GUI for USAXS/SAXS/WAXS data.

Entry point: ``matilda-gui`` (see pyproject.toml).

Usage
-----
From the command line::

    matilda-gui

From Python::

    from matilda.gui.data_reduction import main
    main()
"""

from .main_window import MatildaReductionWindow


def main():
    """Launch the Matilda data-reduction GUI."""
    import sys

    from .._qt import QApplication

    app = QApplication.instance() or QApplication(sys.argv)
    window = MatildaReductionWindow()
    window.show()
    sys.exit(app.exec())


__all__ = ["MatildaReductionWindow", "main"]
