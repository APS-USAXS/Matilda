"""
matilda.gui
===========
GUI subpackage for Matilda — placeholder for future PyQt6/pyqtgraph tools.

Planned tools (not yet implemented)
-------------------------------------
reduction_gui
    Interactive 2-D → 1-D data reduction with user control over masking,
    beam-centre, integration range, and blank selection.  Uses pyqtgraph
    for real-time display of the 2-D detector image and the 1-D result.

survey_tool
    Sample survey / overview: browse all scans in a folder, display
    thumbnails, select scans for manual re-reduction.

analysis_gui
    Launcher / coordinator for pyirena analysis workflows, with live
    parameter editing and result display.

GUI framework
-------------
All GUI code must use PyQt6 (or PySide6 as a drop-in substitute).
All scientific plotting must use pyqtgraph.
matplotlib is NOT to be used in GUI code without explicit approval;
its use is limited to headless file export in matilda.plots.

Dependencies
------------
Declared as optional in pyproject.toml:
    pip install matilda[gui]      # adds PyQt6 + pyqtgraph
Or via conda (already in environment.yml):
    pyqt6
    pyqtgraph
"""
