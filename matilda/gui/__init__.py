"""
matilda.gui
===========
GUI subpackage for Matilda — PySide6/pyqtgraph tools.

Planned tools (not yet implemented)
-------------------------------------
reduction_gui
    Interactive 2-D → 1-D data reduction with user control over masking,
    beam-centre, integration range, and blank selection.  Uses pyqtgraph
    for real-time display of the 2-D detector image and the 1-D result.

sample_plate_setup
    Replaces the Igor Pro "Setup Sample Plates" tool.
    Define sample positions (name, SX, SY, thickness, USAXS/SAXS/WAXS flags),
    load plate geometry templates, click on a plate image to assign positions,
    save/load sets in HDF5, export Bluesky command files (.mac), and
    optionally drive the sample stage via EPICS (Beamline Survey).

survey_tool
    Sample survey / overview: browse all scans in a folder, display
    thumbnails, select scans for manual re-reduction.

analysis_gui
    Launcher / coordinator for pyirena analysis workflows, with live
    parameter editing and result display.

GUI framework
-------------
All GUI code must use PySide6, imported through the ``matilda.gui._qt`` shim —
never ``from PySide6… import`` directly in a panel, and never PyQt6 alongside
it in the same environment (they fight over Qt plugin resolution; on macOS the
symptom is "could not load the Qt platform plugin cocoa").  The shim keeps a
PyQt6 fallback for environments that only have that binding.

All scientific plotting must use pyqtgraph.
matplotlib is NOT to be used in GUI code without explicit approval;
its use is limited to headless file export in matilda.plotData.

Dependencies
------------
Declared as optional in pyproject.toml:
    pip install matilda[gui]      # adds PySide6 + pyqtgraph + pyepics
Or via conda (already in environment.yml):
    pyside6>=6.4,<6.8
    pyqtgraph
"""
