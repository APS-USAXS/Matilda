"""
matilda.reduction
=================
Data-reduction subpackage — placeholder for the planned refactoring of
the flat-file reduction modules into a proper subpackage.

Current state (flat layout in matilda/)
----------------------------------------
The reduction code currently lives as top-level modules inside matilda/:

    convertFlyscan.py   → USAXS flyscan 2-D raw → calibrated I(Q)
    convertUSAXS.py     → USAXS step-scan reduction (same pipeline)
    convertSWAXS.py     → SAXS / WAXS area-detector reduction via pyFAI
    desmearing.py       → Lake/Strobl slit-smearing correction algorithm
    supportFunctions.py → shared numerical helpers (beam-centre, blank, rebin…)
    supportNikaFunctions.py → Nika ↔ Fit2D ↔ pyFAI geometry conversion

Planned future layout (matilda/reduction/)
------------------------------------------
Once the GUI and other tools require importing reduction code as a proper
subpackage, the flat modules above should be moved here and converted to
use relative imports.  Suggested target names:

    flyscan.py          ← convertFlyscan.py
    stepscan.py         ← convertUSAXS.py
    swaxs.py            ← convertSWAXS.py
    desmearing.py       ← desmearing.py  (unchanged)
    support.py          ← supportFunctions.py
    geometry.py         ← supportNikaFunctions.py

Migration note
--------------
The flat modules already use relative imports (``from .hdf5code import …``),
so a move only needs those paths adjusted (e.g. ``from ..hdf5code import …``)
and the console-script entry point in pyproject.toml left pointing at
matilda.matilda:main.  The split is an open follow-up in PLAN.md; nothing
lives in this subpackage yet.
"""
