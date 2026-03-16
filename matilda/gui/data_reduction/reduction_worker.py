"""
matilda.gui.data_reduction.reduction_worker
============================================
QThread worker that runs data reduction on a list of HDF5 files.

Always passes recalculateAllData=True (the user is explicitly requesting
re-reduction with custom parameters, so cached NXcanSAS groups are deleted
and the full pipeline is re-run).

NOTE (future): The converters (convertFlyscan, convertUSAXS, convertSWAXS)
currently do not accept all GUI parameters (thickness, npts, az_min/max,
extrap settings).  When those converters are updated to accept keyword
overrides, the _process_one() method below is the single place to wire them.
The params dict is passed through unchanged so it is ready for that wiring.
"""

import os

try:
    from PySide6.QtCore import QThread, Signal
except ImportError:
    from PyQt6.QtCore import QThread, pyqtSignal as Signal

from .technique_detector import detect_technique

from matilda.convertFlyscan import processFlyscan
from matilda.convertUSAXS import processStepscan
from matilda.convertSWAXS import process2Ddata


class ReductionWorker(QThread):
    """Background thread that processes a list of (path, filename) tuples.

    Signals
    -------
    progress(current, total)
        Emitted after each file attempt (success or failure).
    file_done(filepath, result_dict, technique)
        Emitted when a file is processed successfully.
    file_error(filepath, error_message)
        Emitted when processing raises an exception.
    all_done()
        Emitted after the last file (or on cancellation).
    """

    progress  = Signal(int, int)         # (current, total)
    file_done = Signal(str, dict, str)   # (filepath, result, technique)
    file_error = Signal(str, str)        # (filepath, error_message)
    all_done  = Signal()

    def __init__(
        self,
        file_list: list[tuple[str, str]],
        blanks: dict[str, tuple[str, str] | None],
        params_by_technique: dict[str, dict],
        parent=None,
    ):
        super().__init__(parent)
        self._file_list = file_list
        self._blanks = blanks
        self._params = params_by_technique
        self._cancelled = False

    def cancel(self):
        """Request cancellation; checked between files."""
        self._cancelled = True

    # ── QThread entry point ───────────────────────────────────────────────────

    def run(self):
        total = len(self._file_list)
        for i, (path, filename) in enumerate(self._file_list):
            if self._cancelled:
                break

            filepath  = os.path.join(path, filename)
            technique = detect_technique(path, filename)
            params    = self._params.get(technique, {})
            blank     = self._resolve_blank(technique)

            try:
                result = self._process_one(path, filename, technique, blank, params)
                self.file_done.emit(filepath, result, technique)
            except Exception as exc:
                self.file_error.emit(filepath, str(exc))

            self.progress.emit(i + 1, total)

        self.all_done.emit()

    # ── Private helpers ───────────────────────────────────────────────────────

    def _resolve_blank(
        self, technique: str
    ) -> tuple[str | None, str | None]:
        """Return (blankPath, blankFilename) with cross-technique USAXS fallback."""
        blank = self._blanks.get(technique)
        if blank:
            return blank

        # USAXS blanks are cross-compatible
        if technique == "Flyscan":
            blank = self._blanks.get("StepScan")
        elif technique == "StepScan":
            blank = self._blanks.get("Flyscan")

        if blank:
            return blank
        return (None, None)

    def _process_one(
        self,
        path: str,
        filename: str,
        technique: str,
        blank: tuple[str | None, str | None],
        params: dict,
    ) -> dict:
        bpath, bfile = blank
        recalc = True   # GUI always forces recalculation

        if technique == "Flyscan":
            result = processFlyscan(
                path, filename,
                blankPath=bpath,
                blankFilename=bfile,
                recalculateAllData=recalc,
            )
        elif technique == "StepScan":
            result = processStepscan(
                path, filename,
                blankPath=bpath,
                blankFilename=bfile,
                recalculateAllData=recalc,
            )
        elif technique in ("SAXS", "WAXS"):
            result = process2Ddata(
                path, filename,
                blankPath=bpath,
                blankFilename=bfile,
                recalculateAllData=recalc,
            )
        else:
            raise ValueError(
                f"Cannot process file with technique {technique!r}: "
                f"{os.path.join(path, filename)}"
            )

        # Attempt to overlay the blank's raw 1D curve in the graph by
        # re-running the converter on the blank file with recalculateAllData=False
        # (uses HDF5 cache — fast).  Silently skipped on any failure.
        blank_rd = self._get_blank_reduced_data(technique, blank)
        if blank_rd is not None:
            result["blankReducedData"] = blank_rd

        return result

    def _get_blank_reduced_data(
        self,
        technique: str,
        blank: tuple[str | None, str | None],
    ) -> dict | None:
        """Return the blank file's reducedData dict, or None on any failure."""
        bpath, bfile = blank
        if not bpath or not bfile:
            return None
        try:
            if technique == "Flyscan":
                r = processFlyscan(bpath, bfile, recalculateAllData=False)
            elif technique == "StepScan":
                r = processStepscan(bpath, bfile, recalculateAllData=False)
            elif technique in ("SAXS", "WAXS"):
                r = process2Ddata(bpath, bfile, recalculateAllData=False)
            else:
                return None
            return r.get("reducedData") or None
        except Exception:
            return None
