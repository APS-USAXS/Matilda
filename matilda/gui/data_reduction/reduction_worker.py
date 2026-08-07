"""
matilda.gui.data_reduction.reduction_worker
============================================
QThread worker that runs data reduction on a list of HDF5 files.

Always passes recalculateAllData=True (the user is explicitly requesting
re-reduction with custom parameters, so cached NXcanSAS groups are deleted
and the full pipeline is re-run).

GUI parameters (thickness, npts, desmearing, minQMinFindRatio) are forwarded
from the params dict to the converter functions.
"""

import logging
import os

from .._qt import QThread, Signal

from .technique_detector import detect_technique

from matilda.convertFlyscan import processFlyscan
from matilda.convertUSAXS import processStepscan
from matilda.convertSWAXS import process2Ddata
from matilda.supportFunctions import findProperBlankScan


class _SignalHandler(logging.Handler):
    """Logging handler that emits records via a Qt signal."""

    def __init__(self, signal):
        super().__init__(level=logging.WARNING)
        self._signal = signal

    def emit(self, record):
        try:
            self._signal.emit(self.format(record))
        except RuntimeError:
            pass  # worker already destroyed


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

    progress    = Signal(int, int)         # (current, total)
    file_done   = Signal(str, dict, str)   # (filepath, result, technique)
    file_error  = Signal(str, str)        # (filepath, error_message)
    log_message = Signal(str)             # (formatted log message)
    blank_auto  = Signal(str, str)        # (technique, blank_filename)
    all_done    = Signal()

    def __init__(
        self,
        file_list: list[tuple[str, str]],
        blanks: dict[str, tuple[str, str] | None],
        params_by_technique: dict[str, dict],
        all_files: list[tuple[str, str]] | None = None,
        parent=None,
    ):
        super().__init__(parent)
        self._file_list = file_list
        self._blanks = blanks
        self._params = params_by_technique
        self._all_files = all_files or []
        self._cancelled = False

    def cancel(self):
        """Request cancellation; checked between files."""
        self._cancelled = True

    # ── QThread entry point ───────────────────────────────────────────────────

    def run(self):
        # Route warnings from converter code to the GUI error log
        handler = _SignalHandler(self.log_message)
        handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
        root_logger = logging.getLogger()
        root_logger.addHandler(handler)

        try:
            total = len(self._file_list)
            for i, (path, filename) in enumerate(self._file_list):
                if self._cancelled:
                    break

                filepath  = os.path.join(path, filename)
                technique = detect_technique(path, filename)
                params    = self._params.get(technique, {})
                blank     = self._resolve_blank(technique, path, filename, params)

                try:
                    result = self._process_one(path, filename, technique, blank, params)
                    self.file_done.emit(filepath, result, technique)
                except Exception as exc:
                    self.file_error.emit(filepath, str(exc))

                self.progress.emit(i + 1, total)

            self.all_done.emit()
        finally:
            root_logger.removeHandler(handler)

    # ── Private helpers ───────────────────────────────────────────────────────

    def _resolve_blank(
        self, technique: str, path: str, filename: str, params: dict,
    ) -> tuple[str | None, str | None]:
        """Return (blankPath, blankFilename).

        When blank_mode is "manual", use only the explicitly assigned blank.
        When blank_mode is "auto (nearest preceding)" (default), try the
        manually assigned blank first, then fall back to automatic detection
        using the same nearest-preceding logic as the matilda daemon.
        """
        blank_mode = params.get("blank_mode", "auto (nearest preceding)")

        # ── Manual-only mode: return assigned blank or nothing ────────────
        if "manual" in blank_mode.lower():
            blank = self._blanks.get(technique)
            if blank:
                return blank
            # USAXS blanks are cross-compatible
            if technique == "Flyscan":
                blank = self._blanks.get("StepScan")
            elif technique == "StepScan":
                blank = self._blanks.get("Flyscan")
            return blank if blank else (None, None)

        # ── Auto mode: try assigned blank first, then auto-detect ─────────
        blank = self._blanks.get(technique)
        if blank:
            return blank
        # USAXS cross-compatible fallback
        if technique == "Flyscan":
            blank = self._blanks.get("StepScan")
        elif technique == "StepScan":
            blank = self._blanks.get("Flyscan")
        if blank:
            return blank

        # No manual blank assigned — build a blank list from the loaded
        # files and use findProperBlankScan (same logic as the daemon).
        # USAXS blanks are cross-compatible (Flyscan ↔ StepScan).
        compatible = {technique}
        if technique in ("Flyscan", "StepScan"):
            compatible = {"Flyscan", "StepScan"}
        blank_candidates = [
            (p, f) for p, f in self._all_files
            if "blank" in f.lower()
            and detect_technique(p, f) in compatible
        ]
        if blank_candidates:
            bp, bf = findProperBlankScan(path, filename, blank_candidates)
            if bp is not None and bf is not None:
                logging.info(f"Auto-selected blank for {filename}: {bf}")
                self.blank_auto.emit(technique, bf)
                return (bp, bf)

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
                num_points=params.get("npts", 500),
                desmear_iter=params.get("desmear_iter", 20),
                extrap_method=params.get("extrap_method", "PowerLaw w flat"),
                extrap_qstart=params.get("extrap_qstart", 0.15),  # keep in sync with converter default
                minQMinFindRatio=params.get("minQMinFindRatio", 1.05),
                thickness_override=params.get("thickness"),
                use_mu=params.get("use_mu", False),
                mu=params.get("mu"),
                per_gram=params.get("per_gram", False),
                density=params.get("density"),
                transmission_override=params.get("transmission_override"),
                qmin_override=params.get("qmin_override"),
                desmear_method=params.get("desmear_method", "lake"),
                gp_length_scale=params.get("gp_length_scale", 0.5),
                gp_kernel=params.get("gp_kernel", "matern32"),
            )
        elif technique == "StepScan":
            result = processStepscan(
                path, filename,
                blankPath=bpath,
                blankFilename=bfile,
                recalculateAllData=recalc,
                desmear_iter=params.get("desmear_iter", 20),
                extrap_method=params.get("extrap_method", "PowerLaw w flat"),
                extrap_qstart=params.get("extrap_qstart", 0.15),  # keep in sync with converter default
                minQMinFindRatio=params.get("minQMinFindRatio", 1.05),
                thickness_override=params.get("thickness"),
                use_mu=params.get("use_mu", False),
                mu=params.get("mu"),
                per_gram=params.get("per_gram", False),
                density=params.get("density"),
                transmission_override=params.get("transmission_override"),
                qmin_override=params.get("qmin_override"),
                desmear_method=params.get("desmear_method", "lake"),
                gp_length_scale=params.get("gp_length_scale", 0.5),
                gp_kernel=params.get("gp_kernel", "matern32"),
            )
        elif technique in ("SAXS", "WAXS"):
            result = process2Ddata(
                path, filename,
                blankPath=bpath,
                blankFilename=bfile,
                recalculateAllData=recalc,
                npts=params.get("npts"),
                thickness_override=params.get("thickness"),
                use_mu=params.get("use_mu", False),
                mu=params.get("mu"),
                per_gram=params.get("per_gram", False),
                density=params.get("density"),
                transmission_override=params.get("transmission_override"),
            )
        else:
            raise ValueError(
                f"Cannot process file with technique {technique!r}: "
                f"{os.path.join(path, filename)}"
            )

        # NOTE: blank curve overlay via _get_blank_reduced_data is disabled.
        # Calling the converter on the blank file immediately after writing
        # to the sample HDF5 causes HDF5 file-locking contention on some
        # platforms, hanging the worker thread indefinitely.  Blank overlay
        # will be re-enabled once a direct h5py read approach is implemented.

        return result
