"""
matilda.gui.data_reduction.ascii_exporter
==========================================
Background worker that exports NXcanSAS data from HDF5 files to ASCII .dat files.

For each HDF5 file that contains processed NXcanSAS groups:
  - Desmeared / calibrated data  → <output>/<subfolder>/<stem>.dat
  - Slit-smeared data (USAXS)    → <output>/<subfolder>/<stem>_SMR.dat

Output subfolder is the last component of the source file's directory
(e.g. data_usaxs, data_saxs, data_waxs, data_usaxs_merged).

File format
-----------
# Matilda ASCII export
# Source: filename.h5
# Q units: 1/angstrom
# I units: [cm2/cm3]
# thickness: 0.1234 mm
# Q(1/A)           I(units)           dI(units)
1.23456e-04   4.56789e+02   1.23456e+01
...
"""

import logging
import os

import h5py
import numpy as np

try:
    from PySide6.QtCore import QThread, Signal
except ImportError:
    from PyQt6.QtCore import QThread, pyqtSignal as Signal


_NXCANSAS_ATTRS = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
_NXCANSAS_ITEMS = {'definition': 'NXcanSAS'}


def _find_nxcansas_groups(hdf_file: h5py.File) -> list[str]:
    """Return all NXcanSAS group paths in *hdf_file*."""
    result: list[str] = []

    def _visit(name, obj):
        if not isinstance(obj, h5py.Group):
            return
        if not all(attr in obj.attrs and obj.attrs[attr] == v
                   for attr, v in _NXCANSAS_ATTRS.items()):
            return
        for item, expected in _NXCANSAS_ITEMS.items():
            if item not in obj:
                return
            val = obj[item][()]
            if isinstance(val, (bytes, bytearray)):
                val = val.decode()
            if val != expected:
                return
        result.append(name)

    hdf_file.visititems(_visit)
    return result


def _read_arr(ds) -> np.ndarray | None:
    """Read a dataset safely, return None if missing or not usable."""
    if ds is None:
        return None
    try:
        arr = np.asarray(ds[()], dtype=float).ravel()
        return arr if arr.size > 0 else None
    except Exception:
        return None


def _read_str_attr(ds, key: str, default: str = "") -> str:
    try:
        v = ds.attrs.get(key, default)
        if isinstance(v, (bytes, bytearray)):
            return v.decode()
        return str(v)
    except Exception:
        return default


def _write_dat(output_path: str, Q, I, dI, header_lines: list[str]):
    """Write Q/I/dI columns to *output_path* with comment header."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        for line in header_lines:
            f.write(f"# {line}\n")
        f.write(f"# {'Q(1/A)':<18} {'I':>18} {'dI':>18}\n")
        if dI is not None and len(dI) == len(Q):
            for q, i, e in zip(Q, I, dI):
                f.write(f"  {q:18.6e}   {i:18.6e}   {e:18.6e}\n")
        else:
            for q, i in zip(Q, I):
                f.write(f"  {q:18.6e}   {i:18.6e}   {'0.0':>18}\n")


def export_file(src_path: str, src_filename: str, output_dir: str) -> tuple[int, list[str]]:
    """Export NXcanSAS data from one HDF5 file.

    Returns (n_exported, list_of_output_paths).
    Returns (0, []) if the file has no NXcanSAS data.
    Raises on I/O errors.
    """
    filepath = os.path.join(src_path, src_filename)
    subfolder = os.path.basename(src_path)
    stem = os.path.splitext(src_filename)[0]
    n_exported = 0
    outputs: list[str] = []

    with h5py.File(filepath, "r") as f:
        groups = _find_nxcansas_groups(f)
        if not groups:
            return 0, []

        for grp_path in groups:
            grp = f[grp_path]
            is_smr = "_SMR" in grp_path

            # Locate sasdata sub-group
            sasdata = grp.get("sasdata")
            if sasdata is None:
                continue

            Q  = _read_arr(sasdata.get("Q"))
            I  = _read_arr(sasdata.get("I"))
            dI = _read_arr(sasdata.get("Idev"))

            if Q is None or I is None or len(Q) < 2:
                continue

            # Build header
            i_ds = sasdata.get("I")
            units    = _read_str_attr(i_ds, "units", "[cm2/cm3]")
            blankname = _read_str_attr(i_ds, "blankname", "")
            thickness = i_ds.attrs.get("thickness", "") if i_ds is not None else ""

            header = [
                f"Matilda ASCII export",
                f"Source: {src_filename}",
                f"NXcanSAS group: {grp_path}",
                f"Q units: 1/angstrom",
                f"I units: {units}",
            ]
            if thickness != "":
                header.append(f"thickness: {float(thickness):.4f} mm")
            if blankname:
                header.append(f"blank: {blankname}")

            suffix = "_SMR" if is_smr else ""
            out_name = f"{stem}{suffix}.dat"
            out_path = os.path.join(output_dir, subfolder, out_name)

            _write_dat(out_path, Q, I, dI, header)
            outputs.append(out_path)
            n_exported += 1

    return n_exported, outputs


class AsciiExportWorker(QThread):
    """Background thread that exports all HDF5 files to ASCII .dat files.

    Signals
    -------
    progress(current, total)
    file_done(src_filename, n_datasets)
    file_skipped(src_filename, reason)
    file_error(src_filename, error_message)
    all_done(n_files_exported, n_files_skipped, n_files_errored)
    """

    progress     = Signal(int, int)
    file_done    = Signal(str, int)     # (filename, n_datasets_written)
    file_skipped = Signal(str, str)     # (filename, reason)
    file_error   = Signal(str, str)     # (filename, error)
    all_done     = Signal(int, int, int)

    def __init__(
        self,
        file_list: list[tuple[str, str]],
        output_dir: str,
        parent=None,
    ):
        super().__init__(parent)
        self._file_list = file_list
        self._output_dir = output_dir
        self._cancelled = False

    def cancel(self):
        self._cancelled = True

    def run(self):
        total = len(self._file_list)
        n_exported = n_skipped = n_errored = 0

        for i, (path, filename) in enumerate(self._file_list):
            if self._cancelled:
                break
            try:
                n, _ = export_file(path, filename, self._output_dir)
                if n > 0:
                    self.file_done.emit(filename, n)
                    n_exported += 1
                else:
                    self.file_skipped.emit(filename, "no NXcanSAS data")
                    n_skipped += 1
            except Exception as exc:
                self.file_error.emit(filename, str(exc))
                n_errored += 1
                logging.warning(f"ASCII export failed for {filename}: {exc}")

            self.progress.emit(i + 1, total)

        self.all_done.emit(n_exported, n_skipped, n_errored)
