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


def _hdf5_str(val) -> str:
    """Convert an HDF5 attribute/dataset value to a plain Python string.

    HDF5 values can come back as np.ndarray, np.bytes_, bytes, or str
    depending on how the file was written and the h5py version.
    """
    if isinstance(val, np.ndarray):
        val = val.flat[0] if val.size > 0 else ""
    if isinstance(val, (bytes, bytearray, np.bytes_)):
        return val.decode("utf-8", errors="replace")
    return str(val)


def _find_nxcansas_groups(hdf_file: h5py.File) -> list[str]:
    """Return all NXcanSAS group paths in *hdf_file*."""
    result: list[str] = []

    def _visit(name, obj):
        if not isinstance(obj, h5py.Group):
            return
        if not all(attr in obj.attrs and _hdf5_str(obj.attrs[attr]) == v
                   for attr, v in _NXCANSAS_ATTRS.items()):
            return
        for item, expected in _NXCANSAS_ITEMS.items():
            if item not in obj:
                return
            if _hdf5_str(obj[item][()]) != expected:
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


def _write_dat(output_path: str, Q, I, dI, header_lines: list[str], dQ=None):
    """Write Q/I/dI[/dQ] columns to *output_path* with comment header."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    has_dI = dI is not None and len(dI) == len(Q)
    has_dQ = dQ is not None and len(dQ) == len(Q)
    with open(output_path, "w") as f:
        for line in header_lines:
            f.write(f"# {line}\n")
        if has_dQ:
            f.write(f"# {'Q(1/A)':<18} {'I':>18} {'dI':>18} {'dQ(1/A)':>18}\n")
        else:
            f.write(f"# {'Q(1/A)':<18} {'I':>18} {'dI':>18}\n")
        for idx, (q, i) in enumerate(zip(Q, I)):
            e = dI[idx] if has_dI else 0.0
            if has_dQ:
                f.write(f"  {q:18.6e}   {i:18.6e}   {e:18.6e}   {dQ[idx]:18.6e}\n")
            else:
                f.write(f"  {q:18.6e}   {i:18.6e}   {e:18.6e}\n")


def _read_group_data(grp: h5py.Group) -> dict | None:
    """Read Q, I, dI, dQ from a NXcanSAS group's sasdata subgroup.

    Returns a dict with keys Q, I, dI, dQ, units, blankname, thickness,
    slit_length, or None if the sasdata sub-group is missing or has no data.
    dQ key name differs: desmeared uses 'Qdev', SMR uses 'dQw'.
    slit_length is read from 'dQl' (present in SMR groups only).
    """
    sasdata = grp.get("sasdata")
    if sasdata is None:
        return None

    Q  = _read_arr(sasdata.get("Q"))
    I  = _read_arr(sasdata.get("I"))
    if Q is None or I is None or len(Q) < 2:
        return None

    dI   = _read_arr(sasdata.get("Idev"))
    # dQ: 'Qdev' for desmeared/calibrated, 'dQw' for slit-smeared
    dQ   = _read_arr(sasdata.get("Qdev"))
    if dQ is None:
        dQ = _read_arr(sasdata.get("dQw"))
    slit = _read_arr(sasdata.get("dQl"))   # scalar stored as 1-elem array for SMR

    i_ds = sasdata.get("I")
    raw_thick = i_ds.attrs.get("thickness", "") if i_ds is not None else ""
    try:
        thickness = float(np.asarray(raw_thick).flat[0]) if raw_thick != "" else ""
    except Exception:
        thickness = ""
    return {
        "Q":          Q,
        "I":          I,
        "dI":         dI,
        "dQ":         dQ,
        "units":      _read_str_attr(i_ds, "units", "[cm2/cm3]"),
        "blankname":  _read_str_attr(i_ds, "blankname", ""),
        "thickness":  thickness,
        "slit_length": float(slit[0]) if slit is not None and len(slit) > 0 else None,
    }


def export_file(src_path: str, src_filename: str, output_dir: str) -> tuple[int, list[str]]:
    """Export NXcanSAS data from one HDF5 file.

    Priority: desmeared/calibrated data is always exported when present.
    Slit-smeared (SMR) data is exported only when no desmeared data exists,
    because SMR data is not directly comparable without knowing the slit length.
    When SMR is exported, the slit length is included in the header.
    dQ (Q resolution) is written as a fourth column when available.

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

        main_groups = [g for g in groups if "_SMR" not in g]
        smr_groups  = [g for g in groups if "_SMR" in g]

        # ── Export desmeared / calibrated (preferred) ─────────────────────
        for grp_path in main_groups:
            data = _read_group_data(f[grp_path])
            if data is None:
                continue

            header = [
                "Matilda ASCII export",
                f"Source: {src_filename}",
                f"NXcanSAS group: {grp_path}",
                "Q units: 1/angstrom",
                f"I units: {data['units']}",
            ]
            if data["thickness"] != "":
                header.append(f"thickness: {float(data['thickness']):.4f} mm")
            if data["blankname"]:
                header.append(f"blank: {data['blankname']}")

            out_path = os.path.join(output_dir, subfolder, f"{stem}.dat")
            _write_dat(out_path, data["Q"], data["I"], data["dI"], header, dQ=data["dQ"])
            outputs.append(out_path)
            n_exported += 1

        # ── Export slit-smeared only when no desmeared data was written ───
        if n_exported == 0:
            for grp_path in smr_groups:
                data = _read_group_data(f[grp_path])
                if data is None:
                    continue

                header = [
                    "Matilda ASCII export (slit-smeared — desmeared data not available)",
                    f"Source: {src_filename}",
                    f"NXcanSAS group: {grp_path}",
                    "Q units: 1/angstrom",
                    f"I units: {data['units']}",
                ]
                if data["slit_length"] is not None:
                    header.append(f"slit length: {data['slit_length']:.6f} 1/angstrom")
                if data["thickness"] != "":
                    header.append(f"thickness: {float(data['thickness']):.4f} mm")
                if data["blankname"]:
                    header.append(f"blank: {data['blankname']}")

                out_path = os.path.join(output_dir, subfolder, f"{stem}_SMR.dat")
                _write_dat(out_path, data["Q"], data["I"], data["dI"], header, dQ=data["dQ"])
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
                logging.exception(f"ASCII export failed for {filename}: {exc}")

            self.progress.emit(i + 1, total)

        self.all_done.emit(n_exported, n_skipped, n_errored)
