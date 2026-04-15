"""
matilda.gui.data_reduction.technique_detector
=============================================
Heuristic detection of USAXS/SAXS/WAXS technique from an HDF5 file.

Detection priority:
  1. In-memory cache (avoid re-opening files).
  2. File extension: .hdf → area-detector (SAXS or WAXS); .h5 → USAXS.
  3. HDF5 content: key-presence heuristics.
  4. Parent-folder name suffix (_saxs, _waxs, _usaxs, FLyscans, StepScans).

Returns one of: 'Flyscan', 'StepScan', 'SAXS', 'WAXS', 'Unknown'.
"""

import os

import h5py

# In-memory cache: absolute path → technique string
_CACHE: dict[str, str] = {}

TECHNIQUES = ("Flyscan", "StepScan", "SAXS", "WAXS", "Unknown")


def detect_technique(path: str, filename: str) -> str:
    """Return one of: 'Flyscan', 'StepScan', 'SAXS', 'WAXS', 'Unknown'."""
    full_path = os.path.normpath(os.path.join(path, filename))
    if full_path in _CACHE:
        return _CACHE[full_path]

    technique = _detect(path, filename, full_path)
    _CACHE[full_path] = technique
    return technique


def clear_cache() -> None:
    """Flush the detection cache (e.g. after a file is replaced on disk)."""
    _CACHE.clear()


# ── Private helpers ───────────────────────────────────────────────────────────

def _detect(path: str, filename: str, full_path: str) -> str:
    ext = os.path.splitext(filename)[1].lower()
    if ext in (".hdf", ".hdf5"):
        return _detect_area_detector(full_path, path)
    elif ext in (".h5", ".nxs"):
        return _detect_usaxs(full_path, path)
    else:
        return _detect_from_folder(path)


def _detect_area_detector(full_path: str, folder: str) -> str:
    """Distinguish SAXS vs WAXS by looking for SAXS-specific metadata keys."""
    try:
        with h5py.File(full_path, "r") as f:
            meta = f.get("/entry/instrument/bluesky/metadata")
            if meta is not None:
                if "pin_ccd_tilt_x" in meta or "pin_ccd_center_x_pixel" in meta:
                    return "SAXS"
                if "waxs_ccd_center_x" in meta or "waxs_ccd_center_x_pixel" in meta:
                    return "WAXS"
    except Exception:
        pass

    # Fallback: parent folder name
    folder_lower = os.path.basename(folder).lower()
    if "_saxs" in folder_lower or folder_lower == "saxs":
        return "SAXS"
    if "_waxs" in folder_lower or folder_lower == "waxs":
        return "WAXS"

    return "Unknown"


def _detect_usaxs(full_path: str, folder: str) -> str:
    """Distinguish Flyscan vs StepScan by looking for step-scan-specific data key."""
    try:
        with h5py.File(full_path, "r") as f:
            # StepScan has angular-position array
            if "/entry/data/a_stage_r" in f:
                return "StepScan"

            # Check explicit scan_type metadata
            meta = f.get("/entry/instrument/bluesky/metadata")
            if meta is not None:
                scan_type = meta.get("scan_type")
                if scan_type is not None:
                    val = scan_type[()]
                    if isinstance(val, (bytes, bytearray)):
                        val = val.decode()
                    val = str(val).lower()
                    if "fly" in val:
                        return "Flyscan"
                    if "step" in val or "uascan" in val:
                        return "StepScan"

            # Flyscans often have a continuous-scan marker
            if "/entry/data/AR_start" in f or "/entry/data/upd2" in f:
                return "Flyscan"
    except Exception:
        pass

    # Fallback: parent folder name
    return _detect_from_folder(folder)


def _detect_from_folder(folder: str) -> str:
    folder_lower = os.path.basename(folder).lower()
    if "flyscan" in folder_lower or "_usaxs" in folder_lower:
        return "Flyscan"
    if "stepscan" in folder_lower:
        return "StepScan"
    if "_saxs" in folder_lower or folder_lower == "saxs":
        return "SAXS"
    if "_waxs" in folder_lower or folder_lower == "waxs":
        return "WAXS"
    return "Unknown"
