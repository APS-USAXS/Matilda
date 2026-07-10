"""Unit tests for the GUI technique detector.

The module is loaded directly from its file to avoid importing the GUI
package __init__ (which requires PySide6/PyQt6).
"""

import importlib.util
import os

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

_MODULE_PATH = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "matilda", "gui", "data_reduction", "technique_detector.py",
)


def _load_detector():
    spec = importlib.util.spec_from_file_location("technique_detector", _MODULE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture
def td():
    mod = _load_detector()
    mod.clear_cache()
    return mod


# ── folder-name fallback ─────────────────────────────────────────────────────

@pytest.mark.parametrize("folder,expected", [
    ("/data/exp_usaxs", "Flyscan"),
    ("/data/FLyscans", "Flyscan"),
    ("/data/StepScans", "StepScan"),
    ("/data/exp_saxs", "SAXS"),
    ("/data/exp_waxs", "WAXS"),
    ("/data/other", "Unknown"),
])
def test_detect_from_folder(td, folder, expected):
    assert td._detect_from_folder(folder) == expected


# ── content-based detection (2.10: metadata read from /entry/Metadata) ──────

def test_detect_saxs_from_entry_metadata(td, tmp_path):
    fp = tmp_path / "scan_0001.hdf"
    with h5py.File(fp, "w") as f:
        f["/entry/Metadata/pin_ccd_tilt_x"] = 0.5
        f["/entry/data/data"] = np.zeros((4, 4))
    assert td.detect_technique(str(tmp_path), "scan_0001.hdf") == "SAXS"


def test_detect_waxs_from_entry_metadata(td, tmp_path):
    fp = tmp_path / "scan_0001.hdf"
    with h5py.File(fp, "w") as f:
        f["/entry/Metadata/waxs_ccd_center_x"] = 100.0
        f["/entry/data/data"] = np.zeros((4, 4))
    assert td.detect_technique(str(tmp_path), "scan_0001.hdf") == "WAXS"


def test_detect_stepscan_from_content(td, tmp_path):
    fp = tmp_path / "scan_0001.h5"
    with h5py.File(fp, "w") as f:
        f["/entry/data/a_stage_r"] = np.zeros(10)
    assert td.detect_technique(str(tmp_path), "scan_0001.h5") == "StepScan"


def test_detect_flyscan_fallback_to_folder(td, tmp_path):
    folder = tmp_path / "exp_usaxs"
    folder.mkdir()
    fp = folder / "scan_0001.h5"
    with h5py.File(fp, "w") as f:
        f["/entry/other"] = 1
    assert td.detect_technique(str(folder), "scan_0001.h5") == "Flyscan"


def test_cache(td, tmp_path):
    folder = tmp_path / "exp_saxs"
    folder.mkdir()
    (folder / "a_0001.hdf").touch()          # not a valid HDF5 -> folder fallback
    first = td.detect_technique(str(folder), "a_0001.hdf")
    assert first == "SAXS"
    assert td.detect_technique(str(folder), "a_0001.hdf") == "SAXS"  # cached
