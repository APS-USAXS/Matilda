"""End-to-end smoke tests: reduce real TestData files through the full chain.

The converters MUTATE their input HDF5 files (cached results are written
back), so every test works on copies in tmp_path.

Flyscan and step-scan chains need only h5py/scipy/matplotlib.
SAXS/WAXS (process2Ddata) additionally needs pyFAI and is skipped when absent.
"""

import os
import shutil

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("h5py")
pytest.importorskip("matplotlib")

FLY_SAMPLE = ("TestSet/FLyscans", "PU_25C_orig_0061.h5")
FLY_BLANK = ("TestSet/FLyscans", "HeaterBlank_0060.h5")
STEP_SAMPLE = ("TestSet/StepScans", "SRM3607_300pts_0097.h5")
STEP_BLANK = ("TestSet/StepScans", "AirBlank_300pts_0110.h5")
SAXS_SAMPLE = ("TestSet/SAXS", "PU_25C_orig_0061.hdf")
SAXS_BLANK = ("TestSet/SAXS", "HeaterBlank_0060.hdf")


def _stage(testdata_dir, tmp_path, *entries):
    """Copy TestData files into tmp_path; skip test if any file is missing."""
    names = []
    for subdir, name in entries:
        src = os.path.join(testdata_dir, subdir, name)
        if not os.path.isfile(src):
            pytest.skip(f"test file not available: {src}")
        shutil.copy(src, tmp_path / name)
        names.append(name)
    return names


def _assert_calibrated(sample):
    cd = sample["CalibratedData"]
    assert cd["Intensity"] is not None, "no calibrated data produced"
    assert cd["Q"] is not None
    assert len(cd["Q"]) == len(cd["Intensity"]) == len(cd["Error"])
    assert len(cd["Q"]) > 50
    assert np.all(np.isfinite(cd["Intensity"]))
    assert np.all(np.array(cd["Q"], dtype=float) > 0)
    assert np.all(np.array(cd["Error"], dtype=float) >= 0)


def test_flyscan_end_to_end(testdata_dir, tmp_path):
    from matilda.convertFlyscan import processFlyscan
    sample_name, blank_name = _stage(testdata_dir, tmp_path, FLY_SAMPLE, FLY_BLANK)
    S = processFlyscan(str(tmp_path), sample_name,
                       blankPath=str(tmp_path), blankFilename=blank_name,
                       recalculateAllData=True)
    _assert_calibrated(S)
    # transmission must be physically sensible (regression for 1.3-class bugs)
    T = S["CalibratedData"].get("MeasuredTransmission")
    assert T is not None and 0 < T <= 1.5

    # cached round trip: second call reads back from the file
    S2 = processFlyscan(str(tmp_path), sample_name,
                        blankPath=str(tmp_path), blankFilename=blank_name,
                        recalculateAllData=False)
    assert S2["CalibratedData"]["Intensity"] is not None
    assert len(S2["CalibratedData"]["Q"]) == len(S["CalibratedData"]["Q"])


def test_stepscan_end_to_end(testdata_dir, tmp_path):
    from matilda.convertUSAXS import processStepscan
    sample_name, blank_name = _stage(testdata_dir, tmp_path, STEP_SAMPLE, STEP_BLANK)
    S = processStepscan(str(tmp_path), sample_name,
                        blankPath=str(tmp_path), blankFilename=blank_name,
                        recalculateAllData=True)
    _assert_calibrated(S)
    # step-scan transmission was reciprocal before fix 1.3 — must be < 1
    T = S["CalibratedData"].get("MeasuredTransmission")
    assert T is not None and 0 < T < 1, f"unphysical transmission {T}"


def test_saxs_end_to_end(testdata_dir, tmp_path):
    pytest.importorskip("pyFAI")
    from matilda.convertSWAXS import process2Ddata
    sample_name, blank_name = _stage(testdata_dir, tmp_path, SAXS_SAMPLE, SAXS_BLANK)
    S = process2Ddata(str(tmp_path), sample_name,
                      blankPath=str(tmp_path), blankFilename=blank_name,
                      recalculateAllData=True)
    cd = S["CalibratedData"]
    assert cd["Intensity"] is not None
    assert len(cd["Q"]) > 50
    assert np.all(np.isfinite(cd["Intensity"]))
    T = S["calib2DData"]["transmission"]
    assert 0 < T <= 1.5
