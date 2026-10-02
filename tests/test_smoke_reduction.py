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

# FX4 counting chain (picoamps), from the 2026-09-26 commissioning session.
# See TestData/FX4Set/description.md for the known quirks of these files —
# notably the unphysical transmission, which is why the FX4 tests below check
# the reduction runs and stays finite rather than checking T < 1.
FX4_FLY_SAMPLE = ("FX4Set/flyscan", "GC_SRM3600_0131.h5")
FX4_FLY_BLANK = ("FX4Set/flyscan", "Blank_0130.h5")
FX4_STEP_SAMPLE = ("FX4Set/usaxs", "GC_SRM3600_0044.h5")
FX4_STEP_BLANK = ("FX4Set/usaxs", "Blank_0043.h5")
FX4_SAXS_SAMPLE = ("FX4Set/saxs", "GC_SRM3600_0044.hdf")
FX4_SAXS_BLANK = ("FX4Set/saxs", "Blank_0043.hdf")
FX4_WAXS_SAMPLE = ("FX4Set/waxs", "GC_SRM3600_0044.hdf")
FX4_WAXS_BLANK = ("FX4Set/waxs", "Blank_0043.hdf")


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


# ── FX4 counting chain ───────────────────────────────────────────────────────

def test_fx4_flyscan_end_to_end(testdata_dir, tmp_path):
    from matilda.convertFlyscan import processFlyscan
    sample_name, blank_name = _stage(testdata_dir, tmp_path,
                                     FX4_FLY_SAMPLE, FX4_FLY_BLANK)
    S = processFlyscan(str(tmp_path), sample_name,
                       blankPath=str(tmp_path), blankFilename=blank_name,
                       recalculateAllData=True)
    assert S["RawData"]["chain"] == "FX4"
    _assert_calibrated(S)
    # transmission comes from trans_*_counts alone: this chain records no
    # trans_*_gain, and defaulting those to 1 must not break the double ratio
    T = S["CalibratedData"]["MeasuredTransmission"]
    assert 0.90 < T < 0.97, f"glassy carbon transmission {T} off the measured 0.936"
    # absolute scale: SRM 3600 sits near 30 cm2/cm3 across its plateau
    Q = np.asarray(S["CalibratedData"]["Q"])
    I = np.asarray(S["CalibratedData"]["Intensity"])
    plateau = I[(Q > 0.01) & (Q < 0.1)]
    assert 10 < np.median(plateau) < 100, f"SRM 3600 plateau at {np.median(plateau)}"


def test_fx4_stepscan_end_to_end(testdata_dir, tmp_path):
    from matilda.convertUSAXS import processStepscan
    sample_name, blank_name = _stage(testdata_dir, tmp_path,
                                     FX4_STEP_SAMPLE, FX4_STEP_BLANK)
    S = processStepscan(str(tmp_path), sample_name,
                        blankPath=str(tmp_path), blankFilename=blank_name,
                        recalculateAllData=True)
    assert S["RawData"]["chain"] == "FX4"
    _assert_calibrated(S)
    # the FX4 branch must not reintroduce a gain or dwell divisor
    assert np.allclose(S["reducedData"]["UPD_gains"], 1.0)

    S2 = processStepscan(str(tmp_path), sample_name,
                         blankPath=str(tmp_path), blankFilename=blank_name,
                         recalculateAllData=False)
    assert len(S2["CalibratedData"]["Q"]) == len(S["CalibratedData"]["Q"])


@pytest.mark.parametrize("sample,blank", [
    pytest.param(FX4_SAXS_SAMPLE, FX4_SAXS_BLANK, id="saxs"),
    pytest.param(FX4_WAXS_SAMPLE, FX4_WAXS_BLANK, id="waxs"),
])
def test_fx4_area_detector_end_to_end(testdata_dir, tmp_path, sample, blank):
    pytest.importorskip("pyFAI")
    from matilda.convertSWAXS import process2Ddata
    sample_name, blank_name = _stage(testdata_dir, tmp_path, sample, blank)
    S = process2Ddata(str(tmp_path), sample_name,
                      blankPath=str(tmp_path), blankFilename=blank_name,
                      recalculateAllData=True)
    assert S["RawData"]["chain"] == "FX4"
    cd = S["CalibratedData"]
    assert len(cd["Q"]) > 50
    assert np.all(np.isfinite(cd["Intensity"]))
    assert 0 < S["calib2DData"]["transmission"] <= 1.5


def test_fx4_and_scaler_frames_are_not_mixed(testdata_dir, tmp_path):
    """An FX4 sample against a scaler blank is counts vs picoamps: refuse."""
    pytest.importorskip("pyFAI")
    from matilda.convertSWAXS import process2Ddata
    sample_name, blank_name = _stage(testdata_dir, tmp_path,
                                     FX4_SAXS_SAMPLE, SAXS_BLANK)
    with pytest.raises(ValueError, match="not comparable"):
        process2Ddata(str(tmp_path), sample_name,
                      blankPath=str(tmp_path), blankFilename=blank_name,
                      recalculateAllData=True)
