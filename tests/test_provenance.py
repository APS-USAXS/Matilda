"""Reduction provenance (sasprocess) and override handling without raw-file edits.

* saveNXcanSAS writes Sample["provenance"] to <entry>/sasprocess/sasnote and
  readMyNXcanSAS reads it back, so a cached Sample carries MSAXSCorrection,
  calculated_qmin and calculated_mu like a fresh one.
* A thickness override is recorded in the reduced entry only; the raw
  /entry/sample/thickness stays untouched.
"""

import os
import shutil

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from matilda import hdf5code as hc  # noqa: E402

FLY_SAMPLE = ("TestSet/FLyscans", "PU_25C_orig_0061.h5")
FLY_BLANK = ("TestSet/FLyscans", "HeaterBlank_0060.h5")
SAXS_SAMPLE = ("TestSet/SAXS", "PU_25C_orig_0061.hdf")
SAXS_BLANK = ("TestSet/SAXS", "HeaterBlank_0060.hdf")


def _stage(testdata_dir, tmp_path, *entries):
    names = []
    for subdir, name in entries:
        src = os.path.join(testdata_dir, subdir, name)
        if not os.path.isfile(src):
            pytest.skip(f"test file not available: {src}")
        shutil.copy(src, tmp_path / name)
        names.append(name)
    return names


def _sample_dict(provenance=None):
    S = {
        "CalibratedData": {
            "Intensity": np.ones(5), "Q": np.linspace(0.01, 0.1, 5),
            "Error": np.ones(5) * 0.1, "dQ": np.ones(5) * 0.001,
            "units": "[cm2/cm3]", "Kfactor": None, "OmegaFactor": None,
            "blankname": None, "thickness": 1.0, "MeasuredTransmission": 0.8,
        },
        "RawData": {"filename": "scan", "metadata": {}, "sample": {"name": "testsample"}},
        "reducedData": {"Intensity": np.ones(5), "Q": np.linspace(0.01, 0.1, 5),
                        "Error": np.ones(5) * 0.1, "wavelength": 0.59},
    }
    if provenance is not None:
        S["provenance"] = provenance
    return S


# ── unit level ───────────────────────────────────────────────────────────────

def test_sasprocess_written_and_read_back(tmp_path):
    prov = {"technique": "Flyscan", "counting_chain": "FX4", "desmear_method": "abel",
            "abel_smooth_w": 0.031, "qmin_used": 1.2e-4, "qmin_overridden": True,
            "MSAXSCorrection": 0.97, "calculated_qmin": 1.1e-4, "calculated_mu": 3.5,
            "PeakToPeakTransmission": 0.82, "thickness_override": None}
    hc.saveNXcanSAS(_sample_dict(prov), str(tmp_path), "out.h5")
    with h5py.File(tmp_path / "out.h5") as f:
        proc = f["entry/testsample/sasprocess"]
        assert proc.attrs["canSAS_class"] == "SASprocess"
        assert proc.attrs["NX_class"] == "NXprocess"
        note = proc["sasnote"]
        assert note["abel_smooth_w"][()] == pytest.approx(0.031)
        assert note["qmin_used"].attrs["units"] == "1/angstrom"
        assert note["qmin_overridden"][()] == 1          # bool stored as 0/1
        assert "thickness_override" not in note          # None is skipped
        # existing fields untouched
        assert "entry/testsample/sasdata/I" in f

    with h5py.File(tmp_path / "out.h5", "a") as f:
        for g in ("entry/metadata", "entry/instrument", "entry/sample"):
            f.require_group(g)
    S = hc.readMyNXcanSAS(str(tmp_path), "out.h5", isUSAXS=True)
    assert S["provenance"]["technique"] == "Flyscan"
    assert S["provenance"]["counting_chain"] == "FX4"
    assert S["CalibratedData"]["MSAXSCorrection"] == pytest.approx(0.97)
    assert S["CalibratedData"]["calculated_qmin"] == pytest.approx(1.1e-4)
    assert S["CalibratedData"]["calculated_mu"] == pytest.approx(3.5)
    assert S["reducedData"]["PeakToPeakTransmission"] == pytest.approx(0.82)


def test_old_file_without_sasprocess_still_reads(tmp_path):
    hc.saveNXcanSAS(_sample_dict(), str(tmp_path), "out.h5")
    with h5py.File(tmp_path / "out.h5", "a") as f:
        del f["entry/testsample/sasprocess"]
        for g in ("entry/metadata", "entry/instrument", "entry/sample"):
            f.require_group(g)
    S = hc.readMyNXcanSAS(str(tmp_path), "out.h5", isUSAXS=True)
    assert "provenance" not in S
    assert np.allclose(S["CalibratedData"]["Intensity"], 1.0)


def test_run_identifier_and_tdev(tmp_path):
    hc.saveNXcanSAS(_sample_dict(), str(tmp_path), "scan_0007.h5")
    with h5py.File(tmp_path / "scan_0007.h5") as f:
        assert f["entry/testsample/run"][()] == b"scan_0007"
        assert np.isnan(f["entry/testsample/sastransmission_spectrum/Tdev"][0])

    with h5py.File(tmp_path / "uid.h5", "w") as f:
        f.create_dataset("entry/entry_identifier", data="fbb558ae-0e9a")
    hc.saveNXcanSAS(_sample_dict(), str(tmp_path), "uid.h5")
    with h5py.File(tmp_path / "uid.h5") as f:
        assert f["entry/testsample/run"][()] == b"fbb558ae-0e9a"


def test_desmear_dispatch_reports_abel_width():
    from matilda.desmearing_methods import desmear_dispatch
    q = np.geomspace(1e-4, 0.1, 300)
    smr = 1e3 * q ** -3 + 1.0
    info = {}
    out = desmear_dispatch(q, smr, 0.01 * smr, np.gradient(q), slitLength=0.03,
                           method="abel", abel_num_mc=0, info=info)
    assert out[1] is not None
    assert info["desmear_method"] == "abel"
    assert info["abel_auto_smooth"] is True
    assert info["abel_smooth_w"] > 0


# ── end to end ───────────────────────────────────────────────────────────────

def test_flyscan_overrides_recorded_raw_file_untouched(testdata_dir, tmp_path):
    from matilda.convertFlyscan import processFlyscan
    sample_name, blank_name = _stage(testdata_dir, tmp_path, FLY_SAMPLE, FLY_BLANK)
    with h5py.File(tmp_path / sample_name) as f:
        raw_thickness = float(np.ravel(f["entry/sample/thickness"][()])[0])

    S = processFlyscan(str(tmp_path), sample_name, blankPath=str(tmp_path),
                       blankFilename=blank_name, recalculateAllData=True,
                       thickness_override=2.5, qmin_override=2e-4)
    prov = S["provenance"]
    assert prov["technique"] == "Flyscan"
    assert prov["thickness_used"] == pytest.approx(2.5)
    assert prov["thickness_file"] == pytest.approx(raw_thickness)
    assert prov["thickness_overridden"] is True
    assert prov["qmin_overridden"] is True
    assert prov["qmin_override"] == pytest.approx(2e-4)
    assert prov["transmission_overridden"] is False
    assert prov["desmear_method"] == "abel"
    assert prov["abel_smooth_w"] > 0

    with h5py.File(tmp_path / sample_name) as f:
        assert float(np.ravel(f["entry/sample/thickness"][()])[0]) == raw_thickness
        assert "entry/sample/thickness_original" not in f
        entry = hc._nxcansas_entry_name(f, S["RawData"]["sample"]["name"], sample_name)
        assert f[f"entry/{entry}/sassample/thickness"][()] == pytest.approx(2.5)
        assert f[f"entry/{entry}/sasprocess/sasnote/thickness_used"][()] == pytest.approx(2.5)
        assert f"entry/{entry}_SMR/sasprocess" in f

    # cached read: complete CalibratedData
    S2 = processFlyscan(str(tmp_path), sample_name, blankPath=str(tmp_path),
                        blankFilename=blank_name, recalculateAllData=False)
    cd, cd2 = S["CalibratedData"], S2["CalibratedData"]
    for key in ("MSAXSCorrection", "calculated_qmin", "calculated_mu", "MeasuredTransmission"):
        assert cd2[key] == pytest.approx(cd[key]), key
    assert cd2["thickness"] == pytest.approx(2.5)
    assert S2["reducedData"]["FWHM"] == pytest.approx(S["reducedData"]["FWHM"])


def test_flyscan_override_matches_old_numbers(testdata_dir, tmp_path):
    """Dropping the raw-file write must not change the reduced numbers: an
    override still gives the same I(Q) as a file that holds that thickness."""
    from matilda.convertFlyscan import processFlyscan
    a = tmp_path / "a"
    b = tmp_path / "b"
    a.mkdir()
    b.mkdir()
    sa, ba = _stage(testdata_dir, a, FLY_SAMPLE, FLY_BLANK)
    sb, bb = _stage(testdata_dir, b, FLY_SAMPLE, FLY_BLANK)
    with h5py.File(b / sb, "a") as f:
        del f["entry/sample/thickness"]
        f["entry/sample/thickness"] = 2.5
    Sa = processFlyscan(str(a), sa, blankPath=str(a), blankFilename=ba,
                        recalculateAllData=True, thickness_override=2.5)
    Sb = processFlyscan(str(b), sb, blankPath=str(b), blankFilename=bb,
                        recalculateAllData=True)
    assert np.array_equal(Sa["CalibratedData"]["Intensity"], Sb["CalibratedData"]["Intensity"])
    assert Sa["CalibratedData"]["Kfactor"] == Sb["CalibratedData"]["Kfactor"]


def test_saxs_thickness_override_raw_file_untouched(testdata_dir, tmp_path):
    pytest.importorskip("pyFAI")
    from matilda.convertSWAXS import process2Ddata
    sample_name, blank_name = _stage(testdata_dir, tmp_path, SAXS_SAMPLE, SAXS_BLANK)
    with h5py.File(tmp_path / sample_name) as f:
        raw_thickness = float(np.ravel(f["entry/sample/thickness"][()])[0])
    S = process2Ddata(str(tmp_path), sample_name, blankPath=str(tmp_path),
                      blankFilename=blank_name, recalculateAllData=True,
                      thickness_override=raw_thickness * 2)
    assert S["CalibratedData"]["thickness"] == pytest.approx(raw_thickness * 2)
    prov = S["provenance"]
    assert prov["technique"] == "SAXS"
    assert prov["thickness_overridden"] is True
    assert prov["thickness_file"] == pytest.approx(raw_thickness)
    with h5py.File(tmp_path / sample_name) as f:
        assert float(np.ravel(f["entry/sample/thickness"][()])[0]) == raw_thickness
        assert "entry/sample/thickness_original" not in f
