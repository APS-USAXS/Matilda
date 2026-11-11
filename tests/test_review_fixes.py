"""Regression tests for the fixes from the November 2026 code review.

Each test names the review item it guards.  The converters MUTATE their input
HDF5 files, so every test works on copies in tmp_path (see test_smoke_reduction).
"""

import logging
import os

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("scipy")

from tests.test_smoke_reduction import (  # noqa: E402
    FLY_BLANK,
    FLY_SAMPLE,
    FX4_FLY_BLANK,
    FX4_FLY_SAMPLE,
    FX4_STEP_BLANK,
    FX4_STEP_SAMPLE,
    SAXS_BLANK,
    STEP_BLANK,
    STEP_SAMPLE,
    _stage,
)


def _set_sample_name(filepath, name):
    with h5py.File(filepath, "a") as f:
        del f["entry/sample/name"]
        f["entry/sample/name"] = name


def _reduce_fly(tmp_path, sample, blank, **kw):
    from matilda.convertFlyscan import processFlyscan
    return processFlyscan(str(tmp_path), sample, blankPath=str(tmp_path),
                          blankFilename=blank, recalculateAllData=True, **kw)


# ── 3.1 saveNXcanSAS must never overwrite raw data ───────────────────────────

@pytest.mark.parametrize("name,raw_group,expected_entry", [
    ("sample", "entry/sample", "sample_NXcanSAS"),
    ("metadata", "entry/metadata", "metadata_NXcanSAS"),
    ("", "entry/flyScan", "GC_SRM3600_0131"),
    ("QRS_data", None, "QRS_data_NXcanSAS"),
])
def test_sample_name_colliding_with_raw_group_keeps_raw_data(
        testdata_dir, tmp_path, name, raw_group, expected_entry):
    from matilda.hdf5code import readMyNXcanSAS
    sample, blank = _stage(testdata_dir, tmp_path, FX4_FLY_SAMPLE, FX4_FLY_BLANK)
    fp = tmp_path / sample
    with h5py.File(fp, "r") as f:
        raw_before = {k: f["entry/sample"][k][()] for k in f["entry/sample"]
                      if isinstance(f["entry/sample"][k], h5py.Dataset)}
        thickness_before = f["entry/sample/thickness"][()]
    _set_sample_name(fp, name)

    for _ in range(2):              # a re-reduction must hit the same entry
        S = _reduce_fly(tmp_path, sample, blank)

    with h5py.File(fp, "r") as f:
        # raw data intact
        assert "entry/flyScan/upd_current" in f
        assert f["entry/sample/thickness"][()] == thickness_before
        assert set(raw_before) <= set(f["entry/sample"].keys())
        if raw_group is not None:
            assert raw_group in f
            assert "definition" not in f[raw_group] or raw_group == "entry"
        # reduced data written under the safe name
        assert f"entry/{expected_entry}/sasdata/I" in f
        assert f"entry/{expected_entry}_SMR/sasdata/I" in f
        assert f["entry"].attrs["default"] == expected_entry
        assert "entry/QRS_data/Intensity" in f

    back = readMyNXcanSAS(str(tmp_path), sample, isUSAXS=True)
    assert np.allclose(back["CalibratedData"]["Intensity"], S["CalibratedData"]["Intensity"])


def test_normal_sample_name_location_is_unchanged(testdata_dir, tmp_path):
    """pyirena, DataReporter and Irena read /entry/<sample name>."""
    sample, blank = _stage(testdata_dir, tmp_path, FX4_FLY_SAMPLE, FX4_FLY_BLANK)
    fp = tmp_path / sample
    _set_sample_name(fp, "GlassyCarbon")
    _reduce_fly(tmp_path, sample, blank)
    with h5py.File(fp, "r") as f:
        assert "entry/GlassyCarbon/sasdata/I" in f
        assert "entry/GlassyCarbon_SMR/sasdata/I" in f
        assert f["entry"].attrs["default"] == "GlassyCarbon"


def test_slash_in_sample_name_does_not_nest_groups(testdata_dir, tmp_path):
    sample, blank = _stage(testdata_dir, tmp_path, FX4_FLY_SAMPLE, FX4_FLY_BLANK)
    _set_sample_name(tmp_path / sample, "a/b")
    _reduce_fly(tmp_path, sample, blank)
    with h5py.File(tmp_path / sample, "r") as f:
        assert "entry/a_b/sasdata/I" in f
        assert "entry/a" not in f


def test_missing_timestamp_does_not_fail_the_save(testdata_dir, tmp_path):
    """3.7: timeStamp is informational; its absence must not lose a reduction."""
    from matilda import convertFlyscan
    from matilda.supportFunctions import importFlyscan as real_import

    def import_without_timestamp(path, filename):
        raw = real_import(path, filename)
        raw["metadata"].pop("timeStamp", None)
        return raw

    sample, blank = _stage(testdata_dir, tmp_path, FX4_FLY_SAMPLE, FX4_FLY_BLANK)
    mp = pytest.MonkeyPatch()
    mp.setattr(convertFlyscan, "importFlyscan", import_without_timestamp)
    try:
        S = _reduce_fly(tmp_path, sample, blank)
    finally:
        mp.undo()
    assert S["CalibratedData"]["Intensity"] is not None


# ── 3.2 / 3.3 short fly scans and Qmin override ──────────────────────────────

def test_short_flyscan_after_qmin_override_has_smr_dq(testdata_dir, tmp_path):
    """3.2: 51-800 SMR points used to raise KeyError 'SMR_dQ'."""
    sample, blank = _stage(testdata_dir, tmp_path, FX4_FLY_SAMPLE, FX4_FLY_BLANK)
    S = _reduce_fly(tmp_path, sample, blank, qmin_override=0.2)
    cd = S["CalibratedData"]
    n = len(cd["SMR_Qvec"])
    assert 50 < n <= 800, f"test needs a short SMR curve, got {n} points"
    assert len(cd["SMR_dQ"]) == n
    assert np.all(np.asarray(cd["SMR_dQ"]) > 0)
    assert cd["Intensity"] is not None


def test_qmin_override_when_auto_qmin_fails(testdata_dir, tmp_path):
    """3.3: formatting calculated_qmin=None raised TypeError."""
    sample, blank = _stage(testdata_dir, tmp_path, FLY_SAMPLE, FLY_BLANK)
    S = _reduce_fly(tmp_path, sample, blank, minQMinFindRatio=1e30, qmin_override=1e-3)
    assert S["CalibratedData"]["calculated_qmin"] is None
    assert S["CalibratedData"]["Intensity"] is not None


def test_point_spacing_dq_matches_rebin_convention():
    from matilda.supportFunctions import point_spacing_dQ
    assert np.allclose(point_spacing_dQ([1.0, 2.0, 4.0, 7.0]), [1.0, 1.0, 2.0, 3.0])
    assert np.allclose(point_spacing_dQ([5.0]), [0.0])


# ── 3.4 stale cached reductions are recomputed ───────────────────────────────

def _clear_version_stamps(fp):
    from matilda.hdf5code import REDUCTION_VERSION_ATTR
    with h5py.File(fp, "a") as f:
        def strip(_name, obj):
            if REDUCTION_VERSION_ATTR in obj.attrs:
                del obj.attrs[REDUCTION_VERSION_ATTR]
        f.visititems(strip)


@pytest.mark.parametrize("reducer,pair", [
    ("fly", (FX4_FLY_SAMPLE, FX4_FLY_BLANK)),
    ("step", (FX4_STEP_SAMPLE, FX4_STEP_BLANK)),
])
def test_unstamped_cached_reduction_is_recomputed(testdata_dir, tmp_path, caplog, reducer, pair):
    from matilda.convertFlyscan import processFlyscan
    from matilda.convertUSAXS import processStepscan
    from matilda.hdf5code import REDUCTION_VERSION, REDUCTION_VERSION_ATTR
    process = processFlyscan if reducer == "fly" else processStepscan
    sample, blank = _stage(testdata_dir, tmp_path, *pair)
    kw = dict(blankPath=str(tmp_path), blankFilename=blank)
    process(str(tmp_path), sample, recalculateAllData=True, **kw)

    with caplog.at_level(logging.INFO):
        caplog.clear()
        process(str(tmp_path), sample, recalculateAllData=False, **kw)
    assert "Using existing processed data" in caplog.text      # current: reused

    _clear_version_stamps(tmp_path / sample)                   # pre-fix file
    with caplog.at_level(logging.INFO):
        caplog.clear()
        S = process(str(tmp_path), sample, recalculateAllData=False, **kw)
    assert "Using existing processed data" not in caplog.text
    assert "recalculating" in caplog.text
    assert "MeasuredTransmission" in S["CalibratedData"]        # a fresh reduction
    with h5py.File(tmp_path / sample, "r") as f:
        assert f["entry/QRS_data"].attrs[REDUCTION_VERSION_ATTR] == REDUCTION_VERSION


def test_stale_stepscan_blank_cache_is_recalculated(testdata_dir, tmp_path):
    from matilda.convertUSAXS import getBlankStepscan
    from matilda.supportFunctions import BLANK_CACHE_VERSION
    (blank,) = _stage(testdata_dir, tmp_path, STEP_BLANK)
    first = getBlankStepscan(str(tmp_path), blank, recalculateAllData=True)
    assert int(np.ravel(first["BlankCacheVersion"])[0]) == BLANK_CACHE_VERSION
    with h5py.File(tmp_path / blank, "a") as f:
        del f["entry/blankData/BlankCacheVersion"]
        f["entry/blankData/Intensity"][...] = 0.0          # poison the stale cache
    again = getBlankStepscan(str(tmp_path), blank)
    assert np.allclose(again["Intensity"], first["Intensity"])


# ── 3.5 USAXS sample and blank must share a counting chain ───────────────────

def test_fx4_flyscan_refuses_scaler_blank(testdata_dir, tmp_path):
    sample, blank = _stage(testdata_dir, tmp_path, FX4_FLY_SAMPLE, FLY_BLANK)
    with pytest.raises(ValueError, match="not comparable"):
        _reduce_fly(tmp_path, sample, blank)


def test_scaler_stepscan_refuses_fx4_blank(testdata_dir, tmp_path):
    from matilda.convertUSAXS import processStepscan
    sample, blank = _stage(testdata_dir, tmp_path, STEP_SAMPLE, FX4_STEP_BLANK)
    with pytest.raises(ValueError, match="not comparable"):
        processStepscan(str(tmp_path), sample, blankPath=str(tmp_path),
                        blankFilename=blank, recalculateAllData=True)


# ── 3.7 SAXS/WAXS blanks are cached like any other file ──────────────────────

def test_saxs_blank_is_not_reintegrated_every_cycle(testdata_dir, tmp_path, caplog):
    pytest.importorskip("pyFAI")
    from matilda.convertSWAXS import process2Ddata
    (blank,) = _stage(testdata_dir, tmp_path, SAXS_BLANK)
    # the daemon passes a blank even when the file itself is a blank
    process2Ddata(str(tmp_path), blank, blankPath=str(tmp_path), blankFilename=blank,
                  recalculateAllData=True)
    with caplog.at_level(logging.INFO):
        caplog.clear()
        process2Ddata(str(tmp_path), blank, blankPath=str(tmp_path), blankFilename=blank)
    assert "Using existing processed data" in caplog.text


# ── 3.7 / 5.3 desmeared error bars are reproducible ──────────────────────────

def test_abel_monte_carlo_errors_are_reproducible():
    from matilda.desmearing_methods import desmear_dispatch
    rng = np.random.default_rng(3)
    q = np.geomspace(1e-4, 0.1, 300)
    I_true = 1e3 / (1 + (q / 2e-3) ** 4) + 0.1
    err = 0.03 * I_true
    I = I_true + rng.normal(0, err)
    dq = np.full_like(q, 1e-5)
    runs = [desmear_dispatch(q, I, err, dq, slitLength=0.03, method="abel") for _ in range(2)]
    assert np.array_equal(runs[0][1], runs[1][1])
    assert np.array_equal(runs[0][2], runs[1][2])
    # a different curve gets different draws (seed comes from the data)
    other = desmear_dispatch(q, I * 1.001, err, dq, slitLength=0.03, method="abel")
    assert not np.allclose(other[2] / 1.001, runs[0][2], rtol=1e-6)


# ── 3.7 smooth_r_data survives a bad dwell time ──────────────────────────────

def test_smooth_r_data_with_zero_channel_time():
    from matilda.supportFunctions import smooth_r_data
    n = 200
    q = np.linspace(1e-4, 1e-2, n)
    I = 1e3 * np.exp(-q * 300) + 1.0
    t = np.full(n, 0.01)
    t[150] = 0.0
    t[160] = np.nan
    out = smooth_r_data(I, q, np.zeros(n), 0.01 * I, t, time_base=1.0,
                        fx4_noise_model=True)
    assert np.all(np.isfinite(out["Intensity"]))
    assert out["Intensity"][150] == pytest.approx(I[150])


# ── 3.6 external-tool retry policy ───────────────────────────────────────────

def test_failed_jobs_are_retried_a_bounded_number_of_times():
    from matilda.matilda import _MAX_TOOL_ATTEMPTS, _job_settled, _record_job
    tracked = {}
    key = ("p", "f.h5")
    for _ in range(_MAX_TOOL_ATTEMPTS - 1):
        _record_job(tracked, key, False, "test")
        assert not _job_settled(tracked, key)
    _record_job(tracked, key, False, "test")
    assert _job_settled(tracked, key)

    done = ("p", "g.h5")
    assert not _job_settled(tracked, done)
    _record_job(tracked, done, True, "test")
    assert _job_settled(tracked, done)


def test_merge_waits_for_reduced_partner(tmp_path, monkeypatch):
    import matilda.matilda as m
    usaxs = tmp_path / "exp_usaxs"
    saxs = tmp_path / "exp_saxs"
    usaxs.mkdir()
    saxs.mkdir()
    for p in (usaxs / "S1_x_0001.h5", saxs / "S1_x_0001.hdf"):
        with h5py.File(p, "w") as f:
            f.create_group("entry")
    (tmp_path / m._MERGE_CONFIG_FILENAME).write_text("{}")
    calls = []
    monkeypatch.setattr(m, "_runMergeData", lambda *a: calls.append(a) or True)
    monkeypatch.setattr(m, "_runPyirenaAnalysis", lambda *a: None)
    merged = {}
    m._checkAndRunMerge([(str(usaxs), "S1_x_0001.h5")], "USAXS", merged, {})
    assert calls == [] and merged == {}           # postponed, not recorded

    for p in (usaxs / "S1_x_0001.h5", saxs / "S1_x_0001.hdf"):
        with h5py.File(p, "a") as f:
            g = f.create_group("entry/S1")
            g.attrs["canSAS_class"] = "SASentry"
            g.attrs["NX_class"] = "NXsubentry"
            g["definition"] = "NXcanSAS"
    m._checkAndRunMerge([(str(saxs), "S1_x_0001.hdf")], "SAXS", merged, {})
    assert len(calls) == 1


# ── 5.4 batch reprocessing sends step scans to the step-scan reducer ─────────

def test_process_usaxs_folder_splits_fly_and_step(testdata_dir, tmp_path, monkeypatch):
    import matilda.matilda as m
    usaxs = tmp_path / "exp_usaxs"
    usaxs.mkdir()
    _stage(testdata_dir, usaxs, FX4_FLY_SAMPLE, FX4_FLY_BLANK, FX4_STEP_SAMPLE, FX4_STEP_BLANK)
    seen = {}
    monkeypatch.setattr(m, "processFlyscans",
                        lambda s, b, **kw: seen.setdefault("fly", (s, b)))
    monkeypatch.setattr(m, "processStepscans",
                        lambda s, b, **kw: seen.setdefault("step", (s, b)))
    m.processUSAXSFolder(str(tmp_path))
    names = {k: ([os.path.basename(f) for _, f in s], [f for _, f in b])
             for k, (s, b) in seen.items()}
    assert names["fly"] == (["Blank_0130.h5", "GC_SRM3600_0131.h5"], ["Blank_0130.h5"])
    assert names["step"] == (["Blank_0043.h5", "GC_SRM3600_0044.h5"], ["Blank_0043.h5"])
