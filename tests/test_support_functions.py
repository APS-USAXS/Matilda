"""Unit tests for matilda.supportFunctions pure helpers."""

import numpy as np
import pytest

pytest.importorskip("scipy")
pytest.importorskip("h5py")

from matilda import supportFunctions as sf


# ── find_crossing_index ──────────────────────────────────────────────────────

def test_find_crossing_index_normal():
    assert sf.find_crossing_index(np.array([0.0, 0.5, 1.0, 2.0]), 0.9) == 2


def test_find_crossing_index_never_crossed_returns_int():
    idx = sf.find_crossing_index(np.array([1.0, 2.0, 3.0]), 10.0)
    assert isinstance(idx, int)          # was a float before the fix
    assert 0 <= idx < 3 or idx == 0


# ── _parse_filename_info ─────────────────────────────────────────────────────

def test_parse_filename_info():
    base, num, ext = sf._parse_filename_info("Scan_001.h5")
    assert (base, num, ext) == ("Scan_", 1, ".h5")


def test_parse_filename_info_no_number():
    assert sf._parse_filename_info("plainname.h5") == (None, None, None)


# ── findProperBlankScan ──────────────────────────────────────────────────────

def test_find_proper_blank_nearest_preceding():
    scan_path = "/data/user/exp_usaxs"
    blanks = [
        ("/data/user/exp_usaxs", "Blank_0040.h5"),
        ("/data/user/exp_usaxs", "Blank_0045.h5"),
        ("/data/user/exp_usaxs", "Blank_0060.h5"),   # after the sample: excluded
    ]
    bp, bf = sf.findProperBlankScan(scan_path, "Sample_0050.h5", blanks)
    assert bf == "Blank_0045.h5"


def test_find_proper_blank_none_found():
    bp, bf = sf.findProperBlankScan("/data/user/exp_usaxs", "Sample_0001.h5",
                                    [("/data/user/exp_usaxs", "Blank_0040.h5")])
    assert bp is None and bf is None


# ── subtract_data ────────────────────────────────────────────────────────────

def test_subtract_data_basic():
    X = np.linspace(0.001, 1, 50)
    Y1 = np.full(50, 5.0)
    Y2 = np.full(50, 2.0)
    E = np.full(50, 0.1)
    x, ydiff, esub, ratio = sf.subtract_data(X, Y1, E, X, Y2.copy(), E)
    assert np.allclose(ydiff, 3.0, rtol=1e-6)
    assert np.allclose(ratio, 2.5, rtol=1e-6)
    assert np.allclose(esub, np.sqrt(2) * 0.1, rtol=1e-6)


def test_subtract_data_y2_min_zero_no_inf():
    """Y2 containing an exact 0 must not produce log(0) = -inf (fixed 2.9)."""
    X = np.linspace(0.001, 1, 50)
    Y1 = np.full(50, 2.0)
    Y2 = np.linspace(0, 1, 50)          # min == 0 exactly
    E = np.full(50, 0.1)
    _, ydiff, _, _ = sf.subtract_data(X, Y1, E, X, Y2.copy(), E)
    assert np.all(np.isfinite(ydiff))


# ── rebin_QRSdata ────────────────────────────────────────────────────────────

def test_rebin_qrsdata_normal():
    Wx = np.linspace(1e-5, 0.3, 3000)
    Wy = Wx ** -1
    Ws = 0.01 * Wy
    q, i, e, dq = sf.rebin_QRSdata(Wx, Wy, Ws, 200)
    assert len(q) < 3000
    assert np.all(np.isfinite(e))
    assert len(q) == len(i) == len(e) == len(dq)


def test_rebin_qrsdata_too_few_highq_points_returns_unchanged():
    """<2 points above Q=0.0002 must not raise IndexError (fixed 2.9)."""
    Wx = np.array([1e-5, 5e-5, 1e-4])
    Wy = np.ones(3)
    Ws = np.ones(3)
    q, i, e, dq = sf.rebin_QRSdata(Wx, Wy, Ws, 200)
    assert len(q) == 3


# ── smooth_r_data ────────────────────────────────────────────────────────────

def _synthetic_scan(n=2000):
    q = np.linspace(1e-5, 0.3, n)
    intensity = 1e3 * np.exp(-q / 0.01) + 1.0
    r_error = 0.01 * intensity
    gidx = np.repeat([0.0, 1.0, 2.0, 3.0, 4.0], n // 5)
    gidx[400] = np.nan          # masked points
    gidx[800] = np.nan
    return q, intensity, r_error, gidx


def test_smooth_r_data_with_range_index():
    """smooth_r_data receives the amplifier RANGE INDEX (0-4), fixed 1.1."""
    q, intensity, r_error, gidx = _synthetic_scan()
    meas_time = np.full(len(q), 5e4)    # 0.05 s per point at 1e6 Hz clock
    out = sf.smooth_r_data(intensity.copy(), q, gidx, r_error.copy(),
                           meas_time, replaceNans=True)
    assert np.all(np.isfinite(out["Intensity"]))
    assert len(out["Intensity"]) == len(q)


def test_smooth_r_data_smoothing_branches():
    """0.025 s/point: ranges 0-1 pass through, ranges 2-4 get smoothed."""
    q, intensity, r_error, gidx = _synthetic_scan()
    meas_time = np.full(len(q), 2.5e4)
    out = sf.smooth_r_data(intensity.copy(), q, gidx, r_error.copy(),
                           meas_time, replaceNans=True)
    assert np.all(np.isfinite(out["Intensity"]))


# ── misc helpers ─────────────────────────────────────────────────────────────

def test_filter_nested_dict():
    d = {"keep": 1, "drop": 2, "nested_drop": {"keep": 3}}
    out = sf.filter_nested_dict(d, ["keep"])
    assert out == {"keep": 1}


def test_check_arrays_same_length_raises():
    with pytest.raises(ValueError):
        sf.check_arrays_same_length(np.ones(3), np.ones(4))


def test_gaussian_peak():
    assert sf.gaussian(0.0, 2.0, 0.0, 1.0) == pytest.approx(2.0)


def test_modified_gauss_peak():
    val = sf.modifiedGauss(np.array([0.0]), 2.0, 0.0, 1.0, 2.0)
    assert val[0] == pytest.approx(2.0)


def test_calculate_pd_fly_returns_both_gain_keys():
    """Regression for the duplicate-'UPD_gains'-key bug (1.1)."""
    import inspect
    src = inspect.getsource(sf.calculatePD_Fly)
    assert '"UPD_gainsIndx":GainsIndx' in src
    assert '"UPD_gains":Gains' in src
    assert "AmpReqGain[len(AmpReqGain)-1]" in src   # 1.2 regression guard
