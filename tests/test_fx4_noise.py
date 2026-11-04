"""FX4 fly-scan noise model, smoothing, rebinning and SMR error propagation.

Background and the measurements behind every constant:
docs/fx4-flyscan-uncertainties.md.  The last test is the acceptance test that
the whole chain quotes the right error: it reduces the even and the odd raw
points of a real FX4 fly scan independently and checks that their difference
is as large as their quoted errors say it should be.
"""

import inspect
import os
import shutil

import numpy as np
import pytest

from matilda import fx4support as fx4
import matilda.supportFunctions as sf

pytest.importorskip("h5py")
pytest.importorskip("scipy")


# ── per-interval noise model ─────────────────────────────────────────────────

def test_current_sigma_floor_counts_values_per_read():
    """With no net signal only the floor remains: sigma / sqrt(N * vpr)."""
    sigma = np.full(4, 7.0)
    mean = np.full(4, 2.0)
    total = mean * 190
    out = fx4.current_sigma(sigma, total, mean, 10, np.zeros(4), np.full(4, 0.019))
    assert np.allclose(out, 7.0 / np.sqrt(1900))


def test_current_sigma_shot_term_scales_with_net_current():
    """With no recorded spread only the shot-like term remains."""
    net = np.array([1.0, 100.0, 1e4])
    t = np.full(3, 0.02)
    out = fx4.current_sigma(np.zeros(3), net * 100, net, 10, net, t)
    assert np.allclose(out, np.sqrt(fx4.FX4_SHOT_CHARGE_PC * net / t))


def test_current_sigma_negative_net_gives_no_shot_term():
    out = fx4.current_sigma(np.full(2, 7.0), np.full(2, 200.0), np.full(2, 2.0), 10,
                            np.array([-3.0, -0.1]), np.full(2, 0.02))
    assert np.allclose(out, 7.0 / np.sqrt(1000))


def test_add_relative_floor():
    out = fx4.add_relative_floor(np.array([1e5, 1.0]), np.array([0.0, 0.3]))
    assert out[0] == pytest.approx(1e5 * fx4.FX4_RELATIVE_ERROR_FLOOR)
    assert out[1] == pytest.approx(0.3, rel=1e-5)


# ── smoothing ────────────────────────────────────────────────────────────────

def _flat_noisy_scan(n=3000, level=0.2, noise=1.0, seed=1):
    """Weak signal well inside the noise: many points are negative.

    Smoothing starts at Q = 3e-4 and the lift looks at the last 80% of points,
    so Q starts at 1e-4 to keep that 80% inside the smoothed range, as the
    rocking-curve peak does on a real scan.
    """
    rng = np.random.default_rng(seed)
    q = np.geomspace(1e-4, 0.3, n)
    intensity = level + noise * rng.standard_normal(n)
    return q, intensity, np.full(n, noise), np.full(n, 4.0), np.full(n, 0.019)


def test_fx4_smoothing_reduces_error_by_sqrt_window():
    q, intensity, err, gidx, t = _flat_noisy_scan()
    out = sf.smooth_r_data(intensity, q, gidx, err, t, time_base=1.0, fx4_noise_model=True)
    window = out["SmoothWindow"]
    smoothed = window > 1
    assert smoothed.sum() > 0.8 * len(q)
    assert np.allclose(out["Error"][smoothed], 1.0 / np.sqrt(window[smoothed]))
    assert np.all(window[~smoothed] == 1)


def test_fx4_smoothing_keeps_the_mean_of_weak_data():
    """Linear averaging: no ln(I) bias and no NaN-replaced negative points."""
    q, intensity, err, gidx, t = _flat_noisy_scan(level=0.6, noise=1.0)
    m = q > 1e-3
    assert np.mean(intensity[m] <= 0) > 0.2           # plenty of negative points
    out = sf.smooth_r_data(intensity, q, gidx, err, t, time_base=1.0, fx4_noise_model=True)
    assert np.all(out["Intensity"][m] > 0)             # smoothed: no lift needed
    assert out["Intensity"][m].mean() == pytest.approx(intensity[m].mean(), abs=0.01)


def test_fx4_lift_is_set_by_smoothed_not_raw_minimum():
    """After smoothing, the lift follows the background, not one noise spike."""
    q, intensity, err, gidx, t = _flat_noisy_scan(level=-0.1, noise=1.0)
    m = q > 1e-3
    raw_lift = 1.1 * abs(np.nanmin(intensity[int(0.2 * len(q)):-1]))
    out = sf.smooth_r_data(intensity, q, gidx, err, t, time_base=1.0, fx4_noise_model=True)
    applied = out["Intensity"][m].mean() - intensity[m].mean()
    assert 0 < applied < 0.3 * raw_lift
    assert np.all(out["Intensity"][int(0.2 * len(q)):-1] > 0)


def test_scaler_smoothing_unchanged_apart_from_new_window_key():
    q, intensity, err, gidx, _ = _flat_noisy_scan(level=50.0, noise=1.0)
    meas_time = np.full(len(q), 2.5e4)
    out = sf.smooth_r_data(intensity.copy(), q, gidx, err.copy(), meas_time)
    assert np.array_equal(out["Error"], err)          # error untouched
    assert np.all(out["SmoothWindow"] == 1)


# ── rebinning ────────────────────────────────────────────────────────────────

def test_rebin_counts_only_independent_points():
    q = np.linspace(0.0002, 0.05, 4000)
    i = 1e3 * q ** -2
    e = 0.1 * i
    _, _, e_rms, _ = sf.rebin_QRSdata(q.copy(), i.copy(), e.copy(), 100)
    _, _, e_w1, _ = sf.rebin_QRSdata(q.copy(), i.copy(), e.copy(), 100,
                                     smooth_window=np.ones(4000))
    _, _, e_wide, _ = sf.rebin_QRSdata(q.copy(), i.copy(), e.copy(), 100,
                                       smooth_window=np.full(4000, 1e4))
    # window 1: independent points -> error of the mean, well below the RMS
    assert np.nanmax(e_rms / e_w1) > 3
    # window wider than any bin: the bin holds one independent point -> RMS
    assert np.allclose(e_wide, e_rms)


# ── SMR error does not get multiplied by the peak-to-peak transmission ──────

def test_smr_error_not_scaled_by_peak_to_peak_transmission():
    src = inspect.getsource(sf.calibrateAndSubtractFlyscan)
    code = [line.split("#")[0] for line in src.splitlines()]
    assert not any("SMR_Error" in line and "PeakToPeakTransmission" in line
                   and "=" in line for line in code)


# ── cached blanks from an older reduction are recomputed ─────────────────────

def test_stale_blank_cache_is_recalculated(tmp_path):
    import h5py
    from tests.test_fx4support import _synthetic_fx4_flyscan
    _synthetic_fx4_flyscan(tmp_path / "Blank_0001.h5")
    first = sf.getBlankFlyscan(str(tmp_path), "Blank_0001.h5", recalculateAllData=True)
    assert int(np.ravel(first["BlankCacheVersion"])[0]) == sf.BLANK_CACHE_VERSION
    with h5py.File(tmp_path / "Blank_0001.h5", "a") as f:
        del f["entry/blankData/BlankCacheVersion"]
        f["entry/blankData/Intensity"][...] = 0.0          # poison the stale cache
    again = sf.getBlankFlyscan(str(tmp_path), "Blank_0001.h5")
    assert np.allclose(again["Intensity"], first["Intensity"])
    reused = sf.getBlankFlyscan(str(tmp_path), "Blank_0001.h5")
    assert np.allclose(reused["Intensity"], first["Intensity"])


# ── acceptance: quoted SMR error matches the real noise ──────────────────────

def _halves_ratio(monkeypatch, folder, sample, blank, bands):
    """true/quoted SMR error per Q band from even/odd halves of the raw points.

    Both halves use the full scan's beam centre and peak maximum; otherwise
    each half fits its own peak and that difference masquerades as noise.
    """
    import matilda.convertFlyscan as cf
    from scipy.interpolate import CubicSpline

    orig_import = sf.importFlyscan
    orig_bcc = sf.beamCenterCorrection
    state = {"slice": None, "fixed": {}}

    def importer(path, filename):
        raw = orig_import(path, filename)
        if state["slice"] is not None:
            n = len(raw["UPD_array"])
            for key, value in list(raw.items()):
                if isinstance(value, np.ndarray) and value.ndim == 1 and len(value) == n:
                    raw[key] = value[state["slice"]]
        return raw

    def centre(data, useGauss=1, isBlank=False):
        res = orig_bcc(data, useGauss=useGauss, isBlank=isBlank)
        if state["slice"] is None:
            state["fixed"][isBlank] = (res["Center"], res["Maximum"])
        else:
            x0, peak = state["fixed"][isBlank]
            res["Center"], res["Maximum"] = x0, peak
            res["Q"] = -4 * np.pi * np.sin(np.radians(data["RawData"]["ARangles"] - x0) / 2) / res["wavelength"]
        return res

    for mod in (sf, cf):
        monkeypatch.setattr(mod, "importFlyscan", importer)
        monkeypatch.setattr(mod, "beamCenterCorrection", centre)
    monkeypatch.setattr(cf, "desmear_dispatch", lambda Q, I, E, dQ, **kw: (Q, I, E, dQ))

    def reduce(sl, tag):
        state["slice"] = sl
        work = folder / tag
        work.mkdir()
        for name in (sample, blank):
            shutil.copy(folder / name, work / name)
        out = cf.processFlyscan(str(work), sample, blankPath=str(work),
                                blankFilename=blank, recalculateAllData=True)
        c = out["CalibratedData"]
        return np.asarray(c["SMR_Qvec"]), np.asarray(c["SMR_Int"]), np.asarray(c["SMR_Error"])

    reduce(None, "full")
    qa, ia, ea = reduce(slice(0, None, 2), "even")
    qb, ib, eb = reduce(slice(1, None, 2), "odd")
    ok = ib > 0
    spline = CubicSpline(np.log(qb[ok]), np.log(ib[ok]))
    diff = np.log(ia) - spline(np.log(qa))
    quoted = np.sqrt((ea / ia) ** 2 + np.interp(np.log(qa), np.log(qb), eb / np.abs(ib)) ** 2)
    ratios = []
    for lo, hi in bands:
        m = (qa > lo) & (qa < hi)
        d = diff[m]
        robust_sd = 1.4826 * np.median(np.abs(d - np.median(d)))
        ratios.append(robust_sd / np.median(quoted[m]))
    return ratios


def test_fx4_smr_error_matches_noise_on_real_flyscan(testdata_dir, tmp_path, monkeypatch):
    src = os.path.join(testdata_dir, "FX4Set", "flyscan")
    for name in ("GC_SRM3600_0131.h5", "Blank_0130.h5"):
        shutil.copy(os.path.join(src, name), tmp_path / name)
    bands = [(1e-3, 3e-3), (3e-3, 1e-2), (1e-2, 3e-2), (3e-2, 0.1), (0.1, 0.3)]
    ratios = _halves_ratio(monkeypatch, tmp_path, "GC_SRM3600_0131.h5", "Blank_0130.h5", bands)
    # Before the noise model this pair read 0.6-1.2 only because two errors
    # cancelled; the weak GC_0036 pair read 0.14-0.27.  Now ~0.7-1.1 here.
    for (lo, hi), r in zip(bands, ratios):
        assert 0.55 < r < 1.6, f"true/quoted {r:.2f} for {lo:g}-{hi:g} 1/A"
