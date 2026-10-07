"""Unit tests for matilda.desmearing (Lake/Strobl desmearing)."""

import numpy as np
import pytest

pytest.importorskip("scipy")

from matilda.desmearing import desmearData, extendData


def test_desmear_early_return_is_four_values():
    """Failure path must match the normal 4-value return signature (2.1)."""
    out = desmearData(None, None, None, None, slitLength=0.03)
    assert out == (None, None, None, None)


def test_desmear_empty_input():
    q, i, e, dq = desmearData(np.array([]), np.array([]),
                              np.array([]), np.array([]), slitLength=0.03)
    assert q is None and i is None


def test_desmear_missing_slitlength():
    q, i, e, dq = desmearData(np.ones(10), np.ones(10),
                              np.ones(10), np.ones(10), slitLength=None)
    assert q is None


def test_desmear_synthetic_power_law():
    """Full desmearing run on synthetic slit-smeared power-law data."""
    Q = np.logspace(-4, -0.5, 300)
    Int = 1e6 * Q ** -3 + 100.0
    Err = 0.01 * Int
    dQ = np.gradient(Q)
    DSM_Q, DSM_I, DSM_E, DSM_dQ = desmearData(
        Q, Int.copy(), Err.copy(), dQ,
        slitLength=0.03, MaxNumIter=5, ExtrapMethod="PowerLaw w flat")
    assert DSM_I is not None
    assert len(DSM_I) == len(Q)
    assert np.all(np.isfinite(DSM_I))
    assert np.all(DSM_E >= 0)           # errors are made non-negative
    # Desmearing a decaying power law increases intensity at low Q
    assert DSM_I[0] > Int[0]


@pytest.mark.parametrize("method", ["flat", "Power law", "Porod", "PowerLaw w flat"])
def test_extend_data_methods(method):
    Q = np.logspace(-4, -0.5, 200)
    Int = 1e4 * Q ** -3.5 + 50.0
    Err = 0.01 * Int
    Q2, Int2, Err2, failed = extendData(Q.copy(), Int.copy(), Err.copy(),
                                        0.03, Q[-1] / 2, method)
    assert not failed
    assert len(Q2) > len(Q)                     # data were extended
    assert np.all(np.isfinite(Int2))
    assert np.all(np.diff(Q2[len(Q) - 1:]) > 0)  # extension is increasing in Q


def test_extend_data_rejects_bad_slit():
    with pytest.raises(ValueError):
        extendData(np.ones(10), np.ones(10), np.ones(10), 5.0, 0.5, "flat")


def test_extend_data_rejects_nan_intensity():
    Int = np.ones(10)
    Int[3] = np.nan
    with pytest.raises(ValueError):
        extendData(np.linspace(0.01, 0.1, 10), Int, np.ones(10), 0.03, 0.05, "flat")


# --- truncated-Abel desmearing ------------------------------------------------

SLIT = 0.03


def _model(q):
    """Three-level-ish USAXS model with a flat background."""
    return 1e3 / (1 + (q * 2000) ** 4) + 5.0 / (1 + (q * 120) ** 4) + 1e-6 * q ** -2.0 + 0.05


def _smear(q, slit=SLIT):
    """Accurate slit smearing by direct quadrature: (1/L)*int_0^L I(sqrt(q^2+y^2)) dy."""
    from scipy.integrate import quad
    return np.array([quad(lambda y: _model(np.hypot(qq, y)), 0, slit, limit=200)[0] / slit
                     for qq in q])


def test_abel_recovers_noise_free_model():
    """Noise-free closure test: the inversion must return the true intensity."""
    from matilda.desmearing_methods import desmear_dispatch
    q = np.logspace(-4, np.log10(0.3), 300)
    smr = _smear(q)
    Q, I, E, dQ = desmear_dispatch(q, smr, 1e-3 * smr, np.gradient(q), SLIT,
                                   method="abel", abel_auto_smooth=False,
                                   abel_smooth_w=0.0, abel_num_mc=0)
    ratio = I / _model(Q)
    assert np.sqrt(np.mean((ratio - 1) ** 2)) < 0.01
    assert len(Q) == len(q) and len(E) == len(q) and len(dQ) == len(q)


def test_abel_beats_lake_on_noisy_data():
    """The whole point of the method: it does not amplify noise the way Lake does."""
    from matilda.desmearing_methods import desmear_dispatch
    q = np.logspace(-4, np.log10(0.3), 300)
    smr = _smear(q)
    rng = np.random.default_rng(0)
    noisy = smr * (1 + 0.02 * rng.standard_normal(len(q)))
    err = 0.02 * smr
    _, Ia, Ea, _ = desmear_dispatch(q, noisy, err, np.gradient(q), SLIT,
                                    method="abel", abel_num_mc=5)
    _, Il, _, _ = desmear_dispatch(q, noisy.copy(), err.copy(), np.gradient(q), SLIT,
                                   method="lake")
    rms_abel = np.sqrt(np.mean((Ia / _model(q) - 1) ** 2))
    rms_lake = np.sqrt(np.mean((Il / _model(q) - 1) ** 2))
    assert rms_abel < 0.5 * rms_lake
    assert np.all(np.isfinite(Ea)) and np.all(Ea >= 0)


def test_abel_auto_width_survives_outliers():
    """A few bad points must not collapse the auto smoothing width.

    Under a mean-based chi^2, two or three outliers reach chi^2 = 1 on their own,
    the width drops to ~0 and the result is as noisy as Lake. The criterion is
    therefore the robust (median-based) reduced chi^2.
    """
    from matilda.desmearing_methods import _abel_find_smooth_width
    q = np.logspace(-4, np.log10(0.3), 300)
    smr = _smear(q)
    rng = np.random.default_rng(3)
    y = smr * (1 + 0.02 * rng.standard_normal(len(q)))
    err = 0.02 * smr
    w_clean, _ = _abel_find_smooth_width(q, y, err)
    y[[57, 212, 268]] *= [3.0, 0.3, 4.0]          # three bad points
    w_spiked, _ = _abel_find_smooth_width(q, y, err)
    assert w_clean > 0.03
    assert w_spiked > 0.5 * w_clean


def test_dispatch_lake_matches_desmear_data_exactly():
    """method='lake' must stay bit-for-bit identical to the historical path."""
    from matilda.desmearing_methods import desmear_dispatch
    q = np.logspace(-4, -0.5, 200)
    Int = 1e6 * q ** -3 + 100.0
    err = 0.01 * Int
    dQ = np.gradient(q)
    a = desmearData(q.copy(), Int.copy(), err.copy(), dQ.copy(), slitLength=SLIT,
                    ExtrapMethod="PowerLaw w flat", ExtrapQstart=None, MaxNumIter=20)
    b = desmear_dispatch(q.copy(), Int.copy(), err.copy(), dQ.copy(), SLIT,
                         method="lake", extrap_method="PowerLaw w flat", max_iter=20)
    assert all(np.array_equal(x, y) for x, y in zip(a, b))
