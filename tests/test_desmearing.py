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
