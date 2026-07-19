"""
desmearing_methods.py
=====================
Selectable slit-desmearing methods for USAXS, with a single dispatcher that has
the same signature/return as ``desmearing.desmearData``. Lake stays the default
and is untouched (this module calls it); the new options are Gaussian-process
(GP / "Huang") desmearing, which suppress noise and stay physical (smooth,
positive) on weak / near-zero / over-subtracted data where Lake produces spikes.

Public entry point
------------------
    desmear_dispatch(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength, *,
                     method='lake', length_scale_decades=0.5, kernel='matern32',
                     sigma_log=4.0, extrap_method='PowerLaw w flat',
                     extrap_qstart=None, max_iter=20)
        -> (DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ)

Methods
-------
    'lake'          : the existing Matilda Lake iterative desmearing (unchanged).
    'gp'            : Bayesian GP desmearing (multiplicative model, credibility
                      band -> DSM_Error). kernel 'matern32' (default, robust on
                      featured data) or 'rbf' (very smooth, featureless data).
    'gp_lake_mean'  : GP anchored on a *smoothed* Lake curve (robust low-q edge).

Numpy-only; no new dependencies. Developed and validated in the 'Better
desmearing' sandbox (see that project's PLAN / FINDINGS).
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------------------
# Slit forward operator (linear matrix) on an extended log-q working grid
# ---------------------------------------------------------------------------
class _SlitOperator:
    def __init__(self, q_meas, q_work, M, slit_length):
        self.q_meas = q_meas
        self.q_work = q_work
        self.M = M
        self.slit_length = slit_length

    @classmethod
    def build(cls, q_meas, slit_length, n_extra=40, n_slit=200):
        q_meas = np.asarray(q_meas, float)
        L = float(slit_length)
        qmax = q_meas[-1]
        q_top = np.sqrt(qmax ** 2 + L ** 2)
        if n_extra > 0 and q_top > qmax:
            ext = np.exp(np.linspace(np.log(qmax), np.log(q_top), n_extra + 1)[1:])
            q_work = np.concatenate([q_meas, ext])
        else:
            q_work = q_meas.copy()
        logqw = np.log(q_work)
        u = np.linspace(0.0, L, n_slit)
        w = np.full(n_slit, (u[1] - u[0]))
        w[0] *= 0.5
        w[-1] *= 0.5
        w /= L
        Nm, Nw = len(q_meas), len(q_work)
        M = np.zeros((Nm, Nw))
        for i in range(Nm):
            arg = np.sqrt(q_meas[i] ** 2 + u ** 2)
            arg = np.clip(arg, q_work[0], q_work[-1])
            idx = np.clip(np.searchsorted(logqw, np.log(arg)) - 1, 0, Nw - 2)
            x0, x1 = logqw[idx], logqw[idx + 1]
            t = (np.log(arg) - x0) / (x1 - x0)
            np.add.at(M[i], idx, w * (1.0 - t))
            np.add.at(M[i], idx + 1, w * t)
        return cls(q_meas, q_work, M, L)


def _kernel(x, ell, sigma, kind):
    d = np.abs(x[:, None] - x[None, :]) / ell
    if kind == "rbf":
        K = np.exp(-0.5 * d ** 2)
    elif kind == "matern32":
        s3 = np.sqrt(3.0) * d
        K = (1.0 + s3) * np.exp(-s3)
    else:
        raise ValueError(f"unknown kernel {kind!r}")
    return (sigma ** 2) * K


def _smooth_positive_trend(q, y, sig, ls_decades):
    """Smooth, strictly positive base trend: log-space (geometric-mean) moving
    average after clipping to a noise-level floor. Preserves power laws over many
    decades and tolerates near-zero / negative points."""
    q = np.asarray(q, float); y = np.asarray(y, float); sig = np.asarray(sig, float)
    n = len(q)
    e = sig[np.isfinite(sig) & (sig > 0)]
    floor = float(np.median(e)) if e.size else 1e-6 * np.nanmax(np.abs(y))
    floor = max(floor, 1e-300)
    ly = np.log(np.clip(y, floor, None))
    decades = np.log10(q[-1] / q[0]) if q[-1] > q[0] else 1.0
    ppd = n / max(decades, 1e-6)
    win = max(3, int(round(0.5 * ls_decades * ppd)))
    if win % 2 == 0:
        win += 1
    pad = win // 2
    lyp = np.pad(ly, pad, mode="reflect")
    sm = np.convolve(lyp, np.ones(win) / win, mode="valid")[:n]
    return np.exp(sm)


def _trend_on_work(q, base_measured, q_work, tail_order=2):
    Nm = len(q)
    mu = np.empty(len(q_work))
    mu[:Nm] = base_measured
    if len(q_work) > Nm:
        k = tail_order + 1
        p = np.polyfit(np.log(q[-k:]), np.log(base_measured[-k:]), 1)[0]
        mu[Nm:] = base_measured[-1] * (q_work[Nm:] / q[-1]) ** p
    return np.clip(mu, 1e-300, None)


def _smooth_curve(q, I, decades=0.3):
    """Robust smoothing of a curve over a window of ``decades`` in q, by a
    rolling **median** in linear intensity. Used by the 'lake_smooth' method:
    it keeps Lake's null-space-safe shape (no coherent oscillation), is robust to
    Lake's up/down spikes, and — unlike geometric-mean (log) smoothing — does not
    sit systematically below the data (log-averaging is biased low for noisy
    data; the median is not)."""
    q = np.asarray(q, float); I = np.asarray(I, float)
    n = len(q)
    span = np.log10(q[-1] / q[0]) if q[-1] > q[0] else 1.0
    ppd = n / max(span, 1e-6)
    win = max(3, int(round(decades * ppd)))
    if win % 2 == 0:
        win += 1
    pad = win // 2
    Ip = np.pad(I, pad, mode="reflect")
    out = np.empty(n)
    for i in range(n):
        out[i] = np.median(Ip[i:i + win])
    return out


def _smooth_loglog(q, I, win_frac=0.04):
    """Smooth a (possibly spiky, possibly non-positive) curve in log-log for use
    as a GP prior mean. Drops non-positive points, geometric-mean smooths."""
    q = np.asarray(q, float); I = np.asarray(I, float)
    m = np.isfinite(q) & (q > 0) & np.isfinite(I) & (I > 0)
    q, I = q[m], I[m]
    n = len(q)
    if n < 5:
        return q, I
    win = max(3, int(win_frac * n))
    if win % 2 == 0:
        win += 1
    pad = win // 2
    lp = np.pad(np.log(I), pad, mode="reflect")
    sm = np.convolve(lp, np.ones(win) / win, mode="valid")[:n]
    return q, np.exp(sm)


# ---------------------------------------------------------------------------
# GP desmearing core (multiplicative model, damped Gauss-Newton, credibility)
# ---------------------------------------------------------------------------
def _gp_core(q, y, err, slit_length, length_scale_decades=0.5, sigma_log=4.0,
             kernel="matern32", err_floor_frac=0.01, rel_err_floor=1e-2,
             n_sigma_band=2.0, max_iter=30, tol=1e-4, n_extra=40, n_slit=200,
             jitter=1e-8, resolution_aware=True, res_alpha=1.0,
             low_q_guard=True, prior_mean_q=None, prior_mean_I=None):
    q = np.asarray(q, float); y = np.asarray(y, float); err = np.asarray(err, float)
    m = np.isfinite(q) & (q > 0) & np.isfinite(y)
    q, y = q[m], y[m]
    err = err[m] if err.shape == m.shape else np.full(m.sum(), np.nan)
    L = float(slit_length)

    # Absolute + relative (>=1%, the instrument floor) error weighting.
    sig = err.astype(float)
    good = np.isfinite(sig) & (sig > 0)
    if good.any():
        med = float(np.median(sig[good]))
        sig = np.where(good, sig, med)
    else:
        med = float(np.median(np.abs(np.diff(y))) or (np.median(np.abs(y)) + 1e-30))
        sig = np.full_like(y, med)
    sig = np.maximum(sig, err_floor_frac * med)
    sig = np.maximum(sig, rel_err_floor * np.abs(y))
    Sinv = 1.0 / sig ** 2

    if prior_mean_I is not None:
        pmq = np.asarray(prior_mean_q, float); pmI = np.asarray(prior_mean_I, float)
        base = np.exp(np.interp(np.log(q), np.log(pmq), np.log(np.clip(pmI, 1e-300, None))))
    else:
        base = _smooth_positive_trend(q, y, sig, length_scale_decades)

    op = _SlitOperator.build(q, L, n_extra=n_extra, n_slit=n_slit)
    M = op.M; q_work = op.q_work; Nm, Nw = M.shape
    mu0 = _trend_on_work(q, base, q_work); lnmu0 = np.log(mu0)

    # GP prior on s over x = ln q. Base (stationary) length scale:
    xln = np.log(q_work); ell0 = length_scale_decades * np.log(10.0)
    if resolution_aware:
        # Slit-limited resolution: a point at q is smeared over [q, sqrt(q^2+L^2)],
        # i.e. a half-width in ln q of  W(q) = 0.5*ln(1 + (L/q)^2). At low q
        # (q << L, flat plateau) W is large — the data cannot resolve structure
        # there. Enforce a *local* length scale >= that resolution by warping the
        # coordinate u = integral dx / ell_local, then using a unit-length kernel
        # in u. This gives a non-stationary GP that smooths the unconstrained
        # low-q null space (no ringing) while keeping full resolution at high q.
        W = 0.5 * np.log1p((L / q_work) ** 2)
        ell_local = np.maximum(ell0, res_alpha * W)
        inv = 1.0 / ell_local
        u = np.concatenate([[0.0], np.cumsum(0.5 * (inv[1:] + inv[:-1]) * np.diff(xln))])
        K = _kernel(u, 1.0, sigma_log, kernel) + jitter * np.eye(Nw)
    else:
        K = _kernel(xln, ell0, sigma_log, kernel) + jitter * np.eye(Nw)
    Kinv = np.linalg.inv(K)
    smax = np.log(1e8)

    def objective(s):
        x = np.exp(lnmu0 + s); f = M @ x; r = (y - f) / sig
        return float(r @ r + s @ (Kinv @ s)), x, f

    s = np.zeros(Nw)
    obj, x, f = objective(s)
    for _ in range(max_iter):
        J = M * x[None, :]
        A = (J.T * Sinv) @ J + Kinv
        b = (J.T * Sinv) @ (y - f) - Kinv @ s
        delta = np.linalg.solve(A, b)
        alpha = 1.0; improved = False
        for _ls in range(30):
            s_try = np.clip(s + alpha * delta, -smax, smax)
            obj_try, x_try, f_try = objective(s_try)
            if obj_try < obj:
                improved = True; break
            alpha *= 0.5
        if not improved:
            break
        step = float(np.max(np.abs(s_try - s)))
        s, obj, x, f = s_try, obj_try, x_try, f_try
        if step < tol:
            break

    J = M * x[None, :]
    A = (J.T * Sinv) @ J + Kinv
    var = np.clip(np.diag(np.linalg.inv(A)), 0.0, None)
    sd = np.sqrt(var)
    I = x[:Nm]; sdm = sd[:Nm]
    err_lin = I * sdm
    lower = I * np.exp(-n_sigma_band * sdm)
    upper = I * np.exp(+n_sigma_band * sdm)

    if low_q_guard:
        # Weak, signal-limited flat low-q data (a Guinier plateau) leaves a broad
        # low-q block in the operator's near-null space when slit >> q_min: many
        # low-q smeared points integrate the same range, so the fit is free to
        # ramp the desmeared curve down away from the flat plateau (an artefact —
        # desmearing a flat region should leave it flat). Guard: over the
        # under-constrained low-q block (small total column weight in M), do not
        # let the desmeared intensity fall below the local *smeared* plateau
        # level, which is the correct value for a flat region. A genuine rising
        # low-q power law (resolution-limited, strong signal) is left untouched —
        # there the desmeared exceeds the smeared, so the clamp never binds.
        wcol = M.sum(axis=0)[:Nm]
        thr = 0.5 * np.median(wcol)
        k = 0
        while k < Nm and wcol[k] < thr:
            k += 1
        if k >= 1:
            ref = float(np.median(y[:k]))        # smeared plateau level of the weak block
            I = I.copy()
            i = 0
            while i < Nm and I[i] < ref:          # flatten the contiguous low-q ramp
                sc = ref / max(I[i], 1e-300)
                I[i] = ref
                err_lin[i] *= sc; lower[i] *= sc; upper[i] *= sc
                i += 1
    return q, I, err_lin, lower, upper


def _match_dq(SMR_dQ, q_in, q_out):
    if SMR_dQ is None:
        return None
    a = np.asarray(SMR_dQ, float)
    if a.ndim == 0:
        return a
    if len(a) == len(q_out):
        return a
    return np.interp(np.asarray(q_out, float), np.asarray(q_in, float), a)


# ---------------------------------------------------------------------------
# Dispatcher — same signature/return as desmearData
# ---------------------------------------------------------------------------
def desmear_dispatch(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=None, *,
                     method="lake", length_scale_decades=0.5, kernel="matern32",
                     sigma_log=4.0, extrap_method="PowerLaw w flat",
                     extrap_qstart=None, max_iter=20):
    """Desmear with the selected method. Returns (DSM_Qvec, DSM_Int, DSM_Error,
    DSM_dQ). method='lake' reproduces the current Matilda behaviour exactly."""
    if SMR_Int is None or len(SMR_Int) == 0 or slitLength is None:
        return None, None, None, None

    method = (method or "lake").lower()

    if method in ("lake", "strobl", ""):
        from .desmearing import desmearData
        return desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ,
                           slitLength=slitLength, ExtrapMethod=extrap_method,
                           ExtrapQstart=extrap_qstart, MaxNumIter=max_iter)

    if method in ("lake_smooth", "smoothed_lake"):
        # Lake, then geometric-mean smoothing (window = length_scale_decades).
        # Best for very noisy flat-plateau data, where the GP would ring in the
        # slit operator's low-q null space but Lake's incoherent noise smooths
        # cleanly to the physical plateau.
        from .desmearing import desmearData
        lq, lI, lE, ldQ = desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ,
                                      slitLength=slitLength, ExtrapMethod=extrap_method,
                                      ExtrapQstart=extrap_qstart, MaxNumIter=max_iter)
        win = length_scale_decades if length_scale_decades and length_scale_decades > 0 else 0.3
        return lq, _smooth_curve(lq, lI, decades=win), lE, ldQ

    if method in ("gp", "gp_matern", "gp_rbf", "huang", "huang_gp"):
        k = "rbf" if method == "gp_rbf" else kernel
        q, I, sd, lo, hi = _gp_core(
            SMR_Qvec, SMR_Int, SMR_Error, slitLength,
            length_scale_decades=length_scale_decades, kernel=k,
            sigma_log=sigma_log, max_iter=max(max_iter, 30))
        return q, I, sd, _match_dq(SMR_dQ, SMR_Qvec, q)

    if method in ("gp_lake_mean", "gp_lake"):
        from .desmearing import desmearData
        lq, lI, lE, ldQ = desmearData(
            SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=slitLength,
            ExtrapMethod=extrap_method, ExtrapQstart=extrap_qstart, MaxNumIter=max_iter)
        smq, smI = _smooth_loglog(lq, lI)
        q, I, sd, lo, hi = _gp_core(
            SMR_Qvec, SMR_Int, SMR_Error, slitLength,
            length_scale_decades=length_scale_decades, kernel=kernel,
            sigma_log=sigma_log, max_iter=max(max_iter, 30),
            prior_mean_q=smq, prior_mean_I=smI)
        return q, I, sd, _match_dq(SMR_dQ, SMR_Qvec, q)

    raise ValueError(f"unknown desmear method {method!r} "
                     "(valid: lake, gp, gp_rbf, gp_lake_mean)")


# Map GUI display strings -> (method_key, kernel_key)
GUI_METHOD_MAP = {
    "Lake":                 ("lake", "matern32"),
    "Lake (smoothed)":      ("lake_smooth", "matern32"),
    "Huang GP":             ("gp", "matern32"),
    "Huang GP + Lake mean": ("gp_lake_mean", "matern32"),
}
GUI_KERNEL_MAP = {
    "Matérn-3/2": "matern32",
    "Matern-3/2": "matern32",
    "RBF":        "rbf",
}
