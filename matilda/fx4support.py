"""
fx4support.py
=============
Counting-chain detection and FX4-specific helpers.

The 12-ID-E counting chain changed on 2026-09-26 from
"Femto amplifier + V/F converter + Struck scaler" to **FX4 electrometers**.
All four USAXS data formats (fly scan, step scan, SAXS frame, WAXS frame)
changed with it.  The single consequence that matters for reduction:

    Detector values are gain-independent PICOAMPS, not counts.
    There is no amplifier gain to divide by, no per-range gain table,
    and for the fly scan no dwell-time term either:

        I(q) = upd_current / I0_current

Some field names survived the change with different units and meaning
(``UPD`` and ``I0`` in the step-scan file are the worst offenders), so the
chain can NOT be inferred from field names.  Every format declares
``counting_chain``; that is the only thing to branch on.

Reference: ``bits_usaxs/docs/FX4_data_formats.md``.

Public API
----------
detect_counting_chain(h5obj)  -> 'FX4' or 'scaler'
is_fx4(h5obj)                 -> bool
warn_ring_overflows(...)      -> log a warning when FX4 dropped samples
as_scalar(value)              -> unwrap 1-element arrays / decode bytes
lookup(group, name, default)  -> tolerant scalar read from an h5py group
"""

import logging

import numpy as np

# --- chain names ------------------------------------------------------------
CHAIN_FX4 = "FX4"
CHAIN_SCALER = "scaler"

# Where ``counting_chain`` lives, per format.  Checked in this order; the
# first hit wins.  Entries are (path, attribute-name-or-None).
#
# NOTE on the fly scan: FX4_data_formats.md documents the marker as a dataset
# at /entry/flyScan/counting_chain, but saveFlyData.xml actually writes it as
# an ATTRIBUTE of /entry/program_name (together with config_version).  Both
# are checked so either spelling is recognised.
_CHAIN_LOCATIONS = (
    ("/entry/flyScan/counting_chain", None),            # fly scan (as documented)
    ("/entry/program_name", "counting_chain"),          # fly scan (as written)
    ("/entry/counting_chain", None),                    # SAXS / WAXS frame
    ("/entry/instrument/bluesky/metadata/counting_chain", None),  # step scan
)

# Relative uncertainty assigned to an FX4 current reading when the electrometer
# did not record a per-point sigma (step scans, and fly scans written without
# upd_sigma/I0_sigma).  The FX4 reports a mean over many samples, so the true
# uncertainty is sigma_sample/sqrt(N) — which is not in the file for those
# formats.  This is a placeholder to be replaced once the FX4's own noise has
# been characterised; it sets the error bars, never the intensities.
FX4_RELATIVE_CURRENT_ERROR = 0.01


def as_scalar(value):
    """Unwrap an HDF5 value into a plain Python scalar or str.

    h5py hands back 0-d arrays, 1-element arrays and bytes depending on how
    the writer created the dataset; callers just want the value.
    """
    if isinstance(value, np.ndarray):
        if value.size == 0:
            return None
        value = value.ravel()[0]
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode("utf-8", errors="replace")
    return value


def detect_counting_chain(h5obj):
    """Return ``'FX4'`` or ``'scaler'`` for an open h5py File.

    ``counting_chain`` absent means the old scaler chain — that is the whole
    test, per FX4_data_formats.md §1.  Never branch on ``config_version``:
    the step-scan file's ``program_name@config_version`` is the NeXus
    writer's schema version and did not move with this conversion.
    """
    for path, attr in _CHAIN_LOCATIONS:
        if path not in h5obj:
            continue
        node = h5obj[path]
        if attr is None:
            try:
                raw = node[()]
            except (TypeError, ValueError):
                continue
        else:
            if attr not in node.attrs:
                continue
            raw = node.attrs[attr]
        text = as_scalar(raw)
        if text is None:
            continue
        return CHAIN_FX4 if str(text).strip() == CHAIN_FX4 else CHAIN_SCALER
    return CHAIN_SCALER


def is_fx4(h5obj):
    """True when this file came from the FX4 chain (values in picoamps)."""
    return detect_counting_chain(h5obj) == CHAIN_FX4


def lookup(group, name, default=None):
    """Read a scalar ``name`` from an h5py group, or *default* when absent."""
    if group is None or name not in group:
        return default
    value = as_scalar(group[name][()])
    return default if value is None else value


def lookup_many(groups, name, default=None):
    """Read a scalar ``name`` from the first group that has it.

    The FX4 fly-scan configuration moved several fields out of
    ``/entry/metadata`` and into ``/entry/flyScan``, so the readers try both.
    """
    for group in groups:
        if group is not None and name in group:
            value = as_scalar(group[name][()])
            if value is not None:
                return value
    return default


def warn_ring_overflows(values, filename, what="FX4"):
    """Warn when the FX4 driver's ring buffer overflowed.

    A non-zero count means the oldest samples in an averaging interval were
    discarded, so every mean in the file is biased toward the end of its
    interval.  The data are still readable; the means are not trustworthy.
    """
    for label, value in values.items():
        if value is None:
            continue
        try:
            count = float(as_scalar(value))
        except (TypeError, ValueError):
            continue
        if count:
            logging.warning(
                f"{filename}: {what} {label} = {count:g} (expected 0). "
                "Samples were dropped and the recorded means are biased "
                "toward the end of each interval."
            )


def range_indexed_array(range_index, table, n_points, default=0.0):
    """Expand a per-range lookup *table* to a per-point array.

    Parameters
    ----------
    range_index : array-like or scalar or None
        FX4 range (0-4) at each point.  A scalar is broadcast to every point;
        None yields *default* everywhere.
    table : dict
        Maps int range -> value (for example dark current per range).
    n_points : int
        Length of the output array.
    default : float
        Value used where the range is unknown or missing from *table*.
    """
    out = np.full(n_points, float(default), dtype=float)
    if range_index is None:
        return out
    idx = np.asarray(range_index, dtype=float)
    if idx.ndim == 0:
        idx = np.full(n_points, float(idx))
    elif idx.size < n_points:
        idx = np.pad(idx, (0, n_points - idx.size), constant_values=np.nan)
    else:
        idx = idx[:n_points]
    for i in range(n_points):
        if np.isnan(idx[i]):
            continue
        value = table.get(int(idx[i]))
        if value is not None:
            out[i] = float(value)
    return out


def mean_of_samples(total, mean, fallback=1.0):
    """Number of samples averaged into each FX4 mean: ``N = total / mean``.

    Used for the error of the mean (``sigma/sqrt(N)``) and, in the fly scan,
    for the derived ``channel_time``.  Returns *fallback* where the division
    is undefined, so the result is always safe to take sqrt of.
    """
    if total is None or mean is None:
        return None
    total = np.asarray(total, dtype=float)
    mean = np.asarray(mean, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        samples = np.abs(total / mean)
    samples = np.where(np.isfinite(samples) & (samples >= 1.0), samples, fallback)
    return samples


def ratio_error(numerator, num_error, denominator, den_error):
    """Propagate uncertainty through ``numerator / denominator``.

    Both inputs are currents in pA; the result is the absolute uncertainty of
    the ratio.  Zero or non-finite inputs give a zero contribution rather than
    an inf/NaN that would poison everything downstream.
    """
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    num_error = np.broadcast_to(np.asarray(num_error, dtype=float), numerator.shape)
    den_error = np.broadcast_to(np.asarray(den_error, dtype=float), denominator.shape)

    with np.errstate(divide="ignore", invalid="ignore"):
        rel_num = np.where(numerator != 0, num_error / numerator, 0.0)
        rel_den = np.where(denominator != 0, den_error / denominator, 0.0)
        ratio = np.where(denominator != 0, numerator / denominator, 0.0)
        error = np.abs(ratio) * np.sqrt(rel_num**2 + rel_den**2)
    return np.where(np.isfinite(error), error, 0.0)
