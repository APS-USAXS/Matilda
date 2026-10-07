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


# --- FX4 noise / smoothing treatment ----------------------------------------
# Three switches, each isolating one correction, so they can be turned off
# individually while comparing against Igor.  All three apply to the FX4 chain
# ONLY; the scaler chain is untouched in every case.

# (1) Rebinning returns the error OF THE MEAN, sigma_rms/sqrt(n), instead of the
#     RMS of the n input errors (which is the error of a single point).  The
#     bins are disjoint, so for independent inputs sqrt(n) is simply correct.
#     The inputs have been smoothed and are therefore correlated, which makes
#     this a slight over-reduction; it is left to cancel against (3), where the
#     smoothing itself deliberately does NOT reduce the error.  Net effect is
#     within ~30% of an exact effective-N treatment on commissioning data.
FX4_ERROR_OF_MEAN_ON_REBIN = True

# (2) OFF by default -- see the measurement below before turning it on.
#     The electrometer reports values already averaged over `values_per_read`
#     raw samples behind an analogue filter, so the number of INDEPENDENT
#     samples behind each reading is N/values_per_read, not N = total/mean.
#     Dividing sigma by sqrt(N) therefore understates the PER-RAW-POINT error
#     by sqrt(values_per_read) = 3.2 on the commissioning files, which is the
#     3-5x shortfall recorded in docs/fx4-reduction-guide.md.
#
#     That reasoning is right, but applying it here double-counts, because the
#     pipeline also fails to propagate the noise reduction of smooth_r_data
#     (which averages 6-7 points and leaves the error untouched).  Measured on
#     GC_SRM3600_0131 by reducing two disjoint halves of the raw points
#     independently and comparing them -- the one estimator the smoothing's
#     correlation cannot fool -- true/quoted at the SMR stage is:
#
#         switch 1 on, this off :  1.34 / 0.87 / 0.90   (Q decades 1e-3..1)
#         switch 1 on, this on  :  0.42 / 0.27 / 0.29   -> 3x too conservative
#
#     So the two unpropagated effects currently cancel, and the honest fix is
#     to carry an effective sample count through smoothing and rebinning and
#     then enable this.  Until then, leave it off.
FX4_USE_EFFECTIVE_SAMPLES = False

# (3) The FX4 fly scan has no per-point amplifier range (the reading is
#     gain-independent and `upd_lurange` is written once, after the scan), so
#     every point would get the slowest smoothing time.  Derive a pseudo-range
#     from the signal level instead: on the scaler chain the amplifier range
#     WAS a proxy for signal level, and `rwave_smooth_times` was tuned against
#     it, so a signal-derived index reuses those constants self-consistently.
#     Thresholds are on the raw `upd_current`, descending; a point below the
#     last threshold gets FX4_SMOOTHING_RANGE_DEFAULT (the slowest range).
FX4_SIGNAL_DERIVED_SMOOTHING_RANGE = True
FX4_SMOOTHING_RANGE_THRESHOLDS = ((1e6, 0), (12000.0, 1), (800.0, 2), (30.0, 3))
FX4_SMOOTHING_RANGE_DEFAULT = 4


def signal_derived_range_index(upd_current,
                               thresholds=FX4_SMOOTHING_RANGE_THRESHOLDS,
                               default=FX4_SMOOTHING_RANGE_DEFAULT):
    """Pseudo amplifier-range index (0-4) from the FX4 signal level.

    See FX4_SIGNAL_DERIVED_SMOOTHING_RANGE.  Returns float so the array can
    carry NaN the way the scaler chain's masked points do.
    """
    upd = np.asarray(upd_current, dtype=float)
    index = np.full(upd.shape, float(default))
    for level, range_index in sorted(thresholds, key=lambda t: t[0]):
        index[upd > level] = float(range_index)
    return index


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
