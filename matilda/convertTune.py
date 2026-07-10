"""
convertTune.py
==============
Live tuning-curve extraction for the Matilda beamline service.

The USAXS instrument runs three step-scan tuning macros — ``tune_ar``,
``tune_mr`` and ``tune_a2rp`` — whose data live only in DataBroker and are
reachable through the Tiled server as JSON (they write no HDF5 file).  This
module queries Tiled for the last N tune scans of each type, downloads their
measured arrays, and hands them to plotData.plotTuneResults() so the service
loop can publish live tuning images (tune_ar.jpg, tune_mr.jpg, tune_a2rp.jpg)
next to the existing SAXS/WAXS/USAXS plots.

Axis mapping is derived from the scan metadata, not hardcoded:
    x-axis = motors[0]          (a_stage_r / m_stage_r / a_stage_r2p)
    y-axis = first non-scaler0 detector   (UPD for ar/a2rp, I0 for mr)

Everything here is failure-tolerant: on any network error (dead server, or an
offline installation with no Tiled at all) the functions log and return empty
results, so the main service loop never crashes.

Versions
--------
0.1  2026-07-03  initial version (JIL / Claude)

Public API
----------
getTuneResults(plan_name, NumScans=NumberOfTunesToShow, LastNdays=NumberOfDaysToLookBackTunes)
    Search + download; return an ordered (newest-first) list of tune result
    dicts ready for plotting.

Module parameters
-----------------
NumberOfTunesToShow          default number of tune curves per plot (single knob)
TUNE_PLAN_NAMES              the three tuning plan names handled by the service
NumberOfDaysToLookBackTunes  time window for the Tiled search
"""

import logging
import numpy as np

from .readfromtiled import FindLastTuneScans, tiled_get_primary_data


# --- module parameters (the single place to tune behavior) ---
# Number of most-recent tune curves overlaid on each plot.  Keep <= 10 so the
# tab10 color cycle in plotData.py stays distinct.
NumberOfTunesToShow = 5

# The tuning plans handled by the live service, in display order.
TUNE_PLAN_NAMES = ("tune_ar", "tune_mr", "tune_a2rp")

# How many days back to search for tune scans.
NumberOfDaysToLookBackTunes = 1


def _derive_detector_key(data, motor, detector):
    """Return the data-dict key holding the detector (y-axis) counts.

    Prefers the detector name supplied by the metadata.  Falls back to the
    first array column that is not the motor, its ``_user_setpoint``, or a
    ``scaler0*`` channel.

    Parameters
    ----------
    data : dict
        Flat primary-data dict (column name -> list of values).
    motor : str or None
        Motor / x-axis column name.
    detector : str or None
        Detector name from metadata, if any.

    Returns
    -------
    str or None
        Key present in *data* to use as the y-axis, or None if none suitable.
    """
    if detector and detector in data:
        return detector

    for key in data:
        if key == motor:
            continue
        if key.endswith("_user_setpoint"):
            continue
        if key.startswith("scaler0"):
            continue
        return key
    return None


def downloadTuneData(uid, motor, detector):
    """Download one tune scan and return its (x, y) arrays.

    Parameters
    ----------
    uid : str
        Run uid (Tiled id) of the tune scan.
    motor : str or None
        Motor / x-axis column name (from metadata).
    detector : str or None
        Detector / y-axis column name (from metadata); may be None, in which
        case it is derived from the returned data.

    Returns
    -------
    (numpy.ndarray, numpy.ndarray, str, str) or None
        (x, y, motor_key, detector_key) on success, or None if the data could
        not be downloaded or the required columns are missing.  Never raises.
    """
    data = tiled_get_primary_data(uid)
    if not data:
        return None

    # Resolve the x-axis key: prefer metadata motor, else look for a readback.
    motor_key = motor if (motor and motor in data) else None
    if motor_key is None:
        # last-resort: a column that is neither a setpoint nor a scaler channel
        for key in data:
            if not key.endswith("_user_setpoint") and not key.startswith("scaler0"):
                motor_key = key
                break

    detector_key = _derive_detector_key(data, motor_key, detector)

    if motor_key is None or detector_key is None:
        logging.error(
            f"Tune {uid}: could not resolve axes "
            f"(motor={motor_key}, detector={detector_key}); keys={list(data.keys())}"
        )
        return None

    try:
        x = np.asarray(data[motor_key], dtype=float)
        y = np.asarray(data[detector_key], dtype=float)
    except Exception as e:
        logging.error(f"Tune {uid}: failed to convert arrays: {e}")
        return None

    if x.size == 0 or y.size == 0 or x.size != y.size:
        logging.error(
            f"Tune {uid}: bad array sizes x={x.size} y={y.size}, skipping"
        )
        return None

    return x, y, motor_key, detector_key


def getTuneResults(plan_name, NumScans=NumberOfTunesToShow,
                   LastNdays=NumberOfDaysToLookBackTunes):
    """Search Tiled and download the last N tune scans of one plan type.

    Parameters
    ----------
    plan_name : str
        Tuning plan name ('tune_ar', 'tune_mr', 'tune_a2rp').
    NumScans : int, optional
        Maximum number of tune curves to return.  Default NumberOfTunesToShow.
    LastNdays : int, optional
        Time window for the search.  Default NumberOfDaysToLookBackTunes.

    Returns
    -------
    list of dict
        One dict per successfully downloaded scan, newest-first, with keys:
            'uid', 'scan_id', 'time', 'motor', 'detector',
            'x' (numpy.ndarray), 'y' (numpy.ndarray)
        Scans that fail to download are skipped.  Returns an empty list on
        network failure or when the server is unreachable (offline installs).
    """
    meta = FindLastTuneScans(plan_name, NumScans, LastNdays)
    if not meta:
        return []

    results = []
    for m in meta:
        downloaded = downloadTuneData(m["uid"], m.get("motor"), m.get("detector"))
        if downloaded is None:
            continue
        x, y, motor_key, detector_key = downloaded
        results.append({
            "uid": m["uid"],
            "scan_id": m.get("scan_id"),
            "time": m.get("time"),
            "motor": motor_key,
            "detector": detector_key,
            "x": x,
            "y": y,
        })
    logging.info(f"Tune {plan_name}: downloaded {len(results)}/{len(meta)} curves")
    return results


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    for pn in TUNE_PLAN_NAMES:
        res = getTuneResults(pn, NumberOfTunesToShow, NumberOfDaysToLookBackTunes)
        print(f"{pn}: {[(d['scan_id'], d['motor'], d['detector'], len(d['x'])) for d in res]}")
