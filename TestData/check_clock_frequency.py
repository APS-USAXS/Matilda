#!/usr/bin/env python3
"""
check_clock_frequency.py
========================
Diagnostic: determine the V-to-F / scaler clock frequency from data alone.

The idea: sum(TimePerPoint)/frequency must equal the real (wall-clock)
duration of the measurement.  Testing candidate frequencies (1e6, 1e7)
against the known scan duration tells you which constant is correct —
no hardware trace needed.

Usage:
    python check_clock_frequency.py scanfile1.h5 [scanfile2.h5 ...]

Works on both file types (auto-detected):
    * USAXS flyscan  — sums /entry/flyScan/mca1
      (matilda flyscan code currently assumes 1e6)
    * USAXS step scan — sums /entry/data/seconds
      (matilda step-scan code currently assumes 1e7; this also reveals
       whether 'seconds' is raw scaler counts or already seconds)

The script prints the total measurement time under each interpretation and,
when start/end timestamps are present in the file, compares against the
wall-clock duration and prints a verdict.  Note the summed count time is a
LOWER bound on wall-clock time (fly return strokes, motor moves and readout
overhead add to wall-clock), so the right candidate gives a total <= wall
clock and within roughly a factor of ~2 of it; wrong candidates are off by
a decade.

Requires only h5py + numpy (no matilda imports).
"""

import sys
import datetime
import h5py
import numpy as np

CANDIDATE_FREQUENCIES = (1e6, 1e7)


def _read_scalar(f, path):
    """Return dataset value at path (decoded if bytes), or None."""
    if path not in f:
        return None
    val = f[path][()]
    if isinstance(val, np.ndarray) and val.size == 1:
        val = val.item()
    if isinstance(val, bytes):
        val = val.decode('utf-8', errors='replace')
    return val


def _parse_time(val):
    """Parse an ISO-8601 timestamp string to datetime, or return None."""
    if val is None:
        return None
    try:
        return datetime.datetime.fromisoformat(str(val).strip())
    except ValueError:
        return None


def _wall_clock_seconds(f):
    """Best-effort wall-clock duration from timestamps in the file."""
    # explicit duration dataset, if the writer stored one
    dur = _read_scalar(f, '/entry/duration')
    if dur is not None:
        try:
            return float(dur), 'from /entry/duration'
        except (TypeError, ValueError):
            pass
    start = _parse_time(_read_scalar(f, '/entry/start_time'))
    end = _parse_time(_read_scalar(f, '/entry/end_time'))
    if start is not None and end is not None:
        return (end - start).total_seconds(), 'from start_time/end_time'
    return None, None


def _verdict(total_counts, wall_seconds):
    """Print per-candidate totals and a verdict line."""
    best = None
    for freq in CANDIDATE_FREQUENCIES:
        total_s = total_counts / freq
        line = f"    assuming clock {freq:8.0e} Hz -> total measurement time {total_s:12.3f} s"
        if wall_seconds is not None:
            ratio = total_s / wall_seconds if wall_seconds > 0 else float('inf')
            line += f"   ({ratio:6.2%} of wall clock)"
            # right answer: total <= wall clock, same order of magnitude
            if 0.2 <= ratio <= 1.05:
                line += "   <== PLAUSIBLE"
                best = freq
        print(line)
    if wall_seconds is not None:
        print(f"    wall-clock duration: {wall_seconds:.1f} s")
        if best is not None:
            print(f"    VERDICT: clock frequency is {best:.0e} Hz")
        else:
            print("    VERDICT: no candidate fits — compare totals with the "
                  "known scan duration manually (also consider that "
                  "'seconds' may already be in seconds, see below).")
    else:
        print("    (no start/end timestamps found — compare the totals above "
              "with the scan duration you know from the acquisition setup)")


def check_flyscan(f, filename):
    mca1 = np.ravel(np.array(f['/entry/flyScan/mca1']))
    total_counts = float(np.sum(mca1))
    print(f"  FLYSCAN  {filename}")
    print(f"    points: {mca1.size},  sum(mca1) = {total_counts:.6e} counts")
    stored = _read_scalar(f, '/entry/flyScan/mca_clock_frequency')
    if stored is not None:
        print(f"    mca_clock_frequency stored in file: {stored} "
              "(matilda code says this value is wrong; shown for reference)")
    wall, src = _wall_clock_seconds(f)
    if src:
        print(f"    wall clock {src}")
    _verdict(total_counts, wall)
    print("    matilda flyscan code currently assumes 1e6 Hz.")


def check_stepscan(f, filename):
    seconds = np.ravel(np.array(f['/entry/data/seconds']))
    total = float(np.sum(seconds))
    print(f"  STEP SCAN  {filename}")
    print(f"    points: {seconds.size},  sum(seconds dataset) = {total:.6e}")
    print(f"    per-point values: min {np.min(seconds):.4g}, "
          f"median {np.median(seconds):.4g}, max {np.max(seconds):.4g}")
    wall, src = _wall_clock_seconds(f)
    if src:
        print(f"    wall clock {src}")
    # Interpretation 1: dataset is already in seconds
    print(f"    if 'seconds' is ALREADY seconds -> total {total:12.3f} s"
          + (f"   ({total / wall:6.2%} of wall clock)" if wall else ""))
    # Interpretations 2+: dataset is scaler counts of a candidate clock
    _verdict(total, wall)
    print("    matilda step-scan code currently divides by 1e7 Hz.")
    print("    NOTE: typical step-scan count times are ~0.1-10 s per point —")
    print("    the per-point values above under the right interpretation")
    print("    must land in that range.")


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 1
    for filepath in argv[1:]:
        print("=" * 70)
        try:
            with h5py.File(filepath, 'r') as f:
                if '/entry/flyScan/mca1' in f:
                    check_flyscan(f, filepath)
                elif '/entry/data/seconds' in f:
                    check_stepscan(f, filepath)
                else:
                    print(f"  {filepath}: neither /entry/flyScan/mca1 nor "
                          "/entry/data/seconds found — not a USAXS scan file?")
        except (OSError, KeyError) as exc:
            print(f"  {filepath}: cannot read ({exc})")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
