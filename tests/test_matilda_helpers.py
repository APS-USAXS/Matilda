"""Unit tests for pure helpers in matilda.matilda (service orchestrator).

Importing matilda.matilda pulls in the full converter chain (pyFAI,
tifffile, requests) — skip cleanly when those heavy deps are absent.
"""

import os

import pytest

pytest.importorskip("pyFAI")
pytest.importorskip("requests")

from matilda.matilda import (
    _extract_sample_key,
    _find_matching_partner,
    _remember_file,
    extract_number_from_filename,
)


# ── extract_number_from_filename ─────────────────────────────────────────────

@pytest.mark.parametrize("filename,expected", [
    ("Sample_55C_10min_1234.h5", 1234),
    ("Sample_0044.hdf", 44),
    ("Sample_0044.hdf5", 44),
    ("blank_0006.H5", 6),           # case-insensitive extension
    ("scan_12.nxs", 12),
    ("noNumber.h5", 0),
    ("wrong_123.txt", 0),
])
def test_extract_number_from_filename(filename, expected):
    assert extract_number_from_filename(filename) == expected


# ── _remember_file (bounded FIFO tracker) ────────────────────────────────────

def test_remember_file_fifo_eviction():
    tracked = {}
    for k in range(105):
        _remember_file(tracked, ("p", k), maxlen=100)
    assert len(tracked) == 100
    assert ("p", 104) in tracked        # newest kept
    assert ("p", 0) not in tracked      # oldest evicted


def test_remember_file_membership():
    tracked = {}
    _remember_file(tracked, ("path", "file.h5"))
    assert ("path", "file.h5") in tracked


# ── _extract_sample_key ──────────────────────────────────────────────────────

def test_extract_sample_key():
    assert _extract_sample_key("Sample1_55C_10min_1234.h5") == ("Sample1", "1234")
    assert _extract_sample_key("plain.h5") is None


# ── _find_matching_partner ───────────────────────────────────────────────────

def _make_pair(tmp_path, usaxs_name, saxs_names):
    usaxs = tmp_path / "exp_usaxs"
    saxs = tmp_path / "exp_saxs"
    usaxs.mkdir()
    saxs.mkdir()
    (usaxs / usaxs_name).touch()
    for name in saxs_names:
        (saxs / name).touch()
    return str(usaxs), str(saxs)


def test_find_matching_partner_basic(tmp_path):
    usaxs, saxs = _make_pair(tmp_path, "SampleA_10min_0044.h5",
                             ["SampleA_10min_0044.hdf"])
    partner = _find_matching_partner(usaxs, "SampleA_10min_0044.h5", "_saxs")
    assert partner is not None
    assert partner[1] == "SampleA_10min_0044.hdf"
    assert os.path.basename(partner[0]) == "exp_saxs"


def test_find_matching_partner_disambiguates_by_stem(tmp_path):
    """With multiple candidates the full-stem match wins (Phase 3)."""
    usaxs, saxs = _make_pair(tmp_path, "SampleA_10min_0044.h5",
                             ["SampleA_20min_0044.hdf",
                              "SampleA_10min_0044.hdf"])
    partner = _find_matching_partner(usaxs, "SampleA_10min_0044.h5", "_saxs")
    assert partner[1] == "SampleA_10min_0044.hdf"


def test_find_matching_partner_missing_folder(tmp_path):
    usaxs = tmp_path / "exp_usaxs"
    usaxs.mkdir()
    (usaxs / "S_0001.h5").touch()
    assert _find_matching_partner(str(usaxs), "S_0001.h5", "_saxs") is None
