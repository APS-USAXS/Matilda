"""Unit tests for matilda.hdf5code (HDF5 / NXcanSAS I/O)."""

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")

from matilda import hdf5code as hc


# ── save_dict_to_hdf5 / load_dict_from_hdf5 ─────────────────────────────────

def test_save_dict_skips_none_and_overwrites(tmp_path):
    """None values must be skipped and re-saving must not raise (2.8)."""
    fp = tmp_path / "t.h5"
    d = {"a": np.arange(5), "b": None, "nested": {"c": 1.5, "d": None}}
    with h5py.File(fp, "w") as f:
        hc.save_dict_to_hdf5(d, "entry/test/", f)
        hc.save_dict_to_hdf5(d, "entry/test/", f)      # overwrite path
        loaded = hc.load_dict_from_hdf5(f, "entry/test/")
    assert "b" not in loaded
    assert "d" not in loaded["nested"]
    assert loaded["nested"]["c"] == 1.5
    assert np.array_equal(loaded["a"], np.arange(5))


def test_load_dict_decodes_bytes(tmp_path):
    fp = tmp_path / "t.h5"
    with h5py.File(fp, "w") as f:
        f["entry/name"] = b"hello"
        loaded = hc.load_dict_from_hdf5(f, "entry/")
    assert loaded["name"] == "hello"


# ── readMyNXcanSAS ───────────────────────────────────────────────────────────

def _minimal_usaxs_file(fp, with_qrs=True):
    with h5py.File(fp, "w") as f:
        if with_qrs:
            g = f.create_group("entry/QRS_data")
            g["Intensity"] = np.ones(10)
            g["Q"] = np.linspace(0.001, 0.1, 10)
            g["Error"] = np.ones(10) * 0.1
        f.create_group("entry/metadata")
        f.create_group("entry/instrument")
        f.create_group("entry/sample")


def test_readmynxcansas_smr_keys_are_none_not_tuples(tmp_path):
    """Regression for the (None,) trailing-comma bug (1.5)."""
    fp = tmp_path / "scan.h5"
    _minimal_usaxs_file(fp)
    S = hc.readMyNXcanSAS(str(tmp_path), "scan.h5", isUSAXS=True)
    for key in ("SMR_Qvec", "SMR_Int", "SMR_Error", "SMR_dQ", "slitLength"):
        assert S["CalibratedData"][key] is None, f"{key} is not None"


# ── saveNXcanSAS ─────────────────────────────────────────────────────────────

def _sample_dict():
    return {
        "CalibratedData": {
            "Intensity": np.ones(5), "Q": np.linspace(0.01, 0.1, 5),
            "Error": np.ones(5) * 0.1, "dQ": np.ones(5) * 0.001,
            "units": "[cm2/cm3]", "Kfactor": None, "OmegaFactor": None,
            "blankname": None, "thickness": None,
        },
        "RawData": {
            "filename": "scan",
            "metadata": {"timeStamp": "2026-07-08"},
            "sample": {"name": "testsample"},
        },
        "reducedData": {
            "Intensity": np.ones(5), "Q": np.linspace(0.01, 0.1, 5),
            "Error": np.ones(5) * 0.1,
        },
    }


def test_savenxcansas_none_attributes(tmp_path):
    """None blankname/thickness must not raise TypeError (2.9)."""
    hc.saveNXcanSAS(_sample_dict(), str(tmp_path), "out.h5")
    with h5py.File(tmp_path / "out.h5") as f:
        assert "entry/testsample/sasdata/I" in f
        assert "entry/QRS_data/Intensity" in f


def test_savenxcansas_roundtrip(tmp_path):
    sample = _sample_dict()
    sample["CalibratedData"]["blankname"] = "blank_0001.h5"
    sample["CalibratedData"]["thickness"] = 1.2
    hc.saveNXcanSAS(sample, str(tmp_path), "out.h5")
    # add minimal raw groups so readMyNXcanSAS(isUSAXS=True) can complete
    with h5py.File(tmp_path / "out.h5", "a") as f:
        f.require_group("entry/metadata")
        f.require_group("entry/instrument")
        f.require_group("entry/sample")
    S = hc.readMyNXcanSAS(str(tmp_path), "out.h5", isUSAXS=True)
    assert np.allclose(S["CalibratedData"]["Intensity"], 1.0)
    assert S["CalibratedData"]["blankname"] == "blank_0001.h5"


# ── readGenericNXcanSAS ──────────────────────────────────────────────────────

def test_readgeneric_malformed_returns_none(tmp_path):
    """Malformed NXcanSAS entry must return None, not NameError (2.9)."""
    fp = tmp_path / "bad.h5"
    with h5py.File(fp, "w") as f:
        g = f.create_group("entry/bad")
        g.attrs["canSAS_class"] = "SASentry"
        g.attrs["NX_class"] = "NXsubentry"
        g["definition"] = "NXcanSAS"
    assert hc.readGenericNXcanSAS(str(tmp_path), "bad.h5") is None


def test_readgeneric_no_entries(tmp_path):
    fp = tmp_path / "empty.h5"
    with h5py.File(fp, "w") as f:
        f.create_group("entry")
    assert hc.readGenericNXcanSAS(str(tmp_path), "empty.h5") is None


# ── find_matching_groups ─────────────────────────────────────────────────────

def test_find_matching_groups(tmp_path):
    fp = tmp_path / "t.h5"
    with h5py.File(fp, "w") as f:
        g = f.create_group("entry/match")
        g.attrs["canSAS_class"] = "SASentry"
        g.attrs["NX_class"] = "NXsubentry"
        g["definition"] = "NXcanSAS"
        g2 = f.create_group("entry/nomatch")
        g2.attrs["canSAS_class"] = "SASentry"
        found = hc.find_matching_groups(
            f,
            {"canSAS_class": "SASentry", "NX_class": "NXsubentry"},
            {"definition": "NXcanSAS"},
        )
    assert found == ["entry/match"]
