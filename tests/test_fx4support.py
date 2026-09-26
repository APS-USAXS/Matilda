"""Tests for FX4 counting-chain detection and the FX4 reduction branches.

The FX4 electrometers replaced the Femto + V/F + Struck chain on 2026-09-26.
Real FX4 step-scan / SAXS / WAXS files exist; the FX4 fly scan had not been
commissioned when this was written, so that branch is exercised against a
synthetic file built to the deployed ADconfigs/Flyscan_config/saveFlyData.xml
layout — including the case where that configuration writes no
``/entry/metadata`` group at all.
"""

import h5py
import numpy as np
import pytest

from matilda import fx4support as fx4
from matilda.fx4support import CHAIN_FX4, CHAIN_SCALER


# ── chain detection ──────────────────────────────────────────────────────────

def _write(tmp_path, name, build):
    path = tmp_path / name
    with h5py.File(path, "w") as handle:
        build(handle)
    return path


def test_detect_flyscan_attribute(tmp_path):
    """saveFlyData writes counting_chain as an attribute of program_name."""
    def build(f):
        ds = f.create_dataset("/entry/program_name", data=b"saveFlyData.py")
        ds.attrs["config_version"] = "2.0"
        ds.attrs["counting_chain"] = "FX4"
    with h5py.File(_write(tmp_path, "fly.h5", build)) as f:
        assert fx4.detect_counting_chain(f) == CHAIN_FX4


def test_detect_flyscan_dataset(tmp_path):
    """FX4_data_formats.md documents it as /entry/flyScan/counting_chain."""
    def build(f):
        f.create_dataset("/entry/flyScan/counting_chain", data=b"FX4")
    with h5py.File(_write(tmp_path, "fly2.h5", build)) as f:
        assert fx4.detect_counting_chain(f) == CHAIN_FX4


def test_detect_area_detector(tmp_path):
    def build(f):
        f.create_dataset("/entry/counting_chain", data=np.array([b"FX4"]))
    with h5py.File(_write(tmp_path, "saxs.hdf", build)) as f:
        assert fx4.is_fx4(f)


def test_detect_step_scan(tmp_path):
    def build(f):
        f.create_dataset("/entry/instrument/bluesky/metadata/counting_chain",
                         data=b"FX4")
    with h5py.File(_write(tmp_path, "step.h5", build)) as f:
        assert fx4.detect_counting_chain(f) == CHAIN_FX4


def test_absent_marker_means_scaler(tmp_path):
    """Absent counting_chain -> old chain.  That is the whole test."""
    def build(f):
        ds = f.create_dataset("/entry/program_name", data=b"saveFlyData.py")
        ds.attrs["config_version"] = "1.3"
    with h5py.File(_write(tmp_path, "old.h5", build)) as f:
        assert fx4.detect_counting_chain(f) == CHAIN_SCALER


def test_step_scan_writer_version_is_not_a_chain_marker(tmp_path):
    """program_name@config_version = '1.0' on a step scan is the NeXus
    writer's schema version and did not move with the conversion."""
    def build(f):
        ds = f.create_dataset("/entry/program_name", data=b"NXWriterUascan")
        ds.attrs["config_version"] = "1.0"
        f.create_dataset("/entry/instrument/bluesky/metadata/counting_chain",
                         data=b"FX4")
    with h5py.File(_write(tmp_path, "uascan.h5", build)) as f:
        assert fx4.detect_counting_chain(f) == CHAIN_FX4


# ── helpers ──────────────────────────────────────────────────────────────────

def test_range_indexed_array_scalar_and_per_point():
    table = {0: 1.0, 4: 8.35}
    assert np.allclose(fx4.range_indexed_array(4, table, 3), [8.35] * 3)
    assert np.allclose(fx4.range_indexed_array([0, 4, 2], table, 3),
                       [1.0, 8.35, 0.0])          # 2 not in table -> default
    assert np.allclose(fx4.range_indexed_array(None, table, 2), [0.0, 0.0])


def test_mean_of_samples():
    # 1000 samples of 5 pA sum to 5000
    assert np.allclose(fx4.mean_of_samples([5000.0], [5.0]), [1000.0])
    # undefined division falls back rather than producing inf
    assert np.allclose(fx4.mean_of_samples([1.0], [0.0]), [1.0])


def test_ratio_error_matches_hand_propagation():
    err = fx4.ratio_error([100.0], [1.0], [50.0], [2.0])
    expected = (100 / 50) * np.sqrt((1 / 100) ** 2 + (2 / 50) ** 2)
    assert err[0] == pytest.approx(expected)


def test_ratio_error_is_finite_on_zero_denominator():
    assert fx4.ratio_error([1.0], [1.0], [0.0], [1.0])[0] == 0.0


def test_as_scalar_decodes_bytes_and_unwraps_arrays():
    assert fx4.as_scalar(np.array([b"FX4"])) == "FX4"
    assert fx4.as_scalar(np.array(3.5)) == 3.5
    assert fx4.as_scalar(np.array([])) is None


def test_warn_ring_overflows(caplog):
    fx4.warn_ring_overflows({"ring_overflows": 7}, "sample_0001.h5")
    assert "ring_overflows = 7" in caplog.text
    caplog.clear()
    fx4.warn_ring_overflows({"ring_overflows": 0}, "sample_0001.h5")
    assert caplog.text == ""


# ── synthetic FX4 fly scan ───────────────────────────────────────────────────

AR_CENTER = 8.8439
WAVELENGTH = 0.5904
N_POINTS = 2000
CHANNEL_TIME = 0.05        # s per PSO interval
I0_CURRENT = 5.5e4         # pA, steady monitor
UPD_RANGE = 1              # autoranger range the scan ended on
DARK = {0: 0.0, 1: 0.05, 2: 0.5, 3: 1.2, 4: 8.4}          # pA, per FX4 range
DARK_ERR = {0: 0.0, 1: 0.005, 2: 0.02, 3: 0.04, 4: 0.08}
PEAK_R = 1.0e4             # rocking-curve peak, in R = upd/I0 units
PEAK_WIDTH_DEG = 2.0e-4
Q_WIDTH = 2.0e-4           # 1/A, where the slit-smeared wings roll over


def _synthetic_ar_sweep():
    """AR positions of a fly scan: dense at the peak, sparse at high q.

    The real trajectory varies the stage speed, so the pulse positions are
    not uniform; a uniform sweep would put two points across the 5e-4 deg
    rocking curve and leave the peak fit nothing to work with.
    """
    before = AR_CENTER + np.geomspace(5e-2, 1e-5, 300)
    after = AR_CENTER - np.geomspace(1e-5, 5e-2, N_POINTS - 300)
    return np.concatenate((before, after))       # descending, as on the stage


def _synthetic_fx4_flyscan(path, transmission=1.0, excess=0.0,
                           write_metadata_group=True):
    """Write an FX4 fly-scan file matching saveFlyData.xml v2.0.

    *transmission* attenuates the rocking-curve peak (so the sample/blank
    peak-to-peak ratio is meaningful) and *excess* adds sample scattering
    above the blank's slit-smeared wings, which is what the Qmin search in
    calibrateAndSubtractFlyscan looks for.

    *write_metadata_group* False reproduces the deployed configuration, in
    which an over-long XML comment swallows the ``<group name="metadata">``
    tag so every one of those fields lands in ``/entry/flyScan`` instead.
    """
    ar = _synthetic_ar_sweep()
    dtheta_deg = ar - AR_CENTER
    q = np.abs(4 * np.pi * np.sin(np.radians(dtheta_deg) / 2) / WAVELENGTH)

    R = transmission * PEAK_R * np.exp(-(dtheta_deg / PEAK_WIDTH_DEG) ** 2)
    R = R + 2.0 * (1 + (q / Q_WIDTH) ** 2) ** -1.5        # blank wings
    R = R + excess * (1 + (q / Q_WIDTH) ** 2) ** -1.75    # sample scattering
    upd = R * I0_CURRENT + DARK[UPD_RANGE]
    i0 = np.full(N_POINTS, I0_CURRENT)
    n_samples = CHANNEL_TIME / 1e-3

    with h5py.File(path, "w") as f:
        prog = f.create_dataset("/entry/program_name", data=b"saveFlyData.py")
        prog.attrs["config_version"] = "2.0"
        prog.attrs["counting_chain"] = "FX4"

        fly = f.create_group("/entry/flyScan")
        # AR_PulsePositions carries one extra leading point and zero padding
        fly["AR_PulsePositions"] = np.concatenate(([ar[0]], ar, np.zeros(500)))
        fly["upd_current"] = upd
        fly["I0_current"] = i0
        fly["upd_sigma"] = 0.02 * upd
        fly["I0_sigma"] = 0.02 * i0
        fly["upd_total"] = upd * n_samples
        fly["I0_total"] = i0 * n_samples
        fly["channel_time"] = np.full(N_POINTS, CHANNEL_TIME)
        fly["sample_time"] = 1e-3
        fly["values_per_read"] = 100
        fly["n_points"] = N_POINTS
        fly["ring_overflows"] = 0
        fly["ring_overflows_I0"] = 0
        fly["upd_lurange"] = UPD_RANGE
        for i in range(5):
            fly[f"upd_bkg{i}"] = DARK[i]
            fly[f"upd_bkgErr{i}"] = DARK_ERR[i]
            fly[f"upd_amp_change_mask_time{i}"] = 0.01

        md = f.create_group("/entry/metadata") if write_metadata_group else fly
        md["detector_distance"] = 1073.0
        md["UPDsize"] = 4.95
        md["AR_center"] = AR_CENTER
        md["timeStamp"] = b"2026-09-27 09:00:00"
        md["trans_pin_counts"] = 1.2e7
        md["trans_pin_gain"] = 1.0
        md["trans_I0_counts"] = 6.6e4
        md["trans_I0_gain"] = 1.0

        f["/entry/instrument/source/incident_wavelength"] = WAVELENGTH
        f["/entry/instrument/monochromator/wavelength"] = WAVELENGTH
        f["/entry/instrument/monochromator/energy"] = 21.0
        f["/entry/sample/thickness"] = 1.0
        f["/entry/sample/name"] = b"synthetic"


@pytest.fixture(params=[True, False],
                ids=["with-metadata-group", "metadata-fields-in-flyScan"])
def fx4_flyscan_pair(request, tmp_path):
    """A blank and a scattering sample FX4 fly scan in the same folder."""
    _synthetic_fx4_flyscan(tmp_path / "Blank_0001.h5",
                           write_metadata_group=request.param)
    _synthetic_fx4_flyscan(tmp_path / "Sample_0002.h5", transmission=0.75,
                           excess=4.0, write_metadata_group=request.param)
    return tmp_path


def test_import_fx4_flyscan(fx4_flyscan_pair):
    from matilda.supportFunctions import importFlyscan
    raw = importFlyscan(str(fx4_flyscan_pair), "Sample_0002.h5")
    assert raw["chain"] == CHAIN_FX4
    assert len(raw["ARangles"]) == N_POINTS
    assert np.allclose(raw["TimeInSec"], CHANNEL_TIME)
    # gain-independent chain: I0 gain is exactly 1, not the Femto value
    assert raw["metadata"]["I0Gain"] == 1.0
    assert raw["metadata"]["detector_distance"] == 1073.0
    assert raw["RangeIndex"] == UPD_RANGE


def test_fx4_flyscan_prefers_the_fx4_dark_currents(tmp_path):
    """Femto-era upd_bkg* in /entry/metadata must not win over the FX4 ones.

    Both can be present: /entry/metadata gets usxLAX:pd01:seq02:bkg*, which is
    stale on this chain, while /entry/flyScan gets usxFX4:FX4:seq01:bkg*.
    """
    from matilda.supportFunctions import importFlyscan
    path = tmp_path / "Sample_0002.h5"
    _synthetic_fx4_flyscan(path, write_metadata_group=True)
    with h5py.File(path, "a") as f:
        for i in range(5):
            f[f"/entry/metadata/upd_bkg{i}"] = 1.0e6      # a Femto-era count
        # a field with no FX4 counterpart still comes from /entry/metadata
        assert f["/entry/metadata/UPDsize"][()] == 4.95

    raw = importFlyscan(str(tmp_path), "Sample_0002.h5")
    assert raw["metadata"]["upd_bkg1"] == DARK[1]
    assert raw["metadata"]["UPDsize"] == 4.95


def test_fx4_flyscan_intensity_is_a_plain_current_ratio(fx4_flyscan_pair):
    """I(q) = (upd_current - dark) / I0_current: no gain, no dwell term."""
    from matilda.supportFunctions import importFlyscan, calculatePD_Fly
    Sample = {"RawData": importFlyscan(str(fx4_flyscan_pair), "Sample_0002.h5")}
    reduced = calculatePD_Fly(Sample)
    raw = Sample["RawData"]
    expected = (raw["UPD_array"] - DARK[UPD_RANGE]) / raw["Monitor"]
    assert np.allclose(reduced["Intensity"], expected)
    assert np.allclose(reduced["UPD_gains"], 1.0)


def test_fx4_flyscan_errors_use_the_recorded_sigmas(fx4_flyscan_pair):
    from matilda.supportFunctions import (importFlyscan, calculatePD_Fly,
                                          calculatePDErrorFly)
    Sample = {"RawData": importFlyscan(str(fx4_flyscan_pair), "Sample_0002.h5")}
    Sample["reducedData"] = calculatePD_Fly(Sample)
    error = calculatePDErrorFly(Sample)["Error"]
    assert np.all(np.isfinite(error))
    assert np.all(error > 0)
    # 2% spread within each interval over 50 samples is ~0.28% on the mean,
    # and the same again on I0, so the ratio's relative error stays under 1%.
    relative = error / Sample["reducedData"]["Intensity"]
    assert np.nanmedian(relative) < 0.01


def test_fx4_flyscan_full_reduction(fx4_flyscan_pair):
    """End to end, including blank subtraction and desmearing."""
    from matilda.convertFlyscan import processFlyscan
    Sample = processFlyscan(str(fx4_flyscan_pair), "Sample_0002.h5",
                            blankPath=str(fx4_flyscan_pair),
                            blankFilename="Blank_0001.h5",
                            recalculateAllData=True)
    calibrated = Sample["CalibratedData"]
    assert calibrated["Q"] is not None
    assert len(calibrated["Q"]) > 50
    assert np.all(np.isfinite(calibrated["Intensity"]))
    assert calibrated["units"] == "[cm2/cm3]"


# ── FX4 step scan ────────────────────────────────────────────────────────────

def test_fx4_step_scan_has_no_gain_or_dwell_term():
    """Guard the property that distinguishes the chains, on source."""
    import inspect
    from matilda import convertUSAXS as cu
    src = inspect.getsource(cu._CorrectUPDStepFX4)
    assert "(UPD_array - Bckg_corr) / Monitor" in src
    assert "AmpGain" not in src
