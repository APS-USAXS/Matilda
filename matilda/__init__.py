"""
matilda
=======
Live data-processing package for the APS USAXS/SAXS/WAXS instrument at
sector 12-ID-E, Argonne National Laboratory.

Matilda runs as a systemd user service on ``usaxscontrol.xray.aps.anl.gov``.
Every 15 seconds it polls the Tiled server for new scan data, reduces any
new files to calibrated 1-D I(Q) curves, optionally runs pynika for
detector-geometry recalibration, and saves summary JPEG plots to a
web-visible directory for live monitoring.

Service entry-point
-------------------
    serv_matilda.sh  →  python matilda/matilda.py

The script is run *directly* (not as ``python -m matilda``), which adds the
``matilda/`` directory to sys.path.  This allows bare module imports
(``from convertFlyscan import …``) to work.  See the Bug Register below.

Architecture overview
---------------------

Orchestrator
~~~~~~~~~~~~
matilda/matilda.py
    Main 15-second polling loop.  Queries Tiled, compares lists, dispatches
    to pynika calibration if AgBehenateLaB6 calibrant is detected, then
    calls the processXxx() batch functions and plotXxx() for display.

    Batch functions exposed here:
        processFlyscans(ListOfScans, ListOfBlanks, …)
        processStepscans(ListOfScans, ListOfBlanks, …)
        processADscans(ListOfScans, ListOfBlanks, …)
        processUSAXSFolder(path)

    pynika integration (issue #28):
        _checkAndRunPynikaCalibration()  /  _runPynikaCalibration()
        Triggered when filename contains 'AgBehenateLaB6' (case-insensitive).
        Calls pynika via ``conda run -p <env_path> pynika --auto-fit --save-to-pvs``.

Data-reduction modules (current flat layout, future: matilda/reduction/)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matilda/convertFlyscan.py
    USAXS flyscan: HDF5 NXsas → normalised I(Q) → desmeared calibrated I(Q).
    Pipeline: importFlyscan → calculatePD_Fly → beamCenterCorrection →
              getBlankFlyscan → normalizeByTransmission →
              calibrateAndSubtractFlyscan → desmearData → saveNXcanSAS.

matilda/convertUSAXS.py
    USAXS step-scan: same pipeline as flyscan (identical file format).
    Also exports rebinData() used by convertFlyscan.

matilda/convertSWAXS.py
    SAXS / WAXS area-detector: HDF5 → pyFAI azimuthal integration → I(Q).
    Converts Nika geometry (beam-centre, tilts) to pyFAI PONI via Fit2D.
    Pipeline: importADData → calibrateAD2DData → reduceADData → saveNXcanSAS.

matilda/desmearing.py
    Lake / Strobl iterative slit-smearing correction.  Ported from Igor Pro.
    Entry point: desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength, …)

I/O modules (future: matilda/io/)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matilda/hdf5code.py
    HDF5 / NXcanSAS read-write helpers.
    Key functions: saveNXcanSAS, readMyNXcanSAS, readGenericNXcanSAS,
                   save_dict_to_hdf5, load_dict_from_hdf5, find_matching_groups.

matilda/readfromtiled.py
    Tiled server HTTP queries (REST API).
    Key functions: FindLastScanData, FindLastBlankScan, FindScanDataByName.
    Target server: http://usaxscontrol.xray.aps.anl.gov:8000 / catalog: usaxs_MongoDB.

Support modules (future: matilda/support/)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matilda/supportFunctions.py
    Shared numerical helpers: importFlyscan, calculatePD_Fly, beamCenterCorrection,
    smooth_r_data, getBlankFlyscan, normalizeByTransmission,
    calibrateAndSubtractFlyscan, rebinData, findProperBlankScan,
    read_group_to_dict, filter_nested_dict, subtract_data.

matilda/supportNikaFunctions.py
    Nika → Fit2D → pyFAI PONI geometry conversion.
    Key function: convert_Nika_to_Fit2D(*, SSD, pix_size, BCX, BCY, HorTilt, VertTilt, wavelength)

Plotting (future: matilda/plots/)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matilda/plotData.py
    Headless matplotlib plot export (JPEG files to web directory).
    Key functions: plotUSAXSResults, plotSWAXSResults.
    NOTE: matplotlib use is intentional here for headless/server output.
          All interactive GUI plotting must use pyqtgraph (see matilda/gui/).

GUI subpackage (planned — matilda/gui/)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
matilda/gui/__init__.py
    Placeholder.  Future PyQt6 / pyqtgraph GUI tools:
    - reduction_gui  : interactive 2-D → 1-D reduction with parameter control
    - survey_tool    : scan browser / overview
    - analysis_gui   : pyirena analysis launcher

Future integrations
-------------------
pynika (issue #28 — implemented)
    https://github.com/jilavsky/pynika
    Geometry auto-calibration from AgBehenateLaB6 standard patterns.
    Conda env: /home/beams/USAXS/.conda/envs/pynika

pyirena (planned)
    https://github.com/jilavsky/pyirena
    Post-reduction data analysis (model fitting, size distributions, etc.).
    Will be called after data reduction, similar to the pynika pattern.

Data dictionary schema
----------------------
All processXxx() functions return a list of result dicts with this structure::

    Sample = {
        "RawData": {
            "filename": str,
            "samplename": str,          # from HDF5 sample/name
            "metadata": dict,           # filtered subset of /entry/Metadata
            "instrument": dict,         # from /entry/instrument
            "sample": dict,             # from /entry/sample
            "control": dict,            # from /entry/control  (SAXS/WAXS only)
            "data": np.ndarray,         # raw 2-D image        (SAXS/WAXS only)
        },
        "reducedData": {
            "Q": np.ndarray,            # 1/Å
            "Intensity": np.ndarray,    # normalised, arbitrary units
            "Error": np.ndarray,
            "dQ": np.ndarray,           # USAXS only
        },
        "CalibratedData": {
            "Q": np.ndarray,            # 1/Å
            "Intensity": np.ndarray,    # cm²/cm³ (absolute scale)
            "Error": np.ndarray,
            "dQ": np.ndarray,
            "units": str,               # '[cm2/cm3]'
            "blankname": str,
            "thickness": float,         # mm
            "Kfactor": float or None,   # USAXS only
            "OmegaFactor": float or None,
            "SMR_Int": np.ndarray,      # USAXS only: slit-smeared data
            "SMR_Qvec": np.ndarray,
            "SMR_Error": np.ndarray,
            "SMR_dQ": np.ndarray,
            "slitLength": float,
        },
        "BlankData": {                  # only present when blank was subtracted
            "Q": np.ndarray,
            "Intensity": np.ndarray,
            "Error": np.ndarray,
        },
        "calib2DData": {                # SAXS/WAXS only
            "data": np.ndarray,         # calibrated 2-D image before integration
            "blankname": str,
            "transmission": float,
        },
    }

Known bugs / technical debt (do not fix without review)
---------------------------------------------------------
BUG-01  FIXED (feature/30-documentation).
        Bare module imports converted to relative imports (``from .xxx import``).
        main() function added to matilda.py; entry point enabled in pyproject.toml.
        Service launch: update serv_matilda.sh to use ``python -m matilda.matilda``
        (or the installed ``matilda`` console script).

BUG-02  FIXED (feature/30-documentation).
        ``import six`` removed from hdf5code.py; ``six.u()`` replaced with
        direct string references.

BUG-03  FIXED (feature/30-documentation).
        ``print(uri)`` replaced with ``logging.debug(f"{uri=}")`` in all three
        functions in readfromtiled.py.

BUG-04  FIXED (feature/30-documentation).
        Bare ``except:`` replaced with ``except Exception:`` in readfromtiled.py
        and ``except KeyError:`` in supportFunctions.py (key-lookup fallback).

BUG-05  FIXED (feature/30-documentation).
        3-space comment indentation corrected to 4 spaces in plotData.py.
        Unused ``import pprint as pp`` also removed from plotData.py.

BUG-06  FIXED (feature/30-documentation).
        Duplicate ``import pprint`` removed from convertFlyscan.py;
        only ``import pprint as pp`` remains.

BUG-07  FIXED (feature/30-documentation).
        matplotlib imports removed from convertFlyscan.py, convertUSAXS.py,
        convertSWAXS.py, and supportFunctions.py.  Active debug plt calls
        commented out (re-enable by uncommenting and adding the local import).
        supportFunctions.py: ``import matplotlib.pyplot as plt`` moved inside
        the ``if debugme:`` block so it is only loaded when debugging.

BUG-08  FIXED (feature/30-documentation).
        ``matilda/runTests.py`` deleted; tests/manual_test.py is the replacement.

BUG-09  FIXED (feature/30-documentation).
        Log directory now read from ``MATILDA_LOG_DIR`` env variable; default is
        ``~/.local/share/matilda/log/matilda.log`` (dev machines).
        serv_matilda.sh sets ``MATILDA_LOG_DIR=/share1/log/matilda`` for the service.
        Rotating log: 1 MB × 4 files = 4 MB max on disk.

BUG-10  FIXED (feature/30-documentation).
        hdf5code.py now reads version via importlib.metadata.version('Matilda');
        falls back to 'unknown' when the package is not installed.

Public processing API (importable after ``pip install -e .``)
-------------------------------------------------------------
BUG-01 is now fixed, so all imports below work directly::

    from matilda.matilda import processFlyscans, processStepscans, processADscans
    from matilda.matilda import processUSAXSFolder
    from matilda.convertFlyscan  import processFlyscan
    from matilda.convertUSAXS    import processStepscan
    from matilda.convertSWAXS    import process2Ddata
    from matilda.readfromtiled   import FindLastScanData, FindLastBlankScan
    from matilda.hdf5code        import saveNXcanSAS, readMyNXcanSAS
"""
