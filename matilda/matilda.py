#!/usr/bin/env python3
'''
Matilda - USAXS/SAXS/WAXS data processing tool

When run as a script, it will process the data from the last scan and blank scan

User facing functions are defined in this file.

processADscans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False) 
    - SAXS or WAXS data, finding the appropriate blank for each.
    - ListOfScans is a list of tuples with path and filename of the scan data
    - ListOfBlanks is a list of tuples with path and filename of the blank data
    - forceFirstBlank set to true will skip Blank selection and use the first blank in the list
    - recalculateAllData set to true will force reprocessing of the data, otherwise it
    - returns a list of dictionaries with the reduced data
    - No plotting, just processing

processFlyscans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False)
    - Flyscan data, finding the appropriate blank for each.
    - ListOfScans is a list of tuples with path and filename of the scan data
    - ListOfBlanks is a list of tuples with path and filename of the blank data
    - forceFirstBlank set to true will skip Blank selection and use the first blank in the list
    - recalculateAllData set to true will force reprocessing of the data, otherwise it
    - returns a list of dictionaries with the reduced data
    - No plotting, just processing
    
processStepScans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False)
    - Step scan data, finding the appropriate blank for each.
    - ListOfScans is a list of tuples with path and filename of the scan data
    - ListOfBlanks is a list of tuples with path and filename of the blank data
    - forceFirstBlank set to true will skip Blank selection and use the first blank in the list
    - recalculateAllData set to true will force reprocessing of the data, otherwise it
    - returns a list of dictionaries with the reduced data
    - No plotting, just processing

processUSAXSFolder(path)
    - will process (with forced reprocessing) all the scans in the given path (set recalculateAllData=True)
    - assumes USAXS data are in _usaxs folder
    - assumes SAXS data are in _saxs folder
    - assumes WAXS data are in _waxs folder
    - will process each data scan in the folder
    - no plotting, just processing

When run as main, it will process the data from the last 10 scans and blank scan
for each USAXSstep, Flyscan, SAXS, WAXS data types.


TODO: add to each processXYZ option to force Blank if Blanks is only one and avoid checking on order number. 
        This is to enable use for Igor or elsewhere to use different blank than was measured prior experiment. 

'''

import glob as _glob
import pprint as pp
import numpy as np
import socket
import re
import time
import logging
import datetime
from logging.handlers import RotatingFileHandler
import os
import subprocess


from .readfromtiled import FindLastScanData, FindLastBlankScan, FindLastTuneScans
from .convertFlyscan import processFlyscan
from .convertUSAXS import processStepscan
from .convertSWAXS import process2Ddata
from .convertTune import getTuneResults, TUNE_PLAN_NAMES, NumberOfTunesToShow, NumberOfDaysToLookBackTunes
from .supportFunctions import findProperBlankScan
from .plotData import plotUSAXSResults, plotSWAXSResults, plotTuneResults


#set ImagePath to None to prevent saving images
#imagePath = None  # Do not save images
imagePath = '/home/joule/WEBUSAXS/www_live/'  # Path to save images
#imagePath = '/home/parallels/Desktop/'  # Path to save images

# Conda setup for external tools (pynika, pyirena).
# Full path to the conda executable used at this beamline.
CONDA_EXECUTABLE = '/APSshare/miniconda/x86_64/bin/conda'
# Full path to pynika's conda environment (mirrors the pynika-gui launch script).
PYNIKA_CONDA_ENV_PATH = '/home/beams/USAXS/.conda/envs/pynika'
# Full path to pyirena's conda environment.
PYIRENA_CONDA_ENV_PATH = '/home/beams/USAXS/.conda/envs/pyirena'

# Name of the per-folder JSON config that triggers pyirena analysis.
_PYIRENA_CONFIG_FILENAME = 'pyirena_config.json'
# Name of the per-experiment JSON config that triggers USAXS+SAXS merging.
_MERGE_CONFIG_FILENAME = 'merge_config.json'
# Subfolder name for merged USAXS+SAXS output files.
_MERGED_FOLDER_NAME = 'data_usaxs_merged'
# Maximum number of (path, filename) entries kept in each per-technique
# analyzed-files set.  Large enough to cover user transitions without growing
# without bound.
_MAX_PYIRENA_ANALYZED = 100

# Regex pattern to detect AgBehenateLaB6 calibrant files (any capitalisation)
_CALIBRANT_PATTERN = re.compile(r'agbehenatelab6', re.IGNORECASE)

NumberOfDaysToLookBack = 1  # Number of days to look back for scans
NumberOfDaysToLookBackBlanks = 5  # Number of days to look back for blanks
NumberOfImagesInGraphs = 10  # Number of images to show in the graphs

#recalculateAllData = False  # Set to True to recalculate all data, False to use existing data
    
# Configure logging
# Log directory: use MATILDA_LOG_DIR env variable if set, otherwise default to
# ~/.local/share/matilda/log (standard XDG user data location).
# On the beamline service, serv_matilda.sh sets MATILDA_LOG_DIR=/share1/log/matilda.
_default_log_dir = os.path.join(os.path.expanduser("~"), ".local", "share", "matilda", "log")
log_dir = os.environ.get("MATILDA_LOG_DIR", _default_log_dir)
os.makedirs(log_dir, exist_ok=True)
log_file = os.path.join(log_dir, 'matilda.log')
# Rotating log: 1 MB per file, 3 backups → max 4 MB total on disk.
handler = RotatingFileHandler(log_file, maxBytes=1_000_000, backupCount=3)
logging.basicConfig(
    handlers=[handler],
    level=logging.INFO,
    #level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)


# --- pynika auto-calibration helpers ---

def _runPynikaCalibration(path, filename, instrument_type, calibrated_set):
    """
    Run pynika on a single AgBehenateLaB6 calibrant file to recalibrate the
    instrument geometry and push results to EPICS PVs.

    Parameters
    ----------
    path : str
        Directory containing the HDF5 file.
    filename : str
        HDF5 filename of the calibrant scan.
    instrument_type : str
        'SAXS' or 'WAXS' — passed to pynika --instrument flag.
    calibrated_set : set
        Mutable set of (path, filename) tuples already processed this session.
        Updated in-place on success to prevent re-running the same file.
    """
    file_key = (path, filename)
    if file_key in calibrated_set:
        logging.info(f"Pynika calibration already run for {filename}, skipping.")
        return

    filepath = os.path.join(path, filename)
    if not os.path.exists(filepath):
        logging.error(f"Calibrant file not found: {filepath}")
        return

    logging.info(f"Running pynika calibration: {filepath} [{instrument_type}]")
    cmd = [CONDA_EXECUTABLE, 'run', '--no-capture-output', '-p', PYNIKA_CONDA_ENV_PATH,
           'pynika', '--file', filepath, '--instrument', instrument_type, '--auto-fit', '--save-to-pvs']
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=300,
        )
        if result.returncode == 0:
            logging.info(f"Pynika calibration succeeded for {filename}")
            calibrated_set.add(file_key)
            # Keep the set bounded to match the Tiled look-back window
            while len(calibrated_set) > NumberOfImagesInGraphs:
                calibrated_set.pop()
        else:
            logging.error(
                f"Pynika calibration failed for {filename} "
                f"(exit {result.returncode}): {result.stderr.strip()}"
            )
    except subprocess.TimeoutExpired:
        logging.error(f"Pynika calibration timed out for {filename}")
    except Exception as e:
        logging.error(f"Pynika calibration error for {filename}: {e}", exc_info=True)


def _checkAndRunPynikaCalibration(ListOfScans, instrument_type, calibrated_set):
    """
    Scan a list of (path, filename) tuples.  For any file whose name contains
    'AgBehenateLaB6' (case-insensitive) that has not been calibrated yet this
    session, run pynika calibration before data reduction proceeds.
    """
    for scan_path, scan_filename in ListOfScans:
        if _CALIBRANT_PATTERN.search(scan_filename):
            _runPynikaCalibration(scan_path, scan_filename, instrument_type, calibrated_set)


# --- pyirena auto-analysis helpers ---

def _runPyirenaAnalysis(path, filename, analyzed_set):
    """
    Run pyirena fit_pyirena() on a single reduced data file if
    pyirena_config.json is present in the same folder.

    Parameters
    ----------
    path : str
        Directory containing the reduced HDF5 file.
    filename : str
        HDF5 filename of the reduced scan.
    analyzed_set : set
        Mutable set of (path, filename) tuples already analyzed this session.
        Updated in-place on success to prevent re-running the same file.
        Bounded to _MAX_PYIRENA_ANALYZED entries to limit memory use.
    """
    # Blanks are not analyzed.
    if 'blank' in filename.lower():
        return

    file_key = (path, filename)
    if file_key in analyzed_set:
        return

    config_file = os.path.join(path, _PYIRENA_CONFIG_FILENAME)
    if not os.path.isfile(config_file):
        return

    data_file = os.path.join(path, filename)
    if not os.path.exists(data_file):
        logging.error(f"Pyirena: data file not found: {data_file}")
        return

    logging.info(f"Running pyirena analysis: {data_file}")
    python_snippet = (
        "from pyirena.batch import fit_pyirena; "
        f"fit_pyirena({data_file!r}, {config_file!r}, "
        "save_to_nexus=True, with_uncertainty=False, n_mc_runs=10)"
    )
    cmd = [CONDA_EXECUTABLE, 'run', '--no-capture-output', '-p', PYIRENA_CONDA_ENV_PATH,
           'python', '-c', python_snippet]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            logging.info(f"Pyirena analysis succeeded for {filename}")
            analyzed_set.add(file_key)
            while len(analyzed_set) > _MAX_PYIRENA_ANALYZED:
                analyzed_set.pop()
        else:
            logging.error(
                f"Pyirena analysis failed for {filename} "
                f"(exit {result.returncode}): {result.stderr.strip()}"
            )
    except subprocess.TimeoutExpired:
        logging.error(f"Pyirena analysis timed out for {filename}")
    except Exception as e:
        logging.error(f"Pyirena analysis error for {filename}: {e}", exc_info=True)


def _checkAndRunPyirenaAnalysis(ListOfScans, analyzed_set):
    """
    For each (path, filename) in ListOfScans, run pyirena analysis if
    pyirena_config.json is present in that folder and the file has not
    been analyzed yet this session.  Blank files are skipped automatically.
    """
    for scan_path, scan_filename in ListOfScans:
        _runPyirenaAnalysis(scan_path, scan_filename, analyzed_set)


# --- USAXS+SAXS auto-merge helpers ---

def _extract_sample_key(filename):
    """
    Extract a (prefix, scan_number) key from a data filename.

    The prefix is everything before the first '_' and the scan number is
    the token between the last '_' and the file extension.  Two files
    represent the same sample when both parts match.

    Example: 'Sample1_55C_10min_1234.h5' → ('Sample1', '1234')

    Returns None if the filename does not contain at least one '_'.
    """
    stem = os.path.splitext(filename)[0]  # strip extension
    parts = stem.split('_')
    if len(parts) < 2:
        return None
    return (parts[0], parts[-1])


def _find_matching_partner(path, filename, partner_suffix):
    """
    Given a reduced data file, find the matching file in the partner
    technique folder.

    Parameters
    ----------
    path : str
        Directory of the current file (e.g. '.../data_usaxs').
    filename : str
        Filename of the current reduced scan.
    partner_suffix : str
        Folder suffix of the partner technique ('_usaxs' or '_saxs').

    Returns
    -------
    (partner_path, partner_filename) or None
    """
    key = _extract_sample_key(filename)
    if key is None:
        return None
    prefix, scan_number = key

    parent_dir = os.path.dirname(path)
    # Build partner folder: replace the last _xxx suffix with partner_suffix.
    current_folder_name = os.path.basename(path)
    # Strip the technique suffix (e.g. '_usaxs', '_saxs') from the folder name.
    # The folder name is like 'data_usaxs' or '05_02_UserName_saxs' — the
    # technique suffix is always at the end.
    for suffix in ('_usaxs', '_saxs'):
        if current_folder_name.endswith(suffix):
            base_name = current_folder_name[:-len(suffix)]
            break
    else:
        return None  # folder doesn't follow expected naming

    partner_folder = os.path.join(parent_dir, base_name + partner_suffix)
    if not os.path.isdir(partner_folder):
        return None

    # Glob for HDF5 files matching the same prefix and scan number.
    # Only consider .h5, .hdf, .hdf5 extensions — exclude images (.jpg, .tiff, etc.).
    matches = []
    for ext in ('h5', 'hdf', 'hdf5'):
        pattern = os.path.join(partner_folder, f"{prefix}_*_{scan_number}.{ext}")
        matches.extend(_glob.glob(pattern))
    if not matches:
        return None

    partner_file = os.path.basename(matches[0])
    return (partner_folder, partner_file)


def _runMergeData(usaxs_path, usaxs_filename, saxs_path, saxs_filename,
                  config_file, merged_set):
    """
    Merge a USAXS+SAXS file pair using pyirena merge_data().

    USAXS is always file1 (lower Q, absolute scale), SAXS is file2
    (higher Q).  On both success and failure the merge key is added to
    *merged_set* so the pair is not retried this session.

    Parameters
    ----------
    usaxs_path, usaxs_filename : str
        Path and filename of the reduced USAXS file.
    saxs_path, saxs_filename : str
        Path and filename of the reduced SAXS file.
    config_file : str
        Full path to the merge_config.json file.
    merged_set : set
        Mutable set of merge keys already processed this session.

    Returns
    -------
    bool
        True if the merge succeeded, False otherwise.
    """
    usaxs_file = os.path.join(usaxs_path, usaxs_filename)
    saxs_file = os.path.join(saxs_path, saxs_filename)
    merge_key = (usaxs_file, saxs_file)

    parent_dir = os.path.dirname(usaxs_path)
    output_folder = os.path.join(parent_dir, _MERGED_FOLDER_NAME)
    os.makedirs(output_folder, exist_ok=True)

    logging.info(f"Merging USAXS+SAXS: {usaxs_filename} + {saxs_filename}")
    python_snippet = (
        "from pyirena.batch import merge_data; "
        f"merge_data({usaxs_file!r}, {saxs_file!r}, "
        f"config_file={config_file!r}, "
        f"output_folder={output_folder!r}, "
        "save_to_nexus=True, verbose=True)"
    )
    cmd = [CONDA_EXECUTABLE, 'run', '--no-capture-output', '-p',
           PYIRENA_CONDA_ENV_PATH, 'python', '-c', python_snippet]
    success = False
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
        if result.returncode == 0:
            logging.info(f"Merge succeeded: {usaxs_filename} + {saxs_filename}")
            success = True
        else:
            logging.error(
                f"Merge failed for {usaxs_filename} + {saxs_filename} "
                f"(exit {result.returncode}): {result.stderr.strip()}"
            )
    except subprocess.TimeoutExpired:
        logging.error(
            f"Merge timed out for {usaxs_filename} + {saxs_filename}"
        )
    except Exception as e:
        logging.error(
            f"Merge error for {usaxs_filename} + {saxs_filename}: {e}",
            exc_info=True,
        )

    # Always record the pair so it is not retried this session.
    merged_set.add(merge_key)
    while len(merged_set) > _MAX_PYIRENA_ANALYZED:
        merged_set.pop()

    return success


def _checkAndRunMerge(ListOfScans, technique, merged_set, analyzed_merged_set):
    """
    For each scan in *ListOfScans*, check whether the matching partner
    file (USAXS↔SAXS) exists and a merge_config.json is present.  If
    so, merge the pair and optionally run pyirena analysis on the output.

    Parameters
    ----------
    ListOfScans : list of (path, filename)
        Scans just processed in the current cycle.
    technique : str
        'USAXS' or 'SAXS' — indicates which technique *ListOfScans*
        belongs to.
    merged_set : set
        Tracks (usaxs_file, saxs_file) pairs already merged this session.
    analyzed_merged_set : set
        Tracks merged output files already sent to pyirena this session.
    """
    if technique == "USAXS":
        partner_suffix = '_saxs'
    elif technique == "SAXS":
        partner_suffix = '_usaxs'
    else:
        return

    for scan_path, scan_filename in ListOfScans:
        if 'blank' in scan_filename.lower():
            continue

        partner = _find_matching_partner(scan_path, scan_filename,
                                         partner_suffix)
        if partner is None:
            continue
        partner_path, partner_filename = partner

        # Determine which is USAXS and which is SAXS.
        if technique == "USAXS":
            usaxs_path, usaxs_filename = scan_path, scan_filename
            saxs_path, saxs_filename = partner_path, partner_filename
        else:
            saxs_path, saxs_filename = scan_path, scan_filename
            usaxs_path, usaxs_filename = partner_path, partner_filename

        usaxs_file = os.path.join(usaxs_path, usaxs_filename)
        saxs_file = os.path.join(saxs_path, saxs_filename)
        merge_key = (usaxs_file, saxs_file)
        if merge_key in merged_set:
            continue

        # Check for merge_config.json in the parent data folder.
        parent_dir = os.path.dirname(usaxs_path)
        config_file = os.path.join(parent_dir, _MERGE_CONFIG_FILENAME)
        if not os.path.isfile(config_file):
            continue

        success = _runMergeData(usaxs_path, usaxs_filename,
                                saxs_path, saxs_filename,
                                config_file, merged_set)

        if success:
            # Run pyirena analysis on the merged output if config exists.
            merged_folder = os.path.join(parent_dir, _MERGED_FOLDER_NAME)
            # The merged file is named after the USAXS (file1) input.
            merged_filename = usaxs_filename
            _runPyirenaAnalysis(merged_folder, merged_filename,
                                analyzed_merged_set)


# --- user facing functions ---

def processUSAXSFolder(path):
    """
    Batch-reprocess all USAXS/SAXS/WAXS scans found under a single date folder.

    Intended for offline re-reduction (e.g. from Igor or a notebook) where you
    want to force-recalculate every file regardless of cached results.  Expects
    the standard APS folder layout produced by the Bluesky acquisition system:

        <path>/
            <name>_usaxs/   — USAXS flyscan files (.h5)
            <name>_saxs/    — SAXS area-detector files (.hdf)
            <name>_waxs/    — WAXS area-detector files (.hdf)

    All three sub-folder types must be present.  Files are sorted by the
    trailing scan number embedded in the filename (see extract_number_from_filename).

    Parameters
    ----------
    path : str
        Absolute path to the date/experiment folder that contains the three
        instrument sub-folders.

    Notes
    -----
    * recalculateAllData is forced True — existing cached results are deleted
      and recomputed from raw HDF5 data.
    * No plotting is performed; results are only written back to the HDF5 files.
    * The function silently returns if path does not exist.

    TODO: pynika auto-calibration is not called here (only in the live loop).
    """
    # Get the list of folders in the path folder
    if not os.path.exists(path):
        logging.error(f"The path {path} does not exist.")
        return
    # Get the list of folders in the path folder
    folders = os.listdir(path)
    # Filter the list to include only the folders that end with _usaxs
    usaxs_folder = [folder for folder in folders if folder.endswith('_usaxs')][0]
    
    saxs_folder = [folder for folder in folders if folder.endswith('_saxs')][0]
    
    waxs_folder = [folder for folder in folders if folder.endswith('_waxs')][0]
    
    #get list of files in the usaxs folder
    usaxs_files = []
    files = os.listdir(os.path.join(path, usaxs_folder))
    # Filter the list to include only the files with the extension .h5
    # Sort by the number before .hdf
    h5_files = [file for file in files if file.endswith('.h5')]
    usaxs_files = sorted(h5_files, key=extract_number_from_filename)
    usaxs_blanks = [file for file in usaxs_files if 'blank' in file.lower()]

    #get list of files in the saxs folder
    saxs_files = []
    files = os.listdir(os.path.join(path, saxs_folder))
    # Filter the list to include only the files with the extension .hdf
    # Sort by the number before .hdf
    h5_files = [file for file in files if file.endswith('.hdf')]
    saxs_files = sorted(h5_files, key=extract_number_from_filename)
    saxs_blanks = [file for file in saxs_files if 'blank' in file.lower()]

    
    #get list of files in the waxs folder
    waxs_files = []
    files = os.listdir(os.path.join(path, waxs_folder))
    # Filter the list to include only the files with the extension .hdf
    # Sort by the number before .hdf
    h5_files = [file for file in files if file.endswith('.hdf')]
    waxs_files = sorted(h5_files, key=extract_number_from_filename)
    waxs_blanks = [file for file in waxs_files if 'blank' in file.lower()]
           
    #process all flyscans
    logging.info("Processing USAXS Flyscans")
    ListOfScans = [(os.path.join(path, usaxs_folder), file) for file in usaxs_files]
    ListOfBlanks = [(os.path.join(path, usaxs_folder), file) for file in usaxs_blanks]
    data = processFlyscans(ListOfScans, ListOfBlanks, recalculateAllData=True)

    logging.info("Processing SAXS data")
    ListOfScans = [(os.path.join(path, saxs_folder), file) for file in saxs_files]
    ListOfBlanks = [(os.path.join(path, saxs_folder), file) for file in saxs_blanks]
    data = processADscans(ListOfScans, ListOfBlanks, recalculateAllData=True)


    logging.info("Processing WAXS data")
    ListOfScans = [(os.path.join(path, waxs_folder), file) for file in waxs_files]
    ListOfBlanks = [(os.path.join(path, waxs_folder), file) for file in waxs_blanks]
    data = processADscans(ListOfScans, ListOfBlanks, recalculateAllData=True)
    



#process list of flyscans, finding the appropriate blank for each.
def processFlyscans(ListOfScans, ListOfBlanks, recalculateAllData=False,forceFirstBlank=False):
    """
    Processes a list of flyscans, finding the appropriate blank for each.
    The correct blank is in the same path and has the number _XYZ
    that is smaller than the sample's number.
    forceFirstBlank=True will skip Blank selection and use the first blank in the list
    recalculateAllData=True will force reprocessing of the data, otherwise it
    """
    results=[]
    logging.info("Processing of flyscans with blanks")    
    for scan_path, scan_filename in ListOfScans:
        try:
            if forceFirstBlank:
                # Use the first blank in the list
                selected_blank_path = ListOfBlanks[0][0]
                selected_blank_filename = ListOfBlanks[0][1]
                #logging.info(f"User forced use of blank path: {selected_blank_path}, filename: {selected_blank_filename}")
            else:
                selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
                 #logging.info(f"Automatically selected blank path: {selected_blank_path}, filename: {selected_blank_filename}")
           
            logging.info(f"Flyscan data processing: path {scan_path}, filename {scan_filename}, blank path {selected_blank_path}, blank filename {selected_blank_filename}")
            result = processFlyscan(scan_path, scan_filename,
                                        blankPath=selected_blank_path,
                                        blankFilename=selected_blank_filename,
                                        recalculateAllData=recalculateAllData)          # recalculateAllData will be True for reprocessing
            results.append(result)
        except Exception as e:
            logging.error(f"Error processing scan {scan_filename}: {e}", exc_info=True)
    return results


# Process the step scan data files
def processStepscans(ListOfScans, ListOfBlanks, recalculateAllData=False, forceFirstBlank=False):
    """
    Process a list of USAXS step-scan files, pairing each with an appropriate blank.

    Mirrors processFlyscans() but calls processStepscan() (convertUSAXS module)
    instead of processFlyscan().  The blank-selection logic is identical.

    Parameters
    ----------
    ListOfScans : list of (str, str)
        (path, filename) tuples for step-scan HDF5 files to process.
    ListOfBlanks : list of (str, str)
        (path, filename) tuples for available blank/background measurements.
    recalculateAllData : bool, optional
        When True, delete any cached reduced data and recompute from raw.
        Default False (reuse cached results if present).
    forceFirstBlank : bool, optional
        When True, always pair each scan with ListOfBlanks[0] instead of
        searching for the nearest-preceding blank by scan number.
        Default False.

    Returns
    -------
    list of dict
        One result dictionary per successfully processed scan (same structure
        as returned by processStepscan).  Failed scans are logged and skipped.
    """
    results=[]
    logging.info("Processing of Step scans with blanks")    
    for scan_path, scan_filename in ListOfScans:
        try:
            if forceFirstBlank:
                # Use the first blank in the list
                selected_blank_path = ListOfBlanks[0][0]
                selected_blank_filename = ListOfBlanks[0][1]
                #logging.info(f"User forced use of blank path: {selected_blank_path}, filename: {selected_blank_filename}")
            else:
                selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
                 #logging.info(f"Automatically selected blank path: {selected_blank_path}, filename: {selected_blank_filename}")
           
            logging.info(f"Stepscan data processing: path {scan_path}, filename {scan_filename}, blank path {selected_blank_path}, blank filename {selected_blank_filename}")
            result = processStepscan(scan_path, scan_filename,
                                        blankPath=selected_blank_path,
                                        blankFilename=selected_blank_filename,
                                        recalculateAllData=recalculateAllData)          # recalculateAllData will be True for reprocessing
            results.append(result)
        except Exception as e:
            logging.error(f"Error processing scan {scan_filename}: {e}", exc_info=True)
    return results



#process SAXS/WAXS data, finding the appropriate blank for each.
def processADscans(ListOfScans, ListOfBlanks,recalculateAllData=False,forceFirstBlank=False):
    """
    Processes a list of SAXS/WAXS data, finding the appropriate blank for each.
    The correct blank is in the same path and has the largest number XYZ
    that is smaller than the sample's number.
    """
    results=[]
    #logging.info("Processing of SAXS/WAXS with blanks")    
    logging.info("Processing of SAXS/WAXS with blanks")    
    for scan_path, scan_filename in ListOfScans:
        try:
            if forceFirstBlank:
                # Use the first blank in the list
                selected_blank_path = ListOfBlanks[0][0]
                selected_blank_filename = ListOfBlanks[0][1]
                #logging.info(f"User forced use of blank path: {selected_blank_path}, filename: {selected_blank_filename}")
            else:
                selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
                #logging.info(f"Automatically selected blank path: {selected_blank_path}, filename: {selected_blank_filename}")

            logging.info(f"2D data processing: path {scan_path}, filename {scan_filename}, blank path {selected_blank_path}, blank filename {selected_blank_filename}")
            result = process2Ddata(scan_path, scan_filename,
                                        blankPath=selected_blank_path,
                                        blankFilename=selected_blank_filename,
                                        recalculateAllData=recalculateAllData)          # recalculateAllData will be True for reprocessing
            results.append(result)
        except Exception as e:
            logging.error(f"Error processing AD scan {scan_filename}: {e}", exc_info=True)
    return results


def extract_number_from_filename(filename):
    """Extract number from filename, return 0 if no number found"""
    match = re.search(r'_(\d+)\.hdf', filename)
    return int(match.group(1)) if match else 0



def main():
    """Entry point for the Matilda service (15-second polling loop).

    Invoked by the ``matilda`` console script installed by pyproject.toml,
    or directly via ``python -m matilda.matilda``.
    """
    try:
        listofFlyscansOld=dict()
        listofStepScansOld=dict()
        listofSAXSOld=dict()
        listOfWAXSOld=dict()
        # Track AgBehenateLaB6 calibrant files already processed by pynika this
        # session.  Separate sets for SAXS and WAXS because files share the same
        # name across detectors but live in different folders.
        calibratedSAXSFiles = set()
        calibratedWAXSFiles = set()
        # Track files already submitted to pyirena analysis this session.
        # Separate sets per technique because each folder has its own config.
        analyzedFlyscanFiles = set()
        analyzedStepFiles = set()
        analyzedSAXSFiles = set()
        analyzedWAXSFiles = set()
        # Track USAXS+SAXS pairs already merged and merged files already
        # analyzed this session.  Bounded like the sets above.
        mergedUSAXSSAXSFiles = set()
        analyzedMergedFiles = set()
        # Track the uid list of the last-plotted tune scans per plan type so we
        # only re-download and re-plot when a new tune scan appears.
        listOfTunesOld = {pn: [] for pn in TUNE_PLAN_NAMES}
        while True:
            logging.info("New round of processing started at : %s", datetime.datetime.now()) 

            logging.info("Processing the Flyscans")
            ListOfScans = FindLastScanData("Flyscan",NumberOfImagesInGraphs,NumberOfDaysToLookBack)
            #need to check, if we have same list as before and if yes, skip everything
            if ListOfScans == listofFlyscansOld or len(ListOfScans) == 0:
                logging.info('No new Flyscan data found')
            else:
                #path, filename = ListOfScans[-1]
                listOfBlanks = FindLastBlankScan("Flyscan",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBackBlanks) 
                #listOfBlanks = FindLastBlankScan("Flyscan",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBack) 
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                results = processFlyscans(ListOfScans, listOfBlanks)
                plotUSAXSResults(results, imagePath, isFlyscan=True)
                _checkAndRunPyirenaAnalysis(ListOfScans, analyzedFlyscanFiles)
                _checkAndRunMerge(ListOfScans, "USAXS", mergedUSAXSSAXSFiles, analyzedMergedFiles)
                listofFlyscansOld = ListOfScans
            

            logging.info("Processing the Step scans")
            ListOfScans = FindLastScanData("uascan",NumberOfImagesInGraphs,NumberOfDaysToLookBack)
            #need to check, if we have same list as before and if yes, skip everything
            if ListOfScans == listofStepScansOld or len(ListOfScans) == 0:
                logging.info('No new Step scans data found')
            else:
                #path, filename = ListOfScans[-1]
                listOfBlanks = FindLastBlankScan("uascan",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBackBlanks)
                #listOfBlanks = FindLastBlankScan("uascan",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBack)
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                results = processStepscans(ListOfScans, listOfBlanks)
                plotUSAXSResults(results, imagePath, isFlyscan=False)
                _checkAndRunPyirenaAnalysis(ListOfScans, analyzedStepFiles)
                _checkAndRunMerge(ListOfScans, "USAXS", mergedUSAXSSAXSFiles, analyzedMergedFiles)
                listofStepScansOld = ListOfScans
    
            #process SAXS and WAXS data
            logging.info("Processing the SAXS")
            ListOfScans = FindLastScanData("SAXS",NumberOfImagesInGraphs,NumberOfDaysToLookBack)
            if ListOfScans == listofSAXSOld or len(ListOfScans) == 0:
                logging.info('No new SAXS data found')
            else:
                listOfBlanks = FindLastBlankScan("SAXS",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBackBlanks)
                #listOfBlanks = FindLastBlankScan("SAXS",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBack)
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                _checkAndRunPynikaCalibration(ListOfScans, "SAXS", calibratedSAXSFiles)
                results = processADscans(ListOfScans, listOfBlanks)
                plotSWAXSResults(results, imagePath, isSAXS=True)
                _checkAndRunPyirenaAnalysis(ListOfScans, analyzedSAXSFiles)
                _checkAndRunMerge(ListOfScans, "SAXS", mergedUSAXSSAXSFiles, analyzedMergedFiles)
                listofSAXSOld = ListOfScans

            logging.info("Processing the WAXS")
            ListOfScans = FindLastScanData("WAXS",NumberOfImagesInGraphs,NumberOfDaysToLookBack)
            if ListOfScans == listOfWAXSOld or len(ListOfScans) == 0:
                 #if len(ListOfScans) == 0:
                logging.info('No new WAXS data found')
            else:
                listOfBlanks = FindLastBlankScan("WAXS",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBackBlanks)
                #listOfBlanks = FindLastBlankScan("WAXS",path=None, NumScans=NumberOfImagesInGraphs, lastNdays=NumberOfDaysToLookBack)
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                _checkAndRunPynikaCalibration(ListOfScans, "WAXS", calibratedWAXSFiles)
                results = processADscans(ListOfScans, listOfBlanks)
                plotSWAXSResults(results, imagePath, isSAXS=False)
                _checkAndRunPyirenaAnalysis(ListOfScans, analyzedWAXSFiles)
                listOfWAXSOld = ListOfScans

            # Process live tuning scans (tune_ar, tune_mr, tune_a2rp).
            # These have no HDF5 file; data comes from Tiled only.  Safe against
            # a missing/offline server: FindLastTuneScans returns [] and we skip.
            logging.info("Processing the Tune scans")
            for pn in TUNE_PLAN_NAMES:
                tuneMeta = FindLastTuneScans(pn, NumberOfTunesToShow, NumberOfDaysToLookBackTunes)
                uidList = [m["uid"] for m in tuneMeta]
                if uidList == listOfTunesOld[pn] or len(uidList) == 0:
                    logging.info(f'No new {pn} tune data found')
                    continue
                results = getTuneResults(pn, NumberOfTunesToShow, NumberOfDaysToLookBackTunes)
                plotTuneResults(results, imagePath, pn)
                listOfTunesOld[pn] = uidList

            logging.info('Sleeping for 15 seconds')
            time.sleep(15)
    except KeyboardInterrupt:
            logging.info('Keyboard interrupt')


if __name__ == "__main__":
    main()
