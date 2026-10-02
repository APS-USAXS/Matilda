"""
convertUSAXS.py
===============
Reduce USAXS step-scan HDF5 files to calibrated 1-D I(Q) data.

Main entry point
----------------
processStepscan(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False)

Also exports
------------
rebinData(Q, Intensity, Error, numBins)   — used by convertFlyscan.py

Data flow
---------
importFlyscan() from supportFunctions    — read raw arrays (step and fly share format)
calculatePD_Fly() / beamCenterCorrection()
calibrateAndSubtractFlyscan()
desmearData()
saveNXcanSAS() / readMyNXcanSAS()         — cache in HDF5

Notes
-----
* Despite the module name, the reduction pipeline is nearly identical to
  convertFlyscan.py because step-scan and fly-scan share the same NXsas
  file format at USAXS.
* matplotlib is imported; plt.show() calls are all commented out.
  Target for removal when GUI work begins.
* rebinData() is also called by convertFlyscan — it lives here due to
  historical code organisation but logically belongs in supportFunctions.

TODO: fix calibration and saving (original TODO from module creation).
TODO: reduceStepScanToQR and reduceFlyscanToQR mentioned in original header
      may be legacy names — verify against current function names.
"""

import os
import re
import h5py
import numpy as np
import logging
from .fx4support import (CHAIN_FX4, FX4_RELATIVE_CURRENT_ERROR, as_scalar,
                         detect_counting_chain, range_indexed_array, ratio_error)
from .supportFunctions import read_group_to_dict, filter_nested_dict, check_arrays_same_length
from .supportFunctions import beamCenterCorrection
from .supportFunctions import calibrateAndSubtractFlyscan, load_dict_from_hdf5, save_dict_to_hdf5
from .hdf5code import saveNXcanSAS, readMyNXcanSAS
from .hdf5code import clearAndCheckCachedReduction, writeThicknessOverride
from .supportFunctions import empty_calibrated_data
from .supportFunctions import normalizeByTransmission, transmissionTerms
from .desmearing import desmearData
from .plotData import plotUSAXSResults


# CONFIRMED 2026-07-08 (JIL): 1e7 Hz = Joerger scaler internal clock, correct
# for step scans (flyscans use a 1e6 Hz MCA clock instead — see
# supportFunctions.FLYSCAN_MCA_CLOCK_HZ).  Both are correct for their
# geometry; do not unify.  There is no clock on the FX4 chain: the
# electrometer reports a mean current and importStepScan supplies seconds.
JOERGER_CLOCK_HZ = 1e7


# This code first reduces data to QR and if provided with Blank, it will do proper data calibration, subtraction, and even desmearing
# It will check if QR/NXcanSAS data exist and if not, it will create properly calibrated NXcanSAS in the Nexus file
# If exist and recalculateAllData is False, it will reuse old ones. This is done for plotting.
def processStepscan(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False,
                     desmear_iter=20, extrap_method='PowerLaw w flat',
                     extrap_qstart=0.15, minQMinFindRatio=1.05, thickness_override=None,
                     use_mu=False, mu=None, per_gram=False, density=None,
                     transmission_override=None, qmin_override=None):
    """Reduce a single USAXS step-scan HDF5 file to calibrated 1-D I(Q).

    Structurally identical to processFlyscan() (convertFlyscan module) —
    see that function's docstring for parameter and return-value details.
    The distinction is that step scans use point-by-point detector readout
    rather than continuous motion, which affects the raw-data structure read
    by importFlyscan().

    Parameters
    ----------
    path : str
        Directory containing the sample HDF5 file.
    filename : str
        Filename of the sample HDF5 file (.h5).
    blankPath : str or None, optional
        Directory of blank HDF5 file.  None → only raw QR data, no calibration.
    blankFilename : str or None, optional
        Filename of blank HDF5 file.  None → only raw QR data.
    recalculateAllData : bool, optional
        True → delete cached NXcanSAS data and recompute.  Default False.
    desmear_iter : int, optional
        Maximum Lake/Strobl desmearing iterations.  Default 20.
    extrap_method : str, optional
        High-Q extrapolation method for desmearing.  Default 'PowerLaw w flat'.
    extrap_qstart : float, optional
        Q value above which extrapolation is applied.  Default 0.15 Å⁻¹.
    minQMinFindRatio : float, optional
        Threshold for Q-minimum selection after blank subtraction.  Default 1.05.
    thickness_override : float or None, optional
        When not None, use this value (mm) instead of the HDF5 thickness.

    Returns
    -------
    dict
        Same structure as processFlyscan(): RawData, reducedData,
        CalibratedData (when blank provided).
    """
    # Open the HDF5 file in read/write mode
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r+') as hdf_file:
        # Cache bookkeeping (shared with processFlyscan): with a blank we
        # require the full desmeared NXcanSAS entry, otherwise QRS_data is enough.
        requireCalibrated = (blankFilename is not None and blankPath is not None
                             and "blank" not in filename.lower())
        if clearAndCheckCachedReduction(hdf_file, filename, recalculateAllData, requireCalibrated):
            # exists, so lets reuse the data from the file
            Sample = readMyNXcanSAS(path, filename, isUSAXS=True)
            logging.info(f"Using existing processed data from file {filename}.")
            return Sample

        else:
            Sample = dict()
            if thickness_override is not None:
                writeThicknessOverride(hdf_file,
                                       '/entry/instrument/bluesky/metadata/sample_thickness_mm',
                                       thickness_override, filename)
            Sample["RawData"]=importStepScan(path, filename)                #import data
            Sample["reducedData"]=(createUPDGainsAndBkgErrArrays(Sample))
            Sample["reducedData"].update(CorrectUPDGainsStep(Sample))    # Correct UPD gains=CorrectUPDGainsStep(Sample)    # Correct UPD gains, this is the first step in data reduction
            Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0))
            Sample["reducedData"].update(calculatePDErrorStep(Sample))          # Calculate UPD error, mostly the same as in Igor                
            # Sample["reducedData"].update(smooth_r_data(Sample["reducedData"]["Intensity"],     #smooth data data
            #                                         Sample["reducedData"]["Q"], 
            #                                         Sample["reducedData"]["UPD_gains"], 
            #                                         Sample["reducedData"]["Error"], 
            #                                         Sample["RawData"]["TimePerPoint"],
            #                                         replaceNans=True))                 

            if (
                blankPath is not None
                and blankFilename is not None
                and blankFilename != filename
                and "blank" not in filename.lower()
            ):
                # pass recalculateAllData through so a forced reprocess also
                # invalidates the cached blank (was hardcoded False before)
                Sample["BlankData"]=getBlankStepscan(blankPath, blankFilename,recalculateAllData=recalculateAllData)
                Sample["reducedData"].update(normalizeByTransmission(Sample))          # Normalize sample by dividing by transmission for subtraction
                Sample["CalibratedData"]=(calibrateAndSubtractFlyscan(Sample, minQMinFindRatio=minQMinFindRatio, thickness_override=thickness_override, use_mu=use_mu, mu=mu, per_gram=per_gram, density=density, transmission_override=transmission_override, qmin_override=qmin_override))
                Sample["CalibratedData"].update(calculatedQStep(Sample))
                SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                if len(SMR_Qvec) > 50:  # some data were found. Call this success? 
                    # NOTE: no binning down done for step scans, if we collected many points, we need them. 
                    slitLength=Sample["CalibratedData"]["slitLength"]
                    SMR_Int =Sample["CalibratedData"]["SMR_Int"]
                    SMR_Error =Sample["CalibratedData"]["SMR_Error"]
                    SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                    SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]
                    DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ = desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=slitLength,ExtrapMethod=extrap_method,ExtrapQstart=extrap_qstart, MaxNumIter=desmear_iter)
                    desmearedData={
                        "Intensity":DSM_Int,
                        "Q":DSM_Qvec,
                        "Error":DSM_Error,
                        "dQ":DSM_dQ,
                        "units":Sample["CalibratedData"]["units"],
                        }
                    Sample["CalibratedData"].update(desmearedData)
                else:
                    logging.warning(f"Not enough data points in SMR_Qvec ({len(SMR_Qvec)}) to proceed with desmearing or rebinning. "
                                    "Skipping desmearing and rebinning steps. ")
                    #set calibrated data in the structure to None
                    Sample["CalibratedData"] = empty_calibrated_data()

            else:
                #set calibrated data in the structure to None
                Sample["CalibratedData"] = empty_calibrated_data()
        # Ensure all changes are written and close the HDF5 file
        hdf_file.flush()
    # The 'with' statement will automatically close the file when the block ends
    saveNXcanSAS(Sample,path, filename)
    return Sample

def getBlankStepscan(blankPath, blankFilename, recalculateAllData=False):
      # will reduce the blank linked as input into Igor BL_R_Int 
      # after reducing this first time, data are saved in Nexus file for subsequent use. 
      # We get the BL_QRS and calibration data as result.
    # Open the HDF5 file in read/write mode
    location = 'entry/blankData/'
    Filepath = os.path.join(blankPath, blankFilename)
    with h5py.File(Filepath, 'r+') as hdf_file:
            # Check if the group 'location' exists, if yes, either delete if asked for or use. 
            if recalculateAllData:
                if location in hdf_file:
                    # Delete the group is exists and requested
                    del hdf_file[location]
                    logging.info(f"Deleted existing group 'entry/blankData' in {blankFilename}.")

            if location in hdf_file:
                # exists, so lets reuse the data from the file
                Blank = dict()
                Blank = load_dict_from_hdf5(hdf_file, location)
                logging.info(f"Used existing Blank data from {blankFilename}")
                return Blank
            else:
                Blank = dict()
                Blank["RawData"]=importStepScan(blankPath, blankFilename)         #import data
                (BlTransCounts, BlTransGain,
                 BlI0Counts, BlI0Gain) = transmissionTerms(Blank['RawData'])
                Blank["BlankData"]= (createUPDGainsAndBkgErrArrays(Blank))  
                Blank["BlankData"].update(CorrectUPDGainsStep(Blank))       # Creates Intensity with corrected gains and background subtraction
                Blank["BlankData"].update(calculatePDErrorStep(Blank, isBlank=True))          # Calculate UPD error, mostly the same as in Igor                
                Blank["BlankData"].update(beamCenterCorrection(Blank,useGauss=0, isBlank=True)) #Beam center correction
                Blank["BlankData"].update({"blankname":blankFilename})      # add the name of the blank file
                Blank["BlankData"].update({"BlTransCounts":BlTransCounts})  # add the BlTransCounts
                Blank["BlankData"].update({"BlTransGain":BlTransGain})      # add the BlTransGain
                Blank["BlankData"].update({"BlI0Counts":BlI0Counts})        # add the BlI0Counts
                Blank["BlankData"].update({"BlI0Gain":BlI0Gain})            # add the BlTransGain
                # Blank["BlankData"].update(smooth_r_data(Blank["BlankData"]["Intensity"],     #smooth data data
                #                                         Blank["BlankData"]["Q"], 
                #                                         Blank["BlankData"]["UPD_gains"], 
                #                                         Blank["BlankData"]["Error"], 
                #                                         Blank["RawData"]["TimePerPoint"],
                #                                         replaceNans=True )) 
                # we need to return just the BlankData part 
                BlankData=dict()
                BlankData=Blank["BlankData"]
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                save_dict_to_hdf5(BlankData, location, hdf_file)
                logging.info(f"Appended new Blank data to 'entry/blankData' in {blankFilename}.")
                return BlankData


def createUPDGainsAndBkgErrArrays(Sample):
    """Per-point amplifier gain and dark-current error for a step scan."""
    if Sample["RawData"].get("chain") == CHAIN_FX4:
        # Gain-independent picoamps: the gain is 1 everywhere and the dark
        # current's error is the FX4 sequence program's bkgErr for the range
        # in use at that point — already in pA, with no dwell-time factor.
        n_points = len(Sample["RawData"]["UPD_array"])
        return {"UPD_gains": np.ones(n_points),
                "UPD_bkgErr": range_indexed_array(Sample["RawData"].get("RangeIndex"),
                                                  Sample["RawData"].get("Bkg_err_map") or {},
                                                  n_points)}
    # Create UPD_gains and UPD_bkgErr arrays based on the AmpGain values
    AmpGain = Sample["RawData"]["AmpGain"]
    Bkg_map = Sample["RawData"]["Bkg_map"]
    TimePerPoint = Sample["RawData"]["TimeInSec"]
    UPD_gains = np.zeros_like(AmpGain, dtype=float)
    UPD_bkgErr = np.zeros_like(AmpGain, dtype=float)

    # Assign values based on AmpGain.  Match with tolerance (EPICS-sourced
    # floats may not be bit-exact) and warn on unknown gains, which would
    # otherwise silently leave UPD_gains 0 (division by zero downstream).
    known_gains = {1e4: "1e4", 1e6: "1e6", 1e8: "1e8", 1e10: "1e10", 1e12: "1e12"}
    unknown_gains = set()
    for i, gain in enumerate(AmpGain):
        matched_key = None
        for gval, gkey in known_gains.items():
            if np.isclose(gain, gval, rtol=1e-3):
                matched_key = gkey
                UPD_gains[i] = gval
                break
        if matched_key is not None:
            UPD_bkgErr[i] = Bkg_map[matched_key] * TimePerPoint[i]
        else:
            unknown_gains.add(float(gain))
    if unknown_gains:
        logging.warning(f"Unknown UPD amplifier gain values {sorted(unknown_gains)}; "
                        "gain/background left at 0 for those points.")

    result = dict()
    result["UPD_gains"] = UPD_gains
    result["UPD_bkgErr"] = UPD_bkgErr
    return result

def  calculatedQStep(Sample):
    # Calculate Q from the Sample data
    #TODO: add function adding dQ
    #	variable InstrumentQresolution = 2*pi*sin(BlankWidth/3600*pi/180)/Wavelength
    # dQ=InstrumentQresolution/2    # This is done by calculating Q from the ARangles and Wavelength
    Wavelength = Sample["reducedData"]["wavelength"]
    BlankWidth = Sample["BlankData"]["FWHM"]        # This is the width of the beam in degrees
    Qvec=Sample["CalibratedData"]["SMR_Qvec"]  # This is the Q vector from the calibrated data
    # make a copy of Qvec to avoid modifying the original
    BlankWidth_rad = np.radians(BlankWidth)
    InstrumentQresolution = (2*np.pi*np.sin(BlankWidth_rad)/Wavelength)/2   # Calculate dQ
    #dQ is copy of Qvec filled with InstrumentQresolution
    dQ = np.full_like(Qvec, InstrumentQresolution, dtype=float)  # Fill dQ with InstrumentQresolution
    # Now we can create the result dictionary
    result = dict()
    result["SMR_dQ"] = dQ
    return result

                

def calculatePDErrorStep(Sample, isBlank=False):
    """Uncertainty of the normalised step-scan signal, for either chain."""
    if Sample["RawData"].get("chain") == CHAIN_FX4:
        return _calculatePDErrorStepFX4(Sample, isBlank=isBlank)
    return _calculatePDErrorStepScaler(Sample, isBlank=isBlank)


def _calculatePDErrorStepFX4(Sample, isBlank=False):
    """FX4 step scan: propagate a relative current error through UPD / I0.

    Unlike the fly scan, the uascan file records no per-point sigma for the
    electrometer currents, so the counting statistics that the scaler chain
    relied on simply are not available.  Each channel is given a relative
    uncertainty of FX4_RELATIVE_CURRENT_ERROR, the detector channel is
    combined in quadrature with the measured dark-current error, and the two
    are propagated through the ratio.  This sets error bars only — the
    intensities are untouched.
    """
    raw = Sample["RawData"]
    UPD_array = np.asarray(raw["UPD_array"], dtype=float)
    Monitor = np.asarray(raw["Monitor"], dtype=float)
    block = Sample["BlankData"] if isBlank else Sample["reducedData"]
    updBkg = np.asarray(block.get("UPD_bkg", 0.0), dtype=float)
    updBkgErr = np.asarray(block.get("UPD_bkgErr", 0.0), dtype=float)

    sigma_upd = np.sqrt((FX4_RELATIVE_CURRENT_ERROR * UPD_array)**2 + updBkgErr**2)
    sigma_i0 = FX4_RELATIVE_CURRENT_ERROR * np.abs(Monitor)
    return {"Error": ratio_error(UPD_array - updBkg, sigma_upd, Monitor, sigma_i0)}


def _calculatePDErrorStepScaler(Sample, isBlank=False):
    # TODO : Igor code uses same for Setp and FLyscan...
    #OK, another incarnation of the error calculations...
    UPD_array = Sample["RawData"]["UPD_array"]
    # USAXS_PD = Sample["reducedData"]["Intensity"]
    MeasTime = Sample["RawData"]["TimeInSec"]    #measurement time in seconds per point
    if isBlank:
        UPD_gains=Sample["BlankData"]["UPD_gains"]
        UPD_bkgErr = Sample["BlankData"]["UPD_bkgErr"]    
    else:
        UPD_gains=Sample["reducedData"]["UPD_gains"]
        UPD_bkgErr = Sample["reducedData"]["UPD_bkgErr"]    

    Monitor = Sample["RawData"]["Monitor"]
    I0Gain=Sample["RawData"]["I0gain"]
    VToFFactor = Sample["RawData"]["VToFFactor"]                    #this is mca1 frequency, hardwired to 1e6 
    SigmaUSAXSPD=np.sqrt(UPD_array*(1+0.0001*UPD_array))		    #this is our USAXS_PD error estimate, Poisson error + 1% of value
    SigmaPDwDC=np.sqrt(SigmaUSAXSPD**2+(MeasTime*UPD_bkgErr)**2)    #This should include now measured error for background
    SigmaPDwDC=SigmaPDwDC/(VToFFactor*UPD_gains)
    A=(UPD_array)/(VToFFactor*UPD_gains)		                    #without dark current subtraction
    SigmaMonitor= np.sqrt(Monitor)		                            #these calculations were done for 10^6 
    ScaledMonitor = Monitor
    A = np.where(np.isnan(A), 0.0, A)
    SigmaMonitor = np.where(np.isnan(SigmaMonitor), 0.0, SigmaMonitor)
    SigmaPDwDC = np.where(np.isnan(SigmaPDwDC), 0.0, SigmaPDwDC)
    ScaledMonitor = np.where(np.isnan(ScaledMonitor), 0.0, ScaledMonitor)
    SigmaRwave=np.sqrt((A**2 * SigmaMonitor**4)+(SigmaPDwDC**2 * ScaledMonitor**4)+((A**2 + SigmaPDwDC**2) * ScaledMonitor**2 * SigmaMonitor**2))
    SigmaRwave=SigmaRwave/(ScaledMonitor*(ScaledMonitor**2-SigmaMonitor**2))
    SigmaRwave=SigmaRwave * I0Gain			#fix for use of I0 gain here, the numbers were too low due to scaling of PD by I0Gain
    Error=SigmaRwave / 5                    # this is the error in the USAXS data, it is not the same as in Igor, but it is close enough for now
    result = {"Error":Error}
    return result


## Stepscan main code here
def importStepScan(path, filename):
    """Read a USAXS step-scan (uascan) NeXus file into the RawData dictionary.

    Handles both counting chains (see fx4support).  The array names did NOT
    change across the conversion — ``/entry/data/UPD`` and ``/entry/data/I0``
    exist in both — only their units and meaning did:

    * ``scaler`` — counts accumulated over ``/entry/data/seconds`` ticks of
      the 1e7 Hz Joerger clock, to be divided by the Femto amplifier gain in
      ``upd_autorange_controls_gain``.
    * ``FX4``    — gain-independent mean picoamps.  There is no gain array
      and no ``seconds`` array; the per-point amplifier range needed for the
      dark-current lookup is ``/entry/data/fx4_autorange_lurange``.

    This is the format most likely to be mis-reduced, which is why the chain
    is taken from ``/entry/instrument/bluesky/metadata/counting_chain`` and
    never guessed from the field names.
    """
    with h5py.File(os.path.join(path, filename), 'r') as file:
        chain = detect_counting_chain(file)
        #read various data sets
        #AR angle
        dataset = file['/entry/data/a_stage_r']
        ARangles = np.ravel(np.array(dataset))
        n_points = len(ARangles)
        RangeIndex = None
        Bkg_err_map = None
        if chain == CHAIN_FX4:
            # No 'seconds' column: the FX4 mean is time-independent, so the
            # writer records none.  Reconstruct the requested count time from
            # the plan arguments — it is used for reporting only.
            TimePerPoint = _fx4StepCountTime(file, n_points)
            # Gain-independent picoamps: the gain the scaler chain divided out
            # is exactly 1 here, for both detector and monitor.
            I0gain = np.ones(n_points)
            AmpGain = np.ones(n_points)
            RangeIndex = np.ravel(np.array(file['/entry/data/fx4_autorange_lurange'])) \
                if '/entry/data/fx4_autorange_lurange' in file else None
            if RangeIndex is None:
                logging.warning(f"{filename}: FX4 step scan without "
                                "fx4_autorange_lurange; no dark-current subtraction.")
        else:
            #time per point
            dataset = file['/entry/data/seconds']
            TimePerPoint = np.ravel(np.array(dataset))
            # I0 gain
            dataset = file['/entry/data/I0_autorange_controls_gain']
            I0gain = np.ravel(np.array(dataset))
            #Arrays for gains during data collection
            dataset = file['/entry/data/upd_autorange_controls_gain']
            AmpGain = np.ravel(np.array(dataset))
        #I0 - Monitor (new name: I0, old name: I0_USAXS)
        dataset = file['/entry/data/I0'] if '/entry/data/I0' in file else file['/entry/data/I0_USAXS']
        Monitor = np.ravel(np.array(dataset))
        #UPD (new name: UPD, old name: PD_USAXS)
        dataset = file['/entry/data/UPD'] if '/entry/data/UPD' in file else file['/entry/data/PD_USAXS']
        UPD_array = np.ravel(np.array(dataset))
        #metadata
        keys_to_keep = ['SAD_mm', 'SDD_mm', 'thickness', 'title', 'useSBUSAXS',
                        'intervals', 'VToFFactor'
                    ]
        metadata_group = file['/entry/instrument/bluesky/metadata']
        metadata_dict = read_group_to_dict(metadata_group)     
        metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
        #add more values to metadata_dict
        data = file['/entry/instrument/bluesky/streams/baseline/terms_USAXS_transmission_I0_counts/value'] 
        USAXSPinT_I0Counts = data[0]
        data = file['/entry/instrument/bluesky/streams/baseline/terms_USAXS_transmission_I0_gain/value']
        USAXSPinT_I0Gain = data[0]
        data = file['/entry/instrument/bluesky/streams/baseline/terms_USAXS_transmission_diode_counts/value']
        USAXSPinT_pinCounts = data[0]
        data = file['/entry/instrument/bluesky/streams/baseline/terms_USAXS_transmission_diode_gain/value']
        USAXSPinT_pinGain = data[0]
        data = file['/entry/instrument/bluesky/streams/baseline/terms_USAXS_transmission_count_time/value']
        USAXSPinT_Time = data[0]
        # Fix for former pin<->I0 swap: the diode stream goes to trans_pin_*,
        # the I0 stream to trans_I0_*.  With the old swapped assignment
        # MeasuredTransmission in calibrateAndSubtractFlyscan computed the
        # RECIPROCAL of the intended value for step scans.
        # ⚗️ validate step-scan transmission/calibration against Igor.
        metadata_dict['trans_pin_counts'] = USAXSPinT_pinCounts
        metadata_dict['trans_pin_gain'] = USAXSPinT_pinGain
        metadata_dict['trans_I0_counts'] = USAXSPinT_I0Counts
        metadata_dict['trans_I0_gain'] = USAXSPinT_I0Gain
        metadata_dict['trans_I0_time'] = USAXSPinT_Time
        data = file['/entry/start_time']
        timeStamp = data[()]
        timeStamp = timeStamp.decode('utf-8')        
        #add some missing or incorrectly named parameters to match FLyscan
        data_sdd = file['/entry/instrument/bluesky/metadata/SDD_mm']
        SDDmm = data_sdd[()]
        data = file['/entry/instrument/bluesky/streams/baseline/terms_USAXS_diode_upd_size/value']
        UPDsize = data[0]
        data = file['/entry/instrument/bluesky/metadata/sample_thickness_mm']
        SampleThickness = data[()] 
        metadata_dict['detector_distance'] = SDDmm
        metadata_dict['UPDsize'] = UPDsize
        metadata_dict['timeStamp'] = timeStamp
        #Instrument
        instrument_group = file['/entry/instrument/monochromator']
        instrument_dict = read_group_to_dict(instrument_group)        
        #Sample
        sample_group = file['/entry/sample']
        sample_dict = read_group_to_dict(sample_group)
        sample_dict['thickness'] = SampleThickness

        # now backgrounds for UPD subtraction later
        baseline = "/entry/instrument/bluesky/streams/baseline/"
        if chain == CHAIN_FX4:
            # FX4 dark currents, in pA, keyed by the range index 0-4 that
            # fx4_autorange_lurange reports at each point.  The Femto
            # upd_autorange_controls_* records still exist in the baseline but
            # are stale on this chain — do not read them.
            Bkg_map = {}
            Bkg_err_map = {}
            for i in range(5):
                bkg = f"{baseline}fx4_autorange_ranges_range{i}_background/value"
                err = f"{baseline}fx4_autorange_ranges_range{i}_background_error/value"
                Bkg_map[i] = float(file[bkg][0]) if bkg in file else 0.0
                Bkg_err_map[i] = float(file[err][0]) if err in file else 0.0
        else:
            # these are the locations of the background values...
            # /entry/instrument/bluesky/streams/baseline/I0_autorange_controls_ranges_gain0_background/value, it is array of start adn end values.
            Bkg0 = file[f"{baseline}upd_autorange_controls_ranges_gain0_background/value"][0]
            Bkg1 = file[f"{baseline}upd_autorange_controls_ranges_gain1_background/value"][0]
            Bkg2 = file[f"{baseline}upd_autorange_controls_ranges_gain2_background/value"][0]
            Bkg3 = file[f"{baseline}upd_autorange_controls_ranges_gain3_background/value"][0]
            Bkg4 = file[f"{baseline}upd_autorange_controls_ranges_gain4_background/value"][0]
            # Create a dictionary to map AmpGain values to their corresponding background values
            Bkg_map = {
                "1e4": Bkg0,
                "1e6": Bkg1,
                "1e8": Bkg2,
                "1e10": Bkg3,
                "1e12": Bkg4
            }
    # Call the function with your arrays
    check_arrays_same_length(ARangles, TimePerPoint, Monitor, UPD_array)
    #Package these results into dictionary
    data_dict = {"filename": os.path.splitext(filename)[0],
                "chain": chain,
                "ARangles":ARangles,
                "TimePerPoint": TimePerPoint,
                "TimeInSec": (TimePerPoint if chain == CHAIN_FX4
                              else TimePerPoint / JOERGER_CLOCK_HZ),
                "Monitor":Monitor,
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "I0gain": I0gain,
                "RangeIndex": RangeIndex,
                # V/F converter factor; meaningless on the FX4 chain, where the
                # electrometer reports current directly.
                "VToFFactor": 1.0 if chain == CHAIN_FX4 else 1e6,
                "sample": sample_dict,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                "Bkg_map": Bkg_map,
                "Bkg_err_map": Bkg_err_map,
                }
    return data_dict


def _fx4StepCountTime(file, n_points):
    """Per-point count time for an FX4 uascan, in seconds.

    The FX4 reading is a mean current, so the writer records no ``seconds``
    column.  ``uascan`` computes the dwell from ``count_time`` and, when
    ``useDynamicTime`` is set, scales it by thirds across the scan (base/3
    over the first third, base over the second, 2*base over the last).  That
    is reproduced here.  The value is informational — nothing in the FX4
    reduction divides by it.
    """
    md_path = "/entry/instrument/bluesky/metadata/"
    count_time = None
    if md_path + "plan_args" in file:
        text = str(as_scalar(file[md_path + "plan_args"][()]))
        match = re.search(r"count_time:\s*([0-9.eE+-]+)", text)
        if match:
            try:
                count_time = float(match.group(1))
            except ValueError:
                count_time = None
    if count_time is None or count_time <= 0:
        logging.warning("FX4 step scan: no usable count_time in plan_args; "
                        "reporting 1 s per point.")
        return np.ones(n_points)

    dynamic = str(as_scalar(file[md_path + "useDynamicTime"][()])).strip().lower() == "true" \
        if md_path + "useDynamicTime" in file else False
    if not dynamic:
        return np.full(n_points, count_time)

    intervals = as_scalar(file[md_path + "intervals"][()]) if md_path + "intervals" in file else None
    intervals = float(intervals) if intervals else float(n_points)
    fraction = np.arange(n_points) / max(intervals, 1.0)
    times = np.full(n_points, count_time)
    times[fraction < 0.33] = count_time / 3
    times[fraction >= 0.66] = count_time * 2
    return times


def CorrectUPDGainsStep(data_dict):
    """Normalised step-scan detector signal, for either counting chain."""
    if data_dict["RawData"].get("chain") == CHAIN_FX4:
        return _CorrectUPDStepFX4(data_dict)
    return _CorrectUPDGainsStepScaler(data_dict)


def _CorrectUPDStepFX4(data_dict):
    """FX4 step scan: I = (UPD - dark) / I0, both in picoamps.

    No gain term and no dwell term — the FX4 reports a gain-independent mean
    current.  The dark current is looked up per point from
    ``fx4_autorange_lurange``; unlike the fly scan, the uascan file does
    record the range at every point.

    Also unlike the scaler chain, there is no one-point gain shift to undo:
    the Femto gain readback lagged its own range change, while the FX4 range
    is read in the same event document as the current.
    """
    raw = data_dict["RawData"]
    UPD_array = np.asarray(raw["UPD_array"], dtype=float)
    Monitor = np.asarray(raw["Monitor"], dtype=float)
    n_points = len(UPD_array)

    Bckg_corr = range_indexed_array(raw.get("RangeIndex"),
                                    raw.get("Bkg_map") or {}, n_points)
    return {"Intensity": (UPD_array - Bckg_corr) / Monitor,
            "UPD_bkg": Bckg_corr}


def _CorrectUPDGainsStepScaler(data_dict):
        # here we will multiply UPD by gain and divide by monitor corrected for its gain.
        # get the needed data from dictionary
    AmpGain = data_dict["RawData"]["AmpGain"]
    UPD_array = data_dict["RawData"]["UPD_array"]
    Monitor = data_dict["RawData"]["Monitor"]
    I0gain = data_dict["RawData"]["I0gain"]
    TimePerPoint= data_dict["RawData"]["TimePerPoint"]
    Bkg_map= data_dict["RawData"]["Bkg_map"]
        # for some  reason, the AmpGain is shifted by one value so we need to duplicate the first value and remove end value. 
    first_value = AmpGain[0]
    AmpGain = np.insert(AmpGain, 0, first_value)
    AmpGain = AmpGain[:-1]
                # change gain masking may not any be necessary... 
                # # need to remove points where gain changes
                # # Find indices where the change occurs
                # change_indices = np.where(np.diff(AmpGain) != 0)[0]
                # change_indices = change_indices +1
                # # fix range changes
                # #Correct UPD for gains so we can find max value loaction
                # UPD_temp = (UPD_array*I0gain)/(AmpGain*Monitor)
                # #remove renage chanegs on thsi array
                # UPD_temp[change_indices] = np.nan
                # # now locate location of max value in UPD_array
                # max_index = np.nanargmax(UPD_temp)
                # # we need to limit change_indices to values less than the location of maximum (before peak) = max_index
                # # this removes the range changes only to before the peak location, does nto seem to work, really
                # #change_indices = change_indices[change_indices < max_index]
                # # Create a copy of the array to avoid modifying the original
                # AmpGain_new = AmpGain.astype(float)                 # Ensure the array can hold NaN values
                # # Set the point before each range change to NaN
                # if len(change_indices) > 0:
                #     AmpGain_new[change_indices] = np.nan
                #Correct UPD for gains with points we  want removed set to Nan
    #this neds to be done on UPD_array before we correct for gains.
    # now we need to make a copy of UPD_array and depending on AmpGain value put int BkgX * TimePerPoint
    # Create a new array to store the results
    Bckg_corr = np.zeros_like(AmpGain, dtype=float)
    # Assign background values based on AmpGain and multiply by TimePerPoint
    # Convert keys to floats in Bkg_map for matching with AmpGain
    Bkg_map_float_keys = {float(k): v for k, v in Bkg_map.items()}

    unknown_gains = set()
    for i, gain in enumerate(AmpGain):
        background_value = Bkg_map_float_keys.get(gain)
        if background_value is None:
            # tolerant match — EPICS-sourced floats may not be bit-exact
            for gval, bval in Bkg_map_float_keys.items():
                if np.isclose(gain, gval, rtol=1e-3):
                    background_value = bval
                    break
        if background_value is None:
            background_value = 0    # unknown gain: no background subtraction
            unknown_gains.add(float(gain))
        Bckg_corr[i] = background_value * TimePerPoint[i]/JOERGER_CLOCK_HZ
    if unknown_gains:
        logging.warning(f"CorrectUPDGainsStep: unknown UPD amplifier gain values "
                        f"{sorted(unknown_gains)}; background set to 0 for those points.")
    # 1e7 confirmed correct for step scans (Joerger scaler internal clock);
    # the 1e6 used elsewhere is the flyscan MCA clock — different hardware.
    # Now we can correct UPD_array for background
    # Remove background from UPD_array
    UPD_array_corr = UPD_array - Bckg_corr       
    #Correct UPD for gains and monitor
    UPD_corrected = (UPD_array_corr*I0gain)/(AmpGain*Monitor)
    result = {"Intensity":UPD_corrected}
    return result


# (legacy commented-out reduceStepScanToQR / PlotResults removed 2026-07;
#  see git history if needed again)


if __name__ == "__main__":
    Sample = dict()
    Samples = []

    # "\\Mac\Home\Desktop\TestData\uascan"
    Sample = processStepscan("/Users/ilavsky/Desktop/Test440", "T07_04_440_4min_1025.h5", blankPath="/Users/ilavsky/Desktop/Test440", blankFilename="center_BLANK_1022.h5", recalculateAllData=True)
    Samples.append(Sample)
    plotUSAXSResults(Samples, "/Users/ilavsky/Desktop/Test440", isFlyscan=False)
    #Sample = reduceStepScanToQR("/home/parallels/Github/Matilda/TestData","USAXS_step.h5")
   # Sample = reduceStepScanToQR(r"C:\Users\ilavsky\Documents\GitHub\Matilda\TestData","USAXS_step.h5")
    #Sample["RawData"]=ImportStepScan("/home/parallels/Github/Matilda","USAXS_step.h5")
        #pp.pprint(Sample)
        #Sample["reducedData"]= CorrectUPDGainsStep(Sample)
        #Sample["reducedData"].update(BeamCenterCorrection(Sample))
        #pp.pprint(Sample["reducedData"])
    #PlotResults(Sample)
    #flyscan
    #Sample = dict()
    #Sample = reduceFlyscanToQR("./TestData","USAXS.h5",recalculateAllData=True)
    # Sample["RawData"]=ImportFlyscan("/home/parallels/Github/Matilda","USAXS.h5")
    # #pp.pprint(Sample)
    # Sample["reducedData"]= CorrectUPDGainsFly(Sample)
    # Sample["reducedData"].update(BeamCenterCorrection(Sample))
    # #pp.pprint(Sample["reducedData"])
    #PlotResults(Sample)

  
    
    
