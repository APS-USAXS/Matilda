"""
convertFlyscan.py
=================
Reduce USAXS flyscan HDF5 files to calibrated 1-D I(Q) data.

Main entry point
----------------
processFlyscan(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False)

processFlyscan(samplePath,samplename,blankPath=blankPath,blankFilename=blankFilename,recalculateAllData=False)
        For example of use see: test_matildaLocal() at the end of this file. 
        Does:
        Convert Flyscan USAXS data from the HDF5 format to the 1Ddata
        If Background is None, return only reduced data, no calibration or subtraction.
        If Background is provided, then do calibration and subtraction.
        Store both reduced data and NXcanSAS data in original hdf file, read from file if they exist and skip data reduction. 
        Only some metadata are kept to keep all more reasonable on size

        returns dictionary of this type:
                Sample["reducedData"]["Intensity"],   
                Sample["reducedData"]["Q"], 
                Sample["reducedData"]["UPD_gains"], 
                Sample["reducedData"]["Error"], 
               
                SMR_Int =Sample["CalibratedData"]["SMR_Int"]
                SMR_Error =Sample["CalibratedData"]["SMR_Error"]
                SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]

               Sample["CalibratedData"]={
                     "Intensity":DSM_Int,
                     "Q":DSM_Qvec,
                     "Error":DSM_Error,
                     "dQ":DSM_dQ,
                     "units":"[cm2/cm3]",

Data flow
---------
importFlyscan()                 — read raw arrays from HDF5 NXsas file
calculatePD_Fly()               — compute normalised detector signal
beamCenterCorrection()          — apply beam-centre angle offset
smooth_r_data()                 — optional smoothing
getBlankFlyscan()               — load and process blank scan
normalizeByTransmission()       — apply transmission correction
calibrateAndSubtractFlyscan()   — subtract blank, apply K-factor / Omega
desmear_dispatch()              — slit-smearing correction (Lake/Strobl or GP method)
saveNXcanSAS() / readMyNXcanSAS() — cache results in the original HDF5 file

Notes
-----
* matplotlib is imported but currently only used for optional/debug plots
  (all plt.show() calls are commented out).  Target for removal when GUI
  work begins.
* pprint is imported twice (as pprint and as pp); one import is redundant.
"""
import h5py
import os
import logging
#from scipy.optimize import curve_fit


# rebinData lives in supportFunctions (was imported via convertUSAXS re-export)
from .supportFunctions import rebinData
from .hdf5code import load_dict_from_hdf5, save_dict_to_hdf5, saveNXcanSAS, readMyNXcanSAS
from .hdf5code import clearAndCheckCachedReduction, writeThicknessOverride
from .supportFunctions import importFlyscan, calculatePD_Fly, beamCenterCorrection, smooth_r_data
from .supportFunctions import getBlankFlyscan, normalizeByTransmission,calibrateAndSubtractFlyscan,calculatePDErrorFly
from .supportFunctions import empty_calibrated_data
from .desmearing_methods import desmear_dispatch


# This code first reduces data to QR and if provided with Blank, it will do proper data calibration, subtraction, and even desmearing
# It will check if QR/NXcanSAS data exist and if not, it will create properly calibrated NXcanSAS in the Nexus file
# If exist and recalculateAllData is False, it will reuse old ones. This is done for plotting.
def processFlyscan(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False,
                    num_points=500, desmear_iter=20, extrap_method='PowerLaw w flat',
                    extrap_qstart=0.15, minQMinFindRatio=1.05, thickness_override=None,
                    use_mu=False, mu=None, per_gram=False, density=None,
                    transmission_override=None, qmin_override=None,
                    desmear_method='lake', gp_length_scale=0.5, gp_kernel='matern32'):
    """Reduce a single USAXS flyscan HDF5 file to calibrated 1-D I(Q).

    Results are cached inside the original HDF5 file as NXcanSAS groups so
    that subsequent calls with recalculateAllData=False are fast (data are
    read from file rather than recomputed).

    Parameters
    ----------
    path : str
        Directory containing the sample HDF5 file.
    filename : str
        Filename of the sample HDF5 file (.h5).
    blankPath : str or None, optional
        Directory containing the blank HDF5 file.  If None, only raw QR
        data are produced (no calibration or blank subtraction).
    blankFilename : str or None, optional
        Filename of the blank HDF5 file.  If None, only raw QR data.
    recalculateAllData : bool, optional
        When True, delete cached NXcanSAS groups and recompute everything.
        Default False.
    num_points : int, optional
        Number of output points after rebinning.  Default 500.
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
        Sample dictionary with keys:
        * RawData       — raw detector arrays and metadata
        * reducedData   — Q, Intensity, Error, UPD_gains (normalised, slit-smeared)
        * CalibratedData — Q, Intensity, Error, dQ, units (desmeared, blank-subtracted)
                          Present only when a blank is provided.
        * SMR data stored under CalibratedData['SMR_*'] keys.
    """
    # Open the HDF5 file in read/write mode
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r+') as hdf_file:
        # Cache bookkeeping (shared with processStepscan): with a blank we
        # require the full desmeared NXcanSAS entry, otherwise QRS_data is enough.
        requireCalibrated = (blankFilename is not None and blankPath is not None
                             and "blank" not in filename.lower())
        if clearAndCheckCachedReduction(hdf_file, filename, recalculateAllData, requireCalibrated):
            # exists, so lets reuse the data from the file
            Sample = readMyNXcanSAS(path, filename, isUSAXS = True)
            logging.info(f"Using existing processed data from file {filename}.")
            return Sample

        else:
            Sample = dict()
            if thickness_override is not None:
                writeThicknessOverride(hdf_file, '/entry/sample/thickness',
                                       thickness_override, filename)
            Sample["RawData"]=importFlyscan(path, filename)                         # import data
            Sample["reducedData"]= calculatePD_Fly(Sample)                          # Creates PD_Intensity with corrected gains and background subtraction
            Sample["reducedData"].update(calculatePDErrorFly(Sample))               # Calculate UPD error, mostly the same as in Igor                
            Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0))   # Beam center correction
            Sample["reducedData"].update(smooth_r_data(Sample["reducedData"]["Intensity"],     #smooth data data
                                                    Sample["reducedData"]["Q"],
                                                    Sample["reducedData"]["UPD_gainsIndx"],    # range INDEX (0-4), not gain values
                                                    Sample["reducedData"]["Error"],
                                                    Sample["RawData"]["TimePerPoint"],
                                                    replaceNans=True))

            if (
                blankPath is not None
                and blankFilename is not None
                and blankFilename != filename
                and "blank" not in filename.lower()
            ):
                Sample["BlankData"]=getBlankFlyscan(blankPath, blankFilename,recalculateAllData=recalculateAllData)
                Sample["reducedData"].update(normalizeByTransmission(Sample))          # Normalize sample by dividing by transmission for subtraction
                Sample["CalibratedData"]=(calibrateAndSubtractFlyscan(Sample, minQMinFindRatio=minQMinFindRatio, thickness_override=thickness_override, use_mu=use_mu, mu=mu, per_gram=per_gram, density=density, transmission_override=transmission_override, qmin_override=qmin_override))
                SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                if len(SMR_Qvec) > 50:  # some data were found. Call this success? 
                    if len(SMR_Qvec) > 800:  # if we have enough data, then rebin and desmear
                        Sample["CalibratedData"].update(rebinData(Sample, num_points=num_points, isSMRData=True))         #Rebin data
                    slitLength=Sample["CalibratedData"]["slitLength"]
                    #DesmearNumberOfIterations = 10
                    SMR_Int =Sample["CalibratedData"]["SMR_Int"]
                    SMR_Error =Sample["CalibratedData"]["SMR_Error"]
                    SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                    SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]
                    DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ = desmear_dispatch(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=slitLength, method=desmear_method, length_scale_decades=gp_length_scale, kernel=gp_kernel, extrap_method=extrap_method, extrap_qstart=extrap_qstart, max_iter=desmear_iter)
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


def reduceFlyscanToQR(path, filename, recalculateAllData=False):
    # Open the HDF5 file in read/write mode
    location = 'entry/displayData/'
    with h5py.File(os.path.join(path, filename), 'r+') as hdf_file:
            # Check if the group 'displayData' exists
            if recalculateAllData:
                if location in hdf_file:
                    # Delete the group only if it exists (first run has none)
                    del hdf_file[location]
                    logging.info("Deleted existing group 'entry/displayData'.")

            if location in hdf_file:
                # exists, so lets reuse the data from the file
                Sample = dict()
                Sample = load_dict_from_hdf5(hdf_file, location)
                logging.info(f"Used existing QR data from {filename}")
                return Sample
            else:
                Sample = dict()
                Sample["RawData"]=importFlyscan(path, filename)         #import data
                Sample["reducedData"]= calculatePD_Fly(Sample)       # Correct gains
                Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0)) #Beam center correction
                Sample["reducedData"].update(rebinData(Sample))         #Rebin data
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                # these are not fully reduced data, this is for web plot purpose. 
                save_dict_to_hdf5(Sample, location, hdf_file)
                logging.info(f"Appended new QR data to 'entry/displayData' in {filename}.")
                return Sample


def test_matildaLocal():

    #does the file exists?
    # e = os.path.isfile("C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/USAXS.h5")
    # if not e:
    #     print("File not found")
    #     return
    # else:
    #     print("File found")
    #open the file
    #samplePath = "C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/TestSet/02_21_Megan_usaxs"
    #samplePath = "Z:/Experiments/USAXS_data/2025/2025-07/07_20_Tifani/07_20_Tifani_usaxs"
    samplePath = "/home/parallels/Desktop/testdata/hematite/hematite_usaxs"
    samplename="hematite_48C_98min_0050.h5"
    blankPath=samplePath 
    blankFilename="CapillaryBlank_0006.h5"
    # assign to `Sample` if you re-enable the debug plotting below
    processFlyscan(samplePath,samplename,blankPath=blankPath,blankFilename=blankFilename,recalculateAllData=True)
    #Sample = processFlyscan(samplePath,blankFilename,blankPath=blankPath,blankFilename=blankFilename,recalculateAllData=False)    
    
    # # this is for testing save/restore from Nexus file... 
    # testme=False 

    # if (testme):
    #     # Specify the path and filename
    #     #file_path = 'C:/Users/ilavsky/Desktop/TestNexus.hdf'  # Replace with your actual file path
    #     file_path = r'\\Mac\Home\Desktop\Data\set1\TestNexus.hdf'  # Replace with your actual file path
    #     # Check if the file exists before attempting to delete it
    #     if os.path.exists(file_path):
    #         try:
    #             # Delete the file
    #             os.remove(file_path)
    #             print(f"File '{file_path}' has been deleted successfully.")
    #         except Exception as e:
    #             print(f"An error occurred while trying to delete the file: {e}")
    #     else:
    #         print(f"The file '{file_path}' does not exist.")
    #     #removed file
    #     saveNXcanSAS(Sample,r"\\Mac\Home\Desktop\Data\set1", "TestNexus.hdf")

    #Sample = {}
    #Sample = readMyNXcanSAS(r"\\Mac\Home\Desktop\Data\set1", samplename)
    #pprint.pprint(Data)
    #Sample['CalibratedData']=Data
    # Q = Sample["reducedData"]["Q"]
    # UPD = Sample["reducedData"]["Intensity"]
    # Error = Sample["reducedData"]["Error"]
    # plt.figure(figsize=(6, 12))
    # plt.plot(Q, UPD, linestyle='-')  # You can customize the marker and linestyle
    # #plt.plot(Q, Intensity, linestyle='-')  # You can customize the marker and linestyle
    # plt.title('Plot of Intensity vs. Q')
    # plt.xlabel('log(Q) [1/A]')
    # plt.ylabel('Intensity')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.grid(True)
    # plt.show() 
    # SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"] 
    # SMR_Int =Sample["CalibratedData"]["SMR_Int"] 
    # #SMR_Error =Sample["CalibratedData"]["SMR_Error"] 
    # DSM_Qvec =Sample["CalibratedData"]["Q"]        # used by debug plot below
    # DSM_Int =Sample["CalibratedData"]["Intensity"]
    #DSM_Error =Sample["CalibratedData"]["Error"]
    # Debug plot (requires matplotlib; commented out for production):
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(6, 12))
    # plt.plot(DSM_Qvec, DSM_Int, linestyle='-')
    # plt.title('Plot of Intensity vs. Q')
    # plt.xlabel('log(Q) [1/A]')
    # plt.ylabel('Intensity')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.grid(True)
    # plt.show()




if __name__ == "__main__":
    #test_matilda()
    test_matildaLocal()