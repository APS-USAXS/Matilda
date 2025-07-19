''' 
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
'''
import h5py
import numpy as np
import pprint
import os
import matplotlib.pyplot as plt
import pprint as pp
import logging
#from scipy.optimize import curve_fit


from supportFunctions import subtract_data 
from convertUSAXS import rebinData
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5, saveNXcanSAS, readMyNXcanSAS, find_matching_groups
from supportFunctions import importFlyscan, calculatePD_Fly, beamCenterCorrection, smooth_r_data
from supportFunctions import getBlankFlyscan, normalizeByTransmission,calibrateAndSubtractFlyscan,calculatePDErrorFly
from desmearing import desmearData


# Thos code first reduces data to QR and if provided with Blank, it will do proper data calibration, subtraction, and even desmearing
# It will check if QR/NXcanSAS data exist and if not, it will create properly calibrated NXcanSAS in teh Nexus file
# If exist and recalculateAllData is False, it will reuse old ones. This is doen for plotting.
def processFlyscan(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False):
    # Open the HDF5 file in read/write mode
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r+') as hdf_file:
        # Check if the group 'location' exists, if yes, bail out as this is all needed. 
        required_attributes = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
        required_items = {'definition': 'NXcanSAS'}
        SASentries =  find_matching_groups(hdf_file, required_attributes, required_items)
        if recalculateAllData:
            # Delete the groups which may have een created by previously run saveNXcanSAS
            location = 'entry/QRS_data/'
            if location is not None and location in hdf_file:
                # Delete the group
                del hdf_file[location]
                logging.info(f"Deleted existing group 'entry/QRS_data' for file {filename}. ")
            location = next((entry + '/' for entry in SASentries if '_SMR' in entry), None)
            if location is not None and location in hdf_file:
                # Delete the group
                del hdf_file[location]
                logging.info(f"Deleted existing group with SMR_data for file {filename}. ")
            location = next((entry + '/' for entry in SASentries if '_SMR' not in entry), None)
            if location is not None and location in hdf_file:
                # Delete the group
                del hdf_file[location]
                logging.info(f"Deleted existing NXcanSAS group for file {filename}. ")


        #Now, we will read the data from the file, if the exist. 
        # More checks... if we have blankname, full NXcanSAS need to exist or recalculate
        # if blankname=None, then we just need the QRS_data group.   

        NXcanSASentry = next((entry + '/' for entry in SASentries if '_SMR' not in entry), None)
        location = None
        if blankFilename is not None and blankPath is not None and "blank" not in filename.lower():
            location = NXcanSASentry        # require we have desmeared data
        else:
            location = 'entry/QRS_data/'            # all we want here are QRS data
        
        if location is not None and location in hdf_file:
            # exists, so lets reuse the data from the file
            Sample = dict()
            Sample = readMyNXcanSAS(path, filename, isUSAXS = True)
            logging.info(f"Using existing processed data from file {filename}.")
            return Sample
        
        else:
            Sample = dict()
            Sample["RawData"]=importFlyscan(path, filename)                         # import data
            Sample["reducedData"]= calculatePD_Fly(Sample)                          # Creates PD_Intensity with corrected gains and background subtraction
            Sample["reducedData"].update(calculatePDErrorFly(Sample))               # Calculate UPD error, mostly the same as in Igor                
            Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0))   # Beam center correction
            Sample["reducedData"].update(smooth_r_data(Sample["reducedData"]["Intensity"],     #smooth data data
                                                    Sample["reducedData"]["Q"], 
                                                    Sample["reducedData"]["UPD_gains"], 
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
                Sample["CalibratedData"]=(calibrateAndSubtractFlyscan(Sample))
                SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                if len(SMR_Qvec) > 50:  # some data were found. Call this success? 
                    if len(SMR_Qvec) > 800:  # if we have enough data, then rebin and desmear
                        Sample["CalibratedData"].update(rebinData(Sample, num_points=500, isSMRData=True))         #Rebin data
                    slitLength=Sample["CalibratedData"]["slitLength"]
                    #DesmearNumberOfIterations = 10
                    SMR_Int =Sample["CalibratedData"]["SMR_Int"]
                    SMR_Error =Sample["CalibratedData"]["SMR_Error"]
                    SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
                    SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]
                    DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ = desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=slitLength,ExtrapMethod='PowerLaw w flat',ExtrapQstart=0.1, MaxNumIter = 20)
                    desmearedData=list()
                    desmearedData={
                        "Intensity":DSM_Int,
                        "Q":DSM_Qvec,
                        "Error":DSM_Error,
                        "dQ":DSM_dQ,
                        "units":"[cm2/cm3]",
                        }
                    Sample["CalibratedData"].update(desmearedData)
                else:
                    logging.warning(f"Not enough data points in SMR_Qvec ({len(SMR_Qvec)}) to proceed with desmearing or rebinning. "
                                    "Skipping desmearing and rebinning steps. ")
                    #set calibrated data in the structure to None 
                    Sample["CalibratedData"] = {"SMR_Qvec":None,
                                                "SMR_Int":None,
                                                "SMR_Error":None,
                                                "SMR_dQ":None,
                                                "Kfactor":None,
                                                "OmegaFactor":None,
                                                "blankname":None,
                                                "thickness":None,
                                                "units":"[cm2/cm3]",
                                                "Intensity":None,
                                                "Q":None,
                                                "Error":None,
                                                "dQ":None,
                                                "slitLength":None,
                                                }                    
            
            else:
                #set calibrated data in the structure to None 
                Sample["CalibratedData"] = {"SMR_Qvec":None,
                                            "SMR_Int":None,
                                            "SMR_Error":None,
                                            "SMR_dQ":None,
                                            "Kfactor":None,
                                            "OmegaFactor":None,
                                            "blankname":None,
                                            "thickness":None,
                                            "units":"[cm2/cm3]",
                                            "Intensity":None,
                                            "Q":None,
                                            "Error":None,
                                            "dQ":None,
                                            "slitLength":None,
                                            }
        # Ensure all changes are written and close the HDF5 file
        hdf_file.flush()
    # The 'with' statement will automatically close the file when the block ends
    saveNXcanSAS(Sample,path, filename)
    return Sample


def reduceFlyscanToQR(path, filename, recalculateAllData=False):
    # Open the HDF5 file in read/write mode
    location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
            # Check if the group 'displayData' exists
            if recalculateAllData:
                # Delete the group
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

    Sample = dict()
    #does the file exists?
    # e = os.path.isfile("C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/USAXS.h5")
    # if not e:
    #     print("File not found")
    #     return
    # else:
    #     print("File found")
    #open the file
    #samplePath = "C:/Users/ilavsky/Documents/GitHub/Matilda/TestData/TestSet/02_21_Megan_usaxs"
    samplePath = "/home/parallels/Desktop/Test"
    samplename="MSI07.h5"
    blankPath="/home/parallels/Desktop/Test" 
    blankFilename="AirBlank_1967.h5"
    Sample = processFlyscan(samplePath,samplename,blankPath=blankPath,blankFilename=blankFilename,recalculateAllData=False)    
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
    DSM_Qvec =Sample["CalibratedData"]["Q"] 
    DSM_Int =Sample["CalibratedData"]["Intensity"] 
    #DSM_Error =Sample["CalibratedData"]["Error"] 
    plt.figure(figsize=(6, 12))
    #plt.plot(SMR_Qvec, SMR_Int, linestyle='-')  # You can customize the marker and linestyle
    plt.plot(DSM_Qvec, DSM_Int, linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of Intensity vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('Intensity')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.show() 




if __name__ == "__main__":
    #test_matilda()
    test_matildaLocal()