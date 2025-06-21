'''
Here we develop new code which then moves to proper package
TODO:
    convertFlyscancalibrated.py
New, calibrated Flyscan code. 
    use: 
processFlyscan(samplePath,sampleName,blankPath=blankPath,blankFilename=blankFilename,deleteExisting=False)
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
                Sample["reducedData"]["PD_range"], 
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


from supportFunctions import subtract_data #read_group_to_dict, filter_nested_dict, check_arrays_same_length
from convertUSAXS import rebinData
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5, saveNXcanSAS, readMyNXcanSAS, find_matching_groups
from supportFunctions import importFlyscan, calculatePD_Fly, beamCenterCorrection, smooth_r_data
from desmearing import desmearData

recalculateAllData = True

MinQMinFindRatio = 1.05

# Thos code first reduces data to QR and if provided with Blank, it will do proper data calibration, subtraction, and even desmearing
# It will check if QR/NXcanSAS data exist and if not, it will create properly calibrated NXcanSAS in teh Nexus file
# If exist and recalculateAllData is False, it will reuse old ones. This is doen for plotting.
def processFlyscan(path, filename, blankPath=None, blankFilename=None, deleteExisting=recalculateAllData):
    # Open the HDF5 file in read/write mode
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r+') as hdf_file:
        # Check if the group 'location' exists, if yes, bail out as this is all needed. 
        required_attributes = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
        required_items = {'definition': 'NXcanSAS'}
        SASentries =  find_matching_groups(hdf_file, required_attributes, required_items)
        if deleteExisting:
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
        # More checks... if we have BlankName, full NXcanSAS need to exist or recalculate
        # if BlankName=None, then we just need the QRS_data group.   

        NXcanSASentry = next((entry + '/' for entry in SASentries if '_SMR' not in entry), None)
        location = None
        if blankFilename is not None and blankPath is not None:
            location = NXcanSASentry        # require we have desmeared data
        else:
            location = 'entry/QRS_data/'            # all we want here are QRS data
        
        if location is not None and location in hdf_file:
            # exists, so lets reuse the data from the file
            Sample = dict()
            Sample = readMyNXcanSAS(path, filename)
            logging.info(f"Using existing data for file {filename}. ")
            print(f"Using existing data for file {filename}. ")
            return Sample
        
        else:
            Sample = dict()
            Sample["RawData"]=importFlyscan(path, filename)                 #import data
            Sample["reducedData"]= calculatePD_Fly(Sample)                  # Creates PD_Intensity with corrected gains and background subtraction
            Sample["reducedData"].update(calculatePDError(Sample))          # Calculate UPD error, mostly the same as in Igor                
            Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=0)) #Beam center correction
            Sample["reducedData"].update(smooth_r_data(Sample["reducedData"]["Intensity"],     #smooth data data
                                                    Sample["reducedData"]["Q"], 
                                                    Sample["reducedData"]["PD_range"], 
                                                    Sample["reducedData"]["Error"], 
                                                    Sample["RawData"]["TimePerPoint"],
                                                    replaceNans=True))                 

            if blankPath is not None and blankFilename is not None and blankFilename != filename:               
                Sample["BlankData"]=getBlankFlyscan(blankPath, blankFilename)
                Sample["reducedData"].update(normalizeByTransmission(Sample))          # Normalize sample by dividing by transmission for subtraction
                Sample["CalibratedData"]=(calibrateAndSubtractFlyscan(Sample))
                #pp.pprint(Sample)
                # TODO: fix rebinning for 3 input waves returning 4 waves with dQ
                SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
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
                # save_dict_to_hdf5(Sample, location, hdf_file)
                # print("Appended new data to 'entry/displayData'.")
            
            else:
                #set calibrated data in the structure to None 
                Sample["CalibratedData"] = {"SMR_Qvec":None,
                                            "SMR_Int":None,
                                            "SMR_Error":None,
                                            "SMR_dQ":None,
                                            "Kfactor":None,
                                            "OmegaFactor":None,
                                            "BlankName":None,
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

def normalizeByTransmission(Sample):
    # This is a simple normalization of the Sample Intensity by transmission. 
    # It will be used for background subtraction.
    PeakIntensitySample = Sample["reducedData"]["Maximum"]
    PeakIntensityBlank = Sample["BlankData"]["Maximum"]
    PeakToPeakTransmission = PeakIntensitySample/PeakIntensityBlank
    #BlankIntensity = Sample["BlankData"]["Intensity"]
    Intensity = Sample["reducedData"]["Intensity"]
    Intensity = Intensity / PeakToPeakTransmission
    Error = Sample["reducedData"]["Error"]
    Error = Error / PeakToPeakTransmission
    result = {"Intensity":Intensity,
            "Error":Error,
            "PeakToPeakTransmission":PeakToPeakTransmission
            }
    return result
    
# TODO: remove deleteExisting=True for operations
def reduceFlyscanToQR(path, filename, deleteExisting=True):
    # Open the HDF5 file in read/write mode
    location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
            # Check if the group 'displayData' exists
            if deleteExisting:
                # Delete the group
                del hdf_file[location]
                #print("Deleted existing group 'entry/displayData'.")

            if location in hdf_file:
                # exists, so lets reuse the data from the file
                Sample = dict()
                Sample = load_dict_from_hdf5(hdf_file, location)
                #print("Used existing data")
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
                #print("Appended new data to 'entry/displayData'.")
                return Sample

def calibrateAndSubtractFlyscan(Sample):
    # This is a step where we subtract and calibrate the sample and Blank. 
    Intensity = Sample["reducedData"]["Intensity"]
    BL_Intensity = Sample["BlankData"]["Intensity"]
    Error = Sample["reducedData"]["Error"]
    BL_Error = Sample["BlankData"]["Error"]
    Q = Sample["reducedData"]["Q"]
    BL_Q = Sample["BlankData"]["Q"]
        
    SMR_Qvec, SMR_Int, SMR_Error, IntRatio = subtract_data(Q, Intensity,Error, BL_Q, BL_Intensity, BL_Error)
    # find Qmin as the first point where we get above 5% of the background avleu and larger than instrument resolution
    # IntRatio = Intensity / BL_Intensity, calculated using interpolation in subtract_data function
    # find point where the IntRatio is larger than 1.05 = MinQMinFindRatio, after Q dependent correction
    FWHMSample = Sample["reducedData"]["FWHM"]
    FWHMBlank = Sample["BlankData"]["FWHM"]
    wavelength =  Sample["reducedData"]["wavelength"]
    PeakToPeakTransmission =  Sample["reducedData"]["PeakToPeakTransmission"]
    SaTransCounts = Sample['RawData']['metadata']['trans_pin_counts']
    SaTransGain = Sample['RawData']['metadata']['trans_pin_gain']
    SaI0Counts = Sample['RawData']['metadata']['trans_I0_counts']
    SaI0Gain = Sample['RawData']['metadata']['trans_I0_gain'] 
    BlTransCounts = Sample['BlankData']['BlTransCounts']
    BlTransGain = Sample['BlankData']['BlTransGain']
    BlI0Counts = Sample['BlankData']['BlI0Counts']
    BlI0Gain = Sample['BlankData']['BlI0Gain']
    MeasuredTransmission = ((SaTransCounts/SaTransGain)/(SaI0Counts/SaI0Gain))/((BlTransCounts /BlTransGain )/(BlI0Counts/BlI0Gain))
    MSAXSCorrection = MeasuredTransmission / PeakToPeakTransmission
    QminSample = 4*np.pi*np.sin(np.radians(FWHMSample)/2)/wavelength
    QminBlank = 4*np.pi*np.sin(np.radians(FWHMBlank)/2)/wavelength
    indexSample = np.searchsorted(Q, QminSample)+1
    indexBlank = np.searchsorted(Q, QminBlank)+1
    # now we need to reproduce the Q correction from Igor. 
    MaxCorrection = 1
    PowerCorrection = 3
    QCorrection =  1 + MaxCorrection*(abs(QminBlank/SMR_Qvec))**PowerCorrection
    QCorrection = np.where(QCorrection < (MaxCorrection + 1), QCorrection, MaxCorrection + 1)
    IntRatio = IntRatio/ QCorrection  # apply the Q correction to the IntRatio
    #we need to set IntRatio values for points up to indexSample to 1:
    IntRatio[:indexSample] = 1
    #plot the IntRatio
    # plt.plot(SMR_Qvec, IntRatio, linestyle='-')  # You can customize the marker and linestyle
    # plt.title('Plot of IntRatio vs. Q')
    # plt.xlabel('log(Q) [1/A]')
    # plt.xlim((1e-5, 0.0005))
    # plt.show()    
    # # Find the first index where IntRatio > MinQMinFindRatio
    # np.argmax returns the first index of True. If all are False, it returns 0.
    potential_first_index = np.argmax(IntRatio > MinQMinFindRatio)
    if IntRatio[potential_first_index] > MinQMinFindRatio:
        indexRatio = potential_first_index
        #print(f"First index where IntRatio > {MinQMinFindRatio} is {indexRatio} with value {IntRatio[indexRatio]}.")
    else:
        # This means no element in IntRatio was > MinQMinFindRatio (argmax returned 0 and IntRatio[0] was not > 1.03)
        indexRatio = len(IntRatio) # Default to end of array if no such point is found
        logging.warning(f"No points found where IntRatio > {MinQMinFindRatio}. Defaulting indexRatio to end of array ({indexRatio}).")
        
    largest_value = max(indexSample, indexBlank, indexRatio)
    # Ensure the start_index is within bounds and there's data to slice
    if largest_value < len(SMR_Qvec):
        SMR_Qvec = SMR_Qvec[largest_value:]    
        SMR_Int = SMR_Int[largest_value:]    
        SMR_Error = SMR_Error[largest_value:]
    else:
        # This case means the calculated start index is at or beyond the end of the array
        logging.warning(f"Calculated start_index_for_data ({largest_value}) is at or beyond array length ({len(SMR_Qvec)}). "
                        "Resulting SMR arrays will be empty. Check Qmin, blank, and IntRatio > 1.03 criteria.")
        SMR_Qvec = np.array([])
        SMR_Int = np.array([])
        SMR_Error = np.array([])
    # now calibration... 
    SDD = Sample["RawData"]["metadata"]['detector_distance']
    UPDSize =  Sample["RawData"]["metadata"]['UPDsize']
    thickness = Sample["RawData"]["sample"]['thickness']
    BLPeakMax = Sample["BlankData"]["Maximum"]
    BlankName = Sample["BlankData"]["BlankName"]
    #Igor:	variable SlitLength=0.5*((4*pi)/wavelength)*sin(PhotoDiodeSize/(2*SDDistance))
    slitLength = 0.5*((4*np.pi)//wavelength)*np.sin(UPDSize/(2*SDD))
    OmegaFactor= (UPDSize/SDD)*np.radians(FWHMBlank)
    Kfactor=BLPeakMax*OmegaFactor*thickness * 0.1 
    #apply calibration
    SMR_Int =  SMR_Int / (Kfactor*MSAXSCorrection) 
    SMR_Error = SMR_Error/ (Kfactor*MSAXSCorrection) 
    SMR_Error = SMR_Error * PeakToPeakTransmission  #this is Igor correction from 2014 which fixes issues with high absowrption well scattering samples. 
    return {"SMR_Qvec":SMR_Qvec,
            "SMR_Int":SMR_Int,
            "SMR_Error":SMR_Error,
            "Kfactor":Kfactor,
            "OmegaFactor":OmegaFactor,
            "BlankName":BlankName,
            "thickness":thickness,
            "slitLength":slitLength,
            "units":"[cm2/cm3]"
            }

def getBlankFlyscan(blankPath, blankFilename, deleteExisting=False):
      # will reduce the blank linked as input into Igor BL_R_Int 
      # after reducing this first time, data are saved in Nexus file for subsequent use. 
      # We get the BL_QRS and calibration data as result.
    # Open the HDF5 file in read/write mode
    location = 'entry/blankData/'
    Filepath = os.path.join(blankPath, blankFilename)
    with h5py.File(Filepath, 'r+') as hdf_file:
            # Check if the group 'location' exists, if yes, either delete if asked for or use. 
            if deleteExisting:
                if location in hdf_file:
                    # Delete the group is exists and requested
                    del hdf_file[location]
                    #print("Deleted existing group 'entry/blankData'.")

            if location in hdf_file:
                # exists, so lets reuse the data from the file
                Blank = dict()
                Blank = load_dict_from_hdf5(hdf_file, location)
                #print("Used existing Blank data")
                return Blank
            else:
                Blank = dict()
                Blank["RawData"]=importFlyscan(blankPath, blankFilename)         #import data
                BlTransCounts = Blank['RawData']['metadata']['trans_pin_counts']
                BlTransGain = Blank['RawData']['metadata']['trans_pin_gain']
                BlI0Counts = Blank['RawData']['metadata']['trans_I0_counts']
                BlI0Gain = Blank['RawData']['metadata']['trans_I0_gain']
                Blank["BlankData"]= calculatePD_Fly(Blank)                  # Creates Intensity with corrected gains and background subtraction
                Blank["BlankData"].update({"BlankName":blankFilename})      # add the name of the blank file
                Blank["BlankData"].update({"BlTransCounts":BlTransCounts})  # add the BlTransCounts
                Blank["BlankData"].update({"BlTransGain":BlTransGain})      # add the BlTransGain
                Blank["BlankData"].update({"BlI0Counts":BlI0Counts})        # add the BlI0Counts
                Blank["BlankData"].update({"BlI0Gain":BlI0Gain})            # add the BlTransGain
                Blank["BlankData"].update(calculatePDError(Blank, isBlank=True))          # Calculate UPD error, mostly the same as in Igor                
                Blank["BlankData"].update(beamCenterCorrection(Blank,useGauss=0, isBlank=True)) #Beam center correction
                Blank["BlankData"].update(smooth_r_data(Blank["BlankData"]["Intensity"],     #smooth data data
                                                        Blank["BlankData"]["Q"], 
                                                        Blank["BlankData"]["PD_range"], 
                                                        Blank["BlankData"]["Error"], 
                                                        Blank["RawData"]["TimePerPoint"],
                                                        replaceNans=True )) 
                # we need to return just the BlankData part 
                BlankData=dict()
                BlankData=Blank["BlankData"]
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                save_dict_to_hdf5(BlankData, location, hdf_file)
                #print("Appended new Blank data to 'entry/blankData'.")
                return BlankData



def calculatePDError(Sample, isBlank=False):
    #OK, another incarnation of the error calculations...
    UPD_array = Sample["RawData"]["UPD_array"]
    # USAXS_PD = Sample["reducedData"]["Intensity"]
    MeasTimeCts = Sample["RawData"]["TimePerPoint"]
    Frequency=1e6   #this is frequency of clock fed into mca1
    MeasTime = MeasTimeCts/Frequency    #measurement time in seconds per point
    if isBlank:
        UPD_gains=Sample["BlankData"]["UPD_gains"]
        UPD_bkgErr = Sample["BlankData"]["UPD_bkgErr"]    
    else:
        UPD_gains=Sample["reducedData"]["UPD_gains"]
        UPD_bkgErr = Sample["reducedData"]["UPD_bkgErr"]    

    Monitor = Sample["RawData"]["Monitor"]
    I0AmpGain=Sample["RawData"]["metadata"]["I0AmpGain"]
    VToFFactor = Sample["RawData"]["VToFFactor"]/10      #this is mca1 frequency, HDF5 writer 1.3 and above needs /10 
    SigmaUSAXSPD=np.sqrt(UPD_array*(1+0.0001*UPD_array))		#this is our USAXS_PD error estimate, Poisson error + 1% of value
    SigmaPDwDC=np.sqrt(SigmaUSAXSPD**2+(MeasTime*UPD_bkgErr)**2) #This should include now measured error for background
    SigmaPDwDC=SigmaUSAXSPD/(Frequency*UPD_gains)
    A=(UPD_array)/(VToFFactor[0]*UPD_gains)		#without dark current subtraction
    SigmaMonitor= np.sqrt(Monitor)		            #these calculations were done for 10^6 
    ScaledMonitor = Monitor
    A = np.where(np.isnan(A), 0.0, A)
    SigmaMonitor = np.where(np.isnan(SigmaMonitor), 0.0, SigmaMonitor)
    SigmaPDwDC = np.where(np.isnan(SigmaPDwDC), 0.0, SigmaPDwDC)
    ScaledMonitor = np.where(np.isnan(ScaledMonitor), 0.0, ScaledMonitor)
    SigmaRwave=np.sqrt((A**2 * SigmaMonitor**4)+(SigmaPDwDC**2 * ScaledMonitor**4)+((A**2 + SigmaPDwDC**2) * ScaledMonitor**2 * SigmaMonitor**2))
    SigmaRwave=SigmaRwave/(ScaledMonitor*(ScaledMonitor**2-SigmaMonitor**2))
    SigmaRwave=SigmaRwave * I0AmpGain			#fix for use of I0 gain here, the numbers were too low due to scaling of PD by I0AmpGain
    Error=SigmaRwave		                    # this is the error in the USAXS data, it is not the same as in Igor, but it is close enough for now
    result = {"Error":Error}
    return result



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
    samplePath = r"\\Mac\Home\Desktop\Data\set1"
    sampleName="SA_R_0325.h5"
    blankPath=r"\\Mac\Home\Desktop\Data\set1" 
    blankFilename="TapeBlank_R_0317.h5"
    Sample = processFlyscan(samplePath,sampleName,blankPath=blankPath,blankFilename=blankFilename,deleteExisting=False)    
    
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

    Sample = {}
    Sample = readMyNXcanSAS(r"\\Mac\Home\Desktop\Data\set1", sampleName)
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