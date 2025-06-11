'''
Here we develop new code which then moves to proper package
TODO:
    convertFlyscancalibrated.py
New, calibrated Flyscan code. 
    use: 
processFlyscan(samplePath,sampleName,blankPath=blankPath,blankFilename=blankFilename,deleteExisting=True)
        For example of use see: test_matildaLocal() at the end of this file. 

    returns dictionary of this type:
            result["SampleName"]=sampleName
            result["BlankName"]=blankName
            result["reducedData"] =  {"Intensity":np.ravel(intensity), 
                              "Q":np.ravel(q),
                              "Error":np.ravel(error)}
            result["CalibratedData"] = {"Intensity":np.ravel(intcalib),
                                    "Q":np.ravel(qcalib),
                                    "Error":np.ravel(errcalib),
                                   }  

            this is when read from Nexus file. TODO: check and make sure it all works...  
                    {'BlankName': 'TapeBlank_R_0317.h5',
                'Error': array([5.41821571e+08, 3.61214381e+08, 1.80607190e+08, 1.45182650e+08,
                    8.21066769e+07, 5.28638233e+07, 9.31836113e+07, 3.44934716e+07,
                        ...
                    3.70688211e-02, 3.70688211e-02, 3.70688211e-02]),
                'Error_Attributes': <Attributes of closed HDF5 object>,
                'Int_attributes': <Attributes of closed HDF5 object>,
                'Intensity': array([4.52126415e+09, 3.89536685e+09, 3.08490075e+09, 2.58092468e+09,
                    2.16199554e+09, 1.83504586e+09, 1.49810441e+09, 1.31218746e+09,
                        ...
                    3.24382536e-01, 3.68735659e-01, 2.65436440e-01]),
                'Kfactor': np.float64(1.5194820557062783e-12),
                'OmegaFactor': np.float64(4.171381472542379e-08),
                'Q': array([1.05916862e-04, 1.09962784e-04, 1.14017257e-04, 1.18080210e-04,
                    1.22151570e-04, 1.26231272e-04, 1.30319251e-04, 1.34415446e-04,
                        ...
                    2.93847771e-01, 2.97917894e-01, 3.00031668e-01]),
                'Q_attributes': <Attributes of closed HDF5 object>,
                'dQ': array([4.04592188e-06, 4.04592188e-06, 4.05447384e-06, 4.06295227e-06,
                    4.07136056e-06, 4.07970183e-06, 4.08797899e-06, 4.09619476e-06,
                        ...
                    8.06617937e-03, 8.17837508e-03, 8.17837508e-03]),
                'dQ_Attributes': <Attributes of closed HDF5 object>,
                'label': 'AB3_R_0318',
                'thickness': np.float64(1.0),
                'units': '1/cm'}




   Does:
    Convert Flyscan USAXS data from the HDF5 format to the 1Ddata
    Decide if we need to do desmearing, not ready yet.   
    TODO: 
    Add background path and name as input, set to None as default, if set as input then do subtarcting and calibration.
        If not provided, stop with reduction before subtraction and return None for calibrated data. 
    Store both reduced data and NXcanSAS data in original hdf file, read from file if they exist and skip data reduction. 
    Only some metadata are kept to keep all more reasonable on size
'''
import h5py
import numpy as np
import pprint
#from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pprint as pp
from supportFunctions import subtract_data #read_group_to_dict, filter_nested_dict, check_arrays_same_length
import os
from convertUSAXS import rebinData
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5, saveNXcanSAS, readNXcanSAS
from convertUSAXS import importFlyscan, calculatePD_Fly, beamCenterCorrection, smooth_r_data
from desmearing import desmearData



#develop calibrated flyscan and step scan data routines
# this will check if NXcanSAS data exist and if not, it will create properly calibrated NXcanSAS
# data. It needs Blank and Sample data.
def processFlyscan(path, filename,blankPath=None, blankFilename=None, deleteExisting=False):
    # Open the HDF5 file in read/write mode
    #location = 'entry/displayData/'
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r+') as hdf_file:
        # Check if the group 'location' exists, if yes, bail out as this is all needed. 
        # if deleteExisting:
        #     # Delete the group
        #     del hdf_file[location]
        #     print("Deleted existing group 'entry/displayData'.")

        # if location in hdf_file:
        #     # exists, so lets reuse the data from the file
        #     Sample = dict()
        #     Sample = load_dict_from_hdf5(hdf_file, location)
        #     print("Used existing data")
        #     return Sample
        # else:
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

            if blankPath is not None and blankFilename is not None:               
                Sample["BlankData"]=getBlankFlyscan(blankPath, blankFilename)
                Sample["reducedData"].update(normalizeBlank(Sample))          # Normalize sample by dividing by transmission for subtraction
                Sample["CalibratedData"]=(calibrateAndSubtractFlyscan(Sample))
                #pp.pprint(Sample)
                # TODO: check calibration
                # TODO: fix rebinning for 3 input waves returning 4 waves with dQ
                Sample["CalibratedData"].update(rebinData(Sample, num_points=500, isSMRData=True))         #Rebin data
                # TODO: desmearing here
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
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
                #save_dict_to_hdf5(Sample, location, hdf_file)
                #print("Appended new data to 'entry/displayData'.")
            else:
                #set calibrated data in the structure to None 
                Sample["CalibratedData"] = {"SMR_Qvec":None,
                                            "SMR_Int":None,
                                            "SMR_Error":None,
                                            "Kfactor":None,
                                            "OmegaFactor":None,
                                            "BlankName":None,
                                            "thickness":None,
                                            "units":"[cm2/cm3]",
                                            "Intensity":None,
                                            "Q":None,
                                            "Error":None,
                                            "dQ":None,
                                            }
            return Sample

def normalizeBlank(Sample):
    # This is a simple normalization of the blank data to the sample data. 
    # It will be used for background subtraction.
    PeakIntensitySample = Sample["reducedData"]["Maximum"]
    PeakIntensityBlank = Sample["BlankData"]["Maximum"]
    PeakToPeakTransmission = PeakIntensitySample/PeakIntensityBlank


    Intensity = Sample["reducedData"]["Intensity"]
    Intensity = Intensity / PeakToPeakTransmission
    Error = Sample["reducedData"]["Error"]
    Error = Error / PeakToPeakTransmission
    result = {"Intensity":Intensity,
            "Error":Error,
            "PeakToPeakTransmission":PeakToPeakTransmission
            }
    return result
    

def calibrateAndSubtractFlyscan(Sample):
    # This is a step wehre we subtract and calibrate the sample and Blank. 
    Intensity = Sample["reducedData"]["Intensity"]
    BL_Intensity = Sample["BlankData"]["Intensity"]
    Error = Sample["reducedData"]["Error"]
    BL_Error = Sample["BlankData"]["Error"]
    Q = Sample["reducedData"]["Q"]
    BL_Q = Sample["BlankData"]["Q"]

    SMR_Qvec, SMR_Int, SMR_Error = subtract_data(Q, Intensity,Error, BL_Q, BL_Intensity, BL_Error)
    # TODO: trim, calibrate, 
    # find Qmin as the first point where we get above 3% of the background avleu and larger than instrument resolution
    IntRatio = Intensity / BL_Intensity
    # find point where the IntRatio is larger than 1.03
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
    indexSample = np.searchsorted(Q, QminSample)
    indexBlank = np.searchsorted(Q, QminBlank)
    indexRatio =  np.searchsorted(IntRatio, 1.03)
    largest_value = max(indexSample, indexBlank, indexRatio)
    SMR_Qvec = SMR_Qvec[largest_value-1 : ]    
    SMR_Int = SMR_Int[largest_value-1 : ]    
    SMR_Error = SMR_Error[largest_value-1 : ]
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
                    print("Deleted existing group 'entry/blankData'.")

            if location in hdf_file:
                # exists, so lets reuse the data from the file
                Blank = dict()
                Blank = load_dict_from_hdf5(hdf_file, location)
                print("Used existing Blank data")
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
                print("Appended new Blank data to 'entry/blankData'.")
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
    sampleName="AB3_R_0318.h5"
    blankPath=r"\\Mac\Home\Desktop\Data\set1" 
    blankFilename="TapeBlank_R_0317.h5"
    Sample = processFlyscan(samplePath,sampleName,blankPath=blankPath,blankFilename=blankFilename,deleteExisting=True)    
    
    # this is for testing save/restore from Nexus file... 
    testme=False 

    if (testme):
        # Specify the path and filename
        #file_path = 'C:/Users/ilavsky/Desktop/TestNexus.hdf'  # Replace with your actual file path
        file_path = r'\\Mac\Home\Desktop\Data\set1\TestNexus.hdf'  # Replace with your actual file path
        # Check if the file exists before attempting to delete it
        if os.path.exists(file_path):
            try:
                # Delete the file
                os.remove(file_path)
                print(f"File '{file_path}' has been deleted successfully.")
            except Exception as e:
                print(f"An error occurred while trying to delete the file: {e}")
        else:
            print(f"The file '{file_path}' does not exist.")
        #removed file
        saveNXcanSAS(Sample,r"\\Mac\Home\Desktop\Data\set1", "TestNexus.hdf")

        Data = readNXcanSAS(r"\\Mac\Home\Desktop\Data\set1", "TestNexus.hdf")
        pprint.pprint(Data)
        Sample = {}
        Sample['CalibratedData']=Data
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