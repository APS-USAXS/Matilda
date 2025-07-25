# these are support functions for Matilda

import h5py
import numpy as np
import logging
import os
import copy
import re
import pprint as pp

from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5, saveNXcanSAS, readMyNXcanSAS, find_matching_groups

#recalculateAllData = False

MinQMinFindRatio = 1.05

## support stuff here


## importFlyscan loads data from flyscan NX file. It should be same for QR pass as well as for calibrated data processing. 
def importFlyscan(path, filename):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    with h5py.File(path+"/"+filename, 'r') as file:
        #read various data sets
        #figure out how many points are in AR angles, this has 1 more point that mca data, usually
        dataset = file['/entry/flyScan/AR_PulsePositions'] 
        ARangles =  np.ravel(np.array(dataset))
        num_elements = ARangles.size - 1 
        ARangles= ARangles[-num_elements:]
        #time per point
        dataset = file['/entry/flyScan/mca1'] 
        TimePerPoint = np.ravel(np.array(dataset))  [-num_elements:]
        #I0 - Monitor
        dataset = file['/entry/flyScan/mca2'] 
        Monitor = np.ravel(np.array(dataset))   [-num_elements:]
        #UPD
        dataset = file['/entry/flyScan/mca3'] 
        UPD_array = np.ravel(np.array(dataset)) [-num_elements:]
        #Arrays for UPD gain changes
        dataset = file['/entry/flyScan/changes_DDPCA300_ampGain'] 
        AmpGain = np.ravel(np.array(dataset))
        dataset = file['/entry/flyScan/changes_DDPCA300_ampReqGain'] 
        AmpReqGain = np.ravel(np.array(dataset))
        dataset = file['/entry/flyScan/changes_DDPCA300_mcsChan'] 
        Channel = np.ravel(np.array(dataset))            
        dataset = file['/entry/flyScan/mca_clock_frequency'] 
        vTof = np.full_like(np.ravel(np.array(dataset)), 1e6, dtype=float)
        #vTof = 1e6  overwrite, the mca_clock_frequency (5e7) value is simply wrong.     
        #metadata
        keys_to_keep = ['AR_center', 'ARenc_0', 'DCM_energy', 'DCM_theta', 'I0Gain','detector_distance',
                        'timeStamp','I0AmpGain',
                        'trans_pin_counts','trans_pin_gain','trans_pin_time','trans_I0_counts','trans_I0_gain',
                        'UPDsize', 'trans_I0_counts', 'trans_I0_gain', 'upd_bkg0', 'upd_bkg1','upd_bkg2','upd_bkg3',
                        'upd_bkgErr0','upd_bkgErr1','upd_bkgErr2','upd_bkgErr3','upd_bkgErr4','upd_bkg_err0',
                        'upd_bkg4','DDPCA300_gain0','DDPCA300_gain1','DDPCA300_gain2','DDPCA300_gain3','DDPCA300_gain4',
                        'upd_amp_change_mask_time0','upd_amp_change_mask_time1','upd_amp_change_mask_time2','upd_amp_change_mask_time3','upd_amp_change_mask_time4',
                    ]
        metadata_group = file['/entry/metadata']
        metadata_dict = read_group_to_dict(metadata_group)
        metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
        # we need this key to be there also... COpy of the other one. 
        I0Gain=metadata_dict["I0AmpGain"]   
        metadata_dict["I0Gain"]=I0Gain
        #Instrument
        keys_to_keep = ['monochromator', 'energy', 'wavelength']
        instrument_group = file['/entry/instrument']
        instrument_dict = read_group_to_dict(instrument_group)
        instrument_dict = filter_nested_dict(instrument_dict, keys_to_keep)
        # sample
        sample_group = file['/entry/sample']
        sample_dict = read_group_to_dict(sample_group)

    # Call the function with your arrays
    check_arrays_same_length(ARangles, TimePerPoint, Monitor, UPD_array)
    #Package these results into dictionary
    data_dict = {"filename": os.path.splitext(filename)[0],
                "ARangles":ARangles, 
                "TimePerPoint": TimePerPoint, 
                "Monitor":Monitor, 
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "Channel": Channel,
                "VToFFactor": vTof,
                "AmpReqGain": AmpReqGain,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                "sample": sample_dict,
                }
    
    return data_dict

# this finds the best blank scan for any scan.
# rules: same folder and order number lower than the sample scan., the closest one. 
def findProperBlankScan(scan_path, scan_filename, ListOfBlanks):
    logging.info(f"Looking for proper Blanks for to process scan: {scan_filename} in path: {scan_path}")

    sample_base, sample_num, sample_ext = _parse_filename_info(scan_filename)

    selected_blank_path = None
    selected_blank_filename = None

    if sample_base is None or sample_num is None:
        logging.warning(f"Could not parse sample filename: {scan_filename} to find its number. Proceeding without blank.")
    else:
        logging.debug(f"Parsed sample {scan_filename}: base='{sample_base}', num={sample_num}, ext='{sample_ext}'")
        candidate_blanks = []
        for bl_path, bl_filename in ListOfBlanks:
            if bl_path == scan_path:  # Blank must be in the same path
                blank_base, blank_num, blank_ext = _parse_filename_info(bl_filename)
                if blank_base and blank_num is not None and blank_ext == sample_ext:
                    if blank_num < sample_num:
                        candidate_blanks.append({
                            'path': bl_path,
                            'filename': bl_filename,
                            'number': blank_num
                        })
        
        if candidate_blanks:
            # Sort candidates by their number in descending order to get the largest number < sample_num
            candidate_blanks.sort(key=lambda b: b['number'], reverse=True)
            best_blank = candidate_blanks[0]
            selected_blank_path = best_blank['path']
            selected_blank_filename = best_blank['filename']
            logging.info(f"Found blank for {scan_filename}: {selected_blank_filename} (num: {best_blank['number']}) in path {selected_blank_path}")
        else:
            logging.warning(f"No suitable blank found for sample {scan_filename} (num: {sample_num}) in path '{scan_path}'.")

        # try:
        #     #    processFlyscan(samplePath,samplename,blankPath=blankPath,blankFilename=blankFilename,recalculateAllData=True)
        #     results.append(processFlyscan(path, filename, blankPath=selected_blank_path,blankFilename=selected_blank_filename,recalculateAllData=True))
        # except:
        #     pass
        #     #print("Done processing the Flyscans")
        # Call the imported processFlyscan (singular) which does the actual data reduction
        logging.info(f"Found blank for {scan_filename} , blank name : {selected_blank_filename}")
    return selected_blank_path, selected_blank_filename   

def getBlankFlyscan(blankPath, blankFilename, recalculateAllData=False):
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
                Blank["RawData"]=importFlyscan(blankPath, blankFilename)         #import data
                BlTransCounts = Blank['RawData']['metadata']['trans_pin_counts']
                BlTransGain = Blank['RawData']['metadata']['trans_pin_gain']
                BlI0Counts = Blank['RawData']['metadata']['trans_I0_counts']
                BlI0Gain = Blank['RawData']['metadata']['trans_I0_gain']
                Blank["BlankData"]= calculatePD_Fly(Blank)                  # Creates Intensity with corrected gains and background subtraction
                Blank["BlankData"].update({"blankname":blankFilename})      # add the name of the blank file
                Blank["BlankData"].update({"BlTransCounts":BlTransCounts})  # add the BlTransCounts
                Blank["BlankData"].update({"BlTransGain":BlTransGain})      # add the BlTransGain
                Blank["BlankData"].update({"BlI0Counts":BlI0Counts})        # add the BlI0Counts
                Blank["BlankData"].update({"BlI0Gain":BlI0Gain})            # add the BlTransGain
                Blank["BlankData"].update(calculatePDErrorFly(Blank, isBlank=True))          # Calculate UPD error, mostly the same as in Igor                
                Blank["BlankData"].update(beamCenterCorrection(Blank,useGauss=0, isBlank=True)) #Beam center correction
                Blank["BlankData"].update(smooth_r_data(Blank["BlankData"]["Intensity"],     #smooth data data
                                                        Blank["BlankData"]["Q"], 
                                                        Blank["BlankData"]["UPD_gains"], 
                                                        Blank["BlankData"]["Error"], 
                                                        Blank["RawData"]["TimePerPoint"],
                                                        replaceNans=True )) 
                # we need to return just the BlankData part 
                BlankData=dict()
                BlankData=Blank["BlankData"]
                # Create the group and dataset for the new data inside the hdf5 file for future use. 
                save_dict_to_hdf5(BlankData, location, hdf_file)
                logging.info(f"Appended new Blank data to 'entry/blankData' in {blankFilename}.")
                return BlankData

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
    
def calibrateAndSubtractFlyscan(Sample):
    # This is a step where we subtract and calibrate the sample and Blank. 
    Intensity = Sample["reducedData"]["Intensity"]
    BL_Intensity = Sample["BlankData"]["Intensity"]
    Error = Sample["reducedData"]["Error"]
    BL_Error = Sample["BlankData"]["Error"]
    Q = Sample["reducedData"]["Q"]
    BL_Q = Sample["BlankData"]["Q"]

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

    #TODO: remove nans from both sets of arrays here. 
    
    #Intensity and Error are corrected for transmission in normalize by transmission above. 
    SMR_Qvec, SMR_Int, SMR_Error, IntRatio = subtract_data(Q, Intensity,Error, BL_Q, BL_Intensity, BL_Error)
    # we need to fix negative intensities as Igor does in IN3_FixNegativeIntensities
    MaxSMR_Int = np.max(SMR_Int)
    #find min value in SMR_Int for points from half to end
    MinSMR_Int = np.min(SMR_Int[int(len(SMR_Int)/2):])
    #if MinSMR_Int < 0, then we need add to SMR_Int enough to lift it above zero
    ScaleByBackground = 1.1         #this is from Igor code IN3_FixNegativeIntensities
    ScaleByIntMax = 3e-11           #this is from Igor code IN3_FixNegativeIntensities
    if MinSMR_Int < 0:
        AddToSMR_Int = abs(MinSMR_Int) * ScaleByBackground + ScaleByIntMax * MaxSMR_Int
        SMR_Int = SMR_Int + AddToSMR_Int
    # now we needd to select Qmin which makes sense for the instrument resolution and sample scattering
    # find Qmin as the first point where we get above 5% of the background level and larger than instrument resolution
    # IntRatio = Intensity / BL_Intensity, calculated using interpolation in subtract_data function
    # find point where the IntRatio is larger than 1.05 = MinQMinFindRatio, after Q dependent correction

    MeasuredTransmission = ((SaTransCounts/SaTransGain)/(SaI0Counts/SaI0Gain))/((BlTransCounts /BlTransGain )/(BlI0Counts/BlI0Gain))
    MSAXSCorrection = MeasuredTransmission / PeakToPeakTransmission
    QminSample = 4*np.pi*np.sin(np.radians(FWHMSample)/2)/wavelength
    QminBlank = 4*np.pi*np.sin(np.radians(FWHMBlank)/2)/wavelength
    indexSample = np.searchsorted(Q, QminSample)+1
    indexBlank = np.searchsorted(Q, QminBlank)+1
    # now we need to reproduce the Q correction for 1.05 from Igor. These are values from there... 
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
    blankname = Sample["BlankData"]["blankname"]
    #Igor:	variable SlitLength=0.5*((4*pi)/wavelength)*sin(PhotoDiodeSize/(2*SDDistance))
    slitLength = 0.5*((4*np.pi)/wavelength)*np.sin(UPDSize/(2*SDD))
    OmegaFactor= (UPDSize/SDD)*np.radians(FWHMBlank)
    Kfactor=BLPeakMax*OmegaFactor*thickness * 0.1 
    #apply calibration
    SMR_Int =  SMR_Int / (Kfactor*MSAXSCorrection) 
    SMR_Error = SMR_Error/ (Kfactor*MSAXSCorrection) 
    SMR_Error = SMR_Error * PeakToPeakTransmission  #this is Igor correction from 2014 which fixes issues with high absorption well scattering samples. 
    return {"SMR_Qvec":SMR_Qvec,
            "SMR_Int":SMR_Int,
            "SMR_Error":SMR_Error,
            "Kfactor":Kfactor,
            "OmegaFactor":OmegaFactor,
            "blankname":blankname,
            "thickness":thickness,
            "slitLength":slitLength,
            "units":"[cm2/cm3]"
            }


def calculatePDErrorFly(Sample, isBlank=False):
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
    I0Gain=Sample["RawData"]["metadata"]["I0Gain"]
    VToFFactor = Sample["RawData"]["VToFFactor"]                    #this is mca1 frequency, hardwired to 1e6 
    SigmaUSAXSPD=np.sqrt(UPD_array*(1+0.0001*UPD_array))		    #this is our USAXS_PD error estimate, Poisson error + 1% of value
    SigmaPDwDC=np.sqrt(SigmaUSAXSPD**2+(MeasTime*UPD_bkgErr)**2)    #This should include now measured error for background
    SigmaPDwDC=SigmaPDwDC/(VToFFactor*UPD_gains)
    A=(UPD_array)/(VToFFactor[0]*UPD_gains)		                    #without dark current subtraction
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




def _parse_filename_info(filename):
    """
    Parses a filename like 'name_XYZ.ext' to extract base, number, and extension.
    Returns (base_prefix, number_int, extension) or (None, None, None) if parsing fails.
    Example: "Scan_001.h5" -> ("Scan_", 1, ".h5")
    """
    name_part, ext_part = os.path.splitext(filename)
    # Regex to find a base name ending with an underscore, followed by digits
    match = re.match(r'(.*_)([0-9]+)$', name_part)
    if match:
        base_prefix = match.group(1)
        number_str = match.group(2)
        try:
            number_int = int(number_str)
            return base_prefix, number_int, ext_part
        except ValueError:
            # This case should ideally not happen if regex matches digits
            return None, None, None
    return None, None, None

 
def calculatePD_Fly(data_dict):
        # create the gains array and corrects UPD for it.
        # Masks deadtimes and range changes
        # get the needed data from dictionary
    ARangles = data_dict["RawData"]["ARangles"]
    AmpGain = data_dict["RawData"]["AmpGain"]
    AmpReqGain = data_dict["RawData"]["AmpReqGain"]
    Channel = data_dict["RawData"]["Channel"]
    metadata_dict = data_dict["RawData"]["metadata"]
    instrument_dict = data_dict["RawData"]["instrument"]
    UPD_array = data_dict["RawData"]["UPD_array"]
    TimePerPoint = data_dict["RawData"]["TimePerPoint"]
    Monitor = data_dict["RawData"]["Monitor"]
    VToFFactor = data_dict["RawData"]["VToFFactor"]

    
        # Create Gains arrays - one for requested and one for real
    I0Gain = metadata_dict["I0Gain"]
    num_elements = UPD_array.size 
    AmpGain_array = np.full(num_elements, AmpGain[len(AmpGain)-1])
    AmpGainReq_array = np.full(num_elements,AmpGain[len(AmpReqGain)-1])

        # Iterate over the Channel array to get index pairs
    for i in range(0, len(Channel)-2, 1):
        start_index = int(Channel[i])
        end_index = int(Channel[i + 1])
            
        # Ensure indices are within bounds
        if start_index < 0 or end_index > len(AmpGain_array) or start_index > end_index:
            raise ValueError("Invalid index range in Channel array.")
        
        # Fill the new array with AmpGain and AmpGainReq values between the two indices
        if(start_index!=end_index):
            AmpGain_array[start_index:end_index] = AmpGain[i]    
            AmpGainReq_array[start_index:end_index] = AmpReqGain[i]
        else:
            AmpGain_array[start_index] = AmpGain[i]    
            AmpGainReq_array[start_index] = AmpReqGain[i]

        # Create a new array res with the same shape as AmpGain_array, initialized with NaN
    GainsIndx = np.full(AmpGain_array.shape, np.nan)
    Gains = np.full(AmpGain_array.shape, np.nan)
    updBkg = np.full(AmpGain_array.shape, np.nan)
    updBkgErr = np.full(AmpGain_array.shape, np.nan)

        # Use a boolean mask to find where two arrays agree
    mask = AmpGain_array == AmpGainReq_array

        # Set the values in res where two agree
    GainsIndx[mask] = AmpGain_array[mask]
        #set to Nan also points in channel array that are not in mask
    for i in range(0, len(Channel)-2, 1):
        s = int(Channel[i])
        GainsIndx[s] = np.nan

        #next, replace the values in Gains array with values looked up from metadata dictionary
        #for now, lets look only for DDPCA300_gain+"gainNumber"
        # also need to create gain measureds background 'upd_bkg0', 'upd_bkg1','upd_bkg2','upd_bkg3', 'upd_bkg4'
    for i in range(0, len(GainsIndx)-1, 1):
        if np.isnan(GainsIndx[i]):
            continue
        gainName = 'DDPCA300_gain'+str(int(GainsIndx[i]))
        Gains[i] = metadata_dict[gainName]
        updBkgName = 'upd_bkg'+str(int(GainsIndx[i]))
        updBkg[i] =  metadata_dict[updBkgName]
        try:
            updBkgErrName = 'upd_bkgErr'+str(int(GainsIndx[i]))
            updBkgErr[i] =  metadata_dict[updBkgErrName]
        except:
            updBkgErrName = 'upd_bkg_err'+str(int(GainsIndx[i]))     #typo in Flyscan schema below 1.3 (before June 2025) 
            updBkgErr[i] =  metadata_dict[updBkgErrName]

        #mask amplifier dead times. This is done by comparing table fo deadtimes from metadata with times after range change. 
    Frequency= 1e6      #VToFFactor[0]/10   #this is frequency of clock fed into mca1/10 for HDF5 writer 1.3 and higher
    TimeInSec = TimePerPoint/Frequency
    #print("Exp. time :", sum(TimeInSec))
    for i in range(0, len(Channel)-1, 1):
        startPnt=Channel[i]
        deadtimeName = 'upd_amp_change_mask_time'+str(int(AmpReqGain[i]))
        deadtime = metadata_dict[deadtimeName]
        elapsed = 0
        indx = int(startPnt)
        while elapsed < deadtime:
            elapsed+=TimeInSec[indx]
            Gains[indx]=np.nan
            #print("Index is:", indx)
            #print("elapsed time is:",elapsed) 
            indx += 1

        #Correct UPD for gains and monitor counts and amplifier gain. 
        # Frequency=1e6  #this is to keep in sync with Igor code. 
    PD_Intensity = ((UPD_array-TimeInSec*updBkg)/(Frequency*Gains)) / (Monitor/I0Gain)  
    PD_error = 0.01*PD_Intensity   #this is fake error for QR conversion,recalculated in calculatePDerror
    result = {"Intensity":PD_Intensity,
              "Error":PD_error,
              "UPD_gains":GainsIndx,
              "UPD_gains":Gains,
              "UPD_bkgErr":updBkgErr}
    return result


def rebinData(data_dict,num_points=200, isSMRData=False):
    # Rebin data to 200+peak area points.
    if isSMRData:
        SMR_Int = data_dict["CalibratedData"]["SMR_Int"]
        SMR_Qvec = data_dict["CalibratedData"]["SMR_Qvec"]
        SMR_Error = data_dict["CalibratedData"]["SMR_Error"]
        SMR_QvecNew, SMR_IntNew, SMR_ErrorNew, SMR_dQ = rebin_QRSdata(SMR_Qvec, SMR_Int,SMR_Error, num_points)
        results = {"SMR_Qvec":SMR_QvecNew,
                "SMR_Int":SMR_IntNew,
                "SMR_Error":SMR_ErrorNew,
                "SMR_dQ":SMR_dQ,
                }    
    else:
        Q_array = data_dict["reducedData"]["Q"]
        R_array = data_dict["reducedData"]["Intensity"]
        S_array = data_dict["reducedData"]["Error"]
        Q_arrayNew, R_arrayNew, S_arrayNew = rebin_QRSdata(Q_array, R_array,S_array, num_points)
        results = {"Q":Q_arrayNew,
                "Intensity":R_arrayNew,
                "Error":S_arrayNew  
                }
    return results

## Common steps go here
def beamCenterCorrection(data_dict, useGauss=1, isBlank=False):
    # Find Peak center and create Q vector.
        #RawData=data_dict["rawData"]
        #reducedData = data_dict["reducedData"]
    ARangles = data_dict["RawData"]["ARangles"]
    instrument_dict = data_dict["RawData"]["instrument"]
    if isBlank:
        UPD_array = data_dict["BlankData"]["Intensity"]
    else:
        UPD_array = data_dict["reducedData"]["Intensity"]

        #plt.figure(figsize=(6, 12))
        #plt.plot(ARangles, UPD_array, marker='o', linestyle='-')  # You can customize the marker and linestyle
        # Remove NaN values from both xdata and ydata
    nan_mask = ~np.isnan(ARangles) & ~np.isnan(UPD_array)
    xdata_clean = ARangles[nan_mask]
    ydata_clean = UPD_array[nan_mask]
        #plt.plot(xdata_clean, ydata_clean, marker='o', linestyle='-') 
        # Find the threshold for the top ~40% of UPD_array
    threshold = np.max(ydata_clean)/2.3
        #pp.pprint(threshold)
        #print(f"Threshold for the top 40% of UPD_array: {threshold}")
        # Filter the data to only include the top 40%
    mask = ydata_clean >= threshold
    xdata_filtered = xdata_clean[mask]
    ydata_filtered = ydata_clean[mask]
        #plt.plot(xdata_filtered, ydata_filtered, marker='o', linestyle='-')

    
    if useGauss:
        # Initial guess for the parameters: amplitude, mean, and standard deviation, 2 fgor d parameter
        
        initial_guess = [np.max(ydata_filtered), xdata_filtered[np.argmax(ydata_filtered)], 0.0001]
        # Fit the Gaussian function to the filtered data
        popt, _ = curve_fit(gaussian, xdata_filtered, ydata_filtered, p0=initial_guess)
        # Extract the fitted parameters
        amplitude, x0, sigma = popt
        # Calculate the FWHM
        fwhm = 2 * np.abs(np.sqrt(2 * np.log(2)) * sigma)        
        # Calculate the predicted y values using the fitted parameters
        y_pred = gaussian(xdata_filtered, *popt)       
    else:
        initial_guess = [np.max(ydata_filtered), xdata_filtered[np.argmax(ydata_filtered)], 0.0001, 1.98]
        #print(initial_guess)
        popt, _ = curve_fit(modifiedGauss, xdata_filtered, ydata_filtered, p0=initial_guess)

        # Extract the fitted parameters
        amplitude, x0, sigma, exponent = popt
      
        # Calculate the FWHM
        # Calculate the half maximum
        half_max = amplitude / 2

        # but next calculation needs to be done over larger q range
        threshold = amplitude / 5   #seems to be failing occasionally to find low enough calculated intensity to cross 50%
        mask = ydata_clean >= threshold
        xdata_calc = xdata_clean[mask]
        ydata_calc = ydata_clean[mask]
     
        # Calculate the predicted y values using the fitted parameters
        ydata_calc = modifiedGauss(xdata_calc, *popt)

        # Find where the array crosses the half maximum
        crossings = np.where((ydata_calc[:-1] < half_max) & (ydata_calc[1:] >= half_max) |
                     (ydata_calc[:-1] >= half_max) & (ydata_calc[1:] < half_max))[0]

        # Calculate fractional crossing indices using linear interpolation
        indices = []
        for i in crossings:
            y1, y2 = ydata_calc[i], ydata_calc[i + 1]
            if y2 != y1:
                fractional_index = i + (half_max - y1) / (y2 - y1)
                indices.append(fractional_index)
                
        # Calculate the FWHM
        if len(indices) >= 2:
            i = int(np.floor(indices[0]))
            y1 = xdata_calc[i]
            y2 = xdata_calc[i+1]
            xdata_calcl = y1 + (indices[0] - i) * (y2 - y1)
            i = int(np.floor(indices[-1]))
            y1 = xdata_calc[i]
            y2 = xdata_calc[i+1]
            xdata_calch = y1 + (indices[-1] - i) * (y2 - y1)           
            fwhm = np.abs(xdata_calch - xdata_calcl)
        else:
            fwhm = np.nan  # FWHM cannot be determined

        # Calculate the residuals
        y_pred = modifiedGauss(xdata_filtered, *popt)

    #use above calculated y_pred to get residuals
    residuals = ydata_filtered - y_pred

        # Calculate the chi-square
        # If you have measurement errors, replace 1 with the variance of ydata_filtered
    chi_square = np.sum((residuals**2) / 1)

        #Make wave vector
    Q_array = np.full(UPD_array.shape, 0)
        #AR_center = metadata_dict["AR_center"]
    try:
        # Try to get the value using the first key
        wavelength = instrument_dict["monochromator"]["wavelength"]
    except KeyError:
        #print(instrument_dict)
        # If the first key doesn't exist, try the second key
        wavelength = instrument_dict["wavelength"]

    Q_array = -1*(4*np.pi*np.sin(np.radians(ARangles-x0)/2)/wavelength)
        #Q_array = (4*np.pi*np.sin(np.radians(ARangles-x0)/2)/wavelength)
        #Q_array_log = np.sign(Q_array)*np.log(np.abs(Q_array))
        #pp.pprint(Q_array)

    results = {"Q":Q_array,
            "Chi-Square":chi_square,
            "Center":x0,
            "Maximum":amplitude,
            "FWHM":fwhm,
            "wavelength":wavelength,    
            }
    return results


def smooth_r_data(intensity, qvector, UPD_gains, r_error, meas_time, replaceNans=True):
    # Smoothing times for different ranges
    rwave_smooth_times = [0.02, 0.02, 0.03, 0.1, 0.4]   # these are [in sec] values for USAXS on 4/20/2025

    # Logarithm of intensity
    # Replace negative or zero intensities with NaN before taking log
    intensity_safe = np.where(intensity > 0, intensity, np.nan)
    temp_int_log = np.log(intensity_safe)
    # for some cases replace nans with interpolated values
    if replaceNans:
        # Create a mask for NaNs
        nan_mask = np.isnan(temp_int_log)
        # Indices of non-NaN values
        x_non_nan = qvector[~nan_mask]
        y_non_nan = temp_int_log[~nan_mask]
        # Create an interpolation function
        interp_func = interp1d(x_non_nan, y_non_nan, kind='linear', fill_value='extrapolate')
        # Replace NaNs with interpolated values
        temp_int_log[nan_mask] = interp_func(qvector[nan_mask])


    smooth_intensity = np.copy(temp_int_log)
    meas_time_sec = meas_time/1e6       # meas_time is still frequency, need time in seconds. 

    def linear_fit(x, a, b):
        return a + b * x
    
    #we need to find i where qvector > 0.0003
    startIndex = find_crossing_index(qvector, 0.0003)

    for i in range(startIndex, len(intensity)):
        if UPD_gains[i] == 1:
            tmp_time = rwave_smooth_times[0]
        elif UPD_gains[i] == 2:
            tmp_time = rwave_smooth_times[1]
        elif UPD_gains[i] == 3:
            tmp_time = rwave_smooth_times[2]
        elif UPD_gains[i] == 4:
            tmp_time = rwave_smooth_times[3]
        else:
            tmp_time = rwave_smooth_times[4]

        if meas_time_sec[i] > tmp_time:
            smooth_intensity[i] = temp_int_log[i]
        else:
            start_points = int(np.ceil(tmp_time / meas_time_sec[i])) + 1
            end_points = start_points

            if (i - start_points) < 0:
                raise ValueError("Bad data, cannot fix this. Likely Flyscan parameters were wrong")

            if i + end_points > len(intensity) - 1:
                end_points = len(intensity) - 1 - i

            if (UPD_gains[i - start_points] != UPD_gains[i]) or (UPD_gains[i + end_points] != UPD_gains[i]):
                temp_r = temp_int_log[i - start_points:i + end_points]
                temp_q = qvector[i - start_points:i + end_points]

                if len(temp_r) > np.isnan(temp_r).sum() + 5:
                    popt, _ = curve_fit(linear_fit, temp_q, temp_r)
                    smooth_intensity[i] = linear_fit(qvector[i], *popt)
                    r_error[i] /= 3
                else:
                    smooth_intensity[i] = temp_int_log[i]
                    r_error[i] = r_error[i]
            else:
                temp_r = temp_int_log[i - start_points:i + end_points + 1]
                temp_q = qvector[i - start_points:i + end_points + 1]
                start_x = temp_q[0]
                end_x = temp_q[-1]
                area = np.trapezoid(temp_r, temp_q)
                smooth_intensity[i] = area / (end_x - start_x)
                r_error[i] = r_error[i]

    intensity = np.exp(smooth_intensity)
    #occasionally we are running out of number resolution here and we are getting 0 as exp(smooth_value)
    #this is a hack to fix that, the dynamic range o fthe nstrumet is at most 11 decades, so let's fix background to that. 
    MaxInt = np.max(intensity)
    intensity[intensity < MaxInt*1e-11 ] = MaxInt*1e-11
    return  {"Intensity":intensity,
              "Error":r_error} 

def find_crossing_index(array, target_value):
    for index, value in enumerate(array):
        if value >= target_value:
            return index
    return 0.1*len(array)  # Return None if the target value is not crossed

# subtract QRS data
def subtract_data(X1, Y1, E1, X2, Y2, E2):
    """
    Interpolates and subtracts the input data sets to return X1, Ydiff, and Esub.

    Parameters:
    X1 (array-like): X1 data points.
    Y1 (array-like): Y1 data points.
    E1 (array-like): Uncertainty in Y1.
    X2 (array-like): X2 data points.
    Y2 (array-like): Y2 data points.
    E2 (array-like): Uncertainty in Y2.

    Returns:
    tuple: (X1, Ydiff, Esub, IntRatio) where
        - X1 is the input X1 data points.
        - Ydiff is the difference between Y1 and interpolated Y2.
        - Esub is the propagated uncertainty.
        - IntRatio is the ratio of Y1 to Y2.
    """    
    # Step 1: Interpolate log(Y2) vs X2 to values of X1
    # we need to do some gymnastics developed in IN2G_LogInterpolateIntensity
    # this offsets intensity high enough that log does not faill.
    # find Y2 minimum
    Y2_min = np.min(Y2)
    #if Y2_min<1e-30, offset whole Y2 by 3*abs(Y2_min)
    offset=0
    if Y2_min<1e-30:
        offset =  3*abs(Y2_min)
    
    Y2 = Y2 + offset
    logY2 = np.log(Y2)
    logY2_interp_func = interp1d(X2, logY2, kind='linear', fill_value='extrapolate')
    logY2_interp = logY2_interp_func(X1) #this is the interpolation. 

    # Convert interpolated logY2 back to Y2interp
    Y2_interp = np.exp(logY2_interp)

    #fix the offset, if it was added:
    if offset>0:
        Y2_interp = Y2_interp - offset

    # Step 2: Subtract Y1 - Y2interp to obtain Ydiff and divide to get IntRatio
    Ydiff = Y1 - Y2_interp

    IntRatio = np.abs(Y1 / Y2_interp)  # Calculate the intensity ratio
    
    # Step 3: Linearly interpolate E2 vs X2 to have E2interp vs X1
    E2_interp_func = interp1d(X2, E2, kind='linear', fill_value='extrapolate')
    E2_interp = E2_interp_func(X1)

    # Step 4: Propagate uncertainties for subtraction to obtain Esub
    Esub = np.sqrt(E1**2 + E2_interp**2)

    # Return the three data sets: X1, Ydiff, and Esub
    return X1, Ydiff, Esub, IntRatio

# Function to check if all arrays have the same length
def check_arrays_same_length(*arrays):
    lengths = [arr.size for arr in arrays]  # Get the size of each array
    if len(set(lengths)) != 1:  # Check if all lengths are the same
        raise ValueError("Not all arrays have the same length.")
    #else:
        #print("All arrays have the same length.")

# Gaussian function
def gaussian(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

#modified gaussian function, gives better fits to peak profiles. 
def modifiedGauss(xvar, a, x0, sigma, exponent):
    base = np.abs(xvar - np.abs(x0)) / (2 * np.abs(sigma))
    # Diagnostics: check for invalid values
    if np.any(base < 0) or np.any(np.isnan(base)) or np.any(np.isnan(exponent)):
        logging.warning("Warning: Invalid value encountered in base or exponent in modifiedGauss.")
        logging.info(f"Base min: {np.nanmin(base)}, max: {np.nanmax(base)}, sigma: {sigma}, exponent: {exponent}")
    try:
        result =np.abs(a) * np.exp(-np.power(base, np.abs(exponent)))
    except Exception as e:
        logging.warning(f"Exception in modifiedGauss: {e}")
        result = np.full_like(xvar, np.nan)
    return result



# Function to recursively read a group and store its datasets in a dictionary
def read_group_to_dict(group):
    data_dict = {}
    for key, item in group.items():
        if isinstance(item, h5py.Dataset):
            # Read the dataset
            data = item[()]
             # Check if the dataset is bytes
            if isinstance(data, bytes):
                # Decode bytes to string
                data = data.decode('utf-8')
            # Check if the dataset is an array with a single element
            elif hasattr(data, 'size') and data.size == 1:
                # Convert to a scalar (number or string)
                data = data.item()
                if isinstance(data, bytes):
                    # Decode bytes to string, the above does not seem to catch this? 
                    data = data.decode('utf-8')
            data_dict[key] = data
        elif isinstance(item, h5py.Group):
            # If the item is a group, recursively read its contents
            data_dict[key] = read_group_to_dict(item)
    return data_dict


# this should not fail if keys on the list are not present
def filter_nested_dict(d, keys_to_keep):
    if isinstance(d, dict):
        return {k: filter_nested_dict(v, keys_to_keep) for k, v in d.items() if k in keys_to_keep and k in d}
    elif isinstance(d, list):
        return [filter_nested_dict(item, keys_to_keep) for item in d]
    else:
        return d    

# def results_to_dataset(results):
#     results = copy.deepcopy(results)
#     ds = xr.Dataset()
#     ds['USAXS_int'] = ('q',results['reducedData']['UPD'])
#     ds['q'] = results['reducedData']['Q_array']
#     del results['reducedData']['UPD']
#     del results['reducedData']['Q_array']
#     ds.update(results['reducedData'])
#     for our_name,raw_name in [('AR_angle','ARangles'),
#                               ('TimePerPoint','TimePerPoint'),
#                               ('Monitor','Monitor'),
#                               ('UPD','UPD_array'),
#                              ]:
#         ds[our_name] = ('flyscan_bin',results['RawData'][raw_name])
#         del results['RawData'][raw_name]
#     for our_name,raw_name in [('AmpGain','AmpGain'),
#                               ('AmpReqGain','AmpReqGain'),
#                               ('amp_change_channel','Channel')
#                              ]:
#         ds[our_name] = ('amp_change_channel',results['RawData'][raw_name])
#         del results['RawData'][raw_name]
                                      
#     ds.attrs.update(results['RawData']['metadata'])
#     del results['RawData']['metadata']
#     ds.attrs['instrument'] = results['RawData']['instrument']
#     del results['RawData']['instrument']
#     ds.update(results['RawData'])

#     return ds


'''
    Converted by AI from Igor code

    rebin_data() routine will rebin data from Q=0.0002 higher using rebin_log_data(), min step is set to distacne between Q[0] and Q[1]
    
    This routine will rebin data on log scale. It will produce new Wx and Wy with new NumberOfPoints
    If MinStep > 0 it will try to set the values so the minimum step on log scale is MinStep
    optional Wsdev is standard deviation for each Wy value, it will be propagated through - sum(sdev^2)/numpnts in each bin. 
    optional Wxwidth will generate width of each new bin in x. NOTE: the edge is half linear distance between the two points, no log  
    skewing is done for edges. Therefore the width is really half of the distance between p-1 and p+1 points.  
    optional W1-5 will be averaged for each bin , so this is way to propagate other data one may need to. 

    typical min step is tempMinStep= Qvec[1] - Qvec[0], the differecne between the two first points.
    Igor code calling this:
    	tempMinStep=DSM_Qvec[1]-DSM_Qvec[0]
		IN2G_RebinLogData(DSM_Qvec,DSM_Int,FlyScanRebinToPoints,tempMinStep,Wsdev=DSM_Error,Wxwidth=DSM_dQ)
        where DSM input waves are with high number of points, and are returned in place as reduced number of points
'''


def rebin_QRSdata(Wx, Wy, Ws, NumberOfPoints):
    '''
    Rebin data based on the condition Q < 0.0002 and Q > 0.0002.
    This function will split the data into two parts and rebin the second part using logarithmic scaling.
    The first part (Q < 0.0002) will remain unchanged, while the second part (Q > 0.0002) will be rebinned.
    The function will return the merged data from both parts.
    The rebinning is done using the rebin_log_data function, which takes care of the logarithmic scaling.     
    '''
    # create Wdx array based on Wx
    Wdx = np.zeros(len(Wx))
    Wdx[1:] = Wx[1:] - Wx[:-1]
    Wdx[0] = Wdx[1]  # Set the first element to the same value as the second element        
    # Split arrays based on the condition Q < 0.0002
    mask_less = Wx < 0.0002
    Wx_less = Wx[mask_less]
    Wy_less = Wy[mask_less]
    Ws_less = Ws[mask_less]
    Wdx_less = Wdx[mask_less]

    # Split arrays based on the condition Q > 0.0002
    mask_greater = Wx > 0.0002
    Wx_greater = Wx[mask_greater]
    Wy_greater = Wy[mask_greater]
    Ws_greater = Ws[mask_greater]
    Wdx_greater = Wdx[mask_greater]

    MinStep = Wx_greater[1] - Wx_greater[0]

    Wx_greater2, Wy_greater2, W1, W2, W3, W4, W5, Ws_greater2, Wxsdev, Wxwidth = rebin_log_data(Wx_greater, Wy_greater, NumberOfPoints, MinStep, Wsdev=Ws_greater, Wxsdev=None, Wxwidth=Wdx_greater, W1=None, W2=None, W3=None, W4=None, W5=None)


    Q_merged = np.concatenate((Wx_less, Wx_greater2))
    Intensity_merged = np.concatenate((Wy_less, Wy_greater2))      
    Error_merged = np.concatenate((Ws_less, Ws_greater2))    
    dQ_merged = np.concatenate((Wdx_less, Wxwidth)) 
    # remove Nans form errors as that seems to break stuff latrer
    # Create a mask for NaNs
    nan_mask = np.isnan(Error_merged)
    # Indices of non-NaN values
    x_non_nan = Q_merged[~nan_mask]
    y_non_nan = Error_merged[~nan_mask]
    # Create an interpolation function
    interp_func = interp1d(x_non_nan, y_non_nan, kind='linear', fill_value='extrapolate')
    # Replace NaNs with interpolated values
    Error_merged[nan_mask] = interp_func(Q_merged[nan_mask])   


    return Q_merged, Intensity_merged, Error_merged, dQ_merged



def rebin_log_data(Wx, Wy, NumberOfPoints, MinStep, Wsdev=None, Wxsdev=None, Wxwidth=None, W1=None, W2=None, W3=None, W4=None, W5=None):
    # Determine which additional waves need to be calculated
    CalcSdev = Wsdev is not None
    CalcXSdev = Wxsdev is not None
    CalcWidth = Wxwidth is not None
    CalcW1 = W1 is not None
    CalcW2 = W2 is not None
    CalcW3 = W3 is not None
    CalcW4 = W4 is not None
    CalcW5 = W5 is not None

    OldNumPnts = len(Wx)
    if 2 * NumberOfPoints > OldNumPnts:
        logging.warning(f"Rebinning requested to {NumberOfPoints} points, but data has only {OldNumPnts} points. No rebinning will be done.")
        # Return original data if no rebinning is done
        return Wx, Wy, W1, W2, W3, W4, W5, Wsdev, Wxsdev, Wxwidth

    if Wx[0] <= 0:
        Wx[0] = Wx[1] / 2
    CorrectStart = Wx[0]

    if MinStep > 0:
        StartX = find_correct_log_scale_start(Wx[0], Wx[-1], NumberOfPoints, MinStep)
    else:
        StartX = CorrectStart

    EndX = StartX + abs(Wx[-1] - Wx[0])
    isGrowing = Wx[0] < Wx[-1]

    tempNewLogDist = np.zeros(NumberOfPoints)
    tempNewLogDistBinWidth = np.zeros(NumberOfPoints)

    logstartX = np.log10(StartX)
    logendX = np.log10(EndX)
    tempNewLogDist = np.logspace(logstartX, logendX, NumberOfPoints, base=10)
    tempNewLogDist += CorrectStart - StartX

    tempNewLogDistBinWidth[1:-1] = tempNewLogDist[2:] - tempNewLogDist[:-2]
    tempNewLogDistBinWidth[0] = tempNewLogDistBinWidth[1]
    tempNewLogDistBinWidth[-1] = tempNewLogDistBinWidth[-2]

    Rebinned_WvX = np.zeros(NumberOfPoints)
    Rebinned_WvY = np.zeros(NumberOfPoints)
    Rebinned_Wv1 = np.zeros(NumberOfPoints) if CalcW1 else None
    Rebinned_Wv2 = np.zeros(NumberOfPoints) if CalcW2 else None
    Rebinned_Wv3 = np.zeros(NumberOfPoints) if CalcW3 else None
    Rebinned_Wv4 = np.zeros(NumberOfPoints) if CalcW4 else None
    Rebinned_Wv5 = np.zeros(NumberOfPoints) if CalcW5 else None
    Rebinned_Wsdev = np.zeros(NumberOfPoints) if CalcSdev else None
    Rebinned_Wxsdev = np.zeros(NumberOfPoints) if CalcXSdev else None

    j = 0
    for i in range(NumberOfPoints):
        cntPoints = 0
        BinHighEdge = tempNewLogDist[i] + tempNewLogDistBinWidth[i] / 2
        while j < OldNumPnts and ((Wx[j] < BinHighEdge) if isGrowing else (Wx[j] > BinHighEdge)):
            Rebinned_WvX[i] += Wx[j]
            Rebinned_WvY[i] += Wy[j]
            if CalcW1:
                Rebinned_Wv1[i] += W1[j]
            if CalcW2:
                Rebinned_Wv2[i] += W2[j]
            if CalcW3:
                Rebinned_Wv3[i] += W3[j]
            if CalcW4:
                Rebinned_Wv4[i] += W4[j]
            if CalcW5:
                Rebinned_Wv5[i] += W5[j]
            if CalcSdev:
                Rebinned_Wsdev[i] += Wsdev[j] ** 2
            if CalcXSdev:
                Rebinned_Wxsdev[i] += Wxsdev[j] ** 2
            cntPoints += 1
            j += 1

        if cntPoints > 0:
            Rebinned_WvX[i] /= cntPoints
            Rebinned_WvY[i] /= cntPoints
            if CalcW1:
                Rebinned_Wv1[i] /= cntPoints
            if CalcW2:
                Rebinned_Wv2[i] /= cntPoints
            if CalcW3:
                Rebinned_Wv3[i] /= cntPoints
            if CalcW4:
                Rebinned_Wv4[i] /= cntPoints
            if CalcW5:
                Rebinned_Wv5[i] /= cntPoints
            if CalcSdev:
                Rebinned_Wsdev[i] = np.sqrt(Rebinned_Wsdev[i] / cntPoints)
            if CalcXSdev:
                Rebinned_Wxsdev[i] = np.sqrt(Rebinned_Wxsdev[i] / cntPoints)

    Wx = Rebinned_WvX
    Wy = Rebinned_WvY

    if CalcW1:
        W1 = Rebinned_Wv1
    if CalcW2:
        W2 = Rebinned_Wv2
    if CalcW3:
        W3 = Rebinned_Wv3
    if CalcW4:
        W4 = Rebinned_Wv4
    if CalcW5:
        W5 = Rebinned_Wv5
    if CalcSdev:
        Wsdev = Rebinned_Wsdev
    if CalcXSdev:
        Wxsdev = Rebinned_Wxsdev
    if CalcWidth:
        Wxwidth = tempNewLogDistBinWidth

    return Wx, Wy, W1, W2, W3, W4, W5, Wsdev, Wxsdev, Wxwidth

def my_find_start_value_func(x1, w):
    """
    Objective function to find the correct start value for logarithmic scaling.

    Parameters:
    - x1: The start value where we need to start with log stepping
    - w: A list or array containing [totalRange, NumSteps, MinStep]

    Returns:
    - The absolute difference between the calculated last minimum step and the desired minimum step
    """
    total_range, num_steps, min_step = w
    last_min_step = 10**(np.log10(x1) + (np.log10(x1 + total_range) - np.log10(x1)) / num_steps) - 10**(np.log10(x1))
    return abs(last_min_step - min_step)

def find_correct_log_scale_start(StartValue, EndValue, NumPoints, MinStep):
    """
    Finds the correct start value for logarithmic scaling using optimization.

    Parameters:
    - StartValue: The initial start value
    - EndValue: The end value
    - NumPoints: The number of points
    - MinStep: The minimum step size

    Returns:
    - The optimal start value for logarithmic scaling
    """
    # Define the parameters for the optimization
    w = [EndValue - StartValue, NumPoints, MinStep]

    # Initial guess for the start value
    x1_initial = StartValue

    # Use scipy.optimize.minimize to find the optimal start value
    result = minimize(my_find_start_value_func, x1_initial, args=(w,), method='L-BFGS-B', bounds=[(1e-10, None)])

    # The optimal start value is in result.x[0]
    return result.x[0]


# # Example usage
# StartValue = 1.0
# EndValue = 10.0
# NumPoints = 100
# MinStep = 0.1

# optimal_start = find_correct_log_scale_start(StartValue, EndValue, NumPoints, MinStep)
# print("Optimal Start Value:", optimal_start)