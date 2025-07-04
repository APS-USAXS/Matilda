'''
    convertSWAXScalibrated.py
        New, calibrated SAXS/WAXS code. 
    use: 
    process2Ddata(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False)

    returns dictionary of this type:
            result["samplename"]=samplename
            result["blankname"]=blankname
            result["reducedData"] =  {"Intensity":np.ravel(intensity), 
                              "Q":np.ravel(q),
                              "Error":np.ravel(error)}
            result["CalibratedData"] = {"Intensity":np.ravel(intcalib),
                                    "Q":np.ravel(qcalib),
                                    "Error":np.ravel(errcalib),
                                    }  
    Does:
    Convert SAXS and WAXS area detector data from the HDF5 format to the 1Ddata
    Converts Nika parameters to Fit2D format and then uses pyFAI to convert to poni format
    Both SAXS and WAXS data give sufficiently same data as Nika to be considered same.  
'''


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import h5py
import pprint as pp
import socket
import os
import tifffile as tiff
import logging

from supportFunctions import read_group_to_dict, filter_nested_dict, subtract_data
from supportNikaFunctions import convert_Nika_to_Fit2D
from readfromtiled import FindLastBlankScan
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5, saveNXcanSAS, readMyNXcanSAS, find_matching_groups

# TODO: split into multiple steps as needed
# Import images for sample and blank as separate calls and get sample and blank objects
#   calibration step, calculate corrections, apply corrections
#   subtract 2D images to get calibrated image as Nika
#   convert to 1D data
# TODO: test me... 



## main code here
def process2Ddata(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    location = 'entry/reducedData/'    #we need to make sure we have separate NXcanSAS data here. Is it still entry? 
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
            location = SASentries[0] if len(SASentries) > 0 else None
            if location is not None and location in hdf_file:
                # Delete the group
                del hdf_file[location]
                logging.info(f"Deleted existing NXcanSAS group for file {filename}. ")


        NXcanSASentry = SASentries[0] if len(SASentries) > 0 else None
        if blankFilename is not None and blankPath is not None:
            location = NXcanSASentry                 # require we have calibrated data
        else:
            location = 'entry/QRS_data/'            # all we want here are QRS data

        if location is not None and location in hdf_file:
            # exists, so lets reuse the data from the file
            Sample = dict()
            Sample = readMyNXcanSAS(path, filename)
            logging.info(f"Using existing processed data from file {filename}.")
            return Sample
        
        
        else:
            Sample = dict()
            Sample = importADData(path, filename)   #this is for sample path and blank
            if "saxs" in path:
                plan_name="SAXS"
            else:
                plan_name="WAXS"

            Sample["reducedData"] = reduceADData(Sample, useRawData=True)   #this generates Int vs Q for raw data plot
                                        # q = Sample["reducedData"]["Q"]
                                        # intensity = Sample["reducedData"]["Intensity"]
                                        # error = Sample["reducedData"]["Error"]            
                                        # samplename = Sample["RawData"]["samplename"]
            
            if blankPath is not None and blankFilename is not None:               
                blank = importADData(blankPath, blankFilename)               #this is for blank path and blank name
                Sample["BlankData"] = reduceADData(blank, useRawData=True)   #this generates Int vs Q for blank data plot
                                        # qcalib = Sample["BlankData"]["Q"]
                                        # intensity = Sample["BlankData"]["Intensity"]
                                        # error = Sample["BlankData"]["Error"]            
                                        # blankname = Sample["RawData"]["blankname"] 
                Sample["calib2DData"] = calibrateAD2DData(Sample, blank)
                                    #returns 2D calibrated data
                                        # result = {"data":calib2Ddata,
                                        #           "blankname":blankname,
                                        #           "transmission":transmission
                Sample["CalibratedData"] = reduceADData(Sample, useRawData=False)  #this generates Calibrated 1D data.
                                        #returns :   
                                        # qcalib= Sample["CalibratedData"]["Q"]
                                        # dqcalib= Sample["CalibratedData"]["dQ"]
                                        # intcalib= Sample["CalibratedData"]["Intensity"]
                                        # errcalib= Sample["CalibratedData"]["Error"]
                                        # Sample["CalibratedData"]["units"]
                                        # blankname = Sample["calib2DData"]["blankname"]
            else:
                Sample["CalibratedData"]=dict()
                Sample["CalibratedData"]["Q"] = None
                Sample["CalibratedData"]["dQ"] = None
                Sample["CalibratedData"]["Intensity"] = None
                Sample["CalibratedData"]["Error"] = None
                Sample["CalibratedData"]["units"] = None
                Sample["calib2DData"]=dict()
                Sample["calib2DData"]["blankname"] = None

        hdf_file.flush()
        # The 'with' statement will automatically close the file when the block ends
    
    saveNXcanSAS(Sample,path, filename)
    return Sample


def reduceADToQR(path, filename):
        tempFilename= os.path.splitext(filename)[0]
        tempSample = {"RawData":{"filename": tempFilename}}
        # label = data_dict["RawData"]["filename"]
        # Q_array = data_dict["reducedData"]["Q_array"]
        # Intensity = data_dict["reducedData"]["PD_intensity"]
        tempSample["reducedData"]=ImportAndReduceAD(path, filename)
        #pp.pprint(tempSample)
        #pp.pprint(tempSample["RawData"]["filename"])
        return tempSample


## main code here
def ImportAndReduceAD(path, filename, recalculateAllData=False):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
        # Check if the group 'displayData' exists
        if recalculateAllData:
            # Delete the group
            del hdf_file[location]
            logging.info(f"Deleted existing group 'entry/displayData' in {filename}.")

        if location in hdf_file:
            # exists, so lets reuse the data from the file
            Sample = dict()
            Sample = load_dict_from_hdf5(hdf_file, location)
            logging.info(f"Used existing QR data from {filename}")
            q = Sample["reducedData"]["Q_array"]
            intensity = Sample["reducedData"]["Intensity"]
            result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}  
            return result
        
        
        else:
            Sample = dict()
            #read various data sets
            logging.info(f"Read file :{filename}")
            dataset = hdf_file['/entry/data/data'] 
            my2DData = np.array(dataset)
            #sample
            sample_group = hdf_file['/entry/sample']
            sample_dict = read_group_to_dict(sample_group)
            #instrument
            instrument_group = hdf_file['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            del instrument_dict['detector']['data'] #this is original of 2-d data, we do not need to cary it here second time. 
            #metadata
            keys_to_keep = ['I000_cts', 'I00_cts', 'I00_gain', 'I0_cts', 'I0_gated',
                            'I0_gain', 'I_scaling', 'Pin_TrI0', 'Pin_TrI0gain', 'Pin_TrI0gain','Pin_TrPD','Pin_TrPDgain',
                            'PresetTime', 'monoE', 'pin_ccd_center_x_pixel','pin_ccd_center_y_pixel',
                            'pin_ccd_tilt_x', 'pin_ccd_tilt_y', 'wavelength', 'waxs_ccd_center_x', 'waxs_ccd_center_y',
                            'waxs_ccd_tilt_x', 'waxs_ccd_tilt_y', 'waxs_ccd_center_x_pixel', 'waxs_ccd_center_y_pixel',
                            'scaler_freq'                     
                        ]        
            metadata_group = hdf_file['/entry/Metadata']
            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            #need to append these structures to RAW data so we can optionally save them in Igor
            Sample["RawData"]=dict()
            Sample["RawData"]["metadata"]= metadata_dict
            Sample["RawData"]["instrument"]= instrument_dict
            Sample["RawData"]["sample"]= sample_dict
            # wavelength, keep in A for Fit2D
            wavelength = instrument_dict["monochromator"]["wavelength"]
            # pixel_size, keep in mm, converted in convert_Nika_to_Fit2D to micron for Fit2D and then to m for pyFAI... 
            pixel_size1 = instrument_dict["detector"]["x_pixel_size"] #in mm in NIka, will convert to micron for Fit2D later
            # assume pixels are square, therefore size2 is not needed. No idea how to fix this in Nika or pyFAI for that matter anyway. 
            #pixel_size2 = instrument_dict["detector"]["y_pixel_size"] #in mm in NIka, will convert to micron for Fit2D later
            # detector_distance, keep in mm for Fit2D
            detector_distance = instrument_dict["detector"]["distance"] #in mm in Nika, in mm in Fit2D
            #logging.info(f"Read metadata")
            if "pin_ccd_tilt_x" in metadata_dict:                       # this is SAXS
                usingWAXS=0
                BCX= instrument_dict["detector"]["beam_center_x"]       #  This will be swapped later in convert_Nika_to_Fit2D 
                BCY = instrument_dict["detector"]["beam_center_y"]      #  This will be swapped later in convert_Nika_to_Fit2D
                HorTilt = metadata_dict["pin_ccd_tilt_x"]               #   keep in degrees for Fit2D
                VertTilt = metadata_dict["pin_ccd_tilt_y"]              #   keep in degrees for Fit2D
            else:                                                       # and this is WAXS
                usingWAXS=1
                BCX = instrument_dict["detector"]["beam_center_x"]      #  This will be swapped later in convert_Nika_to_Fit2D
                BCY = instrument_dict["detector"]["beam_center_y"]      #  This will be swapped later in convert_Nika_to_Fit2D
                HorTilt = metadata_dict["waxs_ccd_tilt_x"]              #   keep in degrees for Fit2D
                VertTilt = metadata_dict["waxs_ccd_tilt_y"]             #   keep in degrees for Fit2D    

            #logging.info(f"Finished reading metadata")
        
        # poni is geometry file for pyFAI, created by converting first to Fit2D and then calling pyFAI conversion function.
        my_poni = convert_Nika_to_Fit2D(SSD=detector_distance, pix_size=pixel_size1, BCX=BCX, BCY=BCY, HorTilt=HorTilt, VertTilt=VertTilt, wavelength=wavelength)

        #create mask here. Duplicate the my2DData and set all values above 1e7 to NaN for WAXS or for SAXS mask all negative intensities
        # the differecne is due to Pilatus vs Eiger handing bad pixels differently. Dectris issue... 
        if usingWAXS:
            mask = np.copy(my2DData)
            mask = 0*mask   # set all values to zero
            mask[my2DData > 1e7] = 1
            mask[:, 511:516] = 1
            mask[:, 1026:1041] = 1
            mask[:, 1551:1556] = 1
        else:
            mask = np.copy(my2DData)
            mask = 0*mask   # set all values to zero
            mask[my2DData < 0] = 1
            # Set the first 4 rows to 1
            mask[:, :4] = 1
            # Set rows 192 to 195 to 1
            mask[:, 242:245] = 1

        #logging.info(f"Finished creating mask")
        
        #now define integrator... using pyFAI here. 
        # this does not work, they really do not have way to pass whole poni in? 
        #ai = AzimuthalIntegrator(poni=my_poni)
        # but this works fine
        ai = AzimuthalIntegrator(dist=my_poni.dist, poni1=my_poni.poni1, poni2=my_poni.poni2, rot1=my_poni.rot1, rot2=my_poni.rot2,
                            rot3=my_poni.rot3, pixel1=my_poni.detector.pixel1, pixel2=my_poni.detector.pixel2, 
                            wavelength=my_poni.wavelength)
        
        #   You can specify the number of bins for the integration
        #   set npt to larger of dimmension of my2DData
        if usingWAXS:
            npt = max(my2DData.shape)
        else:
            npt=200 
        #npt = 1000  # Number of bins, if should be lower
        # Perform azimuthal integration
        q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
        #logging.info(f"Finished 2d to 1D conversion")
        Sample["reducedData"] = dict()
        Sample["reducedData"]["Q_array"] = q
        Sample["reducedData"]["Intensity"] = intensity
        save_dict_to_hdf5(Sample, location, hdf_file)
        logging.info(f"Appended new QR data to 'entry/displayData' in {filename}.")
        result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}
        return result




def importADData(path, filename):
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as hdf_file:
            Sample = dict()
            #read various data sets
            logging.info(f"Read file :{filename}")
            dataset = hdf_file['/entry/data/data'] 
            my2DData = np.array(dataset)
            #metadata
            instrument_group = hdf_file['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            del instrument_dict['detector']['data']
            #metadata
            keys_to_keep = ['I000_cts', 'I00_cts', 'I00_gain', 'I0_cts', 'I0_cts_gated',
                            'TR_cts_gated','TR_cts','TR_gain','I0_Sample',
                            'I0_gain', 'I_scaling', 'Pin_TrI0', 'Pin_TrI0gain', 'Pin_TrPD','Pin_TrPDgain',
                            'PresetTime', 'monoE', 'pin_ccd_center_x_pixel','pin_ccd_center_y_pixel',
                            'pin_ccd_tilt_x', 'pin_ccd_tilt_y', 'wavelength', 'waxs_ccd_center_x', 'waxs_ccd_center_y',
                            'waxs_ccd_tilt_x', 'waxs_ccd_tilt_y', 'waxs_ccd_center_x_pixel', 'waxs_ccd_center_y_pixel',
                            'scaler_freq', 'StartTime',                     
                        ]        
            metadata_group = hdf_file['/entry/Metadata']
            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            sample_group = hdf_file['entry/sample']
            sample_dict = read_group_to_dict(sample_group)
            control_group = hdf_file['/entry/control']
            control_dict = read_group_to_dict(control_group)
            Sample["RawData"] = dict()
            Sample["RawData"]["data"] = my2DData
            Sample["RawData"]["filename"] = filename
            Sample["RawData"]["samplename"] = sample_dict["name"]
            Sample["RawData"]["instrument"] = instrument_dict
            Sample["RawData"]["metadata"] = metadata_dict
            Sample["RawData"]["sample"] = sample_dict
            Sample["RawData"]["control"] = control_dict
            #path different names between SWAXS and USAXS
            Sample["RawData"]["metadata"]["timeStamp"]=Sample["RawData"]["metadata"]["StartTime"]
            #logging.info(f"Finished reading data")
            #logging.info(f"Read data")
            return Sample

def calibrateAD2DData(Sample, Blank):
    '''
        Here is how we are suppose to process the data:
        Int = Corrfactor / I0 / SampleThickness * (Sa2D/Transm * -  I0/I0Blank * Blank2D)
        SolidAngeCorr - is done by pyFAI later, no need to do here... 
        Here is lookup from Nika:
        SAXS and WAXS are same : 
        SampleThickness = entry:sample:thickness
        SampleI0 = entry:control:integral
        SampleMeasurementTime = entry:control:preset
        Corrfactor = entry:Metadata:I_scaling
    '''
    blankname = Blank["RawData"]["filename"]     
    sampleThickness=Sample["RawData"]["sample"]["thickness"]
    #sampleMeasurementTime=Sample["RawData"]["control"]["preset"]
    corrFactor=Sample["RawData"]["metadata"]["I_scaling"]
    #blankMeasurementTime=Blank["RawData"]["control"]["preset"]
    sample2Ddata=Sample["RawData"]["data"]
    blank2Ddata = Blank["RawData"]["data"]
    metadata_dict = Sample["RawData"]["metadata"]
    #tranimsisions... 
    if "pin_ccd_tilt_x" in metadata_dict:                       # this is SAXS
        sampleI0        = Sample["RawData"]["metadata"]["I0_cts"]
        sampleI0gain    = Sample["RawData"]["metadata"]["I0_gain"]        
        blankI0         = Blank["RawData"]["metadata"]["I0_cts"]
        blankI0gain     = Blank["RawData"]["metadata"]["I0_gain"]
        sampleTRDiode     = Sample["RawData"]["metadata"]["Pin_TrPD"]
        sampleTRDiodeGain = Sample["RawData"]["metadata"]["Pin_TrPDgain"]
        blankTRDiode      = Blank["RawData"]["metadata"]["Pin_TrPD"]
        blankTRDiodeGain  = Blank["RawData"]["metadata"]["Pin_TrPDgain"]
        sampleTRI0     = Sample["RawData"]["metadata"]["Pin_TrI0"]
        sampleTRI0gain = Sample["RawData"]["metadata"]["Pin_TrI0gain"]
        blankTRI0      = Blank["RawData"]["metadata"]["Pin_TrI0"]
        blankTRI0gain  = Blank["RawData"]["metadata"]["Pin_TrI0gain"]
    else:                                                       # and this is WAXS
        sampleI0        = Sample["RawData"]["control"]["integral"]
        sampleI0gain    = Sample["RawData"]["metadata"]["I0_gain"]        
        blankI0         = Blank["RawData"]["control"]["integral"]
        blankI0gain     = Blank["RawData"]["metadata"]["I0_gain"]
        sampleTRDiode     = Sample["RawData"]["metadata"]["TR_cts"]
        sampleTRDiodeGain = Sample["RawData"]["metadata"]["TR_gain"]
        blankTRDiode      = Blank["RawData"]["metadata"]["TR_cts"]
        blankTRDiodeGain  = Blank["RawData"]["metadata"]["TR_gain"]
        sampleTRI0        = sampleI0
        sampleTRI0gain    = sampleI0gain
        blankTRI0         = blankI0
        blankTRI0gain     = blankI0gain
 
    detector_distance = Sample["RawData"]["instrument"]["detector"]["distance"] 
    pixel_size = Sample["RawData"]["instrument"]["detector"]["x_pixel_size"]

    transmission = ((sampleTRDiode / sampleTRDiodeGain) / (sampleTRI0 / sampleTRI0gain)) / ((blankTRDiode / blankTRDiodeGain) / (blankTRI0 / blankTRI0gain))
    #print(f"Transmission: {transmission}")
    I0s = sampleI0 / sampleI0gain
    I0b = blankI0 / blankI0gain
    #nika also divides by this as solid angle correction:
    #			variable solidAngle = PixelSizeX / SampleToCCDDistance * PixelSizeY / SampleToCCDDistance
    solidAngle = pixel_size**2 / detector_distance**2

    preFactor = corrFactor /I0s/(sampleThickness*0.1)/solidAngle          #includes mm to cm conversion
    #print(f"Sample Thickness: {sampleThickness}, CorrFactor: {corrFactor}, Sample I0: {I0s}, Blank I0: {I0b}")
    calib2Ddata =preFactor*((sample2Ddata/transmission) - (I0s/I0b)*blank2Ddata)
    #Int = Corrfactor / (sampleI0 / sampleI0gain) / SampleThickness * (Sa2D/Transm * -  I0/I0Blank * Blank2D)
    #Wreturn the calibrated data, Blank name and may be some parameters? 
    result = {"data":calib2Ddata,
            "blankname":blankname,
            "transmission":transmission
            }
    return result
    

def reduceADData(Sample, useRawData=True):
        '''
        Here we take 2D data from Sample and reduce them to 1D 
        These 2D data in  Sample["RawData"]["data"] can be raw as in read only or normalized or even subtracted and calibrated. 
        '''
        if useRawData:
            my2DData = Sample["RawData"]["data"]
            my2DRAWdata = Sample["RawData"]["data"]
            blankname =  ""
        else:
            my2DData = Sample["calib2DData"]["data"] 
            my2DRAWdata = Sample["RawData"]["data"]
            blankname =  Sample["calib2DData"]["blankname"]

        samplename = Sample["RawData"]["samplename"]
        metadata_dict = Sample["RawData"]["metadata"]
        instrument_dict = Sample["RawData"]["instrument"]
        #extract numbers needed to reduce the data here. 
        # wavelength, keep in A for Fit2D
        wavelength = instrument_dict["monochromator"]["wavelength"]
        # pixel_size, keep in mm, converted in convert_Nika_to_Fit2D to micron for Fit2D and then to m for pyFAI... 
        pixel_size1 = instrument_dict["detector"]["x_pixel_size"] #in mm in NIka, will convert to micron for Fit2D later
        # assume pixels are square, therefore size2 is not needed. No idea how to fix this in Nika or pyFAI for that matter anyway. 
        #pixel_size2 = instrument_dict["detector"]["y_pixel_size"] #in mm in NIka, will convert to micron for Fit2D later
        # detector_distance, keep in mm for Fit2D
        detector_distance = instrument_dict["detector"]["distance"] #in mm in Nika, in mm in Fit2D
        #logging.info(f"Read metadata")
        if "pin_ccd_tilt_x" in metadata_dict:                       # this is SAXS
            usingWAXS=0
            BCX= instrument_dict["detector"]["beam_center_x"]       #  This will be swapped later in convert_Nika_to_Fit2D 
            BCY = instrument_dict["detector"]["beam_center_y"]      #  This will be swapped later in convert_Nika_to_Fit2D
            HorTilt = metadata_dict["pin_ccd_tilt_x"]               #   keep in degrees for Fit2D
            VertTilt = metadata_dict["pin_ccd_tilt_y"]              #   keep in degrees for Fit2D
        else:                                                       # and this is WAXS
            usingWAXS=1
            BCX = instrument_dict["detector"]["beam_center_x"]      #  This will be swapped later in convert_Nika_to_Fit2D
            BCY = instrument_dict["detector"]["beam_center_y"]      #  This will be swapped later in convert_Nika_to_Fit2D
            HorTilt = metadata_dict["waxs_ccd_tilt_x"]              #   keep in degrees for Fit2D
            VertTilt = metadata_dict["waxs_ccd_tilt_y"]             #   keep in degrees for Fit2D    

        #logging.info(f"Finished reading metadata")
    
        # poni is geometry file for pyFAI, created by converting first to Fit2D and then calling pyFAI conversion function.
        my_poni = convert_Nika_to_Fit2D(SSD=detector_distance, pix_size=pixel_size1, BCX=BCX, BCY=BCY, HorTilt=HorTilt, VertTilt=VertTilt, wavelength=wavelength)
        #create mask here. Duplicate the my2DData and set all values above 1e7 to NaN for WAXS or for SAXS mask all negative intensities
        # the differecne is due to Pilatus vs Eiger handing bad pixels differently. Dectris issue... 
        if usingWAXS:
            mask = np.copy(my2DRAWdata)
            mask = 0*mask   # set all values to zero
            mask[my2DRAWdata > 1e7] = 1
            mask[:, 511:516] = 1
            mask[:, 1026:1041] = 1
            mask[:, 1551:1556] = 1
        else:
            mask = np.copy(my2DRAWdata)
            mask = 0*mask   # set all values to zero
            mask[my2DRAWdata < 0] = 1
            # Set the first 4 rows to 1
            mask[:, :4] = 1
            # Set rows 192 to 195 to 1
            mask[:, 242:245] = 1

        #logging.info(f"Finished creating mask")
        
        #now define integrator... using pyFAI here. 
        # this does not work, they really do not have way to pass whole poni in? 
        #ai = AzimuthalIntegrator(poni=my_poni)
        # but this works fine
        ai = AzimuthalIntegrator(dist=my_poni.dist, poni1=my_poni.poni1, poni2=my_poni.poni2, rot1=my_poni.rot1, rot2=my_poni.rot2,
                            rot3=my_poni.rot3, pixel1=my_poni.detector.pixel1, pixel2=my_poni.detector.pixel2, 
                            wavelength=my_poni.wavelength)
        
        #   You can specify the number of bins for the integration
        #   set npt to larger of dimmension of my2DData
        if usingWAXS:
            npt = max(my2DData.shape)
        else:
            npt=200 
        #npt = 1000  # Number of bins, if should be lower
        # Perform azimuthal integration
        q, intensity, sigma = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, error_model="poisson", unit="q_A^-1")
        # fake q resolution as distrance between the subsequent q poiints
        dQ = np.zeros_like(q)
        dQ[1:-1] = 0.5 * (q[2:] - q[:-2])
        dQ[0] = q[1] - q[0]
        dQ[-1] = q[-1] - q[-2]
        #logging.info(f"Finished 2d to 1D conversion")
        result = dict()
        result["Q"] = q
        result["dQ"] = dQ
        result["Intensity"] = intensity
        result["Error"] = sigma
        result["samplename"]=samplename
        result["blankname"]=blankname
        result["units"]="[cm2/cm3]"
        #and these are for compatiblity with USAXS
        result["Kfactor"]=None
        result["OmegaFactor"]=None
        result["thickness"] = Sample["RawData"]["sample"]["thickness"] if "thickness" in Sample["RawData"]["sample"] else None
        return result

# def reduceADToQR(path, filename):
#         tempFilename= os.path.splitext(filename)[0]
#         tempSample = {"RawData":{"filename": tempFilename}}
#         # label = data_dict["RawData"]["filename"]
#         # Q_array = data_dict["reducedData"]["Q_array"]
#         # Intensity = data_dict["reducedData"]["PD_intensity"]
#         tempSample["reducedData"]=ImportAndReduceAD(path, filename)
#         #pp.pprint(tempSample)
#         #pp.pprint(tempSample["RawData"]["filename"])
#         return tempSample


def PlotResults(data_dict):
    # result = {"Int_raw":np.ravel(intensity), 
    #           "Q_raw":np.ravel(q),
    #           "Error_raw":np.ravel(error),
    #           "Intensity":np.ravel(intcalib),
    #           "Q":np.ravel(qcalib),
    #           "Error":np.ravel(errcalib),
    #           }  
    Q_red = data_dict["reducedData"]["Q"]
    Int_red = data_dict["reducedData"]["Intensity"]
    Q = data_dict["CalibratedData"]["Q"]
    Intensity = data_dict["CalibratedData"]["Intensity"]    # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q_red, Int_red, linestyle='-')  # You can customize the marker and linestyle
    plt.plot(Q, Intensity, linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of Intensity vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('Intensity')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.show()

    # # Plot ydata against xdata
    # plt.figure(figsize=(6, 12))
    # plt.plot(Q_array, Intensity, linestyle='-')  # You can customize the marker and linestyle
    # plt.title('Plot of Intensity vs. Q')
    # plt.xlabel('log(Q) [1/A]')
    # plt.ylabel('Intensity')
    # plt.xscale('linear')
    # plt.yscale('linear')
    # plt.grid(True)
    # plt.show()


if __name__ == "__main__":
    Sample = dict()
    ##Sample=process2Ddata("./TestData/TestSet/02_21_Megan_waxs","PU_25C_2_0063.hdf")
    #PlotResults(Sample)
    #Sample["reducedData"]=test("/home/parallels/Github/Matilda/TestData","LaB6_45deg.tif")
    #pp.pprint(Sample)
    #PlotResults(Sample)



                ## test for tilts using LaB6 45 deg tilted detector from GSAXS-II goes here
                # to the best of my undestanding, the images loaded from tiff file are mirrored and the values here are just weird. 
                # def test(path, filename):
                #     # read data from tiff file and read the data 
                #     # tiff files are actually loaded differently than HDF5 files. Looks like they are mirrored. 
                #     my2DData = tiff.imread(path+'/'+filename)
                #     wavelength = 0.10798 # in A
                #     # pixel_size
                #     pixel_size1 = 0.1 # x in Nika, in mm
                #     #pixel_size2 = 0.1 # y in Nika, in mm
                #     # detector_distance, in mm
                #     detector_distance = 1004.91 # in Nika, in mm 
                #     # Nika BCX and BCY in pixels
                #     BCY = 886.7     # this is for hdf5 x in Nika
                #     BCX = 1048.21   # this is for hdf5 y in Nika
                #     # read Nika HorTilt and VertTilt 
                #     VertTilt  = -44.7   # this is negative value for horizontal tilt in Nika
                #     HorTilt = 0.02      # this is value for vertical tilt in Nika, not sure if this shoudl be negative. 
                #     # poni is geometry file for pyFAI, created by converting first to Fit2D and then calling pyFAI conversion function.
                #     my_poni = convert_Nika_to_Fit2D(detector_distance, pixel_size1, BCX, BCY, HorTilt, VertTilt, wavelength)
                #     # setup integrator geometry
                #     ai = AzimuthalIntegrator(dist=my_poni.dist, poni1=my_poni.poni1, poni2=my_poni.poni2, rot1=my_poni.rot1, rot2=my_poni.rot2,
                #                        rot3=my_poni.rot3, pixel1=my_poni.detector.pixel1, pixel2=my_poni.detector.pixel2, 
                #                        wavelength=my_poni.wavelength)
                #     #create mask here. Duplicate the my2DData and set all values to be masked to NaN, not used here. 
                #     mask = np.copy(my2DData)
                #     mask = 0*mask           # set all values to zero
                #     # Perform azimuthal integration
                #     # You can specify the number of bins for the integration
                #     #set npt to larger of dimmension of my2DData  `
                #     npt = max(my2DData.shape)
                #     q, intensity = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, unit="q_A^-1")
                #     result = {"Intensity":np.ravel(intensity), "Q_array":np.ravel(q)}
                #     return result