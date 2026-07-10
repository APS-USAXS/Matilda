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
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import h5py
import os
import logging

from .supportFunctions import read_group_to_dict, filter_nested_dict
from .supportNikaFunctions import convert_Nika_to_Fit2D
from .hdf5code import save_dict_to_hdf5, load_dict_from_hdf5, saveNXcanSAS, readMyNXcanSAS, find_matching_groups
from .hdf5code import writeThicknessOverride

# ── Integrator cache ──────────────────────────────────────────────────────────
# pyFAI builds internal lookup tables on the first integrate1d() call for a
# given geometry.  For files in the same folder the geometry is usually
# identical, so reusing the same AzimuthalIntegrator avoids rebuilding those
# tables on every call (3× per file × N files).
_cached_integrator = None
_cached_geometry_key = None


def _get_integrator(my_poni):
    """Return a cached AzimuthalIntegrator if geometry matches, else create one."""
    global _cached_integrator, _cached_geometry_key
    key = (my_poni.dist, my_poni.poni1, my_poni.poni2,
           my_poni.rot1, my_poni.rot2, my_poni.rot3,
           my_poni.detector.pixel1, my_poni.detector.pixel2,
           my_poni.wavelength)
    if _cached_geometry_key == key and _cached_integrator is not None:
        return _cached_integrator
    ai = AzimuthalIntegrator(
        dist=my_poni.dist, poni1=my_poni.poni1, poni2=my_poni.poni2,
        rot1=my_poni.rot1, rot2=my_poni.rot2, rot3=my_poni.rot3,
        pixel1=my_poni.detector.pixel1, pixel2=my_poni.detector.pixel2,
        wavelength=my_poni.wavelength,
    )
    _cached_integrator = ai
    _cached_geometry_key = key
    logging.info("Created new pyFAI AzimuthalIntegrator (geometry changed or first call)")
    return ai


def _geometry_from_dicts(instrument_dict, metadata_dict):
    """Extract Nika geometry from the file dicts and convert to pyFAI PONI.

    Returns (my_poni, usingWAXS).  SAXS files are recognised by the presence
    of ``pin_ccd_tilt_x`` in the metadata.

    Units follow Nika conventions (wavelength A, pixel size mm, distance mm);
    conversion to Fit2D/pyFAI units happens inside convert_Nika_to_Fit2D.
    Pixels are assumed square (y_pixel_size ignored); BCX/BCY are swapped for
    Fit2D inside convert_Nika_to_Fit2D.
    """
    wavelength = instrument_dict["monochromator"]["wavelength"]
    pixel_size1 = instrument_dict["detector"]["x_pixel_size"]
    detector_distance = instrument_dict["detector"]["distance"]
    BCX = instrument_dict["detector"]["beam_center_x"]
    BCY = instrument_dict["detector"]["beam_center_y"]
    if "pin_ccd_tilt_x" in metadata_dict:               # SAXS
        usingWAXS = 0
        HorTilt = metadata_dict["pin_ccd_tilt_x"]       # degrees
        VertTilt = metadata_dict["pin_ccd_tilt_y"]
    else:                                               # WAXS
        usingWAXS = 1
        HorTilt = metadata_dict["waxs_ccd_tilt_x"]
        VertTilt = metadata_dict["waxs_ccd_tilt_y"]
    my_poni = convert_Nika_to_Fit2D(SSD=detector_distance, pix_size=pixel_size1,
                                    BCX=BCX, BCY=BCY, HorTilt=HorTilt,
                                    VertTilt=VertTilt, wavelength=wavelength)
    return my_poni, usingWAXS


def _build_mask(my2DRAWdata, usingWAXS, metadata_dict):
    """Build the bad-pixel/gap mask for azimuthal integration.

    Single source of truth — the inline copies in ImportAndReduceAD and
    reduceADData had diverged before this was factored out.

    Detector geometry has changed over time: WAXS masks branch by detector
    size, SAXS masks by acquisition year (StartTime metadata).  The >1e7 /
    <0 cuts differ because Pilatus and Eiger mark bad pixels differently.
    numpy 2-D shape is (rows, cols); rows is the short axis on these detectors.
    """
    mask = np.zeros_like(my2DRAWdata)
    if usingWAXS:
        if my2DRAWdata.shape[0] >= 256:
            # current WAXS detector: 512 x 2068 (4-tile, gaps at cols 512/1024/1536)
            mask[my2DRAWdata > 1e7] = 1
            mask[:, 511:516] = 1
            mask[:, 1026:1041] = 1
            mask[:, 1551:1556] = 1
        else:
            # old WAXS detector (195 x 981)
            mask[my2DRAWdata > 1e7] = 1
            mask[my2DRAWdata < 0] = 1
            mask[:, 486:494] = 1
            mask[:, 979:980] = 1
            mask[0:3, :] = 1
            mask[193:194, :] = 1
    else:
        # SAXS: pick mask by acquisition year.
        # StartTime format: "2022-12-15 10:48:51.230765" (string from /entry/Metadata).
        start_time_raw = metadata_dict.get("StartTime", "")
        if isinstance(start_time_raw, bytes):
            start_time_raw = start_time_raw.decode("utf-8", errors="replace")
        try:
            acquisition_year = int(str(start_time_raw).strip()[:4])
        except (ValueError, TypeError):
            acquisition_year = 9999   # unknown -> treat as current
        mask[my2DRAWdata < 0] = 1
        mask[my2DRAWdata > 1e7] = 1
        if acquisition_year >= 2023:
            # current SAXS detector
            mask[:, :4] = 1
            mask[:, 242:245] = 1
        else:
            # pre-2023 SAXS detector
            mask[:, 0:7] = 1
            mask[86, 17] = 1
            mask[58, 112] = 1
    return mask


# TODO: split into multiple steps as needed
# Import images for sample and blank as separate calls and get sample and blank objects
#   calibration step, calculate corrections, apply corrections
#   subtract 2D images to get calibrated image as Nika
#   convert to 1D data
# TODO: test me... 



## main code here
def process2Ddata(path, filename, blankPath=None, blankFilename=None, recalculateAllData=False,
                   npts=None, thickness_override=None,
                   use_mu=False, mu=None, per_gram=False, density=None,
                   transmission_override=None):
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
            if thickness_override is not None:
                writeThicknessOverride(hdf_file, '/entry/sample/thickness',
                                       thickness_override, filename)
            Sample = importADData(path, filename)   #this is for sample path and blank
            Sample["reducedData"] = reduceADData(Sample, useRawData=True, npts=npts)   #this generates Int vs Q for raw data plot
                                        # q = Sample["reducedData"]["Q"]
                                        # intensity = Sample["reducedData"]["Intensity"]
                                        # error = Sample["reducedData"]["Error"]
                                        # samplename = Sample["RawData"]["samplename"]

            if blankPath is not None and blankFilename is not None and "blank" not in filename.lower():
                blank = importADData(blankPath, blankFilename)               #this is for blank path and blank name
                Sample["BlankData"] = reduceADData(blank, useRawData=True, npts=npts)   #this generates Int vs Q for blank data plot
                                        # qcalib = Sample["BlankData"]["Q"]
                                        # intensity = Sample["BlankData"]["Intensity"]
                                        # error = Sample["BlankData"]["Error"]
                                        # blankname = Sample["RawData"]["blankname"]
                Sample["calib2DData"] = calibrateAD2DData(Sample, blank, thickness_override=thickness_override, use_mu=use_mu, mu=mu, transmission_override=transmission_override)
                                    #returns 2D calibrated data
                                        # result = {"data":calib2Ddata,
                                        #           "blankname":blankname,
                                        #           "transmission":transmission
                Sample["CalibratedData"] = reduceADData(Sample, useRawData=False, npts=npts, per_gram=per_gram, density=density)  #this generates Calibrated 1D data.
                # Propagate calculated_mu from 2D calibration into CalibratedData
                # so the GUI can read it from a single, consistent place.
                Sample["CalibratedData"]["calculated_mu"] = Sample["calib2DData"].get("calculated_mu")
                                        #returns :   
                                        # qcalib= Sample["CalibratedData"]["Q"]
                                        # dqcalib= Sample["CalibratedData"]["dQ"]
                                        # intcalib= Sample["CalibratedData"]["Intensity"]
                                        # errcalib= Sample["CalibratedData"]["Error"]
                                        # Sample["CalibratedData"]["units"]
                                        # blankname = Sample["calib2DData"]["blankname"]
                Sample = fixOversubtractedData(Sample)
                                        #this simply finds min for Intensity and if it is negatgive, adds to it 1.1*min value. 
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

def fixOversubtractedData(Sample):
    '''
    here we take 1D data Sample["CalibratedData"]["Intensity"] and fix them if they are negative
    '''
    minIntVal = Sample["CalibratedData"]["Intensity"].min()
    if minIntVal < 0:
        Sample["CalibratedData"]["Intensity"] = Sample["CalibratedData"]["Intensity"] - 1.1*minIntVal
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
    with h5py.File(os.path.join(path, filename), 'r+') as hdf_file:
        # Check if the group 'displayData' exists
        if recalculateAllData:
            # Delete the group
            if location in hdf_file:
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
                            'scaler_freq', 'StartTime',      # StartTime needed by _build_mask (SAXS year branch)
                        ]
            metadata_group = hdf_file['/entry/Metadata']
            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            #need to append these structures to RAW data so we can optionally save them in Igor
            Sample["RawData"]=dict()
            Sample["RawData"]["metadata"]= metadata_dict
            Sample["RawData"]["instrument"]= instrument_dict
            Sample["RawData"]["sample"]= sample_dict
        # geometry + mask via the shared helpers (single source of truth)
        my_poni, usingWAXS = _geometry_from_dicts(instrument_dict, metadata_dict)
        mask = _build_mask(my2DData, usingWAXS, metadata_dict)

        ai = _get_integrator(my_poni)

        #   You can specify the number of bins for the integration
        #   set npt to larger of dimmension of my2DData
        #   error_model= “azimuthal” or “poisson” (variance = I), “azimuthal” (variance = (I-<I>)^2)

        if usingWAXS:
            npt = max(my2DData.shape)
            q, intensity, sigma = ai.integrate1d(my2DData, npt, mask=mask, error_model= "azimuthal", correctSolidAngle=True, unit="q_A^-1")
        else:
            npt=200 
            # using azimuth_range=(-30,30) should limit the range of data to what Nika is using for SAXS. 
            q, intensity, sigma = ai.integrate1d(my2DData, npt, mask=mask,azimuth_range=(-30,30), error_model= "azimuthal", correctSolidAngle=True, unit="q_A^-1")

        Sample["reducedData"] = dict()
        Sample["reducedData"]["Q_array"] = q
        Sample["reducedData"]["Intensity"] = intensity
        Sample["reducedData"]["Error"] = sigma
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

def calibrateAD2DData(Sample, Blank, thickness_override=None, use_mu=False, mu=None,
                       transmission_override=None):
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
    sampleThickness = thickness_override if thickness_override is not None else Sample["RawData"]["sample"]["thickness"]
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
    if transmission_override is not None and transmission_override > 0:
        logging.info(f"Using transmission override {transmission_override:.4f} (measured: {transmission:.4f}).")
        transmission = float(transmission_override)
    I0s = sampleI0 / sampleI0gain
    I0b = blankI0 / blankI0gain
    #nika also divides by this as solid angle correction:
    #			variable solidAngle = PixelSizeX / SampleToCCDDistance * PixelSizeY / SampleToCCDDistance
    solidAngle = pixel_size**2 / detector_distance**2

    if use_mu and mu is not None and mu > 0:
        if transmission > 0 and transmission < 1:
            thickness_cm = -np.log(transmission) / mu
        else:
            logging.warning(f"Transmission {transmission:.4f} out of range (0,1) for mu calculation. Falling back to sample thickness.")
            thickness_cm = sampleThickness * 0.1
    else:
        thickness_cm = sampleThickness * 0.1  # mm to cm
    preFactor = corrFactor /I0s/thickness_cm/solidAngle
    #print(f"Sample Thickness: {sampleThickness}, CorrFactor: {corrFactor}, Sample I0: {I0s}, Blank I0: {I0b}")
    calib2Ddata =preFactor*((sample2Ddata/transmission) - (I0s/I0b)*blank2Ddata)
    #Int = Corrfactor / (sampleI0 / sampleI0gain) / SampleThickness * (Sa2D/Transm * -  I0/I0Blank * Blank2D)
    # Compute mu from the (possibly overridden) transmission and the
    # thickness actually used. Useful for display even when not in mu mode.
    if 0 < transmission < 1 and thickness_cm > 0:
        calculated_mu = float(-np.log(transmission) / thickness_cm)
    else:
        calculated_mu = None
    #Wreturn the calibrated data, Blank name and may be some parameters?
    result = {"data":calib2Ddata,
            "blankname":blankname,
            "transmission":transmission,
            "calculated_mu":calculated_mu,
            }
    return result
    

def reduceADData(Sample, useRawData=True, npts=None, per_gram=False, density=None):
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
        # geometry + mask via the shared helpers (single source of truth)
        my_poni, usingWAXS = _geometry_from_dicts(instrument_dict, metadata_dict)
        mask = _build_mask(my2DRAWdata, usingWAXS, metadata_dict)

        ai = _get_integrator(my_poni)

        #   You can specify the number of bins for the integration
        #   set npt to larger of dimmension of my2DData
        if npts is not None:
            npt = npts
        else:
            npt = max(my2DData.shape)
        #npt = 1000  # Number of bins, if should be lower
        # Perform azimuthal integration
        #   error_model= "azimuthal" or “poisson” (variance = I), “azimuthal” (variance = (I-<I>)^2)
  
        q, intensity, sigma = ai.integrate1d(my2DData, npt, mask=mask, correctSolidAngle=True, error_model="azimuthal", unit="q_A^-1")
        # fake q resolution as distrance between the subsequent q points
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
        # Per-gram normalization (calibrated data only)
        if not useRawData and per_gram and density is not None and density > 0:
            result["Intensity"] = result["Intensity"] / density
            result["Error"] = result["Error"] / density
            result["units"] = "[cm2/g]"
        return result

def PlotResults(data_dict):
    # result = {"Int_raw":np.ravel(intensity), 
    #           "Q_raw":np.ravel(q),
    #           "Error_raw":np.ravel(error),
    #           "Intensity":np.ravel(intcalib),
    #           "Q":np.ravel(qcalib),
    #           "Error":np.ravel(errcalib),
    #           }  
    pass    # debug-only function; body is fully commented out below
    # Q_red = data_dict["reducedData"]["Q"]              # used by debug plot below
    # Int_red = data_dict["reducedData"]["Intensity"]
    # Q = data_dict["CalibratedData"]["Q"]
    # Intensity = data_dict["CalibratedData"]["Intensity"]
    # Debug plot (requires matplotlib; commented out for production):
    # import matplotlib.pyplot as plt
    # plt.figure(figsize=(6, 12))
    # plt.plot(Q_red, Int_red, linestyle='-')
    # plt.plot(Q, Intensity, linestyle='-')
    # plt.title('Plot of Intensity vs. Q')
    # plt.xlabel('log(Q) [1/A]')
    # plt.ylabel('Intensity')
    # plt.xscale('log')
    # plt.yscale('log')
    # plt.grid(True)
    # plt.show()

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
    Sample = ImportAndReduceAD("//Mac/Home/Desktop/Test", "R6016ACT_T4_H_1081.hdf", recalculateAllData=True)
    ##Sample=process2Ddata("./TestData/TestSet/02_21_Megan_waxs","PU_25C_2_0063.hdf")
    #PlotResults(Sample)
    #Sample["reducedData"]=test("/home/parallels/Github/Matilda/TestData","LaB6_45deg.tif")
    #pp.pprint(Sample)
    #PlotResults(Sample)



# (historical tilt-test scaffolding for GSAXS-II LaB6 tiff data removed 2026-07;
#  see git history or CodeFragments/ if needed again)