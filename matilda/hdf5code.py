"""
hdf5code.py
===========
HDF5 / NXcanSAS read-write helpers used throughout Matilda.

Public functions
----------------
save_dict_to_hdf5(data_dict, location, hdf_file)
    Recursively write a Python dict into an open HDF5 file at the given path.

load_dict_from_hdf5(hdf_file, location)
    Recursively read an HDF5 group back into a Python dict.

saveNXcanSAS(Sample, path, filename)
    Write processed I(Q) data as an NXcanSAS-compliant NXsubentry into the
    original scan HDF5 file.

readMyNXcanSAS(path, filename)
    Read back the NXcanSAS entry written by saveNXcanSAS().

readGenericNXcanSAS(path, filename)
    Read any NXcanSAS entry from a Nexus file (not necessarily written by
    Matilda).

find_matching_groups(hdf_file, required_attributes, required_items)
    Search an open HDF5 file for groups matching a set of attribute and
    dataset criteria; returns a list of matching group paths.

Notes
-----
* The 'six' library is imported (line 8) as a Python 2/3 compatibility shim.
  It appears unused in Python 3 code — the import can likely be removed once
  verified.  The original author flagged this with '#what is this for???'.
* All functions expect the HDF5 file to follow the NeXus / NXcanSAS
  conventions used by the APS USAXS Bluesky acquisition system.
"""
import h5py
import os
import numpy as np
import datetime
import logging
from importlib.metadata import version as _pkg_version, PackageNotFoundError as _PkgNotFoundError


def readGenericNXcanSAS(path, filename):
    """
    read data from NXcanSAS data in Nexus file. Ignore NXsas data and anything else
    """
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as f:
        # Start at the root
        # Find the NXcanSAS entries 
        # rootgroup=f['/']
        # SASentries=  find_NXcanSAS_entries(rootgroup)
        required_attributes = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
        required_items = {'definition': 'NXcanSAS'}
        SASentries =  find_matching_groups(f, required_attributes, required_items)
        #print(f"Found {len(SASentries)} NXcanSAS entries in the file:")
        #print(SASentries)
        FirstEntry = SASentries[0] if SASentries else None
        if FirstEntry is None:
            logging.warning(f"No NXcanSAS entries found in the file {filename}.")
            return None

        current_location = FirstEntry
        default_location = f[current_location].attrs.get('default')
        if default_location is not None:
            current_location = f"{current_location}/{default_location}".strip('/')
            if current_location in f:
                default_location = f[current_location].attrs.get('default')
                if 'default' in f[current_location].attrs:
                    current_location = f"{current_location}/{default_location}".strip('/')

        logging.debug(f"Data is located at: {current_location}")
        group_or_dataset = f[current_location]
        # Retrieve and print the list of attributes
        attributes = group_or_dataset.attrs
        logging.debug(f"Attributes at '{current_location}':")
        for attr_name, attr_value in attributes.items():
            logging.debug(f"{attr_name}: {attr_value}")

        data_location= current_location+'/'+attributes['signal']
        if data_location in f:
            # Access the dataset at the specified location
            dataset = f[data_location]
            # Read the data into a NumPy array
            intensity = dataset[()] 
            # Retrieve and print the list of attributes
            Int_attributes = dataset.attrs
            units=Int_attributes['units']
            Kfactor = Int_attributes["Kfactor"]
            OmegaFactor = Int_attributes["OmegaFactor"]
            blankname = Int_attributes["blankname"]
            thickness = Int_attributes["thickness"]
            label = Int_attributes["label"]

        data_location= current_location+'/'+attributes['I_axes']
        if data_location in f:
            # Access the dataset at the specified location
            dataset = f[data_location]
            # Read the data into a NumPy array
            Q = dataset[()] 
            # Retrieve and print the list of attributes
            Q_attributes = dataset.attrs
            #for attr_name, attr_value in Q_attributes.items():
            #    print(f"{attr_name}: {attr_value}")

        data_location= current_location+'/'+Int_attributes['uncertainties']
        if data_location in f:
            # Access the dataset at the specified location
            dataset = f[data_location]
            # Read the data into a NumPy array
            Error = dataset[()] 
            # Retrieve and print the list of attributes
            Error_attributes = dataset.attrs


        data_location= current_location+'/'+Q_attributes['resolutions']
        if data_location in f:
            # Access the dataset at the specified location
            dataset = f[data_location]
            # Read the data into a NumPy array
            dQ = dataset[()] 
            # Retrieve and print the list of attributes
            dQ_attributes = dataset.attrs
        Data = {
            'Intensity':intensity,
            'Q':Q,
            'dQ':dQ,
            'Error':Error,
            'units':units,
            'Int_attributes':Int_attributes,
            'Q_attributes':Q_attributes,
            'Error_Attributes':Error_attributes,
            'dQ_Attributes':dQ_attributes,
            "Kfactor":Kfactor,
            "OmegaFactor":OmegaFactor,
            "blankname":blankname,
            "thickness":thickness,
            'label':label,
        }

        # Optionally read SAStransmission_spectrum if present
        trans_path = FirstEntry + '/sastransmission_spectrum/'
        if trans_path in f:
            T_val = _get_h5_value(f, trans_path + 'T')
            lambda_val = _get_h5_value(f, trans_path + 'lambda')
            if T_val is not None:
                Data['transmission'] = float(T_val[0]) if hasattr(T_val, '__len__') else float(T_val)
            if lambda_val is not None:
                Data['wavelength'] = float(lambda_val[0]) if hasattr(lambda_val, '__len__') else float(lambda_val)

        return Data


# ---------------------------------------------------------------------------
# NXcanSAS metadata helper writers
# ---------------------------------------------------------------------------

def _write_SAStransmission_spectrum(nxDataEntry, transmission, wavelength, name='sample'):
    """Write an NXcanSAS SAStransmission_spectrum group.

    Parameters
    ----------
    nxDataEntry : h5py.Group
        The NXsubentry group to write into.
    transmission : float or None
        Scalar transmission value (T = I/I0). Skipped if None.
    wavelength : float or None
        Wavelength in Angstroms. Skipped if None.
    name : str
        'sample' or 'can' per the NXcanSAS specification.
    """
    if transmission is None or wavelength is None:
        return

    grp_name = 'sastransmission_spectrum'
    if grp_name in nxDataEntry:
        del nxDataEntry[grp_name]

    nxtrans = nxDataEntry.create_group(grp_name)
    nxtrans.attrs['NX_class'] = 'NXdata'
    nxtrans.attrs['canSAS_class'] = 'SAStransmission_spectrum'
    nxtrans.attrs['signal'] = 'T'
    nxtrans.attrs['T_axes'] = 'lambda'
    nxtrans.attrs['name'] = name

    ds = nxtrans.create_dataset('T', data=np.array([float(transmission)]))
    ds.attrs['units'] = 'dimensionless'
    ds.attrs['long_name'] = 'Transmission'

    ds = nxtrans.create_dataset('Tdev', data=np.array([0.0]))
    ds.attrs['units'] = 'dimensionless'
    ds.attrs['long_name'] = 'Transmission uncertainty'

    ds = nxtrans.create_dataset('lambda', data=np.array([float(wavelength)]))
    ds.attrs['units'] = 'angstrom'
    ds.attrs['long_name'] = 'Wavelength'


def _write_SASsample(nxDataEntry, samplename, thickness):
    """Write an NXcanSAS SASsample group.

    Parameters
    ----------
    nxDataEntry : h5py.Group
        The NXsubentry group to write into.
    samplename : str
        Sample name.
    thickness : float or None
        Sample thickness in mm.
    """
    grp_name = 'sassample'
    if grp_name in nxDataEntry:
        del nxDataEntry[grp_name]

    nxsample = nxDataEntry.create_group(grp_name)
    nxsample.attrs['NX_class'] = 'NXsample'
    nxsample.attrs['canSAS_class'] = 'SASsample'

    nxsample.create_dataset('name', data=samplename)
    if thickness is not None:
        ds = nxsample.create_dataset('thickness', data=float(thickness))
        ds.attrs['units'] = 'mm'


def _write_SASinstrument(nxDataEntry, wavelength, detector_distance,
                         fwhm=None, beam_center=None, chi_square=None):
    """Write an NXcanSAS SASinstrument group with source and detector info.

    Parameters
    ----------
    nxDataEntry : h5py.Group
        The NXsubentry group to write into.
    wavelength : float or None
        Wavelength in Angstroms.
    detector_distance : float or None
        Sample-to-detector distance in mm.
    fwhm : float or None
        FWHM of rocking curve in degrees (USAXS only).
    beam_center : float or None
        Beam center angle in degrees (USAXS only).
    chi_square : float or None
        Chi-square of the beam center fit (USAXS only).
    """
    if wavelength is None and detector_distance is None:
        return

    grp_name = 'sasinstrument'
    if grp_name in nxDataEntry:
        del nxDataEntry[grp_name]

    nxinstr = nxDataEntry.create_group(grp_name)
    nxinstr.attrs['NX_class'] = 'NXinstrument'
    nxinstr.attrs['canSAS_class'] = 'SASinstrument'
    nxinstr.create_dataset('name', data='APS 12IDE USAXS/SAXS/WAXS')

    if wavelength is not None:
        nxsource = nxinstr.create_group('sassource')
        nxsource.attrs['NX_class'] = 'NXsource'
        nxsource.attrs['canSAS_class'] = 'SASsource'
        ds = nxsource.create_dataset('wavelength', data=float(wavelength))
        ds.attrs['units'] = 'angstrom'
        nxsource.create_dataset('radiation', data='x-ray synchrotron')

    if detector_distance is not None:
        nxdet = nxinstr.create_group('sasdetector')
        nxdet.attrs['NX_class'] = 'NXdetector'
        nxdet.attrs['canSAS_class'] = 'SASdetector'
        nxdet.create_dataset('name', data='detector')
        ds = nxdet.create_dataset('SDD', data=float(detector_distance))
        ds.attrs['units'] = 'mm'

    if fwhm is not None or beam_center is not None or chi_square is not None:
        nxnote = nxinstr.create_group('sasnote')
        nxnote.attrs['NX_class'] = 'NXnote'
        nxnote.attrs['canSAS_class'] = 'SASnote'
        if fwhm is not None:
            ds = nxnote.create_dataset('FWHM', data=float(fwhm))
            ds.attrs['units'] = 'degrees'
            ds.attrs['long_name'] = 'FWHM of rocking curve'
        if beam_center is not None:
            ds = nxnote.create_dataset('beam_center', data=float(beam_center))
            ds.attrs['units'] = 'degrees'
            ds.attrs['long_name'] = 'Beam center angle'
        if chi_square is not None:
            nxnote.create_dataset('chi_square', data=float(chi_square))


def saveNXcanSAS(Sample,path, filename):
    
    #read stuff from the data dictionary
    Intensity = Sample["CalibratedData"]["Intensity"]
    Q = Sample["CalibratedData"]["Q"]
    Error = Sample["CalibratedData"]["Error"]
    dQ = Sample["CalibratedData"]["dQ"]
    units = Sample["CalibratedData"]["units"]
    Kfactor = Sample["CalibratedData"]["Kfactor"] if "Kfactor" in Sample["CalibratedData"] else None
    OmegaFactor = Sample["CalibratedData"]["OmegaFactor"] if "OmegaFactor" in Sample["CalibratedData"] else None
    blankname = Sample["CalibratedData"]["blankname"] if "blankname" in Sample["CalibratedData"] else None
    thickness = Sample["CalibratedData"]["thickness"] if "thickness" in Sample["CalibratedData"] else None
    label = Sample["RawData"]["filename"]
    timeStamp = Sample["RawData"]["metadata"]["timeStamp"]
    samplename = Sample["RawData"]["sample"]["name"]
    if isinstance(samplename, bytes):
        samplename = samplename.decode('utf-8')


    if "SMR_Int" in Sample["CalibratedData"]:
        SMR_Int =Sample["CalibratedData"]["SMR_Int"]
        SMR_Error =Sample["CalibratedData"]["SMR_Error"]
        SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
        SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]
        slitLength=Sample["CalibratedData"]["slitLength"]
    else:
        SMR_Int = None
        SMR_Error = None
        SMR_Qvec = None
        SMR_dQ = None
        slitLength = None

    R_Int = Sample["reducedData"]["Intensity"]
    R_Qvec = Sample["reducedData"]["Q"]
    R_Error=Sample["reducedData"]["Error"]

    if "BlankData" in Sample:
        BL_R_Int = Sample["BlankData"]["Intensity"]
        BL_Q_vec = Sample["BlankData"]["Q"]
        BL_Error = Sample["BlankData"]["Error"]
    else:
        BL_R_Int = None
        BL_Q_vec = None
        BL_Error = None

    # --- Extract metadata for NXcanSAS groups ---
    # Transmission: USAXS stores in CalibratedData, SAXS/WAXS in calib2DData
    transmission = Sample.get("CalibratedData", {}).get("MeasuredTransmission", None)
    if transmission is None:
        transmission = Sample.get("calib2DData", {}).get("transmission", None)

    # Wavelength: USAXS stores in reducedData, SAXS/WAXS in RawData instrument
    wavelength = Sample.get("reducedData", {}).get("wavelength", None)
    if wavelength is None:
        try:
            wavelength = Sample["RawData"]["instrument"]["monochromator"]["wavelength"]
        except (KeyError, TypeError):
            wavelength = Sample.get("RawData", {}).get("metadata", {}).get("wavelength", None)

    # Sample thickness: prefer CalibratedData, fall back to RawData
    sample_thickness = Sample.get("CalibratedData", {}).get("thickness", None)
    if sample_thickness is None:
        sample_thickness = Sample.get("RawData", {}).get("sample", {}).get("thickness", None)

    # Detector distance: USAXS in metadata, SAXS/WAXS in instrument
    detector_distance = Sample.get("RawData", {}).get("metadata", {}).get("detector_distance", None)
    if detector_distance is None:
        try:
            detector_distance = Sample["RawData"]["instrument"]["detector"]["distance"]
        except (KeyError, TypeError):
            pass

    # USAXS-specific metadata (None for SAXS/WAXS)
    fwhm = Sample.get("reducedData", {}).get("FWHM", None)
    beam_center = Sample.get("reducedData", {}).get("Center", None)
    chi_square = Sample.get("reducedData", {}).get("Chi-Square", None)

    #this is Desmeared USAXS data, SLitSmeared data and plot data, all at once.
    # create the HDF5 NeXus file with same structure as our raw data files have...
    Filepath = os.path.join(path, filename)
    logging.info(f"Saving NXcanSAS data to {Filepath}")
    with h5py.File(Filepath, "a") as f:
        # point to the default data to be plotted
        f.attrs['default']          = 'entry'   #our files have one entry input.
        # these are hopefully optional and useful. 
        f.attrs['file_name']        = filename
        f.attrs['file_time']        = timeStamp 
        f.attrs['instrument']       = '12IDE USAXS'
        f.attrs['creator']          = 'Matilda NeXus writer'
        try:
            _matilda_ver = _pkg_version('Matilda')
        except _PkgNotFoundError:
            _matilda_ver = 'unknown'
        f.attrs['Matilda_version']  = _matilda_ver
        f.attrs['NeXus_version']    = '4.3.0' #2025-5-9 4.3.0 is rc, it is current. 
        f.attrs['HDF5_version']     = h5py.version.hdf5_version
        f.attrs['h5py_version']     = h5py.version.version

        # now create the NXentry group called entry if does not exist
        if 'entry' not in f:
            nxentry = f.create_group('entry')    
        
        nxentry = f['entry']
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['canSAS_class'] = 'SASentry'
        nxentry.attrs['default']  = samplename   #modify with the most reduced data.
        
        #add definition as NXsas - this is location of raw AND reduced data
        # Check if 'definition' dataset exists in the entry group and delete it if present
        if 'definition' in nxentry:
            del nxentry['definition']
        nxentry.create_dataset('definition', data='NXsas')
        # other groups should be here from RAW data, so ignore. 

        if Intensity is not None:
            logging.info(f"Wrote Desmeared NXcanSAS group for file {filename}. ")
            # create the NXsubentry group for Desmeared reduced data. 
            newDataPath = "entry/"+samplename
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]
            
            nxDataEntry = f.create_group(newDataPath)
            nxDataEntry.attrs['NX_class'] = 'NXsubentry'
            nxDataEntry.attrs['canSAS_class'] = 'SASentry'
            nxDataEntry.attrs['default'] = 'sasdata'
            nxDataEntry.attrs['title'] = samplename
            #add definition as NXcanSas
            nxDataEntry.create_dataset('definition', data='NXcanSAS')
            #add title as NXcanSas
            nxDataEntry.create_dataset('title', data=samplename)
            #add run (compulsory)
            nxDataEntry.create_dataset('run', data="run_identifier")

            # create the NXdata group for I(Q) for the avergaed data
            nxdata = nxDataEntry.create_group('sasdata')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'Q'      # X axis of default plot
            #nxdata.attrs['Q_indices'] = [1]    # TODO not sure what this means

            # Y axis data
            ds = nxdata.create_dataset('I', data=Intensity)
            ds.attrs['units'] = '1/cm'
            ds.attrs['uncertainties'] = 'Idev'
            ds.attrs['long_name'] = 'cm2/cm3'    # suggested X axis plot label
            ds.attrs['blankname'] = blankname
            ds.attrs['thickness'] = thickness
            ds.attrs['label'] = label
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            if Kfactor is not None:
                ds.attrs['Kfactor'] = Kfactor
            if OmegaFactor is not None:
                ds.attrs['OmegaFactor'] = OmegaFactor

            # X axis data
            ds = nxdata.create_dataset('Q', data=Q)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label
            ds.attrs['resolutions'] = 'Qdev'
        
            # d X axis data
            ds = nxdata.create_dataset('Qdev', data=dQ)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'   
            # dI axis data
            ds = nxdata.create_dataset('Idev', data=Error)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Uncertainties'

            # NXcanSAS metadata groups for desmeared entry
            _write_SAStransmission_spectrum(nxDataEntry, transmission, wavelength)
            _write_SASsample(nxDataEntry, samplename, sample_thickness)
            _write_SASinstrument(nxDataEntry, wavelength, detector_distance,
                                fwhm=fwhm, beam_center=beam_center, chi_square=chi_square)

        if SMR_Int is not None:
            logging.info(f"Wrote SMR NXcanSAS group for file {filename}. ")
            # add the SMR data
            # create the NXsubentry group for Desmeared reduced data. 
            newDataPath = "entry/"+samplename+"_SMR"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            nxDataEntry.attrs['NX_class'] = 'NXsubentry'
            nxDataEntry.attrs['canSAS_class'] = 'SASentry'
            nxDataEntry.attrs['default'] = 'sasdata'
            nxDataEntry.attrs['title'] = samplename
            #add definition as NXcanSas
            nxDataEntry.create_dataset('definition', data='NXcanSAS')
            #add title as NXcanSas
            nxDataEntry.create_dataset('title', data=samplename)
            #add run (compulsory)
            nxDataEntry.create_dataset('run', data="run_identifier")

            # create the NXdata group for I(Q) for the avergaed data
            nxdata = nxDataEntry.create_group('sasdata')
            nxdata.attrs['NX_class'] = 'NXdata'
            nxdata.attrs['canSAS_class'] = 'SASdata'
            nxdata.attrs['signal'] = 'I'      # Y axis of default plot
            nxdata.attrs['I_axes'] = 'Q'      # X axis of default plot
            #nxdata.attrs['Q_indices'] = [1]    # TODO not sure what this means

            # Y axis data
            ds = nxdata.create_dataset('I', data=SMR_Int)
            ds.attrs['units'] = '1/cm'
            ds.attrs['uncertainties'] = 'Idev'
            ds.attrs['long_name'] = 'Intensity[cm2/cm3]'    # suggested X axis plot label
            ds.attrs['Kfactor'] = Kfactor
            ds.attrs['OmegaFactor'] = OmegaFactor
            ds.attrs['blankname'] = blankname
            ds.attrs['thickness'] = thickness
            ds.attrs['label'] = label

            # X axis data
            ds = nxdata.create_dataset('Q', data=SMR_Qvec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label
            ds.attrs['resolutions'] = 'dQw,dQl'
        
            # d X axis data
            ds = nxdata.create_dataset('dQw', data=SMR_dQ)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'dQw (A^-1)'           
            # slitlength
            ds = nxdata.create_dataset('dQl', data=slitLength)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'dQl (A^-1)'   
            # dI axis data
            ds = nxdata.create_dataset('Idev', data=SMR_Error)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Uncertainties'

            # NXcanSAS metadata groups for SMR entry
            _write_SAStransmission_spectrum(nxDataEntry, transmission, wavelength)
            _write_SASsample(nxDataEntry, samplename, sample_thickness)
            _write_SASinstrument(nxDataEntry, wavelength, detector_distance,
                                fwhm=fwhm, beam_center=beam_center, chi_square=chi_square)

        if R_Int is not None:
            logging.info(f"Wrote QRS group for file {filename}. ")
            newDataPath = "entry/"+"QRS_data"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            # R_Int axis data
            ds = nxDataEntry.create_dataset('Intensity', data=R_Int)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            # R_Qvec axis data
            ds = nxDataEntry.create_dataset('Q', data=R_Qvec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q'    # suggested X axis plot label
            # R_Error axis data
            ds = nxDataEntry.create_dataset('Error', data=R_Error)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Error'    # suggested X axis plot label
            
        if BL_R_Int is not None:
            logging.info(f"Wrote Blank data group for file {filename}. ")
            newDataPath = "entry/"+"Blank_data"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            # R_Int axis data
            ds = nxDataEntry.create_dataset('Intensity', data=BL_R_Int)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            ds.attrs['blankname']=blankname 
            # R_Qvec axis data
            ds = nxDataEntry.create_dataset('Q', data=BL_Q_vec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q'    # suggested X axis plot label
            # R_Error axis data
            ds = nxDataEntry.create_dataset('Error', data=BL_Error)
            ds.attrs['units'] = 'arb'
            ds.attrs['long_name'] = 'Error'    # suggested X axis plot label
 
 

    logging.info(f"Wrote NXcanSAS data to file: {filename}")

def readMyNXcanSAS(path, filename, isUSAXS = False):
    """
    Read My own data from NXcanSAS data in Nexus file.
    
    Parameters:
    path (str): The directory path where the file is located.
    filename (str): The name of the Nexus file to read.

    Returns:
    dict: A dictionary containing the read data.
    """    
    Sample = dict()
    Filepath = os.path.join(path, filename)
    with h5py.File(Filepath, 'r') as f:
        # Start at the root
        # Find the NXcanSAS entries 
        required_attributes = {'canSAS_class': 'SASentry', 'NX_class': 'NXsubentry'}
        required_items = {'definition': 'NXcanSAS'}
        SASentries =  find_matching_groups(f, required_attributes, required_items)
        logging.debug(f"Found {SASentries} entries in the file:{filename}")
              
        location = 'entry/QRS_data/'
        if location in f:
            Sample['reducedData'] = dict()
            dataset = _get_h5_value(f, location + "Intensity")
            if dataset is not None:
                Sample['reducedData']['Intensity'] = dataset
            dataset = _get_h5_value(f, location + "Q")
            if dataset is not None:
                Sample['reducedData']['Q'] = dataset
            dataset = _get_h5_value(f, location + "Error")
            if dataset is not None:
                Sample['reducedData']['Error'] = dataset

        location = 'entry/Blank_data/'
        if location in f:
            Sample['BlankData'] = dict()
            dataset = _get_h5_value(f, location + "Intensity")
            if dataset is not None:
                Sample['BlankData']['Intensity'] = dataset
            dataset = _get_h5_value(f, location + "Q")
            if dataset is not None:
                Sample['BlankData']['Q'] = dataset
            dataset = _get_h5_value(f, location + "Error")
            if dataset is not None:
                Sample['BlankData']['Error'] = dataset
            # BL_R_Int = Sample["BlankData"]["Intensity"]
            # BL_Q_vec = Sample["BlankData"]["Q"]
            # BL_Error = Sample["BlankData"]["Error"]    

        #location = 'entry/'+filename.split('.')[0]+'_SMR/'
        # location is the first of entries from SASentries which contains string _SMR
        location = next((entry + '/' for entry in SASentries if '_SMR' in entry), None)
        logging.debug(f"Found SMR entry at: {location}")
        if 'CalibratedData' not in Sample:
            Sample['CalibratedData'] = dict()
        
        if location is not None and location in f:
            isUSAXS = True      #have SMR data, assume USAXS setup

            dataset = _get_h5_value(f, location + "sasdata/I")
            if dataset is not None:
                Sample['CalibratedData']['SMR_Int'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Q")
            if dataset is not None:
                Sample['CalibratedData']['SMR_Qvec'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Idev")
            if dataset is not None:
                Sample['CalibratedData']['SMR_Error'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/dQw")
            if dataset is not None:
                Sample['CalibratedData']['SMR_dQ'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/dQl")
            if dataset is not None:
                Sample['CalibratedData']['slitLength'] = dataset
        else:
            Sample["CalibratedData"] ["SMR_Qvec"] = None,
            Sample["CalibratedData"] ["SMR_Int"] = None,
            Sample["CalibratedData"] ["SMR_Error"] = None,
            Sample["CalibratedData"] ["SMR_dQ"] = None,
            Sample["CalibratedData"] ["slitLength"] = None,

     
        location = next((entry + '/' for entry in SASentries if '_SMR' not in entry), None)
        logging.debug(f"Found NXcanSAS entry at: {location}")
        if 'RawData' not in Sample:          
            Sample['RawData'] = dict()
            Sample["RawData"]["sample"]=dict()

        if location is not None and location in f:
            dataset = _get_h5_value(f, location + "sasdata/I")
            if dataset is not None:
                Sample['CalibratedData']['Intensity'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Q")
            if dataset is not None:
                Sample['CalibratedData']['Q'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Idev")
            if dataset is not None:
                Sample['CalibratedData']['Error'] = dataset
            dataset = _get_h5_value(f, location + "sasdata/Qdev")
            if dataset is not None:
                Sample['CalibratedData']['dQ'] = dataset
            dataset = _get_h5_value(f, location + "title")
            if dataset is not None:
                Sample["RawData"]["sample"]["name"] = dataset

            location = location+'/sasdata/'
            if "I" in f[location]:
                attributes = f[location + "I"].attrs
                Sample['CalibratedData']['units'] = attributes['units']
                Sample['CalibratedData']['blankname'] = attributes["blankname"]
                Sample['CalibratedData']['thickness'] = attributes["thickness"]
                Sample["RawData"]["filename"] = attributes["label"]
                Sample['CalibratedData']['Kfactor'] = attributes["Kfactor"] if "Kfactor" in attributes else None
                Sample['CalibratedData']['OmegaFactor'] = attributes["OmegaFactor"] if "OmegaFactor" in attributes else None

            # Read SAStransmission_spectrum if present (written by saveNXcanSAS)
            # location was modified above to point to sasdata/, go back to parent
            parent_location = location.rsplit('sasdata/', 1)[0]
            trans_location = parent_location + 'sastransmission_spectrum/'
            if trans_location in f:
                T_val = _get_h5_value(f, trans_location + 'T')
                lambda_val = _get_h5_value(f, trans_location + 'lambda')
                if T_val is not None:
                    Sample['CalibratedData']['transmission'] = float(T_val[0]) if hasattr(T_val, '__len__') else float(T_val)
                if lambda_val is not None:
                    Sample['CalibratedData']['wavelength'] = float(lambda_val[0]) if hasattr(lambda_val, '__len__') else float(lambda_val)

            # Read SASsample if present
            sample_location = parent_location + 'sassample/'
            if sample_location in f:
                thick_val = _get_h5_value(f, sample_location + 'thickness')
                if thick_val is not None:
                    Sample['CalibratedData']['thickness'] = float(thick_val)

            # Read SASinstrument if present
            instr_location = parent_location + 'sasinstrument/'
            if instr_location in f:
                if 'metadata' not in Sample.get('RawData', {}):
                    Sample['RawData']['metadata'] = {}
                source_wl = _get_h5_value(f, instr_location + 'sassource/wavelength')
                if source_wl is not None:
                    Sample['RawData']['metadata']['wavelength'] = float(source_wl)
                sdd_val = _get_h5_value(f, instr_location + 'sasdetector/SDD')
                if sdd_val is not None:
                    Sample['RawData']['metadata']['detector_distance'] = float(sdd_val)

        else:
            Sample["RawData"]["filename"] = filename
            Sample['CalibratedData']['Intensity'] = None
            Sample['CalibratedData']['Q'] = None
            Sample['CalibratedData']['Error'] = None
            Sample['CalibratedData']['Kfactor'] = None
            Sample['CalibratedData']['OmegaFactor'] = None
            Sample['CalibratedData']['blankname'] = None
            Sample['CalibratedData']['thickness'] = None
            Sample['CalibratedData']['units'] = None
            Sample['CalibratedData']['Error'] = None

        #and now we need to read the other groups, which are raw data... 
        if isUSAXS : 
            #metadata
            keys_to_keep = ['AR_center', 'ARenc_0', 'DCM_energy', 'DCM_theta', 'I0Gain','detector_distance',
                            'timeStamp','I0AmpGain',
                            'trans_pin_counts','trans_pin_gain','trans_pin_time','trans_I0_counts','trans_I0_gain',
                            'UPDsize', 'trans_I0_counts', 'trans_I0_gain', 'upd_bkg0', 'upd_bkg1','upd_bkg2','upd_bkg3',
                            'upd_bkgErr0','upd_bkgErr1','upd_bkgErr2','upd_bkgErr3','upd_bkgErr4','upd_bkg_err0',
                            'upd_bkg4','DDPCA300_gain0','DDPCA300_gain1','DDPCA300_gain2','DDPCA300_gain3','DDPCA300_gain4',
                            'SAD_mm', 'SDD_mm', 'thickness', 'title', 'useSBUSAXS',
                            'intervals', 'VToFFactor',
                            'upd_amp_change_mask_time0','upd_amp_change_mask_time1','upd_amp_change_mask_time2','upd_amp_change_mask_time3','upd_amp_change_mask_time4',
                        ]
            # Prefer the "classic" location, but fall back to the Bluesky one.
            metadata_group = None
            for md_path in (
                    "/entry/metadata",
                    "/entry/instrument/bluesky/metadata",
            ):
                if md_path in f:
                    metadata_group = f[md_path]
                    break
            if metadata_group is None:
                raise KeyError(
                    f"Could not find metadata group in {filename}. "
                    "Tried: /entry/metadata and /entry/instrument/bluesky/metadata"
                )

            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            # we need this key to be there also... Copy of the other one.
            # Ensure I0AmpGain exists; default to 1e6 if missing.
            metadata_dict["I0AmpGain"] = metadata_dict.get("I0AmpGain", 1e6)
            metadata_dict["I0Gain"] = metadata_dict["I0AmpGain"]
            #Instrument
            keys_to_keep = ['monochromator', 'energy', 'wavelength']
            instrument_group = f['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            instrument_dict = filter_nested_dict(instrument_dict, keys_to_keep)
            # sample
            sample_group = f['/entry/sample']
            sample_dict = read_group_to_dict(sample_group)

            Sample["RawData"]["metadata"] = metadata_dict
            Sample["RawData"]["instrument"] = instrument_dict
            Sample["RawData"]["sample"].update(sample_dict)
        else:       #this is SWAXS
            #metadata
            instrument_group = f['/entry/instrument']
            instrument_dict = read_group_to_dict(instrument_group)
            #occasionally this fails since 'data' does not exist. 
            # now, why this shoudl tno exists is mystery for me... 
            try:
                del instrument_dict['detector']['data']
            except KeyError:
                pass
            #metadata
            keys_to_keep = ['I000_cts', 'I00_cts', 'I00_gain', 'I0_cts', 'I0_cts_gated',
                            'TR_cts_gated','TR_cts','TR_gain','I0_Sample',
                            'I0_gain', 'I_scaling', 'Pin_TrI0', 'Pin_TrI0gain', 'Pin_TrPD','Pin_TrPDgain',
                            'PresetTime', 'monoE', 'pin_ccd_center_x_pixel','pin_ccd_center_y_pixel',
                            'pin_ccd_tilt_x', 'pin_ccd_tilt_y', 'wavelength', 'waxs_ccd_center_x', 'waxs_ccd_center_y',
                            'waxs_ccd_tilt_x', 'waxs_ccd_tilt_y', 'waxs_ccd_center_x_pixel', 'waxs_ccd_center_y_pixel',
                            'scaler_freq', 'StartTime',                     
                        ]        
            metadata_group = f['/entry/Metadata']
            metadata_dict = read_group_to_dict(metadata_group)
            metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
            sample_group = f['entry/sample']
            sample_dict = read_group_to_dict(sample_group)
            control_group = f['/entry/control']
            control_dict = read_group_to_dict(control_group)
            Sample["RawData"]["instrument"] = instrument_dict
            Sample["RawData"]["metadata"] = metadata_dict
            Sample["RawData"]["sample"].update(sample_dict)
            Sample["RawData"]["control"] = control_dict
 
        return Sample

def _get_h5_value(h5file, path):
    """Return dataset value at path or None if missing."""
    if path in h5file:
        return h5file[path][()]
    return None

def save_dict_to_hdf5(dic, location, h5file):
    """
    Save a dictionary to an HDF5 file.

    Parameters:
    dic (dict): The dictionary to save.
    filename (str): The name of the HDF5 file.
    """
    def recursively_save_dict_contents_to_group(h5file, path, dic):
        for key, item in dic.items():
            if isinstance(item, dict):
                # Create a new group for nested dictionaries
                logging.debug(f"Creating group: {path} + {key}")
                group = h5file.create_group(path + key)
                recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
            else:
                # Save numpy arrays and other data types
                h5file[path + key] = item

    recursively_save_dict_contents_to_group(h5file, location, dic)

# # Example usage
# data_dict = {
#     'array': np.array([1, 2, 3]),
#     'value': 42,
#     'nested': {
#         'string': 'hello',
#         'array2': np.array([4, 5, 6])
#     }
# }

#save_dict_to_hdf5(data_dict, 'data.h5')

def load_dict_from_hdf5(hdf_file, location):
    """
    Load a dictionary from an HDF5 file.

    Parameters:
    filename (str): The name of the HDF5 file.

    Returns:
    dict: The loaded dictionary.

    location (str): The location in the HDF5 file to load from.
    'root:DisplayData/'
    """
    def recursively_load_dict_contents_from_group(h5file, path):
        ans = {}
        for key, item in h5file[path].items():
            if isinstance(item, h5py._hl.group.Group):
                ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
            else:
                tempItem = item[()]
                if isinstance(tempItem, bytes):             # Convert bytes to string
                    tempItem = tempItem.decode('utf-8')
                ans[key] = tempItem
        return ans

    return recursively_load_dict_contents_from_group(hdf_file,location)

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


# def find_NXcanSAS_entries(group, path=''):
#     nxcanSAS_entries = []
    
#     for name, item in group.items():
#         current_path = f"{path}/{name}" if path else name
        
#         # Check if the item is a group
#         if isinstance(item, h5py.Group):
#             # Check if the group has the attribute "NXcanSAS"
#             if 'canSAS_class' in item.attrs:
#                 if(item.attrs['canSAS_class'] == 'SASentry'):
#                     if "definition" in item:
#                         definition_data = item["definition"][()]
#                         # Check if "NXcanSAS" is in the definition data
#                         if isinstance(definition_data, bytes):
#                             definition_data = definition_data.decode('utf-8')
                        
#                         print(f"Definition data: {definition_data}")
#                         if definition_data == 'NXcanSAS':
#                             print(f"Found NXcanSAS entry at: {current_path}")
#                             nxcanSAS_entries.append(current_path)
            
#             # Recursively search within the group
#             nxcanSAS_entries.extend(find_NXcanSAS_entries(item, current_path))
    
#     return nxcanSAS_entries

# this code can find any group which contains listed attributes:values and items:values (strings and variables)
# this is general purpose code for HDF5 - Nexus evaluation
def find_matching_groups(hdf5_file, required_attributes, required_items):
    def check_group(name, obj):
        if isinstance(obj, h5py.Group):
            # Check attributes
            attributes_match = all(
                attr in obj.attrs and obj.attrs[attr] == value
                for attr, value in required_attributes.items()
            )
            
            # Check items
            items_match = True
            for item, expected_value in required_items.items():
                if item in obj:
                    actual_value = obj[item][()]
                    # Decode byte strings to regular strings if necessary
                    if isinstance(actual_value, bytes):
                        actual_value = actual_value.decode('utf-8')
                    if actual_value != expected_value:
                        items_match = False
                        break
                else:
                    items_match = False
                    break
            
            if attributes_match and items_match:
                matching_group_paths.append(name)

    matching_group_paths = []

    hdf5_file.visititems(check_group)

    return matching_group_paths
