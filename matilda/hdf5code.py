'''
    this contains needed hdf5 support for matilda
    used by saving blank BL_QRS data. 
    use code: 
    save_dict_to_hdf5(data_dict, 'data.h5') to save a dictionary to an hdf5 file.
    use code: load_dict_from_hdf5('data.h5') to load a dictionary from an hdf5 file.
    use code: load_dict_from_hdf5('data.h5','root:DisplayData/') to load a dictionary from an hdf5 file.

'''
import h5py
import os
import numpy as np
import six  #what is this for???
import datetime
import logging

def readGenericNXcanSAS(path, filename):
    '''
    read data from NXcanSAS data in Nexus file. Ignore NXsas data and anything else
    '''
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
            print("No NXcanSAS entries found in the file.")
            return None

        current_location = FirstEntry
        default_location = f[current_location].attrs.get('default')
        if default_location is not None:
            current_location = f"{current_location}/{default_location}".strip('/')
            if current_location in f:
                default_location = f[current_location].attrs.get('default')
                if 'default' in f[current_location].attrs:
                    current_location = f"{current_location}/{default_location}".strip('/')

        #print(f"Data is located at: {current_location}")    
        group_or_dataset = f[current_location]
        # Retrieve and print the list of attributes
        attributes = group_or_dataset.attrs
        print(f"Attributes at '{current_location}':")
        for attr_name, attr_value in attributes.items():
            print(f"{attr_name}: {attr_value}")

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
            BlankName = Int_attributes["BlankName"]
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
            "BlankName":BlankName,
            "thickness":thickness,
            'label':label,
        }
        return Data

def saveNXcanSAS(Sample,path, filename):
    
    #read stuff from the data dictionary
    Intensity = Sample["CalibratedData"]["Intensity"]
    Q = Sample["CalibratedData"]["Q"]
    Error = Sample["CalibratedData"]["Error"]
    dQ = Sample["CalibratedData"]["dQ"]
    units = Sample["CalibratedData"]["units"]
    Kfactor = Sample["CalibratedData"]["Kfactor"]
    OmegaFactor = Sample["CalibratedData"]["OmegaFactor"]
    BlankName = Sample["CalibratedData"]["BlankName"]
    thickness = Sample["CalibratedData"]["thickness"]
    label = Sample["RawData"]["Filename"]
    timeStamp = Sample["RawData"]["metadata"]["timeStamp"]
    SampleName = Sample["RawData"]["sample"]["name"]
    SampleName = SampleName.decode('utf-8')

    SMR_Int =Sample["CalibratedData"]["SMR_Int"]
    SMR_Error =Sample["CalibratedData"]["SMR_Error"]
    SMR_Qvec =Sample["CalibratedData"]["SMR_Qvec"]
    SMR_dQ =Sample["CalibratedData"]["SMR_dQ"]
    slitLength=Sample["CalibratedData"]["slitLength"]

    R_Int = Sample["reducedData"]["Intensity"]
    R_Qvec = Sample["reducedData"]["Q"]
    R_Error=Sample["reducedData"]["Error"]

    #this is Desmeared USAXS data, SLitSmeared data and plot data, all at once.
    # create the HDF5 NeXus file with same structure as our raw data files have...
    Filepath = os.path.join(path, filename)
    print(f"Saving NXcanSAS data to {Filepath}")
    with h5py.File(Filepath, "a") as f:
        # point to the default data to be plotted
        f.attrs['default']          = 'entry'   #our files have one entry input.
        # these are hopefullyoptional and useful. 
        f.attrs['file_name']        = filename
        f.attrs['file_time']        = timeStamp 
        f.attrs['instrument']       = '12IDE USAXS'
        f.attrs['creator']          = 'Matilda NeXus writer'
        f.attrs['NeXus_version']    = '4.3.0' #2025-5-9 4.3.0 is rc, it is current. 
        f.attrs['HDF5_version']     = six.u(h5py.version.hdf5_version)
        f.attrs['h5py_version']     = six.u(h5py.version.version)

        # now create the NXentry group called entry if does not exist
        if 'entry' not in f:
            nxentry = f.create_group('entry')    
        
        nxentry = f['entry']
        nxentry.attrs['NX_class'] = 'NXentry'
        nxentry.attrs['default']  = SampleName   #modify with the most reduced data.
        
        #add definition as NXsas - this is location of raw AND reduced data
        # Check if 'definition' dataset exists in the entry group and delete it if present
        if 'definition' in nxentry:
            del nxentry['definition']
        nxentry.create_dataset('definition', data='NXsas')
        # other groups should be here from RAW data, so ignore. 

        if Intensity is not None:
            logging.info(f"Wrote Desmeared NXcanSAS group for file {filename}. ")
            # create the NXsubentry group for Desmeared reduced data. 
            newDataPath = "entry/"+SampleName
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]
            
            nxDataEntry = f.create_group(newDataPath)
            nxDataEntry.attrs['NX_class'] = 'NXsubentry'
            nxDataEntry.attrs['canSAS_class'] = 'SASentry'
            nxDataEntry.attrs['default'] = 'sasdata'
            nxDataEntry.attrs['title'] = SampleName
            #add definition as NXcanSas
            nxDataEntry.create_dataset('definition', data='NXcanSAS')
            #add title as NXcanSas
            nxDataEntry.create_dataset('title', data=SampleName)
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
            ds.attrs['Kfactor'] = Kfactor
            ds.attrs['OmegaFactor'] = OmegaFactor
            ds.attrs['BlankName'] = BlankName
            ds.attrs['thickness'] = thickness
            ds.attrs['label'] = label

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

        if SMR_Int is not None:
            logging.info(f"Wrote SMR NXcanSAS group for file {filename}. ")
            # add the SMR data
            # create the NXsubentry group for Desmeared reduced data. 
            newDataPath = "entry/"+SampleName+"_SMR"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            nxDataEntry.attrs['NX_class'] = 'NXsubentry'
            nxDataEntry.attrs['canSAS_class'] = 'SASentry'
            nxDataEntry.attrs['default'] = 'sasdata'
            nxDataEntry.attrs['title'] = SampleName
            #add definition as NXcanSas
            nxDataEntry.create_dataset('definition', data='NXcanSAS')
            #add title as NXcanSas
            nxDataEntry.create_dataset('title', data=SampleName)
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
            ds.attrs['long_name'] = 'cm2/cm3'    # suggested X axis plot label
            ds.attrs['Kfactor'] = Kfactor
            ds.attrs['OmegaFactor'] = OmegaFactor
            ds.attrs['BlankName'] = BlankName
            ds.attrs['thickness'] = thickness
            ds.attrs['label'] = label

            # X axis data
            ds = nxdata.create_dataset('Q', data=SMR_Qvec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q (A^-1)'    # suggested Y axis plot label
            ds.attrs['resolutions'] = 'dQw, dQI'
        
            # d X axis data
            ds = nxdata.create_dataset('dQw', data=SMR_dQ)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'dQw (A^-1)'           
            # slitlength
            ds = nxdata.create_dataset('dQI', data=slitLength)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'dQI (A^-1)'   
            # dI axis data
            ds = nxdata.create_dataset('Idev', data=SMR_Error)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Uncertainties'  

        if R_Int is not None:
            logging.info(f"Wrote QRS group for file {filename}. ")
            newDataPath = "entry/"+"QRS_data"
            if newDataPath in f:
                logging.warning(f"NXcanSAS group {newDataPath} already exists in file {filename}. Overwriting.")
                del f[newDataPath]

            nxDataEntry = f.create_group(newDataPath)
            # R_Int axis data
            ds = nxDataEntry.create_dataset('Intensity', data=R_Int)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Intensity'    # suggested X axis plot label
            # R_Qvec axis data
            ds = nxDataEntry.create_dataset('Q', data=R_Qvec)
            ds.attrs['units'] = '1/angstrom'
            ds.attrs['long_name'] = 'Q'    # suggested X axis plot label
            # R_Error axis data
            ds = nxDataEntry.create_dataset('Error', data=R_Error)
            ds.attrs['units'] = 'cm2/cm3'
            ds.attrs['long_name'] = 'Error'    # suggested X axis plot label
            
    
    #print("wrote file:", filename)

def readMyNXcanSAS(path, filename):
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
        #print(f"Found {SASentries} entries in the file:{filename}")
              
        location = 'entry/QRS_data/'
        if location in f:
            Sample['reducedData'] = dict()
            dataset = f[location + "Intensity"]
            if dataset is not None:
                Sample['reducedData']['Intensity'] = dataset[()]
            dataset = f[location + "Q"]
            if dataset is not None:
                Sample['reducedData']['Q'] = dataset[()]
            dataset = f[location + "Error"]
            if dataset is not None:
                Sample['reducedData']['Error'] = dataset[()]
     
        #location = 'entry/'+filename.split('.')[0]+'_SMR/'
        # location is the first of entries from SASentries which contains string _SMR
        location = next((entry + '/' for entry in SASentries if '_SMR' in entry), None)
        #print(f"Found SMR entry at: {location}")
        if location is not None and location in f:
            Sample['CalibratedData'] = dict()
            dataset = f[location + "sasdata/I"]
            if dataset is not None:
                Sample['CalibratedData']['SMR_Int'] = dataset[()]
            dataset = f[location + "sasdata/Q"]
            if dataset is not None:
                Sample['CalibratedData']['SMR_Qvec'] = dataset[()]
            dataset = f[location + "sasdata/Idev"]
            if dataset is not None:
                Sample['CalibratedData']['SMR_Error'] = dataset[()]
            dataset = f[location + "sasdata/dQw"]
            if dataset is not None:
                Sample['CalibratedData']['SMR_dQ'] = dataset[()]
            dataset = f[location + "sasdata/dQI"]
            if dataset is not None:
                Sample['CalibratedData']['slitLength'] = dataset[()]
        else:
            Sample['CalibratedData'] = dict()
            Sample["CalibratedData"] ["SMR_Qvec"] = None,
            Sample["CalibratedData"] ["SMR_Int"] = None,
            Sample["CalibratedData"] ["SMR_Error"] = None,
            Sample["CalibratedData"] ["SMR_dQ"] = None,
            Sample["CalibratedData"] ["slitLength"] = None,

     
        location = next((entry + '/' for entry in SASentries if '_SMR' not in entry), None)
        #print(f"Found NXcanSAS entry at: {location}")
        if location is not None and location in f:
            Sample['CalibratedData'] = dict()
            Sample['RawData'] = dict()
            Sample['RawData']['sample'] = dict()
            dataset = f[location + "sasdata/I"]
            if dataset is not None:
                Sample['CalibratedData']['Intensity'] = dataset[()]
            dataset = f[location + "sasdata/Q"]
            if dataset is not None:
                Sample['CalibratedData']['Q'] = dataset[()]
            dataset = f[location + "sasdata/Idev"]
            if dataset is not None:
                Sample['CalibratedData']['Error'] = dataset[()]
            dataset = f[location + "sasdata/Qdev"]
            if dataset is not None:
                Sample['CalibratedData']['dQ'] = dataset[()]            
            dataset = f[location + "title"]
            if dataset is not None:
                Sample["RawData"]["sample"]["name"] = dataset[()]
            attributes = f[location + "sasdata/I"].attrs
            Sample['CalibratedData']['units'] = attributes['units']
            Sample['CalibratedData']['Kfactor'] = attributes["Kfactor"]
            Sample['CalibratedData']['OmegaFactor'] = attributes["OmegaFactor"]
            Sample['CalibratedData']['BlankName'] = attributes["BlankName"]
            Sample['CalibratedData']['thickness'] = attributes["thickness"]
            Sample["RawData"]["Filename"] = attributes["label"]
        else:
            Sample['CalibratedData'] = dict()
            Sample['RawData'] = dict()
            Sample["RawData"]["Filename"] = filename
            Sample["RawData"]["sample"] = dict()
            Sample['CalibratedData']['Intensity'] = None
            Sample['CalibratedData']['Q'] = None
            Sample['CalibratedData']['Error'] = None
            Sample['CalibratedData']['Kfactor'] = None
            Sample['CalibratedData']['OmegaFactor'] = None
            Sample['CalibratedData']['BlankName'] = None
            Sample['CalibratedData']['thickness'] = None
            Sample['CalibratedData']['units'] = None
            Sample['CalibratedData']['Error'] = None

     
        return Sample



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
                print(f"Creating group: {path} + {key}")
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
                ans[key] = item[()]
        return ans

    return recursively_load_dict_contents_from_group(hdf_file,location)


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
