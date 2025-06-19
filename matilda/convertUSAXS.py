''' 

    Import data from flyscan or step scan hdf5 file
    and convert to QRdata

    First check if group root:DisplayData exists, if not, load data, process and save for future use. 
    If it exists, load data from it.
    This is to save time for future use.

    Main routines: 
    reduceStepScanToQR
    reduceFlyscanToQR
'''

import os
import h5py
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import pprint as pp
from supportFunctions import read_group_to_dict, filter_nested_dict, check_arrays_same_length
from supportFunctions import gaussian, modifiedGauss
from supportFunctions import beamCenterCorrection, rebinData
from rebinData import rebin_QRSdata
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5



## Flyscan main code here


## Stepscan main code here
def ImportStepScan(path, filename):
    # Open the HDF5 file and read its content, parse content in numpy arrays and dictionaries
    with h5py.File(path+"/"+filename, 'r') as file:
        #read various data sets
        #AR angle
        dataset = file['/entry/data/a_stage_r'] 
        ARangles = np.ravel(np.array(dataset))         
        #time per point
        dataset = file['/entry/data/seconds'] 
        TimePerPoint = np.ravel(np.array(dataset))         
        # I0 gain
        dataset = file['/entry/data/I0_autorange_controls_gain'] 
        I0gain = np.ravel(np.array(dataset)) 
        #I0 - Monitor
        dataset = file['/entry/data/I0'] 
        Monitor = np.ravel(np.array(dataset))  
        #UPD
        dataset = file['/entry/data/UPD'] 
        UPD_array = np.ravel(np.array(dataset))
        #Arrays for gain changes
        dataset = file['/entry/data/upd_autorange_controls_gain'] 
        AmpGain = np.ravel(np.array(dataset))
            #dataset = file['/entry/data/upd_autorange_controls_reqrange'] 
            #AmpReqGain = np.ravel(np.array(dataset))       #this contains only 0 values, useless... 
        #metadata
        keys_to_keep = ['SAD_mm', 'SDD_mm', 'thickness', 'title', 'useSBUSAXS',
                        'intervals', 
                    ]
        metadata_group = file['/entry/instrument/bluesky/metadata']
        metadata_dict = read_group_to_dict(metadata_group)     
        metadata_dict = filter_nested_dict(metadata_dict, keys_to_keep)
        #Instrument
        instrument_group = file['/entry/instrument/monochromator']
        instrument_dict = read_group_to_dict(instrument_group)        
        #Sample
        sample_group = file['/entry/sample']
        sample_dict = read_group_to_dict(sample_group)
        # now backgrounds for UPD subtraction later
        # these are the locations of the background values... 
        # /entry/instrument/bluesky/streams/baseline/I0_autorange_controls_ranges_gain0_background/value, it is array of start adn end values. 
        Bkg0 = file["/entry/instrument/bluesky/streams/baseline/upd_autorange_controls_ranges_gain0_background/value"][0]
        Bkg1 = file["/entry/instrument/bluesky/streams/baseline/upd_autorange_controls_ranges_gain1_background/value"][0]
        Bkg2 = file["/entry/instrument/bluesky/streams/baseline/upd_autorange_controls_ranges_gain2_background/value"][0]
        Bkg3 = file["/entry/instrument/bluesky/streams/baseline/upd_autorange_controls_ranges_gain3_background/value"][0]
        Bkg4 = file["/entry/instrument/bluesky/streams/baseline/upd_autorange_controls_ranges_gain4_background/value"][0]
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
    data_dict = {"Filename": os.path.splitext(filename)[0],
                "ARangles":ARangles, 
                "TimePerPoint": TimePerPoint, 
                "Monitor":Monitor, 
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "I0gain": I0gain,
                "sample": sample_dict,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                "Bkg_map": Bkg_map,
                }
    
    return data_dict
    
def CorrectUPDGainsStep(data_dict):
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
    #TODO: here we also need to subtract background, not done for now anywhere for step scanning. 
    #this neds to be done on UPD_array before we correct for gains.
    # now we need to make a copy of UPD_array and depending on AmpGain value put int BkgX * TimePerPoint
    # Create a new array to store the results
    Bckg_corr = np.zeros_like(AmpGain, dtype=float)
    # Assign background values based on AmpGain and multiply by TimePerPoint
    # Convert keys to floats in Bkg_map for matching with AmpGain
    Bkg_map_float_keys = {float(k): v for k, v in Bkg_map.items()}

    for i, gain in enumerate(AmpGain):
        background_value = Bkg_map_float_keys.get(gain, 0) # Default to 0 if not found
        Bckg_corr[i] = background_value * TimePerPoint[i]/1e7  # Convert to seconds if needed, here we assume TimePerPoint is in microseconds
    # Now we can correct UPD_array for background
    # Remove background from UPD_array
    UPD_array_corr = UPD_array - Bckg_corr       
    #Correct UPD for gains and monitor
    UPD_corrected = (UPD_array_corr*I0gain)/(AmpGain*Monitor)
    result = {"Intensity":UPD_corrected}
    return result


def PlotResults(data_dict):
        # Plot UPD vs Q.
    Q = data_dict["reducedData"]["Q"]
    Intensity = data_dict["reducedData"]["Intensity"]
    
        # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q, Intensity, marker='o', linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of UPD vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('UPD')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.show()


# TODO: remove deleteExisting=True for operations
def reduceStepScanToQR(path, filename, deleteExisting=True):
  # Open the HDF5 file in read/write mode
    location = 'entry/displayData/'
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
        if deleteExisting:
            # Delete the group
            if location in hdf_file:
                del hdf_file[location]
                print("Deleted existing group 'entry/displayData'.")    
        
        if location in hdf_file:
                # # exists, reuse existing data
                Sample = dict()
                Sample = load_dict_from_hdf5(hdf_file, location)
                return Sample
        else:
                Sample = dict()
                Sample["RawData"]=ImportStepScan(path, filename)
                Sample["reducedData"]= CorrectUPDGainsStep(Sample)
                Sample["reducedData"].update(beamCenterCorrection(Sample,useGauss=1))
                # Create the group and dataset for the new data inside the hdf5 file for future use.
                # these are not fully reduced data, this is for web plot purpose.
                save_dict_to_hdf5(Sample, location, hdf_file)
                return Sample



if __name__ == "__main__":
    Sample = dict()
    #Sample = reduceStepScanToQR("/home/parallels/Github/Matilda/TestData","USAXS_step.h5")
    Sample = reduceStepScanToQR(r"C:\Users\ilavsky\Documents\GitHub\Matilda\TestData","USAXS_step.h5")
    #Sample["RawData"]=ImportStepScan("/home/parallels/Github/Matilda","USAXS_step.h5")
        #pp.pprint(Sample)
        #Sample["reducedData"]= CorrectUPDGainsStep(Sample)
        #Sample["reducedData"].update(BeamCenterCorrection(Sample))
        #pp.pprint(Sample["reducedData"])
    #PlotResults(Sample)
    #flyscan
    #Sample = dict()
    #Sample = reduceFlyscanToQR("./TestData","USAXS.h5",deleteExisting=True)
    # Sample["RawData"]=ImportFlyscan("/home/parallels/Github/Matilda","USAXS.h5")
    # #pp.pprint(Sample)
    # Sample["reducedData"]= CorrectUPDGainsFly(Sample)
    # Sample["reducedData"].update(BeamCenterCorrection(Sample))
    # #pp.pprint(Sample["reducedData"])
    PlotResults(Sample)

  
    
    
