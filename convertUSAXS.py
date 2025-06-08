# this file will import data from flyscan or step scan hdf5 file
# it will convert to QR and optionally plot the results


import h5py
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import pprint as pp
from supportFunctions import read_group_to_dict
import os



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


## main code here
def ImportFlyscan(path, filename):
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
        #metadata
        metadata_group = file['/entry/metadata']
        metadata_dict = read_group_to_dict(metadata_group)
        #Instrument
        instrument_group = file['/entry/instrument']
        instrument_dict = read_group_to_dict(instrument_group)


    # Call the function with your arrays
    check_arrays_same_length(ARangles, TimePerPoint, Monitor, UPD_array)
    #Package these results into dictionary
    data_dict = {"Filename": os.path.splitext(filename)[0],
                "ARangles":ARangles, 
                "TimePerPoint": TimePerPoint, 
                "Monitor":Monitor, 
                "UPD_array": UPD_array,
                "AmpGain": AmpGain,
                "Channel": Channel,
                "AmpReqGain": AmpReqGain,
                "metadata": metadata_dict,
                "instrument": instrument_dict,
                }
    
    return data_dict

## main code here
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
        metadata_group = file['/entry/instrument/bluesky/metadata']
        metadata_dict = read_group_to_dict(metadata_group)     
        #Instrument
        instrument_group = file['/entry/instrument/monochromator']
        instrument_dict = read_group_to_dict(instrument_group)        
        #Sample
        sample_group = file['/entry/sample']
        sample_dict = read_group_to_dict(sample_group)


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
                "metadata": instrument_dict,
                "instrument": instrument_dict,
                }
    
    return data_dict
    
def CorrectUPDGainsStep(data_dict):
        # here we will multiply UPD by gain and divide by monitor corrected for its gain.
        # get the needed data from dictionary
    AmpGain = data_dict["RawData"]["AmpGain"]
    UPD_array = data_dict["RawData"]["UPD_array"]
    Monitor = data_dict["RawData"]["Monitor"]
    I0gain = data_dict["RawData"]["I0gain"]
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

        #Correct UPD for gains and monitor
    UPD_corrected = (UPD_array*I0gain)/(AmpGain*Monitor)
    result = {"UPD":UPD_corrected}
    return result


def CorrectUPDGainsFly(data_dict):
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

    
        # Create Gains arrays - one for requested and one for real
    I0AmpGain = metadata_dict["I0AmpGain"]
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
    for i in range(0, len(GainsIndx)-1, 1):
        if np.isnan(GainsIndx[i]):
            continue
        gainName = 'DDPCA300_gain'+str(int(GainsIndx[i]))
        Gains[i] = metadata_dict[gainName]

        #mask amplifier dead times. This is done by comparing table fo deadtimes from metadata with times after range change. 
    Frequency=1e6   #this is frequency of clock fed ito mca1
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

        #Correct UPD for gains and monitor counts and amplfiier gain. 
    UPD_corrected = UPD_array/(Monitor/I0AmpGain)/(Gains)     
    
    result = {"UPD":UPD_corrected}
    return result





def BeamCenterCorrection(data_dict):
    # Find Peak center and create Q vector.
        #RawData=data_dict["rawData"]
        #ReducedData = data_dict["ReducedData"]
    ARangles = data_dict["RawData"]["ARangles"]
    instrument_dict = data_dict["RawData"]["instrument"]
    UPD_array = data_dict["ReducedData"]["UPD"]
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

        # Initial guess for the parameters: amplitude, mean, and standard deviation
    initial_guess = [np.max(ydata_filtered), xdata_filtered[np.argmax(ydata_filtered)], 0.0001]
        #print(initial_guess)

        # Fit the Gaussian function to the filtered data
    popt, _ = curve_fit(gaussian, xdata_filtered, ydata_filtered, p0=initial_guess)

        # Extract the fitted parameters
    amplitude, x0, sigma = popt
        #print(popt)
        # Calculate the FWHM
    fwhm = 2 * np.abs(np.sqrt(2 * np.log(2)) * sigma)

        # Calculate the predicted y values using the fitted parameters
    y_pred = gaussian(xdata_filtered, *popt)

        # Calculate the residuals
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

    results = {"Q_array":Q_array,
            "Chi-Square":chi_square,
            "Center":x0,
            "Maximum":amplitude,
            "FWHM":fwhm    
            }
    return results


def PlotResults(data_dict):
        # Plot UPD vs Q.
    Q_array = data_dict["ReducedData"]["Q_array"]
    UPD = data_dict["ReducedData"]["UPD"]
    
        # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q_array, UPD, marker='o', linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of UPD vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('UPD')
    plt.xscale('log')
    plt.yscale('log')
    plt.grid(True)
    plt.show()


def reduceFlyscanToQR(path, filename):
    # Open the HDF5 file in read/write mode
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
            # Check if the group 'DisplayData' exists
            #if 'DisplayData' in hdf_file:
            #    print("Group 'root/displayData' already exists.")
            #     # Delete the group
            #     #del hdf_file['root/displayData']
            #     #print("Deleted existing group 'root/displayData'.") 
            #     #label = data_dict["RawData"]["Filename"]
            #     #Q_array = data_dict["ReducedData"]["Q_array"]
            #     #UPD = data_dict["ReducedData"]["UPD"]
            #     group = hdf_file['root/displayData']
            #     Sample = {key: group[key][()] for key in group.keys()}
            #     print("Loaded dictionary:")
            #     print(Sample)
            #     return Sample
            # else:
        Sample = dict()
        Sample["RawData"]=ImportFlyscan(path, filename)
            #pp.pprint(Sample)
        Sample["ReducedData"]= CorrectUPDGainsFly(Sample)
        Sample["ReducedData"].update(BeamCenterCorrection(Sample))
            #pp.pprint(Sample["ReducedData"])
            #PlotResults(Sample)
            # Create the group and dataset for the new data
            # newData=dict()
            # newData["RawData"]={"Filename": filename}
            # newData["ReducedData"]={"Q_array":Sample["ReducedData"]["Q_array"],
            #                         "UPD":Sample["ReducedData"]["UPD"],}
            # group = hdf_file.require_group('root')
            # group = hdf_file.require_group('displayData')
            # # Store each key-value pair as a dataset
            # for key, value in newData.items():
            #     # Convert value to a numpy array if it's not already
            #     if not isinstance(value, np.ndarray):
            #         value = np.array(value)
            #     # Create or overwrite the dataset
            #     if key in group:
            #         del group[key]  # Delete existing dataset if it exists
            #     group.create_dataset(key, data=value)

            # print("Appended new data to 'root/displayData'.")
        return Sample


def reduceStepScanToQR(path, filename):
  # Open the HDF5 file in read/write mode
    with h5py.File(path+'/'+filename, 'r+') as hdf_file:
        Sample = dict()
        Sample["RawData"]=ImportStepScan(path, filename)
        #pp.pprint(Sample)
        Sample["ReducedData"]= CorrectUPDGainsStep(Sample)
        Sample["ReducedData"].update(BeamCenterCorrection(Sample))
            #pp.pprint(Sample["ReducedData"])
            #PlotResults(Sample)
        return Sample


if __name__ == "__main__":
    Sample = dict()
    Sample["RawData"]=ImportStepScan("/home/parallels/Github/Matilda","USAXS_step.h5")
    #pp.pprint(Sample)
    Sample["ReducedData"]= CorrectUPDGainsStep(Sample)
    Sample["ReducedData"].update(BeamCenterCorrection(Sample))
    #pp.pprint(Sample["ReducedData"])
    PlotResults(Sample)
    # Sample = dict()
    # Sample["RawData"]=ImportFlyscan("/home/parallels/Github/Matilda","USAXS.h5")
    # #pp.pprint(Sample)
    # Sample["ReducedData"]= CorrectUPDGainsFly(Sample)
    # Sample["ReducedData"].update(BeamCenterCorrection(Sample))
    # #pp.pprint(Sample["ReducedData"])
    # PlotResults(Sample)

  
    
    