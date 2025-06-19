# these are support functions for Matilda

import h5py
import numpy as np
import logging
import copy
import xarray as xr
from scipy.interpolate import interp1d
## support stuff here


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
    logY2 = np.log(Y2)
    logY2_interp_func = interp1d(X2, logY2, kind='linear', fill_value='extrapolate')
    logY2_interp = logY2_interp_func(X1)

    # Convert interpolated logY2 back to Y2interp
    Y2_interp = np.exp(logY2_interp)

    # Step 2: Subtract Y1 - Y2interp to obtain Ydiff and divide to get IntRatio
    Ydiff = Y1 - Y2_interp

    IntRatio = Y1 / Y2_interp  # Calculate the intensity ratio
    
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
def modifiedGauss(xvar, a, x0, sigma, dparameter):
    return a * np.exp(-(np.abs(xvar - x0) / (2*sigma)) ** dparameter)



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
            data_dict[key] = data
        elif isinstance(item, h5py.Group):
            # If the item is a group, recursively read its contents
            data_dict[key] = read_group_to_dict(item)
    return data_dict

# def filter_nested_dict(d, keys_to_keep):
#     if isinstance(d, dict):
#         return {k: filter_nested_dict(v, keys_to_keep) for k, v in d.items() if k in keys_to_keep}
#     elif isinstance(d, list):
#         return [filter_nested_dict(item, keys_to_keep) for item in d]
#     else:
#         return d


# this should not fail if keys on the list are not present
def filter_nested_dict(d, keys_to_keep):
    if isinstance(d, dict):
        return {k: filter_nested_dict(v, keys_to_keep) for k, v in d.items() if k in keys_to_keep and k in d}
    elif isinstance(d, list):
        return [filter_nested_dict(item, keys_to_keep) for item in d]
    else:
        return d    

def results_to_dataset(results):
    results = copy.deepcopy(results)
    ds = xr.Dataset()
    ds['USAXS_int'] = ('q',results['reducedData']['UPD'])
    ds['q'] = results['reducedData']['Q_array']
    del results['reducedData']['UPD']
    del results['reducedData']['Q_array']
    ds.update(results['reducedData'])
    for our_name,raw_name in [('AR_angle','ARangles'),
                              ('TimePerPoint','TimePerPoint'),
                              ('Monitor','Monitor'),
                              ('UPD','UPD_array'),
                             ]:
        ds[our_name] = ('flyscan_bin',results['RawData'][raw_name])
        del results['RawData'][raw_name]
    for our_name,raw_name in [('AmpGain','AmpGain'),
                              ('AmpReqGain','AmpReqGain'),
                              ('amp_change_channel','Channel')
                             ]:
        ds[our_name] = ('amp_change_channel',results['RawData'][raw_name])
        del results['RawData'][raw_name]
                                      
    ds.attrs.update(results['RawData']['metadata'])
    del results['RawData']['metadata']
    ds.attrs['instrument'] = results['RawData']['instrument']
    del results['RawData']['instrument']
    ds.update(results['RawData'])

    return ds

