'''
    Will convert SAXS and WAXS area detector data from the HDF5 format to the 1Ddata
    we will use pyFai to do the conversion
    we will convert Nika parameters to Fit2D format and then use pyFAI to convert to poni format
    Only some metadata are kept to keep all more reasnobale on size

'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import h5py
from supportFunctions import read_group_to_dict, filter_nested_dict
import pprint as pp
import os
import tifffile as tiff
import logging
from convertNikaTopyFAI import convert_Nika_to_Fit2D
from hdf5code import save_dict_to_hdf5, load_dict_from_hdf5


def PlotResults(data_dict):
    # Find Peak center and create Q vector.
    Q_array = data_dict["reducedData"]["Q_array"]
    Intensity = data_dict["reducedData"]["Intensity"]
    # Plot ydata against xdata
    plt.figure(figsize=(6, 12))
    plt.plot(Q_array, Intensity, linestyle='-')  # You can customize the marker and linestyle
    plt.title('Plot of Intensity vs. Q')
    plt.xlabel('log(Q) [1/A]')
    plt.ylabel('Intensity')
    plt.xscale('linear')
    plt.yscale('linear')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    Sample = dict()
    Sample=reduceADToQR("./TestData/TestTiltData","LaB6_tilt7v_0049.hdf")
    #Sample["reducedData"]=test("/home/parallels/Github/Matilda/TestData","LaB6_45deg.tif")
    #pp.pprint(Sample)
    PlotResults(Sample)



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
