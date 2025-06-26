# Convert Nika data to pyFAI using Fit2D 
'''
    Convert Nika SAS geometry parameters to pyFAI using Fit2D format.
    Needs pyFAI library. 
    Convert Nika SDD, pix size, BCX, BCY, HorTilt, verTilt 
    First convert into Fit2D format
    Then use pyFAI to convert Fit2D to poni format.
    NOTE: this seems to wokr only for data from hdf5 files = Nexus NXsas files. 
    Tiff ffiles seem to be loaded differently between Nika and pyFAI. Geomtry conversion (signs/X/Y changes are needed. )
'''

from pyFAI.geometry.fit2d import Fit2dGeometry, convert_from_Fit2d
from pprint import pprint as pp
import numpy as np

# Fit2D parameter definitions/example
# direct_dist = 200.0  # Distance from sample to detector in mm
# center_x = 1000.0    # Beam center X in pixels
# center_y = 1000.0    # Beam center Y in pixels
# tilt = 0.0           # Tilt angle in degrees
# tilt_plan_rotation = 0.0  # Tilt plan rotation in degrees
# pixel_x = 75         # Pixel size in um (X direction)
# pixel_y = 75         # Pixel size in um (Y direction)
# spline_file = None   # Path to spline file if any
# detector = None      # Detector object if any
# wavelength = 0.59    # Wavelength in Angstroms (this is 21keV)


# the main callable function 
def convert_Nika_to_Fit2D(*, SSD, pix_size, BCX, BCY, HorTilt, VertTilt, wavelength):
    '''
    Convert Nika data to Fit2D format.
    Needs pyFAI
    Convert Nika SDD, pix size, BCX, BCY, HorTilt, verTiilt
    first into Fit2D format
    Then use pyFAI to convert the data to Fit2D format.
    units: 
        SDD in mm, 
        pix_size in mm (convert below to um), 
        BCX, BCY in pixels, 
        HorTilt, VertTilt in degrees, 
        wavelength in A
    '''
    # Check for missing or None parameters
    if SSD is None:
        raise TypeError("Missing required parameter: SSD")
    if pix_size is None:
        raise TypeError("Missing required parameter: pix_size")
    if BCX is None:
        raise TypeError("Missing required parameter: BCX")
    if BCY is None:
        raise TypeError("Missing required parameter: BCY")
    if HorTilt is None:
        raise TypeError("Missing required parameter: HorTilt")
    if VertTilt is None:
        raise TypeError("Missing required parameter: VertTilt")
    if wavelength is None:
        raise TypeError("Missing required parameter: wavelength")


    Fit2D_TiltDir, Fit2D_TiltAngle = convert_Nika_to_Fit2D_angles(HorTilt, VertTilt)     

    # Create a Fit2dGeometry object
    fit2d_geometry = Fit2dGeometry(
        directDist=SSD,
        centerX=BCY,
        centerY=BCX,
        tilt=Fit2D_TiltAngle,
        tiltPlanRotation=Fit2D_TiltDir,
        pixelX=pix_size*1000,
        pixelY=pix_size*1000,
        splineFile=None,
        detector=None,
        wavelength=wavelength
    )
    #pp(fit2d_geometry)
    # Convert to PONI object
    poni_object = convert_from_Fit2d(fit2d_geometry)

    # Print the PONI object
    #pp(poni_object)
    return poni_object


def convert_Nika_to_Fit2D_angles(HorTilt, VertTilt):
    '''
    Convert angles from Nika format to Fit2D format.
    The conversion is based on the assumption that the angles are given in degrees.
    The conversion is done by negating the angles and converting them to radians.
    The function returns the converted angles in degrees.
    '''
    # Convert angles from degrees to radians
    Ax_rad = np.radians(HorTilt)
    Ay_rad = np.radians(VertTilt)
    
    # in any of the two tilts is negative, the tilt direction will be negative.
    sign = 1
    if Ax_rad < 0:
        sign = -1
    if Ay_rad < 0:
        sign = -1

    # Calculate the rotation vector components
    # Assuming small angles, the rotation vector can be approximated as:
    # Rx = [1, 0, 0] rotation around x-axis
    # Ry = [0, 1, 0] rotation around y-axis
    # Combined effect in the xy-plane
    Rx = np.array([0, np.sin(Ax_rad), np.cos(Ax_rad)])
    Ry = np.array([np.sin(Ay_rad), 0, np.cos(Ay_rad)])

    # Combined rotation vector
    R_combined = Rx + Ry

    # Calculate tilt direction angle (theta) in the xy-plane
    theta = sign*np.degrees(np.arctan2(R_combined[1], R_combined[0]))

     # Calculate tilt angle (phi) as the magnitude of the rotation vector
    phi = np.degrees(np.arcsin(np.linalg.norm([np.sin(Ax_rad), np.sin(Ay_rad)])))
    return theta, phi

