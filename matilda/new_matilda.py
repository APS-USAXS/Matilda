# this file will grab last N data files of various types.
# it will process them and generate plots
#!/usr/bin/env python3


import matplotlib.pyplot as plt
import pprint as pp
import numpy as np
import socket
import logging
from logging.handlers import RotatingFileHandler
import time
import os
import re
import datetime

from convertUSAXS import reduceStepScanToQR
from readfromtiled import FindLastScanData, FindLastBlankScan
from convertFlyscanNew import processFlyscan,reduceFlyscanToQR
#from developNewStepScan import processStepScan
from convertSAS import reduceADToQR
from supportFunctions import findProperBlankScan

# define any globals here
default_plt_font_size = 7
#imagePath = '/home/joule/WEBUSAXS/www_live/'  # Path to save images
imagePath = '/home/parallels/Desktop/'  # Path to save images
recalculateAllData = True

    
    
    
# Configure logging
# Get the directory of the current script
script_dir = os.path.dirname(os.path.abspath(__file__))
# Define the log directory path
log_dir = os.path.join(script_dir, 'log')
# Create the log directory if it doesn't exist
os.makedirs(log_dir, exist_ok=True)
# Define the log file path
log_file = os.path.join(log_dir, 'matilda.log')
handler = RotatingFileHandler(log_file, maxBytes=200000, backupCount=1)
logging.basicConfig(
    handlers=[handler],
    level=logging.INFO,        # Set the logging level
    format='%(asctime)s - %(levelname)s - %(message)s',  # Format of the log messages
    datefmt='%Y-%m-%d %H:%M:%S'  # Date format
)

# # Log messages
# logging.debug('This is a debug message')
# logging.info('This is an info message')
# logging.warning('This is a warning message')
# logging.error('This is an error message')
# logging.critical('This is a critical message')



# Here we process different types of scans
# Process the Flyscan data files
def processFlyscans(ListOfScans):
    results=[]
    for scan in ListOfScans:
        path = scan[0]
        filename = scan[1]
        #print(f"Processing file: {filename}")
        try:
            results.append(reduceFlyscanToQR(path, filename))
        except:
            pass
    #print("Done processing the Flyscans")
    return results

# Process the step scan data files
def processStepscans(ListOfScans):
    results=[]
    for scan in ListOfScans:
        path = scan[0]
        filename = scan[1]
        #print(f"Processing file: {filename}")
        try:
            results.append(reduceStepScanToQR(path, filename))
        except:
            pass
    #print("Done processing the Step scans")
    return results

# Process SAXS and WAXS data files
def processSASdata(ListOfScans):
    results=[]
    for scan in ListOfScans:
        path = scan[0]
        filename = scan[1]
        #print(f"Processing file: {filename}")
        try:
            results.append(reduceADToQR(path, filename))
        except:
            pass
    #print("Done processing the AD data")
    return results

def get_usaxs_r_plot_style():
    """
    Returns a dictionary with Matplotlib styling parameters for USAXS plots.
    """
    return {
        "figsize": (6, 6),
        "title": 'Plot of Normalized Intensity vs. Q',
        "xlabel": 'log(Q) [1/A]',
        "ylabel": 'Normalized Intensity',
        "xscale": 'log',
        "yscale": 'log',
        "xlim": (1e-5, 1),
        "grid": True,
        "font_size": default_plt_font_size  # Assumes default_plt_font_size is defined globally
    }

def get_usaxs_cal_plot_style():
    """
    Returns a dictionary with Matplotlib styling parameters for USAXS plots.
    """
    return {
        "figsize": (6, 6),
        "title": 'Plot of Calibrated Intensity vs. Q',
        "xlabel": 'log(Q) [1/A]',
        "ylabel": 'Calibrated Intensity',
        "xscale": 'log',
        "yscale": 'log',
        "xlim": (1e-5, 1),
        "grid": True,
        "font_size": default_plt_font_size  # Assumes default_plt_font_size is defined globally
    }

def plotUSAXSResults(ListOfresults, isFlyscan=True):  
    # Number of data sets
    num_data_sets = len(ListOfresults)
    # Choose a colormap
    cmap = plt.get_cmap('viridis')
    # Generate colors from the colormap
    colors = [cmap(i) for i in np.linspace(0, 1, num_data_sets)]
    print(f'Got {num_data_sets} data sets to plot')

    # Get plot styling
    style = get_usaxs_r_plot_style()
    # Set the font size to specific size
    plt.rcParams['font.size'] = style["font_size"]

    # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, color in zip(range(len(ListOfresults)),colors):
        data_dict = ListOfresults[i]
        label = data_dict["RawData"]["Filename"]
        Q_array = data_dict["reducedData"]["Q"]
        UPD = data_dict["reducedData"]["Intensity"]
        plt.plot(Q_array, UPD, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle

    plt.title(style["title"])
    plt.xlabel(style["xlabel"])
    plt.ylabel(style["ylabel"])
    plt.xscale(style["xscale"])
    plt.yscale(style["yscale"])
    plt.xlim(style["xlim"])
    plt.grid(style["grid"])
    # Add legend
    plt.legend()
    # Save the plot as a JPEG image
    if isFlyscan:
        plt.savefig(os.path.join(imagePath, 'usaxs.jpg'), format='jpg', dpi=300)
    else:
        plt.savefig(os.path.join(imagePath, 'stepusaxs.jpg'), format='jpg', dpi=300) # this step scan
    #plt.show()
    plt.close()

   # Get plot styling
    style = get_usaxs_cal_plot_style()
    # Set the font size to specific size
    plt.rcParams['font.size'] = style["font_size"]

   # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, color in zip(range(len(ListOfresults)),colors):
        data_dict = ListOfresults[i]
        label = data_dict["RawData"]["Filename"]
        Q_array = data_dict["CalibratedData"]["Q"]
        UPD = data_dict["CalibratedData"]["Intensity"]
        plt.plot(Q_array, UPD, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle

    plt.title(style["title"])
    plt.xlabel(style["xlabel"])
    plt.ylabel(style["ylabel"])
    plt.xscale(style["xscale"])
    plt.yscale(style["yscale"])
    plt.xlim(style["xlim"])
    plt.grid(style["grid"])
    # Add legend
    plt.legend()
    # Save the plot as a JPEG image
    if isFlyscan:
        plt.savefig(os.path.join(imagePath, 'usaxs_cal.jpg'), format='jpg', dpi=300)
    else:
        plt.savefig(os.path.join(imagePath, 'stepusaxs_cal.jpg'), format='jpg', dpi=300) # this step scan
    #plt.show()
    plt.close()



# def plotSWAXSResults(ListOfresults, isSAXS = True):  
#     # Number of data sets
#     num_data_sets = len(ListOfresults)
#     # Choose a colormap
#     cmap = plt.get_cmap('viridis')
#     # Generate colors from the colormap
#     colors = [cmap(i) for i in np.linspace(0, 1, num_data_sets)]

#     # Set the font size to specific size
#     plt.rcParams['font.size'] = default_plt_font_size 

#     # Plot ydata against xdata
#     plt.figure(figsize=(6, 6))
#     for i, color in zip(range(len(ListOfresults)),colors):
#         data_dict = ListOfresults[i]
#         label = data_dict["RawData"]["Filename"]
#         Q_array = data_dict["reducedData"]["Q_array"]
#         UPD = data_dict["reducedData"]["Intensity"]
#         plt.plot(Q_array, UPD, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle
#     plt.ylabel('Intensity')   
#     if isSAXS:
#         plt.title('Plot of SAXS Intensity vs. Q')   
#         plt.xlabel('log(Q) [1/A]')
#         plt.xscale('log')
#         plt.yscale('log')
#         #plt.xlim(1e-5, 1)
#         plt.grid(True)
#         # Add legend
#         plt.legend()
#         current_hostname = socket.gethostname()
#         if current_hostname == 'usaxscontrol.xray.aps.anl.gov':
#             plt.savefig('/home/joule/WEBUSAXS/www_live/saxs.jpg', format='jpg', dpi=300)
#             plt.close()
#         else:
#             plt.savefig('saxs.jpg', format='jpg', dpi=300)
#         plt.show()
#     else:
#         plt.title('Plot of WAXS Intensity vs. Q')   
#         plt.xlabel('Q [1/A]')
#         plt.xscale('linear')
#         plt.yscale('linear')        
#         #plt.xlim(1e-5, 1)
#         plt.grid(True)
#         # Add legend
#         plt.legend()
#         # Save the plot as a JPEG image
#         current_hostname = socket.gethostname()
#         if current_hostname == 'usaxscontrol.xray.aps.anl.gov':
#             plt.savefig('/home/joule/WEBUSAXS/www_live/waxs.jpg', format='jpg', dpi=300)
#             plt.close()
#         else:
#             plt.savefig('waxs.jpg', format='jpg', dpi=300)
#         plt.show()

def processFlyscans(ListOfScans, ListOfBlanks):
    """
    Processes a list of flyscans, finding the appropriate blank for each.
    The correct blank is in the same path and has the largest number XYZ
    that is smaller than the sample's number.
    """
    results=[]
    logging.info("Starting processing of flyscans with blanks")    
    for scan_path, scan_filename in ListOfScans:
        selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
        result = processFlyscan(scan_path, scan_filename,
                                    blankPath=selected_blank_path,
                                    blankFilename=selected_blank_filename,
                                    deleteExisting=recalculateAllData)          # deleteExisting will be True for reprocessing
        results.append(result)
        #except Exception as e:
        #    logging.error(f"Error processing scan {scan_filename} with blank {selected_blank_filename}: {e}")
    return results




if __name__ == "__main__":
    #try:
    #    while True:
    #logging.info("New round of processing started at : %s", datetime.datetime.now()) 
    #print("New round of processing started at : ", datetime.datetime.now()) 
    #logging.info('Processing the Flyscans')
    print("Processing the Flyscans")
    #ListOfScans = FindLastScanData("Flyscan",1,0)
    ListOfScans = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'R6016HRC_T4_V_1077.h5'],
    ]
                   #['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   #'AirBlank_1076.h5'],]
    path, filename = ListOfScans[0]
    #listOfBlanks = FindLastBlankScan("Flyscan",path, 1,0)
    listOfBlanks = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'AirBlank_1076.h5']]
    print(f'Got list : {ListOfScans}')
    print(f'Got blank list : {listOfBlanks}')
    
    results = processFlyscans(ListOfScans, listOfBlanks)
    
    plotUSAXSResults(results, isFlyscan=True)  
        
    
    #print(f'Got results : {results}')
    #processFlyscan(samplePath,sampleName,blankPath=blankPath,blankFilename=blankFilename,deleteExisting=True)
    # returns dictionary of this type:
    #         result["SampleName"]=sampleName
    #         result["BlankName"]=blankName
    #         result["reducedData"] =  {"Intensity":np.ravel(intensity), 
    #                           "Q":np.ravel(q),
    #                           "Error":np.ravel(error)}
    #         result["CalibratedData"] = {"Intensity":np.ravel(intcalib),
    #                                 "Q":np.ravel(qcalib),
    #                                 "Error":np.ravel(errcalib),
    #                                }  



    #         ListOfresults = processFlyscans(ListOfScans)
    #         logging.info(f'Got list : {ListOfScans}')
    #         if len(ListOfresults) > 0:
    #             plotUSAXSResults(ListOfresults,isFlyscan=True)
    #         else:
    #             logging.info('No Flyscan data found')

    #         logging.info('Processing the step scans')
    #         #print("Processing the Stepscans")
    #         ListOfScans = GetListOfScans("uascan")
    #         ListOfresults = processStepscans(ListOfScans)
    #         if len(ListOfresults) > 0:
    #             plotUSAXSResults(ListOfresults,isFlyscan=False)
    #         else:
    #             logging.info('No Step scan data found') 

    #         #print("Done processing the Step scans")
    #         #print("Processing the SAXS scans")
    #         logging.info('Processing the SAXS')
    #         ListOfScans = GetListOfScans("SAXS")
    #         ListOfresults = processSASdata(ListOfScans)
    #         if len(ListOfresults) > 0:
    #             plotSWAXSResults(ListOfresults,isSAXS = True)
    #         else:
    #             logging.info('No SAXS data found')
               
    #         #print("Processing the WAXS scans")
    #         logging.info('Processing the WAXS')
    #         ListOfScans = GetListOfScans("WAXS")
    #         ListOfresults = processSASdata(ListOfScans)
    #         if len(ListOfresults) > 0:
    #             plotSWAXSResults(ListOfresults, isSAXS = False)
    #         else:
    #             logging.info('No WAXS data found')

    #         # wait for more data, 30s seems reasonable
    #         time.sleep(30)
    # except KeyboardInterrupt:
    #     print("Keyboard interrupt")
