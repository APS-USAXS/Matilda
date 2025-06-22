# this file will grab last N data files of various types.
# it will process them and generate plots
#!/usr/bin/env python3


import pprint as pp
import numpy as np
import socket
import logging
from logging.handlers import RotatingFileHandler
import os


from convertUSAXS import reduceStepScanToQR
from readfromtiled import FindLastScanData, FindLastBlankScan
from convertFlyscanNew import processFlyscan,reduceFlyscanToQR
from convertSWAXSNew import process2Ddata
#from developNewStepScan import processStepScan
from supportFunctions import findProperBlankScan
from plotData import plotUSAXSResults, plotSWAXSResults


#set ImagePath to None to prevent saving images
#imagePath = None  # Do not save images
#imagePath = '/home/joule/WEBUSAXS/www_live/'  # Path to save images
imagePath = '/home/parallels/Desktop/'  # Path to save images


recalculateAllData = False  # Set to True to recalculate all data, False to use existing data

    
    
    
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
# def processFlyscans(ListOfScans):
#     results=[]
#     for scan in ListOfScans:
#         path = scan[0]
#         filename = scan[1]
#         #print(f"Processing file: {filename}")
#         try:
#             results.append(reduceFlyscanToQR(path, filename))
#         except:
#             pass
#     #print("Done processing the Flyscans")
#     return results

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

def processADscans(ListOfScans, ListOfBlanks):
    """
    Processes a list of SAXS/WAXS data, finding the appropriate blank for each.
    The correct blank is in the same path and has the largest number XYZ
    that is smaller than the sample's number.
    """
    results=[]
    logging.info("Starting processing of SAXS/WAXS with blanks")    
    for scan_path, scan_filename in ListOfScans:
        selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
        result = process2Ddata(scan_path, scan_filename,
                                    blankPath=selected_blank_path,
                                    blankFilename=selected_blank_filename,
                                    deleteExisting=recalculateAllData)          # deleteExisting will be True for reprocessing
        results.append(result)
    return results

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
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'R6016HRC_T4_H_V_1084.h5'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'R6016HRC_RB_H_1085.h5'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'R6016ACT_T4_V_H_1086.h5'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'R6016ACT_T4_H_V_1087.h5'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'R6016HRC_T4_V_H_1083.h5'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'AirBlank_1076.h5'],
                   ]
    path, filename = ListOfScans[0]
    #listOfBlanks = FindLastBlankScan("Flyscan",path, 1,0)
    listOfBlanks = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   'AirBlank_1076.h5']]
    print(f'Got list : {ListOfScans}')
    print(f'Got blank list : {listOfBlanks}')
    
    results = processFlyscans(ListOfScans, listOfBlanks)
    
    plotUSAXSResults(results, imagePath, isFlyscan=True)  
        


    print("Processing the SAXS")
    #ListOfScans = FindLastScanData("Flyscan",1,0)
    ListOfScans = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'R6016HRC_T4_V_1077.hdf'],   
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'R6016HRC_T4_H_V_1084.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'R6016HRC_RB_H_1085.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'R6016ACT_T4_V_H_1086.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'R6016ACT_T4_H_V_1087.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'R6016HRC_T4_V_H_1083.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'AirBlank_1076.hdf'],
                   ]
                   #['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   #'AirBlank_1076.hdf'],]
    path, filename = ListOfScans[0]
    #listOfBlanks = FindLastBlankScan("Flyscan",path, 1,0)
    listOfBlanks = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_saxs',
                   'AirBlank_1076.hdf']]
    print(f'Got list : {ListOfScans}')
    print(f'Got blank list : {listOfBlanks}')
    
    results = processADscans(ListOfScans, listOfBlanks)
    
    plotSWAXSResults(results, imagePath, isSAXS = True)  

    print("Processing the WAXS")
    #ListOfScans = FindLastScanData("Flyscan",1,0)
    ListOfScans = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'R6016HRC_T4_V_1077.hdf'],   
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'R6016HRC_T4_H_V_1084.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'R6016HRC_RB_H_1085.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'R6016ACT_T4_V_H_1086.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'R6016ACT_T4_H_V_1087.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'R6016HRC_T4_V_H_1083.hdf'],
                   ['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'AirBlank_1076.hdf'],
                   ]
                   #['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_usaxs',
                   #'AirBlank_1076.hdf'],]
    path, filename = ListOfScans[0]
    #listOfBlanks = FindLastBlankScan("Flyscan",path, 1,0)
    listOfBlanks = [['/home/parallels/Desktop/06_15_Rakesh/06_15_Rakesh_waxs',
                   'AirBlank_1076.hdf']]
    print(f'Got list : {ListOfScans}')
    print(f'Got blank list : {listOfBlanks}')
    
    results = processADscans(ListOfScans, listOfBlanks)
    
    plotSWAXSResults(results, imagePath, isSAXS = False)  


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
