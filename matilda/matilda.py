#!/usr/bin/env python3
'''
Matilda - USAXS/SAXS/WAXS data processing tool

When run as a script, it will process the data from the last scan and blank scan

User facing functions are defined in this file.

processADscans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False) 
    - SAXS or WAXS data, finding the appropriate blank for each.
    - ListOfScans is a list of tuples with path and filename of the scan data
    - ListOfBlanks is a list of tuples with path and filename of the blank data
    - forceFirstBlank set to true will skip Blank selection and use the first blank in the list
    - recalculateAllData set to true will force reprocessing of the data, otherwise it
    - returns a list of dictionaries with the reduced data
    - No plotting, just processing

processFlyscans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False)
    - Flyscan data, finding the appropriate blank for each.
    - ListOfScans is a list of tuples with path and filename of the scan data
    - ListOfBlanks is a list of tuples with path and filename of the blank data
    - forceFirstBlank set to true will skip Blank selection and use the first blank in the list
    - recalculateAllData set to true will force reprocessing of the data, otherwise it
    - returns a list of dictionaries with the reduced data
    - No plotting, just processing
    
TODO: processStepScans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False)
    - Step scan data, finding the appropriate blank for each.
    - ListOfScans is a list of tuples with path and filename of the scan data
    - ListOfBlanks is a list of tuples with path and filename of the blank data
    - forceFirstBlank set to true will skip Blank selection and use the first blank in the list
    - recalculateAllData set to true will force reprocessing of the data, otherwise it
    - returns a list of dictionaries with the reduced data
    - No plotting, just processing

processUSAXSFolder(path)
    - will process (with forced reprocessing) all the scans in the given path (set recalculateAllData=True)
    - assumes USAXS data are in _usaxs folder
    - assumes SAXS data are in _saxs folder
    - assumes WAXS data are in _waxs folder
    - will process each data scan in the folder
    - no plotting, just processing

When run as main, it will process the data from the last 10 scans and blank scan
for each USAXSstep, Flyscan, SAXS, WAXS data types.


TODO: add to each processXYZ option to force Blank if Blanks is only one and avoid checking on order number. 
        This is to enable use for Igor or elsewhere to use different blank than was measured prior experiment. 

'''

import pprint as pp
import numpy as np
import socket
import re
import time
import logging
import datetime
from logging.handlers import RotatingFileHandler
import os


from convertUSAXS import reduceStepScanToQR
from readfromtiled import FindLastScanData, FindLastBlankScan
from convertFlyscanNew import processFlyscan
from convertSWAXSNew import process2Ddata
#from developNewStepScan import processStepScan
from supportFunctions import findProperBlankScan
from plotData import plotUSAXSResults, plotSWAXSResults


#set ImagePath to None to prevent saving images
#imagePath = None  # Do not save images
imagePath = '/home/joule/WEBUSAXS/www_live/'  # Path to save images
#imagePath = '/home/parallels/Desktop/'  # Path to save images

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


# user facing functions

def processUSAXSFolder(path):
    # Get the list of folders in the path folder
    if not os.path.exists(path):
        logging.error(f"The path {path} does not exist.")
        return
    # Get the list of folders in the path folder
    folders = os.listdir(path)
    # Filter the list to include only the folders that end with _usaxs
    usaxs_folder = [folder for folder in folders if folder.endswith('_usaxs')][0]
    
    saxs_folder = [folder for folder in folders if folder.endswith('_saxs')][0]
    
    waxs_folder = [folder for folder in folders if folder.endswith('_waxs')][0]
    
    #get list of files in the usaxs folder
    usaxs_files = []
    files = os.listdir(os.path.join(path, usaxs_folder))
    # Filter the list to include only the files with the extension .h5
    # Sort by the number before .hdf
    h5_files = [file for file in files if file.endswith('.h5')]
    usaxs_files = sorted(h5_files, key=extract_number_from_filename)
    usaxs_blanks = [file for file in usaxs_files if 'blank' in file.lower()]

    #get list of files in the saxs folder
    saxs_files = []
    files = os.listdir(os.path.join(path, saxs_folder))
    # Filter the list to include only the files with the extension .hdf
    # Sort by the number before .hdf
    h5_files = [file for file in files if file.endswith('.hdf')]
    saxs_files = sorted(h5_files, key=extract_number_from_filename)
    saxs_blanks = [file for file in saxs_files if 'blank' in file.lower()]

    
    #get list of files in the waxs folder
    waxs_files = []
    files = os.listdir(os.path.join(path, waxs_folder))
    # Filter the list to include only the files with the extension .hdf
    # Sort by the number before .hdf
    h5_files = [file for file in files if file.endswith('.hdf')]
    waxs_files = sorted(h5_files, key=extract_number_from_filename)
    waxs_blanks = [file for file in waxs_files if 'blank' in file.lower()]
       
    #force recalculateAllData to True
    recalculateAllData = True
    
    #process all flyscans
    logging.info("Processing USAXS Flyscans")
    ListOfScans = [(os.path.join(path, usaxs_folder), file) for file in usaxs_files]
    ListOfBlanks = [(os.path.join(path, usaxs_folder), file) for file in usaxs_blanks]
    data = processFlyscans(ListOfScans, ListOfBlanks, recalculateAllData=recalculateAllData)

    logging.info("Processing SAXS data")
    ListOfScans = [(os.path.join(path, saxs_folder), file) for file in saxs_files]
    ListOfBlanks = [(os.path.join(path, saxs_folder), file) for file in saxs_blanks]
    data = processADscans(ListOfScans, ListOfBlanks, recalculateAllData=recalculateAllData)


    logging.info("Processing WAXS data")
    ListOfScans = [(os.path.join(path, waxs_folder), file) for file in waxs_files]
    ListOfBlanks = [(os.path.join(path, waxs_folder), file) for file in waxs_blanks]
    data = processADscans(ListOfScans, ListOfBlanks, recalculateAllData=recalculateAllData)
    



#process list of flyscans, finding the appropriate blank for each.
def processFlyscans(ListOfScans, ListOfBlanks, recalculateAllData=recalculateAllData,forceFirstBlank=False):
    """
    Processes a list of flyscans, finding the appropriate blank for each.
    The correct blank is in the same path and has the number _XYZ
    that is smaller than the sample's number.
    forceFirstBlank=True will skip Blank selection and use the first blank in the list
    recalculateAllData=True will force reprocessing of the data, otherwise it
    """
    results=[]
    logging.info("Processing of flyscans with blanks")    
    for scan_path, scan_filename in ListOfScans:
        try:
            if forceFirstBlank:
                # Use the first blank in the list
                selected_blank_path = ListOfBlanks[0][0]
                selected_blank_filename = ListOfBlanks[0][1]
                #logging.info(f"User forced use of blank path: {selected_blank_path}, filename: {selected_blank_filename}")
            else:
                selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
                 #logging.info(f"Automatically selected blank path: {selected_blank_path}, filename: {selected_blank_filename}")
           
            logging.info(f"Flyscan data processing: path {scan_path}, filename {scan_filename}, blank path {selected_blank_path}, blank filename {selected_blank_filename}")
            result = processFlyscan(scan_path, scan_filename,
                                        blankPath=selected_blank_path,
                                        blankFilename=selected_blank_filename,
                                        deleteExisting=recalculateAllData)          # deleteExisting will be True for reprocessing
            results.append(result)
        except Exception as e:
            logging.error(f"Error processing scan {scan_filename}: {e}", exc_info=True)
    return results


# Process the step scan data files
def processStepscans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False):
    results=[] 
    for scan in ListOfScans:
        path = scan[0]
        filename = scan[1]
        logging.info(f"Processing step scan: {filename}")
        try:
            results.append(reduceStepScanToQR(path, filename))
        except Exception as e:
            logging.error(f"Failed to process step scan {filename} in {path}: {e}", exc_info=True)
    return results

#process SAXS/WAXS data, finding the appropriate blank for each.
def processADscans(ListOfScans, ListOfBlanks,recalculateAllData=recalculateAllData,forceFirstBlank=False):
    """
    Processes a list of SAXS/WAXS data, finding the appropriate blank for each.
    The correct blank is in the same path and has the largest number XYZ
    that is smaller than the sample's number.
    """
    results=[]
    logging.info("Processing of SAXS/WAXS with blanks")    
    for scan_path, scan_filename in ListOfScans:
        try:
            if forceFirstBlank:
                # Use the first blank in the list
                selected_blank_path = ListOfBlanks[0][0]
                selected_blank_filename = ListOfBlanks[0][1]
                #logging.info(f"User forced use of blank path: {selected_blank_path}, filename: {selected_blank_filename}")
            else:
                selected_blank_path, selected_blank_filename   = findProperBlankScan(scan_path, scan_filename, ListOfBlanks)
                #logging.info(f"Automatically selected blank path: {selected_blank_path}, filename: {selected_blank_filename}")

            logging.info(f"2D data processing: path {scan_path}, filename {scan_filename}, blank path {selected_blank_path}, blank filename {selected_blank_filename}")
            result = process2Ddata(scan_path, scan_filename,
                                        blankPath=selected_blank_path,
                                        blankFilename=selected_blank_filename,
                                        deleteExisting=recalculateAllData)          # deleteExisting will be True for reprocessing
            results.append(result)
        except Exception as e:
            logging.error(f"Error processing AD scan {scan_filename}: {e}", exc_info=True)
    return results


def extract_number_from_filename(filename):
    """Extract number from filename, return 0 if no number found"""
    match = re.search(r'_(\d+)\.hdf', filename)
    return int(match.group(1)) if match else 0



if __name__ == "__main__":
    try:
        listofFlyscansOld=dict()
        listofStepScansOld=dict()
        listofSAXSOld=()
        listOfWAXSOld=()
        while True:
            logging.info("New round of processing started at : %s", datetime.datetime.now()) 

            logging.info("Processing the Flyscans")
            ListOfScans = FindLastScanData("Flyscan",10,5)
            #need to check, if we have same list as before and if yes, skip everything
            if ListOfScans == listofFlyscansOld or len(ListOfScans) == 0:
                logging.info('No new Flyscan data found')
            else:
                #path, filename = ListOfScans[-1]
                listOfBlanks = FindLastBlankScan("Flyscan",path=None, NumScans=10, lastNdays=5) 
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                results = processFlyscans(ListOfScans, listOfBlanks)
                plotUSAXSResults(results, imagePath, isFlyscan=True) 
                listofFlyscansOld = ListOfScans 
            
        
            # logging.info('Processing the step scans')
            # ListOfScans = FindLastScanData("uascan",10,50) 
            # listOfBlanks = FindLastBlankScan("uascan",path=None, NumScans=10, lastNdays=50)
            # ListOfresults = processStepscans(ListOfScans)
            # if len(ListOfresults) > 0:
            #     plotUSAXSResults(ListOfresults,isFlyscan=False)
            # else:
            #     logging.info('No Step scan data found') 

            logging.info("Processing the SAXS")
            ListOfScans = FindLastScanData("SAXS",10,5)
            if ListOfScans == listofSAXSOld or len(ListOfScans) == 0:
                logging.info('No new SAXS data found')
            else:
                listOfBlanks = FindLastBlankScan("SAXS",path=None, NumScans=10, lastNdays=5)
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                results = processADscans(ListOfScans, listOfBlanks)
                plotSWAXSResults(results, imagePath, isSAXS = True) 
                listofSAXSOld = ListOfScans 

            logging.info("Processing the WAXS")
            ListOfScans = FindLastScanData("WAXS",10,5)
            if ListOfScans == listOfWAXSOld or len(ListOfScans) == 0:
                logging.info('No new WAXS data found')
            else:
                listOfBlanks = FindLastBlankScan("WAXS",path=None, NumScans=10, lastNdays=5)
                logging.info(f'Got list : {ListOfScans}')
                logging.info(f'Got blank list : {listOfBlanks}')
                results = processADscans(ListOfScans, listOfBlanks)      
                plotSWAXSResults(results, imagePath, isSAXS = False)  
                listofWAXSOld = ListOfScans 

            time.sleep(30)
    except KeyboardInterrupt:
            logging.info('Keyboard interrupt') 
