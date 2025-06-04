# readTiled.py 
# version 0.1   2025-03-04


# this code will get last N data sets from tiled database for Flyscan, USAXS, SAXS, and WAXS
# it will recreate what IR3BS_GetJSONScanData() does in Indra package In3_ClueSkyreader.ipf
# it wil luse tiled web interface and read data into json.
# the main purpose is to get list fo filenames and paths for the last N scans
# and then use this list to download the data from the server
#the code will aslo remeber when was the last polling and ask for data only since last time. 
# need to handel gracefully the start
# code will trigger appropriate data recduction routines for different data sets
# at the end the code will generate necessary pictures for live data page.
# the code will be run in background and be checked by a cron job every 5 minutes

# method used buitls on https://github.com/BCDA-APS/bdp-tiled/blob/main/demo_client.ipynb
# and follows Igor code to get the right data sets


# import necessary libraries
import requests
import json
import datetime
import time
import socket
import logging



def iso_to_ts(isotime):
    return datetime.datetime.fromisoformat(isotime).timestamp()

def ts_to_iso(time):
    return datetime.datetime.fromtimestamp(time).isoformat()

current_hostname = socket.gethostname()
if current_hostname == 'usaxscontrol.xray.aps.anl.gov':
    server = "usaxscontrol.xray.aps.anl.gov"
else:
    server = "localhost"
port = 8020
catalog = "usaxs"



def print_results_summary(r):
    """We'll use this a few times."""
    xref = dict(First=0, Last=-1)
    for k, v in dict(First=0, Last=-1).items():
        if server ==  "localhost" :
            md = r["data"][v]["attributes"]["metadata"]["selected"]  #this is for UBUNTU VM, usaxscontrol does not have selected
        else:                                          
            md = r["data"][v]["attributes"]["metadata"]     #this is for usaxscontrol 
        #print(md)
        plan_name = md["plan_name"]
        scan_id = r["data"][v]["id"]
        started = ts_to_iso(md["time"])
        hdf5_file = md["hdf5_file"]
        hdf5_path = md["hdf5_path"]
        print(f"{k:5s} run: {plan_name=} started : {started} path: {hdf5_path=} {hdf5_file=} id: {scan_id}")


def convert_results(r):
    OutputList=[]
    for v in range(len(r["data"])):
        if server ==  "localhost" :
            md = r["data"][v]["attributes"]["metadata"]["selected"]  #this is for UBUNTU VM, usaxscontrol does not have selected
        else:                                          
            md = r["data"][v]["attributes"]["metadata"]["selected"]     #this is for usaxscontrol 
        #print(md)
        #plan_name = md["plan_name"]
        #scan_id = r["data"][v]["id"]
        #started = ts_to_iso(md["time"])
        hdf5_file = md["hdf5_file"]
        hdf5_path = md["hdf5_path"]
        #print(f" path: {hdf5_path=} {hdf5_file=}")
        if hdf5_file is not None:
            OutputList.append([hdf5_path,hdf5_file])
    return OutputList
        
#print(f'Search of {catalog=} has {len(r["data"])} runs.')
#print_results_summary(r)


def FindLastScanData(plan_name,NumScans=10):
    #print (FindLastScanData("Flyscan",10))
    #print (FindLastScanData("uascan",10))
    #print (FindLastScanData("SAXS",10))
    #print (FindLastScanData("WAXS",10))
    #print(f"Search for {plan_name=}")
    # Find all runs in a catalog between these two ISO8601 dates.
    # TODO - manage the times by rembering last call and only asking for data since last time
    #start_time = time.time()    #current time in seconds
    end_time = time.time()
    tz = "US/Central"

    #this filters for specific time AND for specific plan_name
    uri = (
        f"http://{server}:{port}"
        "/api/v1/search"
        f"/{catalog}"
        f"?page[limit]={NumScans}"                                                  # 0: all matching, 10 is 10 scans. Must be >0 value
        "&filter[eq][condition][key]=plan_name"                             # does not work... filter by plan_name
        f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
        f"&filter[time_range][condition][since]={(end_time-86400)}"         # time range, start time - 24 hours from now
        f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
        f"&filter[time_range][condition][timezone]={tz}"                    # time range
        "&sort=-time"                                                       # sort by time, -time gives last scans first
        "&fields=metadata"                                                  # return metadata
        "&omit_links=true"                                                  # no links
        "&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
        hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
        )
    logging.info(f"{uri=}")

    #print(f"{uri=}")
    try:
        r = requests.get(uri).json()
        #print(f'Search of {catalog=} has {len(r["data"])} runs.')
        #print_results_summary(r)
        # this is now a list of Flyscan data sets
        ScanList = convert_results(r)
        #print(ScanList)
        #logging.info('Received expected data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except: 
        # url communication failed, happens and shoudl not crash anything.
        # this is workaround.   
        logging.error('Could not get data from tiled server at  usaxscontrol.xray.aps.anl.gov')
        logging.error(f"Failed {uri=}")
        return []




#  http://10.211.55.7:8020/api/v1/search//usaxs/?page[offset]=0
# &page[limit]=100
# &filter[time_range][condition][since]=1738738800
# &filter[time_range][condition][timezone]=US/Central
# &sort=time
# &fields=metadata
# &omit_links=true
# &select_metadata={detectors:start.detectors,motors:start.motors,plan_name:start.plan_name,time:start.time,
# scan_title:start.plan_args.scan_title,title:start.title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}

#this requests all different records from given time range
# uri =(
#     f"http://{server}:{port}"  # standard prefix
#     "/api/v1/search"    # API command
#     f"/{catalog}"       # catalog
#     "?"                 # begin any command options
#     "page[limit]=10"   # 0: all matching
#     "&"                 # separator between any additional options
#     f"filter[time_range][condition][since]={iso_to_ts(start_time)}"
#     f"&filter[time_range][condition][until]={iso_to_ts(end_time)}"
#     f"&filter[time_range][condition][timezone]={tz}"
#     "&sort=time"
#     "&fields=metadata"
#     "&omit_links=true"
#     "&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
#                         hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"
# )
#print(f"{uri=}")
#r = requests.get(uri).json()
#print(r["data"][0])
