'''
    readTiled.py 
    version
    0.2   2025-04-15
    0.3   2025-06-01

    useful functions are:
    readfromtiled.FindLastScanData(plan_name,NumScans=10, LastNdays=1)
    readfromtiled.FindScanDataByName(plan_name,scan_title,NumScans=1,lastNdays=0)
    readfromtiled.FindLastBlankScan(plan_name,NumScans=1, lastNdays=0)

    method used builds on https://github.com/BCDA-APS/bdp-tiled/blob/main/demo_client.ipynb
    and follows Igor code to get the right data sets
'''

# import necessary libraries
import requests
import json
import datetime
import time
import socket
import logging
from typing import Any, Optional


def iso_to_ts(isotime):
    return datetime.datetime.fromisoformat(isotime).timestamp()

def ts_to_iso(time):
    return datetime.datetime.fromtimestamp(time).isoformat()

current_hostname = socket.gethostname()
if current_hostname == 'usaxscontrol.xray.aps.anl.gov':
    server = "usaxscontrol.xray.aps.anl.gov"
else:
    #server = "localhost"
    server = "usaxscontrol.xray.aps.anl.gov"

port = 8020
catalog = "usaxs"
TILED_TIMEOUT = 10  # seconds

def tiled_get(
        *, 
        router: Optional[str] = "search",  # "search" and "metadata" are common"
        timeout: Optional[float] = TILED_TIMEOUT,
        uid: Optional[str] = None,
        **params: Optional[dict[str, Any]],  # Tiled options after the ? in the URL
    ) -> dict:
    """
    Call 'requests.get()' with our server & catalog details.

    EXAMPLES::

        # information about the last run
        params = {"page[limit]": 3, "sort": "-time"}
        run_info = tiled_get(**params)

        # metadata of a specific run
        run_md_info = tiled_get(router="metadata", uid="f3b1c4e5-7f3a-4a3b-8c9d-3e2f1b4c5d6e")
        md = run_md_info["data"]["attributes"]["metadata"]
    """
    url = f"http://{server}:{port}/api/v1/{router}/{catalog}"
    if isinstance(uid, str):
        url += f"/{uid}"
    response = requests.get(url, params=params, timeout=timeout)
    return response.json()


def successful_run(uid: Optional[str] = None) -> bool:
    """
    Was the Bluesky run with this uid successful?

    When 'uid' is None, then report about the last run in the catalog.

    This involves searching for a key in the 'stop' document.
    The tiled server does not have direct way to query for keys that
    are not in the 'start' document.

    EXAMPLES:

        successful_run()  # last run
        successful_run("f3b1c4e5-7f3a-4a3b-8c9d-3e2f1b4c5d6e")  # specific run
    """
    if uid is None:
        params = {
            "page[offset]": 0,
            "page[limit]": 1,
            "sort": "-time",
        }
        run_info = tiled_get(**params)
        last_run_index = 0  # when sorted by reverse time (params["sort"] value)
        md = run_info["data"][last_run_index]["attributes"]["metadata"]
    else:
        run_info = tiled_get(router="metadata", uid=uid)
        md = run_info["data"]["attributes"]["metadata"]
    success = (md.get("stop") or {}).get("exit_status", "?") == "success"
    return success


def print_results_summary(r):
    """We'll use this a few times."""
    xref = dict(First=0, Last=-1)
    for k, v in dict(First=0, Last=-1).items():
        md = r["data"][v]["attributes"]["metadata"]["selected"]  #From 6-1-2025 ["selected"] is in both VM and usaxscontrol tiled
        #print(md)
        plan_name = md["plan_name"]
        scan_id = r["data"][v]["id"]
        started = ts_to_iso(md["time"])
        hdf5_file = md["hdf5_file"]
        hdf5_path = md["hdf5_path"]
        print(f"{k:5s} run: {plan_name=} started : {started} path: {hdf5_path=} {hdf5_file=} id: {scan_id}")


def convert_results(r):
    OutputList=[]
    logging.info(f"Parsing data")
    for v in range(len(r["data"])):
        md = r["data"][v]["attributes"]["metadata"]
        plan_name = md["plan_name"]
        success = (md.get("stop") or {}).get("exit_status", "?") == "success"
        if not success and plan_name is "Flyscan":
            logging.info(f"Skipping a scan which has not succeeded.")
            continue  # skip this run

        md = r["data"][v]["attributes"]["metadata"]["selected"]  #From 6-1-2025 ["selected"] is in both VM and usaxscontrol tiled. 
        #print(md)
        #plan_name = md["plan_name"]
        #scan_id = r["data"][v]["id"]
        #started = ts_to_iso(md["time"])
        hdf5_file = md["hdf5_file"]
        hdf5_path = md["hdf5_path"]
        #print(f" path: {hdf5_path=} {hdf5_file=}")
        if hdf5_file is not None:
            OutputList.append([hdf5_path,hdf5_file])
    logging.info(f"Succeeded, here is OutputList: {OutputList}")
    return OutputList
        
#print(f'Search of {catalog=} has {len(r["data"])} runs.')
#print_results_summary(r)


def FindScanDataByName(plan_name,scan_title,NumScans=1,lastNdays=0):
    #this filters for specific time AND for specific plan_name
    if lastNdays > 0:
        # if LastNdays is set, then we will ask for data from the last N days
        start_time = time.time() - (lastNdays * 86400)
        end_time = time.time()      # current time in seconds
        end_time = end_time - 10    # this is to fix SAXS/WAXS possibly not being done. 
        tz = "US/Central"
        uri = (
            f"http://{server}:{port}"
            "/api/v1/search"
            f"/{catalog}"
            f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
            "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
            f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
            "&filter[eq][condition][key]=title"                                 # filter by title
            f'&filter[eq][condition][value]="{scan_title}"'                     # filter by title value
            f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
            f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
            f"&filter[time_range][condition][timezone]={tz}"                    # time range
            "&sort=-time"                                                       # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            #"&select_metadata={plan_name:start.plan_name,time:start.time,\
            #scan_title:start.plan_args.scan_title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"
                                                                                # select metadata
            )
    else:
        uri = (
            f"http://{server}:{port}"
            "/api/v1/search"
            f"/{catalog}"
            f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
            "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
            f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
            "&filter[eq][condition][key]=title"                                 # filter by title
            f'&filter[eq][condition][value]="{scan_title}"'                     # filter by title value
            "&sort=-time"                                                       # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            #"&select_metadata={plan_name:start.plan_name,time:start.time,\
            #scan_title:start.plan_args.scan_title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   
                                                                                # select metadata
            )
      
    logging.info(f"{uri=}")
    #additional keywords:
    #[noteq] - not equal
    #[contains] - seems same as eq in use, and cannot be made into case insensitive. Not useful. 
    #[in] - in a list of values
    #[notin] - not in a list of values
    #[comparison] - comparison with lt, gt, le, ge for numerical values
    #working examples:
    #http://10.211.55.7:8020/api/v1/search/usaxs/?page[limit]=10&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)blank&sort=-time
    #returns list of "Blank" samples, not not ist of samples contains "blank" in name
    #http://10.211.55.7:8020/api/v1/search/usaxs/?page[limit]=1&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)water*blank&sort=-time
    #returns last scan which conatins case independent "water blank" in name
    #http://10.211.55.7:8020/api/v1/search/usaxs/?page[limit]=1&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)blank&sort=-time&omit_links=true&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}
    #returns last scan which conatisn case independet "water blank" in name

    try:
        r = requests.get(uri).json()
        logging.info(f"Got json for : {plan_name}")
        ScanList = convert_results(r)
        #logging.info('Received expected data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except: 
        # url communication failed, happens and should not crash anything.
        # this is workaround.   
        logging.error('Could not get data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.error(f"Failed {uri=}")
        return []
    

def FindLastBlankScan(plan_name,path=None, NumScans=1, lastNdays=0):
    #this filters for last collected Blank for specific plan_name
    if path is None:
        if lastNdays > 0:
            # if LastNdays is set, then we will ask for data from the last N days
            start_time = time.time() - (lastNdays * 86400)
            end_time = time.time()    #current time in seconds
            tz = "US/Central"
            uri = (
                f"http://{server}:{port}"
                "/api/v1/search"
                f"/{catalog}"
                f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
                "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
                f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
                f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
                f"&filter[time_range][condition][timezone]={tz}"                    # time range
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                #"&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
                #hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
                )
        else:
            # if LastNdays is not set, then we will ask for all data
            uri = (
                f"http://{server}:{port}"
                "/api/v1/search"
                f"/{catalog}"
                f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
                "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
                f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                #"&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
                #hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
                )
    else:
        if lastNdays > 0:
            # if LastNdays is set, then we will ask for data from the last N days
            start_time = time.time() - (lastNdays * 86400)
            end_time = time.time()    #current time in seconds
            tz = "US/Central"
            uri = (
                f"http://{server}:{port}"
                "/api/v1/search"
                f"/{catalog}"
                f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
                "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
                f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                "&filter[regex][condition][key]=hdf5_path"                          # filter by path
                f'&filter[regex][condition][pattern]={path}'                        # filter by path value, if path is provided
                f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
                f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
                f"&filter[time_range][condition][timezone]={tz}"                    # time range
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                #"&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
                #hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
                )
        else:
            # if LastNdays is not set, then we will ask for all data
            uri = (
                f"http://{server}:{port}"
                "/api/v1/search"
                f"/{catalog}"
                f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
                "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
                f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                "&filter[regex][condition][key]=hdf5_path"                          # filter by path
                f'&filter[regex][condition][pattern]={path}'                        # filter by path value, if path is provided
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                #"&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
                #hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
                )
                   
    logging.info(f"{uri=}")

    
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
 

def FindLastScanData(plan_name,NumScans=10, LastNdays=1):
    #print (FindLastScanData("Flyscan",10,LastNdays=1))
    #print (FindLastScanData("uascan",10,LastNdays=1))
    #print (FindLastScanData("SAXS",10,LastNdays=1))
    #print (FindLastScanData("WAXS",10,LastNdays=1))
    #print(f"Search for {plan_name=}")
    # Find all runs in a catalog between these two ISO8601 dates.
    #start_time = time.time()    #current time in seconds
    # we need to fix file not ready issue. Sometimes the last file is simply not ready 
    # when we are asking for it. Let's try to ask for files at least 30 second before now. 
    offsetTime = 20
    if plan_name == "Flyscan" or plan_name == "uascan":
        offsetTime = 95
    # this shifts the querried time by offsetTime seconds to past, providing file flush out the files. 
    end_time = time.time() - offsetTime
    tz = "US/Central"
    if LastNdays > 0:
        # if LastNdays is set, then we will ask for data from the last N days
        start_time = end_time - (LastNdays * 86400)
    
    #this filters for specific time AND for specific plan_name
    if LastNdays > 0:
        uri = (
            f"http://{server}:{port}"
            "/api/v1/search"
            f"/{catalog}"
            f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
            #"&filter[contains][condition][exit_status]"                         # filter by exist status key present
            #&filter[comparison][condition][operator]=gt&filter[comparison][condition][key]=duration&filter[comparison][condition][value]=0.1
            "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
            f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
            f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
            f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds - offsetTime
            f"&filter[time_range][condition][timezone]={tz}"                    # time range
            "&sort=-time"                                                       # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            #"&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
            #hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
            )
    else:
        # if LastNdays is not set, then we will ask for all data
        uri = (
            f"http://{server}:{port}"
            "/api/v1/search"
            f"/{catalog}"
            f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
            #"&filter[contains][condition][exit_status]"                         # filter by exist status key present
            "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
            f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
            "&sort=-time"                                                       # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            #"&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,\
            #hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}"   # select metadata
            )
          
    logging.info(f"{uri=}")
    #print(f"{uri=}")
    try:
        r = requests.get(uri).json()
        # this is now a list of Flyscan data sets
        ScanList = convert_results(r)
        #print(ScanList)
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except: 
        # url communication failed, happens and shoudl not crash anything.
        # this is workaround.   
        logging.error('Could not get data from tiled server at  usaxscontrol.xray.aps.anl.gov')
        logging.error(f"Failed {uri=}")
        return []


# Example usage of the functions

if __name__ == "__main__":
    # Example usage
    #logging.basicConfig(level=logging.INFO)
    plan_name = "Flyscan"
    scan_title = "SRM3600"
    num_scans = 5

    # Find last scans
    last_scans = FindLastScanData(plan_name, num_scans, 0)
    print(f"Last {num_scans} scans for {plan_name}: {last_scans}")

    # Find specific scan by name
    specific_scan = FindScanDataByName(plan_name, scan_title, num_scans, 0)
    print(f"Specific scan '{scan_title}' for {plan_name}: {specific_scan}")

    # Find last blank scan
    last_blank_scan = FindLastBlankScan(plan_name, None, num_scans, 0)
    print(f"Last blank scan for {plan_name}: {last_blank_scan}")
  