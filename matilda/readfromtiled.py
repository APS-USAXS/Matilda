"""
readfromtiled.py
================
Query the USAXS Tiled server for scan metadata and file paths.

Versions
--------
0.2  2025-04-15
0.3  2025-06-01

Public API
----------
FindLastScanData(plan_name, NumScans=10, LastNdays=1)
    Return the most-recent N completed scans of a given Bluesky plan type
    from the last N days.

FindScanDataByName(plan_name, scan_title, NumScans=1, lastNdays=0)
    Return scans matching both a plan name and a scan_title substring.

FindLastBlankScan(plan_name, path=None, NumScans=1, lastNdays=1)
    Return the most-recent blank/background scans for a given plan type,
    optionally restricted to a specific file-system path.

All three functions return a list of [hdf5_path, hdf5_file] pairs suitable
for passing to the processXxx() functions in matilda.py.  On network failure
they return an empty list rather than raising an exception.

Tiled server
------------
Target: http://usaxscontrol.xray.aps.anl.gov:8000
Catalog: usaxs_MongoDB
Hostname detection: if running on usaxscontrol itself, 'localhost' is used
to avoid proxy/firewall issues.

Method based on:
    https://github.com/BCDA-APS/bdp-tiled/blob/main/demo_client.ipynb
and mirrors the Igor macro logic for scan selection.

TODO: bare except clauses (lines ~226, ~346, ~425) should be narrowed to
      'except Exception' to avoid silently swallowing KeyboardInterrupt.
TODO: debug print() calls should be replaced with logging.debug().
"""

# import necessary libraries
import requests
import json
import datetime
import time
import socket
import logging
from typing import Any, Optional


def iso_to_ts(isotime):
    """Convert an ISO-8601 string to a POSIX timestamp (float seconds)."""
    return datetime.datetime.fromisoformat(isotime).timestamp()

def ts_to_iso(time):
    """Convert a POSIX timestamp (float seconds) to a local ISO-8601 string."""
    return datetime.datetime.fromtimestamp(time).isoformat()

current_hostname = socket.gethostname()
if current_hostname in ('usaxscontrol', 'usaxscontrol.xray.aps.anl.gov'):
    server = "localhost"   # avoid proxy/firewall issues when running on the same machine
else:
    server = "usaxscontrol.xray.aps.anl.gov"

port = 8000
#catalog = "raw"
catalog = "usaxs_MongoDB"
TILED_TIMEOUT = 10  # seconds

select_metadata = ",".join([
    "plan_name:start.plan_name",
    "time:start.time",
    "scan_title:start.plan_args.scan_title",
    "hdf5_file:start.hdf5_file",
    "hdf5_path:start.hdf5_path",
])

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
            "sort": "-metadata.start.time",
            #"sort": "-time",
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
    """Print a one-line summary of the first and last run in a Tiled response dict.

    Parameters
    ----------
    r : dict
        Raw JSON response from tiled_get(); must contain a 'data' list.

    Notes
    -----
    Helper used during interactive debugging / development; not called by the
    main processing loop.
    """
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
    """Convert a raw Tiled search response into a list of [path, filename] pairs.

    Skips any run that did not finish successfully (exit_status != 'success').
    Skips runs where hdf5_file is None (e.g. non-file-writing plans).

    Parameters
    ----------
    r : dict
        Raw JSON response from tiled_get() or requests.get().json().

    Returns
    -------
    list of [str, str]
        Each element is [hdf5_path, hdf5_file] for a successfully completed run.
    """
    OutputList=[]
    for v in range(len(r["data"])):
        uid = r["data"][v]["id"]
        raw_md = r["data"][v]["attributes"]["metadata"]
        md = raw_md.get("selected", raw_md)   # usaxscontrol/VM has ["selected"]; OTZ does not
        success = successful_run(uid)
        #if not success and (md["plan_name"] == "Flyscan"):
        if not success :
            tempPlanName=md["plan_name"]
            logging.info(f"Skipping unfinished/aborted scan: {uid} of type : {tempPlanName}")
            continue

        #started = ts_to_iso(md["time"])
        hdf5_file = md["hdf5_file"]
        hdf5_path = md["hdf5_path"]
        #print(f" path: {hdf5_path=} {hdf5_file=}")
        if hdf5_file is not None:
            OutputList.append([hdf5_path,hdf5_file])
    return OutputList
        
#print(f'Search of {catalog=} has {len(r["data"])} runs.')
#print_results_summary(r)


def FindScanDataByName(plan_name, scan_title, NumScans=1, lastNdays=1):
    """Return scans matching both a Bluesky plan name and a scan title.

    Parameters
    ----------
    plan_name : str
        Bluesky plan name to filter on (e.g. 'Flyscan', 'SAXS', 'WAXS', 'uascan').
    scan_title : str
        Exact scan title string to match (stored in start.plan_args.scan_title).
    NumScans : int, optional
        Maximum number of results to return (Tiled page limit). Default 1.
    lastNdays : int, optional
        Restrict search to the last N days.  0 means no time restriction.
        Default 1.

    Returns
    -------
    list of [str, str]
        [hdf5_path, hdf5_file] pairs, empty list on network failure.

    TODO: the two filter[eq] blocks for plan_name and title share the same
          Tiled filter key — the second silently overwrites the first.
          This is a known Tiled API quirk; the title filter may not work.
    TODO: debug print(uri) on line ~217 should be logging.debug().
    """
    #this filters for specific time AND for specific plan_name
    # select_metadata = ",".join([
    #     "plan_name:start.plan_name",
    #     "time:start.time",
    #     "scan_title:start.plan_args.scan_title",
    #     "hdf5_file:start.hdf5_file",
    #     "hdf5_path:start.hdf5_path",
    # ])
    if lastNdays > 0:
        # if LastNdays is set, then we will ask for data from the last N days
        start_time = time.time() - (lastNdays * 86400)
        end_time = time.time()      # current time in seconds
        tz = "US/Central"
        # server works in UTC, no need to provide timezone, but we need to add offset if needed.
        #offset = 6*60*60  # US/Central is UTC-6
        #start_time += offset
        #end_time += offset
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
            "&sort=-time"                                                        # sort by time, -time gives last scans first
            #f"&filter[comparison][condition][key]=start.time&filter[comparison][condition][operator]=ge&filter[comparison][condition][value]={(start_time)}"             # time range, start time - 24 hours from now
            #f"&filter[comparison2][condition][key]=start.time&filter[comparison2][condition][operator]=le&filter[comparison2][condition][value]={end_time}"                 # time range, current time in seconds
            #"&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            f"&select_metadata={{{select_metadata}}}"                               # select metadata
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
            "&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            f"&select_metadata={{{select_metadata}}}"                               # select metadata
            )
      
    #logging.info(f"{uri=}")
    #additional keywords:
    #[noteq] - not equal
    #[contains] - seems same as eq in use, and cannot be made into case insensitive. Not useful. 
    #[in] - in a list of values
    #[notin] - not in a list of values
    #[comparison] - comparison with lt, gt, le, ge for numerical values
    #working examples:
    #http://10.211.55.7:8000/api/v1/search/usaxs/?page[limit]=10&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)blank&sort=-time
    #returns list of "Blank" samples, not not ist of samples contains "blank" in name
    #http://10.211.55.7:8000/api/v1/search/usaxs/?page[limit]=1&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)water*blank&sort=-time
    #returns last scan which conatins case independent "water blank" in name
    #http://10.211.55.7:8000/api/v1/search/usaxs/?page[limit]=1&filter[eq][condition][key]=plan_name&filter[eq][condition][value]=%22WAXS%22&filter[regex][condition][key]=title&filter[regex][condition][pattern]=(?i)blank&sort=-time&omit_links=true&select_metadata={plan_name:start.plan_name,time:start.time,scan_title:start.plan_args.scan_title,hdf5_file:start.hdf5_file,hdf5_path:start.hdf5_path}
    #returns last scan which conatisn case independet "water blank" in name
    print(uri)
    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        #logging.info(f"Got json for : {plan_name}")        #this does not work for some reason? 
        ScanList = convert_results(r)
        #ScanList=[]
        logging.info('Received expected data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except: 
        # url communication failed, happens and should not crash anything.
        logging.error(f'Could not get data from tiled server at {server}')
        logging.error(f"Failed {uri=}")
        return []
    

def FindLastBlankScan(plan_name, path=None, NumScans=1, lastNdays=1):
    """Return the most-recent blank/background scans for a given plan type.

    Searches the Tiled catalog for scans whose title matches the regex
    ``(?i)blank`` (case-insensitive), i.e. any scan whose title contains
    the word "blank".

    Parameters
    ----------
    plan_name : str
        Bluesky plan name (e.g. 'Flyscan', 'SAXS', 'WAXS', 'uascan').
    path : str or None, optional
        If provided, additionally filter by hdf5_path matching this string
        (regex-matched by Tiled).  Default None (no path restriction).
    NumScans : int, optional
        Maximum number of blank scans to return.  Default 1.
    lastNdays : int, optional
        Restrict search to the last N days.  0 means no time restriction.
        Default 1.

    Returns
    -------
    list of [str, str]
        [hdf5_path, hdf5_file] pairs, empty list on network failure.

    TODO: debug print(uri) should be logging.debug().
    """
    #this filters for last collected Blank for specific plan_name
    if path is None:
        if lastNdays > 0:
            # if LastNdays is set, then we will ask for data from the last N days
            start_time = time.time() - (lastNdays * 86400)
            end_time = time.time()    #current time in seconds
            tz = "US/Central"
            # server works in UTC, no need to provide timezone, but we need to add offset if needed.
            #offset = 6*60*60  # US/Central is UTC-6
            #start_time += offset
            #end_time += offset
            uri = (
                f"http://{server}:{port}"
                "/api/v1/search"
                f"/{catalog}"
                f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
                "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
                f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
                #f'&filter[full_text][condition][text]={plan_name}'                   # filter by plan_name value, full text search, should be faster than eq   
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
                f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
                f"&filter[time_range][condition][timezone]={tz}"                    # time range
                #f"&filter[comparison][condition][key]=start.time&filter[comparison][condition][operator]=ge&filter[comparison][condition][value]={(start_time)}"             # time range, start time - 24 hours from now
                #f"&filter[comparison2][condition][key]=start.time&filter[comparison2][condition][operator]=le&filter[comparison2][condition][value]={end_time}"                 # time range, current time in seconds
                #"&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                f"&select_metadata={{{select_metadata}}}"                               # select metadata
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
                #f'&filter[full_text][condition][text]={plan_name}'                   # filter by plan_name value, full text search, should be faster than eq   
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                "&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                f"&select_metadata={{{select_metadata}}}"                               # select metadata
                )
    else:
        if lastNdays > 0:
            # if LastNdays is set, then we will ask for data from the last N days
            start_time = time.time() - (lastNdays * 86400)
            end_time = time.time()    #current time in seconds
            tz = "US/Central"
            # server works in UTC, no need to provide timezone, but we need to add offset if needed.
            #offset = 6*60*60  # US/Central is UTC-6
            #start_time += offset
            #end_time += offset            
            uri = (
                f"http://{server}:{port}"
                "/api/v1/search"
                f"/{catalog}"
                f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
                "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
                f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
                #f'&filter[full_text][condition][text]={plan_name}'                   # filter by plan_name value, full text search, should be faster than eq   
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                "&filter[regex][condition][key]=hdf5_path"                          # filter by path
                f'&filter[regex][condition][pattern]={path}'                        # filter by path value, if path is provided
                f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
                f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
                f"&filter[time_range][condition][timezone]={tz}"                    # time range
                #f"&filter[comparison][condition][key]=start.time&filter[comparison][condition][operator]=ge&filter[comparison][condition][value]={(start_time)}"             # time range, start time - 24 hours from now
                #f"&filter[comparison2][condition][key]=start.time&filter[comparison2][condition][operator]=le&filter[comparison2][condition][value]={end_time}"                 # time range, current time in seconds
                #"&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                f"&select_metadata={{{select_metadata}}}"                               # select metadata
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
                #f'&filter[full_text][condition][text]={plan_name}'                   # filter by plan_name value, full text search, should be faster than eq   
                "&filter[regex][condition][key]=title"                              # filter by title
                f'&filter[regex][condition][pattern]=(?i)blank'                     # filter by title value
                "&filter[regex][condition][key]=hdf5_path"                          # filter by path
                f'&filter[regex][condition][pattern]={path}'                        # filter by path value, if path is provided
                "&sort=-time"                                                       # sort by time, -time gives last scans first
                #"&sort=-metadata.start.time"                                                       # sort by time, -time gives last scans first
                "&fields=metadata"                                                  # return metadata
                "&omit_links=true"                                                  # no links
                f"&select_metadata={{{select_metadata}}}"                               # select metadata
                )
                   
    #logging.info(f"{uri=}")
    print(uri)

    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        ScanList = convert_results(r)
        #logging.info('Received expected data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except: 
        # url communication failed, happens and shoudl not crash anything.
        logging.error(f'Could not get data from tiled server at  {server}')
        logging.error(f"Failed {uri=}")
        return []
 

def FindLastScanData(plan_name, NumScans=10, LastNdays=1):
    """Return the most-recent completed scans for a given Bluesky plan type.

    This is the primary entry point called by the matilda main loop every 15 s.
    Only runs with exit_status == 'success' are included (via convert_results).

    Parameters
    ----------
    plan_name : str
        Bluesky plan name to filter on.  Known values used by matilda:
        'Flyscan', 'uascan', 'SAXS', 'WAXS'.
    NumScans : int, optional
        Maximum number of scans to return (Tiled page[limit]).  Default 10.
    LastNdays : int, optional
        Restrict search to the last N days.  0 means all time.  Default 1.

    Returns
    -------
    list of [str, str]
        [hdf5_path, hdf5_file] pairs sorted newest-first by Tiled (sort=-time).
        Returns an empty list on network failure.

    Notes
    -----
    * A small time offset that was previously applied to work around file-flush
      latency has been removed (see commented-out offsetTime code).
    * The filter[contains][condition][exit_status] clause is a Tiled-specific
      filter that checks whether the 'exit_status' key exists in the stop doc.
    TODO: debug print(f"{uri=}") should be logging.debug().
    """
    #print (FindLastScanData("Flyscan",10,LastNdays=1))
    #print (FindLastScanData("uascan",10,LastNdays=1))
    #print (FindLastScanData("SAXS",10,LastNdays=1))
    #print (FindLastScanData("WAXS",10,LastNdays=1))
    #print(f"Search for {plan_name=}")
    # Find all runs in a catalog between these two ISO8601 dates.
    start_time = 0
    # we need to fix file not ready issue. Sometimes the last file is simply not ready 
    # when we are asking for it. Let's try to ask for files at least 30 second before now. 
    #offsetTime = 20
    ##if plan_name == "Flyscan" or plan_name == "uascan":
     #   offsetTime = 95
    # this shifts the querried time by offsetTime seconds to past, providing file flush out the files. 
    end_time = time.time() #- offsetTime
    tz = "US/Central"
    # server works in UTC, no need to provide timezone, but we need to add offset if needed.
    #offset = 6*60*60  # US/Central is UTC-6
    #end_time += offset    
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
            "&filter[contains][condition][exit_status]"                         # filter by exist status key present
            #&filter[comparison][condition][operator]=gt&filter[comparison][condition][key]=duration&filter[comparison][condition][value]=0.1
            "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
            f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
            #f'&filter[full_text][condition][text]={plan_name}'                   # filter by plan_name value, full text search, should be faster than eq   
            f"&filter[time_range][condition][since]={(start_time)}"             # time range, start time - 24 hours from now
            f"&filter[time_range][condition][until]={end_time}"                 # time range, current time in seconds
            f"&filter[time_range][condition][timezone]={tz}"                    # time range
            #f"&filter[comparison][condition][key]=start.time&filter[comparison][condition][operator]=ge&filter[comparison][condition][value]={(start_time)}"             # time range, start time - 24 hours from now
            #f"&filter[comparison2][condition][key]=start.time&filter[comparison2][condition][operator]=le&filter[comparison2][condition][value]={end_time}"                 # time range, current time in seconds
            "&sort=-time"                                                      # sort by time, -time gives last scans first
            #"&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            f"&select_metadata={{{select_metadata}}}"                               # select metadata
            )
    else:
        # if LastNdays is not set, then we will ask for all data
        uri = (
            f"http://{server}:{port}"
            "/api/v1/search"
            f"/{catalog}"
            f"?page[limit]={NumScans}"                                          # 0: all matching, 10 is 10 scans. Must be >0 value
            "&filter[contains][condition][exit_status]"                         # filter by exist status key present
            "&filter[eq][condition][key]=plan_name"                             # filter by plan_name
            f'&filter[eq][condition][value]="{plan_name}"'                      # filter by plan_name value
            #f'&filter[full_text][condition][text]={plan_name}'                   # filter by plan_name value, full text search, should be faster than eq   
            #"&sort=-metadata.start.time"                                        # sort by time, -time gives last scans first
            "&sort=-time"                                        # sort by time, -time gives last scans first
            "&fields=metadata"                                                  # return metadata
            "&omit_links=true"                                                  # no links
            f"&select_metadata={{{select_metadata}}}"                               # select metadata
            )
          
    #logging.info(f"{uri=}")
    print(f"{uri=}")
    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        # this is now a list of Flyscan data sets
        ScanList = convert_results(r)
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except: 
        # url communication failed, happens and shoudl not crash anything.
        logging.error(f'Could not get data from tiled server at  {server}')
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
  