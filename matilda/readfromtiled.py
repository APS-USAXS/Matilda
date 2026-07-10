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

The three functions above return a list of [hdf5_path, hdf5_file] pairs suitable
for passing to the processXxx() functions in matilda.py.  On network failure
they return an empty list rather than raising an exception.

FindLastTuneScans(plan_name, NumScans=5, LastNdays=1)
    Return metadata dicts for the most-recent tuning scans (tune_ar/tune_mr/
    tune_a2rp).  These carry no HDF5 file; the measured arrays are fetched
    separately via tiled_get_primary_data().

tiled_get_primary_data(uid)
    Download the primary data stream (flat dict of column arrays) for one run.

Both tune helpers return empty/None on network failure so the service loop and
offline installations never crash.

Tiled server
------------
Target: http://usaxscontrol.xray.aps.anl.gov:8000
Catalog: usaxs_MongoDB
Hostname detection: if running on usaxscontrol itself, 'localhost' is used
to avoid proxy/firewall issues.

Method based on:
    https://github.com/BCDA-APS/bdp-tiled/blob/main/demo_client.ipynb
and mirrors the Igor macro logic for scan selection.

Note on duplicate filter keys
-----------------------------
Tiled accepts only ONE condition per filter type (filter[eq], filter[regex])
in a query string; a second block with the same key is ignored.  Where two
conditions of the same type are needed (plan_name + title, or title + path),
the second condition is applied client-side after the response arrives.
"""

# import necessary libraries
import requests
import datetime
import re
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
    # exit_status lets convert_results() check run success without an extra
    # per-uid metadata request (was an N+1 HTTP pattern before).
    "exit_status:stop.exit_status",
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


def _build_search_uri(NumScans, plan_name=None, title_regex=None, lastNdays=0,
                      require_exit_status=False, select_md=None):
    """Build a Tiled /search URI (single place for the query syntax).

    Parameters
    ----------
    NumScans : int
        page[limit]; must be > 0.
    plan_name : str or None
        Adds a filter[eq] condition on plan_name.
    title_regex : str or None
        Adds a filter[regex] condition on title (e.g. '(?i)blank').
    lastNdays : int
        > 0 adds a time_range filter covering the last N days; 0 = all time.
    require_exit_status : bool
        Adds filter[contains][condition][exit_status] — only runs whose stop
        document exists (i.e. finished runs).
    select_md : str or None
        select_metadata template; defaults to the module-level select_metadata.

    Notes
    -----
    * Tiled accepts only ONE condition per filter type ([eq], [regex]) in a
      query string; additional same-type conditions must be applied
      client-side by the caller (see FindScanDataByName / FindLastBlankScan).
    * Sort is always '-time' (newest first); verified to work on the current
      server both with and without a time_range filter.
    * Useful Tiled filter keywords (for reference): [noteq], [contains],
      [in], [notin], [comparison] (lt/gt/le/ge).  Working example:
      /api/v1/search/<catalog>/?page[limit]=10&filter[eq][condition][key]=plan_name
      &filter[eq][condition][value]="WAXS"&filter[regex][condition][key]=title
      &filter[regex][condition][pattern]=(?i)blank&sort=-time
    """
    parts = [
        f"http://{server}:{port}/api/v1/search/{catalog}",
        f"?page[limit]={NumScans}",
    ]
    if require_exit_status:
        parts.append("&filter[contains][condition][exit_status]")
    if plan_name is not None:
        parts.append("&filter[eq][condition][key]=plan_name")
        parts.append(f'&filter[eq][condition][value]="{plan_name}"')
    if title_regex is not None:
        parts.append("&filter[regex][condition][key]=title")
        parts.append(f"&filter[regex][condition][pattern]={title_regex}")
    if lastNdays > 0:
        end_time = time.time()
        start_time = end_time - (lastNdays * 86400)
        # server works in UTC; timezone parameter provided for completeness
        parts.append(f"&filter[time_range][condition][since]={start_time}")
        parts.append(f"&filter[time_range][condition][until]={end_time}")
        parts.append("&filter[time_range][condition][timezone]=US/Central")
    parts.append("&sort=-time")
    parts.append("&fields=metadata")
    parts.append("&omit_links=true")
    parts.append(f"&select_metadata={{{select_md if select_md is not None else select_metadata}}}")
    return "".join(parts)


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
        # Prefer exit_status delivered with the search results (via
        # select_metadata); fall back to a per-uid request only when absent.
        exit_status = md.get("exit_status")
        if exit_status is not None:
            success = (exit_status == "success")
        else:
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

    Notes
    -----
    The title match is applied CLIENT-SIDE: Tiled ignores a second
    filter[eq] block with the same condition key, so only the plan_name
    filter is sent to the server and the title is matched here.  Because
    page[limit] applies before the client-side filter, fewer than NumScans
    results may be returned when many scans share the plan but not the title.
    """
    # plan_name filtered server-side; title matched client-side below
    # (Tiled query syntax reference lives in _build_search_uri)
    uri = _build_search_uri(NumScans, plan_name=plan_name, lastNdays=lastNdays)
    logging.debug(f"{uri=}")
    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        # Client-side title filter (see Notes in docstring): keep only runs
        # whose scan_title matches exactly.
        filtered = []
        for entry in r.get("data", []):
            raw_md = entry["attributes"]["metadata"]
            md = raw_md.get("selected", raw_md)
            if md.get("scan_title") == scan_title:
                filtered.append(entry)
        r["data"] = filtered
        ScanList = convert_results(r)
        logging.info('Received expected data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except Exception:
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
        (regex, applied CLIENT-SIDE — Tiled ignores a second filter[regex]
        block with the same condition key, so the title regex is sent to the
        server and the path is matched here).  Default None (no restriction).
        Note: page[limit] applies before the client-side filter, so fewer
        than NumScans results may be returned.
    NumScans : int, optional
        Maximum number of blank scans to return.  Default 1.
    lastNdays : int, optional
        Restrict search to the last N days.  0 means no time restriction.
        Default 1.

    Returns
    -------
    list of [str, str]
        [hdf5_path, hdf5_file] pairs, empty list on network failure.
    """
    # plan_name (eq) and blank-title (regex) filtered server-side;
    # hdf5_path (when given) matched client-side below.
    uri = _build_search_uri(NumScans, plan_name=plan_name,
                            title_regex="(?i)blank", lastNdays=lastNdays)
    logging.debug(f"{uri=}")
    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        if path is not None:
            # Client-side hdf5_path filter (see docstring): Tiled cannot take
            # two filter[regex] blocks in one query.
            filtered = []
            for entry in r.get("data", []):
                raw_md = entry["attributes"]["metadata"]
                md = raw_md.get("selected", raw_md)
                if re.search(path, md.get("hdf5_path") or ""):
                    filtered.append(entry)
            r["data"] = filtered
        ScanList = convert_results(r)
        #logging.info('Received expected data from tiled server at usaxscontrol.xray.aps.anl.gov')
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except Exception:
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
    """
    # Only finished runs (exit_status present in stop doc) for this plan_name.
    uri = _build_search_uri(NumScans, plan_name=plan_name, lastNdays=LastNdays,
                            require_exit_status=True)
    logging.debug(f"{uri=}")
    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        # this is now a list of Flyscan data sets
        ScanList = convert_results(r)
        logging.info(f"Plan name: {plan_name}, list of scans:{ScanList}")
        return ScanList
    except Exception:
        # url communication failed, happens and shoudl not crash anything.
        logging.error(f'Could not get data from tiled server at  {server}')
        logging.error(f"Failed {uri=}")
        return []


def FindLastTuneScans(plan_name, NumScans=5, LastNdays=1):
    """Return metadata for the most-recent tuning scans of a given plan type.

    Used by the Matilda service to plot live tuning curves (tune_ar, tune_mr,
    tune_a2rp).  Unlike FindLastScanData, tune scans do not write HDF5 files;
    their data lives only in the Tiled/DataBroker catalog and is downloaded
    separately via tiled_get_primary_data().

    Parameters
    ----------
    plan_name : str
        Tuning plan name, e.g. 'tune_ar', 'tune_mr', 'tune_a2rp'.
    NumScans : int, optional
        Maximum number of scans to return (Tiled page[limit]).  Default 5.
    LastNdays : int, optional
        Restrict search to the last N days.  0 means all time.  Default 1.

    Returns
    -------
    list of dict
        One dict per scan, newest-first, with keys:
            'uid'      : str   — run uid (== Tiled id), used to fetch arrays
            'scan_id'  : int   — Bluesky scan_id (may be None)
            'time'     : float — POSIX start time (may be None)
            'motor'    : str   — x-axis motor name (motors[0], may be None)
            'detector' : str   — y-axis detector name (first non-scaler0, may be None)
        Returns an empty list on network failure or if the server is absent
        (offline installations), so the caller never crashes.
    """
    # Tune scans carry no hdf5 file; we need scan_id, time, motors and detectors
    # so the downloader can derive the x-axis (motor) and y-axis (detector) keys.
    tune_select_metadata = ",".join([
        "plan_name:start.plan_name",
        "time:start.time",
        "scan_id:start.scan_id",
        "uid:start.uid",
        "motors:start.motors",
        "detectors:start.detectors",
    ])

    uri = _build_search_uri(NumScans, plan_name=plan_name, lastNdays=LastNdays,
                            select_md=tune_select_metadata)
    logging.debug(f"{uri=}")
    try:
        r = requests.get(uri, timeout=TILED_TIMEOUT).json()
        results = []
        for v in range(len(r["data"])):
            uid = r["data"][v]["id"]
            raw_md = r["data"][v]["attributes"]["metadata"]
            # usaxscontrol/VM may wrap in ["selected"]; OTZ/current does not.
            md = raw_md.get("selected", raw_md)
            motors = md.get("motors") or []
            detectors = md.get("detectors") or []
            motor = motors[0] if motors else None
            # y-axis is the first detector that is not the scaler channel.
            detector = None
            for det in detectors:
                if det and not det.startswith("scaler0"):
                    detector = det
                    break
            results.append({
                "uid": md.get("uid", uid),
                "scan_id": md.get("scan_id"),
                "time": md.get("time"),
                "motor": motor,
                "detector": detector,
            })
        logging.info(f"Plan name: {plan_name}, found {len(results)} tune scans")
        return results
    except Exception:
        # url communication failed; must not crash the service loop.
        logging.error(f'Could not get tune data from tiled server at {server}')
        logging.error(f"Failed {uri=}")
        return []


def tiled_get_primary_data(uid, timeout=TILED_TIMEOUT):
    """Download the primary data stream (flat dict of arrays) for one run.

    Hits the Tiled 'node/full' endpoint, which returns the actual measured
    arrays as JSON:

        GET /api/v1/node/full/{catalog}/{uid}/primary/data?format=json

    For a tune scan this returns keys such as
    {'UPD': [...], 'scaler0_time': [...], 'a_stage_r': [...],
     'a_stage_r_user_setpoint': [...]}.

    Parameters
    ----------
    uid : str
        Run uid (Tiled id) of the scan.
    timeout : float, optional
        Request timeout in seconds.  Default TILED_TIMEOUT.

    Returns
    -------
    dict or None
        Mapping of column name to list of values, or None on failure.  Never
        raises, so callers on the service loop are safe against a dead server.
    """
    uri = (
        f"http://{server}:{port}"
        f"/api/v1/node/full/{catalog}/{uid}/primary/data?format=json"
    )
    logging.debug(f"{uri=}")
    try:
        r = requests.get(uri, timeout=timeout)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, dict):
            return data
        logging.error(f"Unexpected primary-data payload for {uid}: not a dict")
        return None
    except Exception:
        logging.error(f'Could not get primary data for {uid} from tiled at {server}')
        logging.error(f"Failed {uri=}")
        return None


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
  