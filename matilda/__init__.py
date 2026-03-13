"""
matilda
=======
Live data-processing package for the APS USAXS/SAXS/WAXS instrument.

Typical use: run matilda/matilda.py as a script.  It polls the Tiled server
every 15 s and reduces any new scans to calibrated 1-D I(Q) data.

Public processing functions (importable for use from notebooks or Igor)
-----------------------------------------------------------------------
from matilda.matilda import processFlyscans, processStepscans, processADscans
from matilda.matilda import processUSAXSFolder

from matilda.convertFlyscan  import processFlyscan
from matilda.convertUSAXS    import processStepscan
from matilda.convertSWAXS    import process2Ddata

from matilda.readfromtiled   import FindLastScanData, FindLastBlankScan
from matilda.hdf5code        import saveNXcanSAS, readMyNXcanSAS
"""
