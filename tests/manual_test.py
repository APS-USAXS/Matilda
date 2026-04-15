"""
tests/manual_test.py
====================
Manual integration test script — run by hand to verify the full data
reduction pipeline against local test data files.

This is NOT an automated pytest test.  It requires real HDF5 data files
on the local machine and must be run explicitly:

    conda activate matilda
    python tests/manual_test.py

Originally: matilda/runTests.py (kept there for backwards compatibility
with anyone who calls it directly from the matilda/ directory).

TODO: Replace with proper pytest fixtures that use the bundled TestData/
      files so CI can run automatically without manual file paths.
"""

from matilda.matilda import processFlyscans, processStepscans, processADscans

# Adjust paths to your local copies of the test HDF5 files.
# Linux example (parallels VM):
#   ListOfScans = [['/home/parallels/Desktop/AJA','LEWOHtCO2wPos5_25C_464min_1443.h5']]
#   ListOfBlanks = [['/home/parallels/Desktop/AJA','Blank2_40C_8min_0123.h5']]

ListOfScans = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_usaxs',
                'GlassyCarbonM4_B_0035.h5']]
ListOfBlanks = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_usaxs',
                 'AirBlank_0012.h5']]
ListOfScans2D = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_saxs',
                  'GlassyCarbonM4_B_0035.hdf']]
ListOfBlanks2D = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_saxs',
                   'AirBlank_0012.hdf']]

if __name__ == "__main__":
    results = processFlyscans(ListOfScans, ListOfBlanks, recalculateAllData=False, forceFirstBlank=True)
    results = processADscans(ListOfScans2D, ListOfBlanks2D, recalculateAllData=False, forceFirstBlank=True)
    print(results)
