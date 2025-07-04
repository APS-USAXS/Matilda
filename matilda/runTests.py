'''
These are various tests that can be run to check the functionality of the Matilda package.

'''

from matilda import processFlyscans, processStepscans, plotUSAXSResults, processADscans

# THIS MUST run on Linux!  
#ListOfScans = [['/home/parallels/Desktop/AJA','LEWOHtCO2wPos5_25C_464min_1443.h5']]
#ListOfBlanks = [['/home/parallels/Desktop/AJA','Blank2_40C_8min_0123.h5']]
ListOfScans = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_usaxs','GlassyCarbonM4_B_0035.h5']]
ListOfBlanks = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_usaxs','AirBlank_0012.h5']]
ListOfScans2D = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_saxs','GlassyCarbonM4_B_0035.hdf']]
ListOfBlanks2D = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_saxs','AirBlank_0012.hdf']]
imagePath = '/home/parallels/Desktop/'  # Path to save images

if __name__ == "__main__":
    results = processFlyscans(ListOfScans, ListOfBlanks,recalculateAllData=False,forceFirstBlank=True)   
    #results = processStepscans(ListOfScans, ListOfBlanks,recalculateAllData=True,forceFirstBlank=True)
    results = processADscans(ListOfScans2D,ListOfBlanks2D, recalculateAllData=False, forceFirstBlank=True)
    print(results)
    #plotUSAXSResults(results, imagePath, isFlyscan=True) 