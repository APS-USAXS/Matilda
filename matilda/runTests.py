'''
These are various tests that can be run to check the functionality of the Matilda package.

'''

from matilda import processFlyscans, processStepscans, plotUSAXSResults, processADscans

# THIS MUST run on Linux!  
#ListOfScans = [['/home/parallels/Desktop/AJA','LEWOHtCO2wPos5_25C_464min_1443.h5']]
#ListOfBlanks = [['/home/parallels/Desktop/AJA','Blank2_40C_8min_0123.h5']]
ListOfScans = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_usaxs','GlassyCarbonM4_B_0035.h5']]
ListOfBlanks = [['C:/Users/ilavsky/Desktop/06_12_GlassyCarbon/06_12_GlassyCarbon_usaxs','AirBlank_0012.h5']]
#ListOfScans = [[r'C:\Users\ilavsky\Desktop\06_12_GlassyCarbon\06_12_GlassyCarbon_usaxs','GlassyCarbonM4_B_0035.h5']]
#ListOfBlanks = [[r'C:\Users\ilavsky\Desktop\06_12_GlassyCarbon\06_12_GlassyCarbon_usaxs','AirBlank_0012.h5']]
imagePath = '/home/parallels/Desktop/'  # Path to save images

if __name__ == "__main__":
    results = processFlyscans(ListOfScans, ListOfBlanks,recalculateAllData=False,forceFirstBlank=True)   
    #results = processStepscans(ListOfScans, ListOfBlanks,recalculateAllData=True,forceFirstBlank=True)
    #results = processADscans([['//Mac/Home/Desktop/06_15_Rakesh/06_15_Rakesh_saxs','R6016HRC_RB_V_1079.hdf'],], [['//Mac/Home/Desktop/06_15_Rakesh/06_15_Rakesh_saxs','AirBlank_1082.hdf'],], recalculateAllData=True, forceFirstBlank=True)
    print(results)
    #plotUSAXSResults(results, imagePath, isFlyscan=True) 