'''
These are various tests that can be run to check the functionality of the Matilda package.

'''

from matilda import processFlyscans, processStepscans, plotUSAXSResults

# THIS MUST run on Linux!  
#ListOfScans = [['/home/parallels/Desktop/AJA','LEWOHtCO2wPos5_25C_464min_1443.h5']]
#ListOfBlanks = [['/home/parallels/Desktop/AJA','Blank2_40C_8min_0123.h5']]
ListOfScans = [['/home/parallels/Desktop/step','X15_S_CaP_Sotiriou_0351.h5']]
ListOfBlanks = [['/home/parallels/Desktop/step','Blank_Yavitt_0343.h5']]
imagePath = '/home/parallels/Desktop/'  # Path to save images

if __name__ == "__main__":
    #results = processFlyscans(ListOfScans, ListOfBlanks,recalculateAllData=True,forceFirstBlank=True)   
    results = processStepscans(ListOfScans, ListOfBlanks,recalculateAllData=True,forceFirstBlank=True)
    #print(results)
    plotUSAXSResults(results, imagePath, isFlyscan=True) 