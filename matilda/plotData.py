#matilda/plotData.py
'''
These are plots for Matylda data analysis.

'''
import matplotlib.pyplot as plt
import pprint as pp
import logging
import numpy as np
import os



# define any globals here
default_plt_font_size = 7



def plotUSAXSResults(ListOfresults, imagePath, isFlyscan=True):  

    if imagePath is None:
        logging.warning("Image path is None, skipping plotting.")
        return
    
    # Number of data sets
    num_data_sets = len(ListOfresults)
    # Choose a colormap
    cmap = plt.get_cmap('viridis')
    # Generate colors from the colormap
    colors = [cmap(i) for i in np.linspace(0, 1, num_data_sets)]
    print(f'Got {num_data_sets} data sets to plot')

    # Get plot styling
    style = get_usaxs_r_plot_style()
    # Set the font size to specific size
    plt.rcParams['font.size'] = style["font_size"]

    # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, color in zip(range(len(ListOfresults)),colors):
        data_dict = ListOfresults[i]
        label = data_dict["RawData"]["Filename"]
        Q_array = data_dict["reducedData"]["Q"]
        UPD = data_dict["reducedData"]["Intensity"]
        plt.plot(Q_array, UPD, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle

    plt.title(style["title"])
    plt.xlabel(style["xlabel"])
    plt.ylabel(style["ylabel"])
    plt.xscale(style["xscale"])
    plt.yscale(style["yscale"])
    plt.xlim(style["xlim"])
    plt.grid(style["grid"])
    # Add legend
    plt.legend()
    # Save the plot as a JPEG image
    if isFlyscan:
        plt.savefig(os.path.join(imagePath, 'usaxs.jpg'), format='jpg', dpi=300)
    else:
        plt.savefig(os.path.join(imagePath, 'stepusaxs.jpg'), format='jpg', dpi=300) # this step scan
    #plt.show()
    plt.close()

   # Get plot styling
    style = get_usaxs_cal_plot_style()
    # Set the font size to specific size
    plt.rcParams['font.size'] = style["font_size"]

   # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, color in zip(range(len(ListOfresults)),colors):
        data_dict = ListOfresults[i]
        if data_dict["CalibratedData"]["Intensity"] is not None:
            label = data_dict["RawData"]["Filename"]
            Q_array = data_dict["CalibratedData"]["Q"]
            UPD = data_dict["CalibratedData"]["Intensity"]
            plt.plot(Q_array, UPD, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle

    plt.title(style["title"])
    plt.xlabel(style["xlabel"])
    plt.ylabel(style["ylabel"])
    plt.xscale(style["xscale"])
    plt.yscale(style["yscale"])
    plt.xlim(style["xlim"])
    plt.grid(style["grid"])
    # Add legend
    plt.legend()
    # Save the plot as a JPEG image
    if isFlyscan:
        plt.savefig(os.path.join(imagePath, 'usaxs_cal.jpg'), format='jpg', dpi=300)
    else:
        plt.savefig(os.path.join(imagePath, 'stepusaxs_cal.jpg'), format='jpg', dpi=300) # this step scan
    #plt.show()
    plt.close()



def plotSWAXSResults(ListOfresults, imagePath, isSAXS = True):  
    
    if imagePath is None:
        logging.warning("Image path is None, skipping plotting.")
        return
    
    # Number of data sets
    num_data_sets = len(ListOfresults)
    # Choose a colormap
    cmap = plt.get_cmap('viridis')
    # Generate colors from the colormap
    colors = [cmap(i) for i in np.linspace(0, 1, num_data_sets)]

    # Set the font size to specific size
    plt.rcParams['font.size'] = default_plt_font_size 

    # Plot ydata against xdata
    plt.figure(figsize=(6, 6))
    for i, color in zip(range(len(ListOfresults)),colors):
        data_dict = ListOfresults[i]
        label = data_dict["RawData"]["Filename"]
        Q_array = data_dict["reducedData"]["Q"]
        UPD = data_dict["reducedData"]["Intensity"]
        plt.plot(Q_array, UPD, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle
    plt.ylabel('Intensity')   
    if isSAXS:
        plt.title('Plot of SAXS Intensity vs. Q')   
        plt.xlabel('log(Q) [1/A]')
        plt.xscale('log')
        plt.yscale('log')
        #plt.xlim(1e-5, 1)
        plt.grid(True)
        # Add legend
        plt.legend()
        # Save the plot as a JPEG image
        plt.savefig(os.path.join(imagePath, 'saxs.jpg'), format='jpg', dpi=300)
        #plt.show()
        plt.close()

    else:       #this is WAXS data
        plt.title('Plot of WAXS Intensity vs. Q')   
        plt.xlabel('Q [1/A]')
        plt.xscale('linear')
        plt.yscale('linear')        
        #plt.xlim(1e-5, 1)
        plt.grid(True)
        # Add legend
        plt.legend()
        # Save the plot as a JPEG image
        plt.savefig(os.path.join(imagePath, 'waxs_cal.jpg'), format='jpg', dpi=300)
        #plt.show()
        plt.close()

 
    # Calibrated data plotting. 
    # Get plot styling
    style = get_usaxs_cal_plot_style()
   # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, color in zip(range(len(ListOfresults)),colors):
        data_dict = ListOfresults[i]
        if data_dict["CalibratedData"]["Intensity"] is not None:
            label = data_dict["RawData"]["Filename"]
            Q_array = data_dict["CalibratedData"]["Q"]
            Intensity = data_dict["CalibratedData"]["Intensity"]
            plt.plot(Q_array, Intensity, color=color, linestyle='-', label=label)  # You can customize the marker and linestyle

    plt.ylabel('Intensity')   
    if isSAXS:
        plt.title('Plot of SAXS Calibrated Intensity vs. Q')   
        plt.xlabel('log(Q) [1/A]')
        plt.xscale('log')
        plt.yscale('log')
        #plt.xlim(1e-5, 1)
        plt.grid(True)
        # Add legend
        plt.legend()
        # Save the plot as a JPEG image
        plt.savefig(os.path.join(imagePath, 'saxs_calib.jpg'), format='jpg', dpi=300)
        #plt.show()
        plt.close()

    else:       #this is WAXS data
        plt.title('Plot of WAXS Calibrated Intensity vs. Q')   
        plt.xlabel('Q [1/A]')
        plt.xscale('linear')
        plt.yscale('linear')        
        #plt.xlim(1e-5, 1)
        plt.grid(True)
        # Add legend
        plt.legend()
        # Save the plot as a JPEG image
        plt.savefig(os.path.join(imagePath, 'waxs_calib.jpg'), format='jpg', dpi=300)
        #plt.show()
        plt.close()




def get_usaxs_r_plot_style():
    """
    Returns a dictionary with Matplotlib styling parameters for USAXS plots.
    """
    return {
        "figsize": (6, 6),
        "title": 'Plot of Normalized Intensity vs. Q',
        "xlabel": 'log(Q) [1/A]',
        "ylabel": 'Normalized Intensity',
        "xscale": 'log',
        "yscale": 'log',
        "xlim": (1e-5, 1),
        "grid": True,
        "font_size": default_plt_font_size  # Assumes default_plt_font_size is defined globally
    }
    




def get_usaxs_cal_plot_style():
    """
    Returns a dictionary with Matplotlib styling parameters for USAXS plots.
    """
    return {
        "figsize": (6, 6),
        "title": 'Plot of Calibrated Intensity vs. Q',
        "xlabel": 'log(Q) [1/A]',
        "ylabel": 'Calibrated Intensity',
        "xscale": 'log',
        "yscale": 'log',
        "xlim": (1e-5, 1),
        "grid": True,
        "font_size": default_plt_font_size  # Assumes default_plt_font_size is defined globally
    }



