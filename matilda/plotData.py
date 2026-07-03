"""
plotData.py
===========
Matplotlib-based plotting routines for Matilda data analysis.

Generates JPEG summary plots and writes them to a web-visible directory
(imagePath) for live monitoring.  Called after each processing cycle by
the main loop in matilda.py.

Current outputs
---------------
usaxs.jpg / stepusaxs.jpg         — raw USAXS/step-scan I vs Q
usaxs_cal.jpg / stepusaxs_cal.jpg — calibrated USAXS/step-scan I vs Q
saxs.jpg                          — raw SAXS I vs Q (log-log)
saxs_cal.jpg                      — calibrated SAXS I vs Q (log-log)
waxs.jpg                          — raw WAXS I vs Q (linear)
waxs_cal.jpg                      — calibrated WAXS I vs Q (linear)
tune_ar.jpg / tune_mr.jpg / tune_a2rp.jpg — live tuning curves (detector vs motor)

GUI transition note
-------------------
This module uses matplotlib for headless (file-only) output.  Future GUI
work should use pyqtgraph for interactive display.  matplotlib may be kept
for file export or replaced entirely depending on requirements.

TODO: plotUSAXSResults has an off-by-one indentation on the second plot
      block (lines starting with '   # Get plot styling' after first plt.close()).
      Functionally correct but visually misleading.
"""
import matplotlib.pyplot as plt
import logging
import os
import datetime



# define any globals here
default_plt_font_size = 7

# Up to 10 datasets per plot (enforced upstream).
# tab10 gives 10 maximally-distinct colors; line styles add a second
# visual channel so colorblind users can still tell datasets apart.
PLOT_COLORS = [plt.get_cmap('tab10')(i / 10) for i in range(10)]
PLOT_LINESTYLES = ['-', '--', '-.', ':']



def plotUSAXSResults(ListOfresults, imagePath, isFlyscan=True):
    """Save USAXS / step-scan summary plots to JPEG files.

    Produces two JPEG files per call:
    * Raw normalised I vs Q  (usaxs.jpg or stepusaxs.jpg)
    * Calibrated I vs Q      (usaxs_cal.jpg or stepusaxs_cal.jpg)

    Data sets use the tab10 colormap (10 distinct colors) with cycling line
    styles (solid, dashed, dash-dot, dotted) for colorblind accessibility.
    Y-axis is clamped to at most 14 decades below the maximum to avoid empty
    log plots from outlier points.

    Parameters
    ----------
    ListOfresults : list of dict
        Each dict is the result of processFlyscan() or processStepscan().
        Required keys: RawData.filename, reducedData.Q, reducedData.Intensity,
        CalibratedData.Q, CalibratedData.Intensity (None if no blank).
    imagePath : str or None
        Directory to write JPEG files into.  If None, plotting is skipped.
    isFlyscan : bool, optional
        True  → save as usaxs*.jpg (flyscan).
        False → save as stepusaxs*.jpg (step scan).
        Default True.
    """
    if imagePath is None:
        logging.warning("Image path is None, skipping plotting.")
        return
    
    # Number of data sets
    num_data_sets = len(ListOfresults)
    logging.info(f'Got {num_data_sets} USAXS data sets to plot')

    # Get plot styling
    style = get_usaxs_r_plot_style()
    # Set the font size to specific size
    plt.rcParams['font.size'] = style["font_size"]

    # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, data_dict in enumerate(ListOfresults):
        label = data_dict["RawData"]["filename"]
        Q_array = data_dict["reducedData"]["Q"]
        UPD = data_dict["reducedData"]["Intensity"]
        plt.plot(Q_array, UPD, color=PLOT_COLORS[i % 10], linestyle=PLOT_LINESTYLES[i % 4], label=label)

    plt.title(style["title"])
    plt.xlabel(style["xlabel"])
    plt.ylabel(style["ylabel"])
    plt.xscale(style["xscale"])
    plt.yscale(style["yscale"])
    #plt.xlim(style["xlim"])
    plt.xlim()          #autoscales
    plt.ylim()          #autoscales
    plt.grid(style["grid"])
    # Add legend
    plt.legend()
    #limit the maximum range of the plot to at most 10 decades from max:
    # Get the current limits
    xlim = plt.xlim()   #returns a tuple of the form (xmin, xmax)
    ylim = plt.ylim()   #returns a tuple of the form (ymin, ymax)
    # Calculate the new limits,for y we want to keep the max and limit the min to larger of existing value and max/1e12
    # for x we want to keep min or force 1e-5 (440 resolution) is min is smaller. 
    new_xlim = (max(xlim[0], 1e-5), xlim[1])
    new_ylim = (max(ylim[0], ylim[1] / 10**14), ylim[1])
    # Set the new limits
    plt.xlim(new_xlim)
    plt.ylim(new_ylim)
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
    for i, data_dict in enumerate(ListOfresults):
        if data_dict["CalibratedData"]["Intensity"] is not None:
            label = data_dict["RawData"]["filename"]
            Q_array = data_dict["CalibratedData"]["Q"]
            UPD = data_dict["CalibratedData"]["Intensity"]
            plt.plot(Q_array, UPD, color=PLOT_COLORS[i % 10], linestyle=PLOT_LINESTYLES[i % 4], label=label)

    plt.title(style["title"])
    plt.xlabel(style["xlabel"])
    plt.ylabel(style["ylabel"])
    plt.xscale(style["xscale"])
    plt.yscale(style["yscale"])
    plt.xlim()          #autoscales
    plt.ylim()          #autoscales
    plt.grid(style["grid"])
    # Add legend
    plt.legend()
    #limit the maximum range of the plot to at most 10 decades from max:
    # Get the current limits
    xlim = plt.xlim()   #returns a tuple of the form (xmin, xmax)
    ylim = plt.ylim()   #returns a tuple of the form (ymin, ymax)
    # Calculate the new limits,for y we want to keep the max and limit the min to larger of existing value and max/1e12
    # for x we want to keep min or force 1e-5 (440 resolution) is min is smaller. 
    new_xlim = (max(xlim[0], 1e-5), xlim[1])
    new_ylim = (max(ylim[0], ylim[1] / 10**14), ylim[1])    # Set the new limits
    plt.xlim(new_xlim)
    plt.ylim(new_ylim)
    # Save the plot as a JPEG image
    if isFlyscan:
        plt.savefig(os.path.join(imagePath, 'usaxs_cal.jpg'), format='jpg', dpi=300)
    else:
        plt.savefig(os.path.join(imagePath, 'stepusaxs_cal.jpg'), format='jpg', dpi=300) # this step scan
    #plt.show()
    plt.close()



def plotSWAXSResults(ListOfresults, imagePath, isSAXS=True):
    """Save SAXS or WAXS summary plots to JPEG files.

    Produces two JPEG files per call:
    * Raw I vs Q           (saxs.jpg or waxs.jpg)
    * Calibrated I vs Q    (saxs_cal.jpg or waxs_cal.jpg)

    SAXS plots use log-log axes and limit the display range to 4 decades
    below the maximum.  WAXS plots use linear axes with no range clamping.
    Calibrated plots that have no data (CalibratedData.Intensity is None)
    are silently skipped per scan.

    Parameters
    ----------
    ListOfresults : list of dict
        Each dict is the result of process2Ddata().  Required keys:
        RawData.filename, reducedData.Q, reducedData.Intensity,
        CalibratedData.Q, CalibratedData.Intensity (None if no blank).
    imagePath : str or None
        Directory to write JPEG files into.  If None, plotting is skipped.
    isSAXS : bool, optional
        True  → SAXS (log-log, saves saxs*.jpg).
        False → WAXS (linear, saves waxs*.jpg).
        Default True.
    """
    if imagePath is None:
        logging.warning("Image path is None, skipping plotting.")
        return
    
    # Set the font size to specific size
    plt.rcParams['font.size'] = default_plt_font_size 

    # Plot ydata against xdata
    plt.figure(figsize=(6, 6))
    for i, data_dict in enumerate(ListOfresults):
        label = data_dict["RawData"]["filename"]
        Q_array = data_dict["reducedData"]["Q"]
        UPD = data_dict["reducedData"]["Intensity"]
        plt.plot(Q_array, UPD, color=PLOT_COLORS[i % 10], linestyle=PLOT_LINESTYLES[i % 4], label=label)
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
        #limit the maximum range of the plot to at most 10 decades from max:
        # Get the current limits
        xlim = plt.xlim()   #returns a tuple of the form (xmin, xmax)
        ylim = plt.ylim()   #returns a tuple of the form (ymin, ymax)
        # Calculate the new limits, we want to keep tha max and limit the min to larger of existing value and max/1e10
        new_xlim = (max(xlim[0], xlim[1] / 10**4), xlim[1])
        new_ylim = (max(ylim[0], ylim[1] / 10**4), ylim[1])
        # Set the new limits
        plt.xlim(new_xlim)
        plt.ylim(new_ylim)
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
        plt.savefig(os.path.join(imagePath, 'waxs.jpg'), format='jpg', dpi=300)
        #plt.show()
        plt.close()

 
    # Calibrated data plotting. 
    # Get plot styling
    style = get_usaxs_cal_plot_style()
   # Plot ydata against xdata
    plt.figure(figsize=style["figsize"])
    for i, data_dict in enumerate(ListOfresults):
        if data_dict["CalibratedData"]["Intensity"] is not None:
            label = data_dict["RawData"]["filename"]
            Q_array = data_dict["CalibratedData"]["Q"]
            Intensity = data_dict["CalibratedData"]["Intensity"]
            plt.plot(Q_array, Intensity, color=PLOT_COLORS[i % 10], linestyle=PLOT_LINESTYLES[i % 4], label=label)

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
        #limit the maximum range of the plot to at most 10 decades from max:
        # Get the current limits
        xlim = plt.xlim()   #returns a tuple of the form (xmin, xmax)
        ylim = plt.ylim()   #returns a tuple of the form (ymin, ymax)
        # Calculate the new limits, we want to keep tha max and limit the min to larger of existing value and max/1e10
        new_xlim = (max(xlim[0], xlim[1] / 10**5), xlim[1])
        new_ylim = (max(ylim[0], ylim[1] / 10**5), ylim[1])
        # Set the new limits
        plt.xlim(new_xlim)
        plt.ylim(new_ylim)
        # Save the plot as a JPEG image
        plt.savefig(os.path.join(imagePath, 'saxs_cal.jpg'), format='jpg', dpi=300)
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
        plt.savefig(os.path.join(imagePath, 'waxs_cal.jpg'), format='jpg', dpi=300)
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


def plotTuneResults(ListOfTuneResults, imagePath, plan_name):
    """Save a live tuning-curve plot (detector counts vs motor position) to JPEG.

    Overlays the last N tune scans of a single tuning plan on one linear plot,
    newest-first, using the same tab10 colors and cycling line styles as the
    other Matilda plots.  Written as ``<plan_name>.jpg`` (tune_ar.jpg,
    tune_mr.jpg, tune_a2rp.jpg) into *imagePath*.

    Parameters
    ----------
    ListOfTuneResults : list of dict
        Result dicts from convertTune.getTuneResults(). Required keys:
        'x', 'y' (numpy arrays), 'motor', 'detector', 'scan_id', 'time'.
    imagePath : str or None
        Directory to write the JPEG into.  If None, plotting is skipped.
    plan_name : str
        Tuning plan name; also the output filename stem (e.g. 'tune_ar').
    """
    if imagePath is None:
        logging.warning("Image path is None, skipping tune plotting.")
        return

    if not ListOfTuneResults:
        logging.info(f"No tune data to plot for {plan_name}")
        return

    logging.info(f"Got {len(ListOfTuneResults)} {plan_name} tune curves to plot")

    # Set the font size to specific size
    plt.rcParams['font.size'] = default_plt_font_size

    plt.figure(figsize=(6, 6))
    motor_label = plan_name          # sensible fallbacks if metadata is sparse
    detector_label = 'Counts'
    for i, data_dict in enumerate(ListOfTuneResults):
        x = data_dict["x"]
        y = data_dict["y"]
        motor_label = data_dict.get("motor") or motor_label
        detector_label = data_dict.get("detector") or detector_label
        # Legend label: scan_id + local HH:MM:SS when available.
        scan_id = data_dict.get("scan_id")
        tstamp = data_dict.get("time")
        if tstamp is not None:
            timestr = datetime.datetime.fromtimestamp(tstamp).strftime('%H:%M:%S')
        else:
            timestr = ""
        label = f"#{scan_id} {timestr}".strip()
        plt.plot(x, y, color=PLOT_COLORS[i % 10], linestyle=PLOT_LINESTYLES[i % 4], label=label)

    plt.title(f'Tune: {plan_name}')
    plt.xlabel(motor_label)
    plt.ylabel(detector_label)
    plt.xscale('linear')
    plt.yscale('linear')
    plt.grid(True)
    plt.legend()
    # Save the plot as a JPEG image named after the plan.
    plt.savefig(os.path.join(imagePath, f'{plan_name}.jpg'), format='jpg', dpi=300)
    plt.close()
