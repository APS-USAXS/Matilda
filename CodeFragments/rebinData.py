
'''
    Converted by AI from Igor code

    rebin_data() routine will rebin data from Q=0.0002 higher using rebin_log_data(), min step is set to distacne between Q[0] and Q[1]
    
    This routine will rebin data on log scale. It will produce new Wx and Wy with new NumberOfPoints
    If MinStep > 0 it will try to set the values so the minimum step on log scale is MinStep
    optional Wsdev is standard deviation for each Wy value, it will be propagated through - sum(sdev^2)/numpnts in each bin. 
    optional Wxwidth will generate width of each new bin in x. NOTE: the edge is half linear distance between the two points, no log  
    skewing is done for edges. Therefore the width is really half of the distance between p-1 and p+1 points.  
    optional W1-5 will be averaged for each bin , so this is way to propagate other data one may need to. 

    typical min step is tempMinStep= Qvec[1] - Qvec[0], the differecne between the two first points.
    Igor code calling this:
    	tempMinStep=DSM_Qvec[1]-DSM_Qvec[0]
		IN2G_RebinLogData(DSM_Qvec,DSM_Int,FlyScanRebinToPoints,tempMinStep,Wsdev=DSM_Error,Wxwidth=DSM_dQ)
        where DSM input waves are with high number of points, and are returned in place as reduced number of points
'''

import numpy as np
from scipy.optimize import minimize
from scipy.interpolate import interp1d



def rebin_QRSdata(Wx, Wy, Ws, NumberOfPoints):
    '''
    Rebin data based on the condition Q < 0.0002 and Q > 0.0002.
    This function will split the data into two parts and rebin the second part using logarithmic scaling.
    The first part (Q < 0.0002) will remain unchanged, while the second part (Q > 0.0002) will be rebinned.
    The function will return the merged data from both parts.
    The rebinning is done using the rebin_log_data function, which takes care of the logarithmic scaling.     
    '''
    # create Wdx array based on Wx
    Wdx = np.zeros(len(Wx))
    Wdx[1:] = Wx[1:] - Wx[:-1]
    Wdx[0] = Wdx[1]  # Set the first element to the same value as the second element        
    # Split arrays based on the condition Q < 0.0002
    mask_less = Wx < 0.0002
    Wx_less = Wx[mask_less]
    Wy_less = Wy[mask_less]
    Ws_less = Ws[mask_less]
    Wdx_less = Wdx[mask_less]

    # Split arrays based on the condition Q > 0.0002
    mask_greater = Wx > 0.0002
    Wx_greater = Wx[mask_greater]
    Wy_greater = Wy[mask_greater]
    Ws_greater = Ws[mask_greater]
    Wdx_greater = Wdx[mask_greater]

    MinStep = Wx_greater[1] - Wx_greater[0]

    Wx_greater2, Wy_greater2, W1, W2, W3, W4, W5, Ws_greater2, Wxsdev, Wxwidth = rebin_log_data(Wx_greater, Wy_greater, NumberOfPoints, MinStep, Wsdev=Ws_greater, Wxsdev=None, Wxwidth=Wdx_greater, W1=None, W2=None, W3=None, W4=None, W5=None)


    Q_merged = np.concatenate((Wx_less, Wx_greater2))
    Intensity_merged = np.concatenate((Wy_less, Wy_greater2))      
    Error_merged = np.concatenate((Ws_less, Ws_greater2))    
    dQ_merged = np.concatenate((Wdx_less, Wxwidth)) 
    # remove Nans form errors as that seems to break stuff latrer
    # Create a mask for NaNs
    nan_mask = np.isnan(Error_merged)
    # Indices of non-NaN values
    x_non_nan = Q_merged[~nan_mask]
    y_non_nan = Error_merged[~nan_mask]
    # Create an interpolation function
    interp_func = interp1d(x_non_nan, y_non_nan, kind='linear', fill_value='extrapolate')
    # Replace NaNs with interpolated values
    Error_merged[nan_mask] = interp_func(Q_merged[nan_mask])   


    return Q_merged, Intensity_merged, Error_merged, dQ_merged



def rebin_log_data(Wx, Wy, NumberOfPoints, MinStep, Wsdev=None, Wxsdev=None, Wxwidth=None, W1=None, W2=None, W3=None, W4=None, W5=None):
    # Determine which additional waves need to be calculated
    CalcSdev = Wsdev is not None
    CalcXSdev = Wxsdev is not None
    CalcWidth = Wxwidth is not None
    CalcW1 = W1 is not None
    CalcW2 = W2 is not None
    CalcW3 = W3 is not None
    CalcW4 = W4 is not None
    CalcW5 = W5 is not None

    OldNumPnts = len(Wx)
    if 2 * NumberOfPoints > OldNumPnts:
        print("User requested rebinning of data, but old number of points is less than 2*requested number of points, no rebinning done")
        return None

    if Wx[0] <= 0:
        Wx[0] = Wx[1] / 2
    CorrectStart = Wx[0]

    if MinStep > 0:
        StartX = find_correct_log_scale_start(Wx[0], Wx[-1], NumberOfPoints, MinStep)
    else:
        StartX = CorrectStart

    EndX = StartX + abs(Wx[-1] - Wx[0])
    isGrowing = Wx[0] < Wx[-1]

    tempNewLogDist = np.zeros(NumberOfPoints)
    tempNewLogDistBinWidth = np.zeros(NumberOfPoints)

    logstartX = np.log10(StartX)
    logendX = np.log10(EndX)
    tempNewLogDist = np.logspace(logstartX, logendX, NumberOfPoints, base=10)
    tempNewLogDist += CorrectStart - StartX

    tempNewLogDistBinWidth[1:-1] = tempNewLogDist[2:] - tempNewLogDist[:-2]
    tempNewLogDistBinWidth[0] = tempNewLogDistBinWidth[1]
    tempNewLogDistBinWidth[-1] = tempNewLogDistBinWidth[-2]

    Rebinned_WvX = np.zeros(NumberOfPoints)
    Rebinned_WvY = np.zeros(NumberOfPoints)
    Rebinned_Wv1 = np.zeros(NumberOfPoints) if CalcW1 else None
    Rebinned_Wv2 = np.zeros(NumberOfPoints) if CalcW2 else None
    Rebinned_Wv3 = np.zeros(NumberOfPoints) if CalcW3 else None
    Rebinned_Wv4 = np.zeros(NumberOfPoints) if CalcW4 else None
    Rebinned_Wv5 = np.zeros(NumberOfPoints) if CalcW5 else None
    Rebinned_Wsdev = np.zeros(NumberOfPoints) if CalcSdev else None
    Rebinned_Wxsdev = np.zeros(NumberOfPoints) if CalcXSdev else None

    j = 0
    for i in range(NumberOfPoints):
        cntPoints = 0
        BinHighEdge = tempNewLogDist[i] + tempNewLogDistBinWidth[i] / 2
        while j < OldNumPnts and ((Wx[j] < BinHighEdge) if isGrowing else (Wx[j] > BinHighEdge)):
            Rebinned_WvX[i] += Wx[j]
            Rebinned_WvY[i] += Wy[j]
            if CalcW1:
                Rebinned_Wv1[i] += W1[j]
            if CalcW2:
                Rebinned_Wv2[i] += W2[j]
            if CalcW3:
                Rebinned_Wv3[i] += W3[j]
            if CalcW4:
                Rebinned_Wv4[i] += W4[j]
            if CalcW5:
                Rebinned_Wv5[i] += W5[j]
            if CalcSdev:
                Rebinned_Wsdev[i] += Wsdev[j] ** 2
            if CalcXSdev:
                Rebinned_Wxsdev[i] += Wxsdev[j] ** 2
            cntPoints += 1
            j += 1

        if cntPoints > 0:
            Rebinned_WvX[i] /= cntPoints
            Rebinned_WvY[i] /= cntPoints
            if CalcW1:
                Rebinned_Wv1[i] /= cntPoints
            if CalcW2:
                Rebinned_Wv2[i] /= cntPoints
            if CalcW3:
                Rebinned_Wv3[i] /= cntPoints
            if CalcW4:
                Rebinned_Wv4[i] /= cntPoints
            if CalcW5:
                Rebinned_Wv5[i] /= cntPoints
            if CalcSdev:
                Rebinned_Wsdev[i] = np.sqrt(Rebinned_Wsdev[i] / cntPoints)
            if CalcXSdev:
                Rebinned_Wxsdev[i] = np.sqrt(Rebinned_Wxsdev[i] / cntPoints)

    Wx = Rebinned_WvX
    Wy = Rebinned_WvY

    if CalcW1:
        W1 = Rebinned_Wv1
    if CalcW2:
        W2 = Rebinned_Wv2
    if CalcW3:
        W3 = Rebinned_Wv3
    if CalcW4:
        W4 = Rebinned_Wv4
    if CalcW5:
        W5 = Rebinned_Wv5
    if CalcSdev:
        Wsdev = Rebinned_Wsdev
    if CalcXSdev:
        Wxsdev = Rebinned_Wxsdev
    if CalcWidth:
        Wxwidth = tempNewLogDistBinWidth

    return Wx, Wy, W1, W2, W3, W4, W5, Wsdev, Wxsdev, Wxwidth

def my_find_start_value_func(x1, w):
    """
    Objective function to find the correct start value for logarithmic scaling.

    Parameters:
    - x1: The start value where we need to start with log stepping
    - w: A list or array containing [totalRange, NumSteps, MinStep]

    Returns:
    - The absolute difference between the calculated last minimum step and the desired minimum step
    """
    total_range, num_steps, min_step = w
    last_min_step = 10**(np.log10(x1) + (np.log10(x1 + total_range) - np.log10(x1)) / num_steps) - 10**(np.log10(x1))
    return abs(last_min_step - min_step)

def find_correct_log_scale_start(StartValue, EndValue, NumPoints, MinStep):
    """
    Finds the correct start value for logarithmic scaling using optimization.

    Parameters:
    - StartValue: The initial start value
    - EndValue: The end value
    - NumPoints: The number of points
    - MinStep: The minimum step size

    Returns:
    - The optimal start value for logarithmic scaling
    """
    # Define the parameters for the optimization
    w = [EndValue - StartValue, NumPoints, MinStep]

    # Initial guess for the start value
    x1_initial = StartValue

    # Use scipy.optimize.minimize to find the optimal start value
    result = minimize(my_find_start_value_func, x1_initial, args=(w,), method='L-BFGS-B', bounds=[(1e-10, None)])

    # The optimal start value is in result.x[0]
    return result.x[0]


# # Example usage
# StartValue = 1.0
# EndValue = 10.0
# NumPoints = 100
# MinStep = 0.1

# optimal_start = find_correct_log_scale_start(StartValue, EndValue, NumPoints, MinStep)
# print("Optimal Start Value:", optimal_start)