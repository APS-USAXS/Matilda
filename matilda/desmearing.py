'''
    Converted desmearing code from Igor Pro
    This routine will calculate desmeard data using the Lake method and used by USAXS instrument during automatic data reduction 

    final routine is : desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=SlitLength, MaxNumIter = None, ExtrapMethod='flat',ExtrapQstart=None)
    Returns:
            Desmeared: Q, Intensity, Errors, dQ arrays
              
    
    slit length is in q units [1/A] and must be provided
    MaxNumIter = None - if not provided, set to 20 and also using automatic method
    ExtrapMethod='flat', other options are "Porod", "Power law" and "Power law with flat"
    Note that if extrapolation with more complciated functions fails, flat is automatically applied. Untested some of the extrapolations. 
    ExtrapQstart=None - if not present, set to Qmax/1.5 
    Assumes Qmax is significantly larger than slit length, typically Sl=0.03 and Qmax=0.3
    Igor code is in IN3_Calculations with functions with similar names.
    
'''
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import simpson, trapezoid



def extendData(Q_vct, Int_wave, Err_wave, slitLength, Qstart, SelectedFunction):
    if not isinstance(slitLength, (int, float)):
        raise ValueError("Slit length error")
    if slitLength < 0.0001 or slitLength > 1:
        raise ValueError("Weird value for Slit length, please check")
    
    # Check for NaNs or INFs in Int_wave
    if np.isnan(Int_wave).any() or np.isinf(Int_wave).any():
        raise ValueError("Int_wave contains NaNs or INFs")
    
    ProblemsWithQ = ""
    ProblemWithFit = ""
    ProblemsWithInt = ""
    
    DataLengths = len(Q_vct) - 1
    Qstep = ((Q_vct[DataLengths] / Q_vct[DataLengths - 1]) - 1) * Q_vct[DataLengths]
    ExtendByQ = np.sqrt(Q_vct[DataLengths]**2 + (1.5 * slitLength)**2) - Q_vct[DataLengths]
    if ExtendByQ < 2.1 * Qstep:
        ExtendByQ = 2.1 * Qstep
    NumNewPoints = int(np.floor(ExtendByQ / Qstep))
    if NumNewPoints < 1:
        NumNewPoints = 1
    OriginalNumPnts = len(Int_wave)
    if NumNewPoints > OriginalNumPnts:
        NumNewPoints = OriginalNumPnts
    newLength = len(Q_vct) + NumNewPoints
    
    # Find the index to start fitting
    FitFrom = np.searchsorted(Q_vct, Qstart)
    if FitFrom <= 0:
        FitFrom = DataLengths - 10
        ProblemsWithQ = "I did reset Fitting Q range for you..."
    
    # Extend the arrays
    Int_wave = np.resize(Int_wave, newLength)
    Q_vct = np.resize(Q_vct, newLength)
    Err_wave = np.resize(Err_wave, newLength)
    
    # Define fitting functions
    def flat_fnct(x, K0):
        return K0
    
    def power_law_fnct(x, K0, K1, K2):
        return K0 + K1 * x**K2
    
    def porod_fnct(x, K0, K1):
        return K0 + K1 / x**4
    
    # Fit and extend based on the selected function
    if SelectedFunction == "flat":
        try:
            popt, _ = curve_fit(flat_fnct, Q_vct[FitFrom:DataLengths], Int_wave[FitFrom:DataLengths], sigma=Err_wave[FitFrom:DataLengths])
            for i in range(1, NumNewPoints + 1):
                Q_vct[DataLengths + i] = Q_vct[DataLengths] + (ExtendByQ) * (i / NumNewPoints)
                Int_wave[DataLengths + i] = popt[0]
        except RuntimeError:
            ProblemWithFit = "Linear fit function did not converge properly, change function or Q range"
    
    elif SelectedFunction == "Power law":
        initial_guess = [np.min(Int_wave),(Int_wave[FitFrom]-np.min(Int_wave))/(Q_vct[FitFrom]**(-3.5))]
        try:
            popt, _ = curve_fit(power_law_fnct, Q_vct[FitFrom:DataLengths], Int_wave[FitFrom:DataLengths], sigma=Err_wave[FitFrom:DataLengths],p0=initial_guess)
            for i in range(1, NumNewPoints + 1):
                Q_vct[DataLengths + i] = Q_vct[DataLengths] + (ExtendByQ) * (i / NumNewPoints)
                Int_wave[DataLengths + i] = popt[0] + popt[1] * Q_vct[DataLengths + i]**popt[2]
        except RuntimeError:
            ProblemWithFit = "Power law fit function did not converge properly, change function or Q range"
    
    elif SelectedFunction == "Porod":
        initial_guess = [np.min(Int_wave),(Int_wave[FitFrom]-np.min(Int_wave))/(Q_vct[FitFrom]**(-3.5))]
        try:
            popt, _ = curve_fit(porod_fnct, Q_vct[FitFrom:DataLengths], Int_wave[FitFrom:DataLengths], sigma=Err_wave[FitFrom:DataLengths],p0=initial_guess)
            for i in range(1, NumNewPoints + 1):
                Q_vct[DataLengths + i] = Q_vct[DataLengths] + (ExtendByQ) * (i / NumNewPoints)
                Int_wave[DataLengths + i] = popt[0] + popt[1] / Q_vct[DataLengths + i]**4
        except RuntimeError:
            ProblemWithFit = "Porod fit function did not converge properly, change function or Q range"
    
    elif SelectedFunction == "PowerLaw w flat":
        initial_guess = [np.min(Int_wave),(Int_wave[FitFrom]-np.min(Int_wave))/(Q_vct[FitFrom]**(-3.5)), -3.5]
        try:
            popt, _ = curve_fit(power_law_fnct, Q_vct[FitFrom:DataLengths], Int_wave[FitFrom:DataLengths], sigma=Err_wave[FitFrom:DataLengths],p0=initial_guess)
            for i in range(1, NumNewPoints + 1):
                Q_vct[DataLengths + i] = Q_vct[DataLengths] + (ExtendByQ) * (i / NumNewPoints)
                Int_wave[DataLengths + i] = popt[0] + popt[1] * Q_vct[DataLengths + i]**popt[2]
        except RuntimeError:
            ProblemWithFit = "Power Law with flat fit function did not converge properly, change function or Q range"
    
    ExtensionFailed = False
    ErrorMessages = ""
    if ProblemsWithQ:
        ErrorMessages += ProblemsWithQ + "\n"
    if ProblemsWithInt:
        ErrorMessages += ProblemsWithInt + "\n"
    if ProblemWithFit:
        ErrorMessages += ProblemWithFit
    
    if ErrorMessages:
        ExtensionFailed = True
        #print(f"Extending data by average intensity (aka:flat)")
        AveInt = np.mean(Int_wave[FitFrom:DataLengths])
        for i in range(1, NumNewPoints + 1):
            Q_vct[DataLengths + i] = Q_vct[DataLengths] + (ExtendByQ) * (i / NumNewPoints)
            Int_wave[DataLengths + i] = AveInt
        ExtensionFailed = False
    
    return Q_vct, Int_wave, Err_wave, ExtensionFailed

def smearDataOverSlit(Q_vec_sm2, Smear_Q, Smear_Q2, tempQ_vec_sm, tempInt_to_smear):
    """
    Smear data fast function.

    Parameters:
    Q_vec_sm2 : float
        The squared value of the Q vector to be smeared.
    Smear_Q : array-like
        The Q values for smearing. This is the slit Q values. 
    Smear_Q2 : array-like
        The squared Q values for smearing.
    tempQ_vec_sm : array-like
        The temporary Q vector for smearing.
    tempInt_to_smear : array-like
        The temporary intensities to smear.

    Returns:
    float
        The integrated intensity over the slit.
    """
    # Calculate the interpolated Q values
    InterSmear_Q = np.sqrt(Q_vec_sm2 + Smear_Q2)

    # Interpolate the intensity values
    #np.interp(xnew, x, y)
    Smear_Int = np.interp(InterSmear_Q, tempQ_vec_sm, tempInt_to_smear)
    # interp_func = interp1d(tempQ_vec_sm, tempInt_to_smear, kind='linear', fill_value="extrapolate")
    # Smear_Int = interp_func(InterSmear_Q)

    # Integrate the intensity over the slit
    #integrated_intensity = simpson(Smear_Int, Smear_Q)
    integrated_intensity = trapezoid(Smear_Int, Smear_Q)

    return integrated_intensity

def smearIntensityArray(Int_to_smear, Q_vec_sm, slitLength):
    """
    Smear data function.

    Parameters:
    Int_to_smear : array-like
        The intensities to smear.
    Q_vec_sm : array-like
        The Q vector for smearing.
    slitLength : float
        The length of the slit.

    Returns:
    array-like
        The smeared intensities.
    """
    oldNumPnts = len(Q_vec_sm)
    
    # Determine the number of new points
    if oldNumPnts < 300:
        newNumPoints = 2 * oldNumPnts
    else:
        newNumPoints = oldNumPnts + 300
    
    # Extend the Int_to_smear and Q_vec_sm arrays
    tempInt_to_smear = np.resize(Int_to_smear, newNumPoints)
    tempQ_vec_sm = np.resize(Q_vec_sm, newNumPoints)
    
    tempQ_vec_sm[oldNumPnts:] = tempQ_vec_sm[oldNumPnts - 1] + 20 * (tempQ_vec_sm[:newNumPoints - oldNumPnts])
    tempInt_to_smear[oldNumPnts:] = tempInt_to_smear[oldNumPnts - 1] * (1 - (tempQ_vec_sm[oldNumPnts:] - tempQ_vec_sm[oldNumPnts]) / (20 * tempQ_vec_sm[oldNumPnts - 1]))
    
    # Create Smear_Q and Smear_Int arrays
    Smear_Q = np.zeros(oldNumPnts)
    Smeared_int = np.zeros(oldNumPnts)
    
    DataLengths = len(Q_vec_sm)
    
    # Create distribution of points in Smear_Q
    Smear_Q = slitLength * (Q_vec_sm - Q_vec_sm[0]) / (Q_vec_sm[DataLengths - 1] - Q_vec_sm[0])
    
    # Calculate squared values
    Q_vec_sm2 = Q_vec_sm ** 2
    Smear_Q2 = Smear_Q ** 2
    
    # Calculate smeared intensities using the fast function
    Smeared_int = np.array([smearDataOverSlit(Q_vec_sm2[i], Smear_Q, Smear_Q2, tempQ_vec_sm, tempInt_to_smear) for i in range(len(Q_vec_sm2))])
    # Using a loop for better readability
    # for i in range(DataLengths):
    #     # Calculate the smeared intensity for each point
    #     Smeared_int[i] = smearDataOverSlit(Q_vec_sm2[i], Smear_Q, Smear_Q2, tempQ_vec_sm, tempInt_to_smear)

    # Normalize the smeared intensities
    Smeared_int *= 1 / slitLength
    
    return Smeared_int

def calculateLineAverage(WaveY, WaveX, ivalue):
    """
    Calculate line average without doing a line fit.

    Parameters:
    WaveY : array-like
        The Y values of the wave.
    WaveX : array-like
        The X values of the wave.
    ivalue : int
        The index for which to calculate the line average.

    Returns:
    float
        The calculated line average value.
    """
    if ivalue > 1:
        sumx = WaveX[ivalue - 2] + WaveX[ivalue - 1] + WaveX[ivalue] + WaveX[ivalue + 1]+ WaveX[ivalue + 2]
        sumx2 = WaveX[ivalue - 2]**2 +WaveX[ivalue - 1]**2 + WaveX[ivalue]**2 + WaveX[ivalue + 1]**2+ WaveX[ivalue + 2]**2
        sumy = WaveY[ivalue - 2] + WaveY[ivalue - 1] + WaveY[ivalue] + WaveY[ivalue + 1]+ WaveY[ivalue + 2]
        sumxy = (WaveX[ivalue - 2] * WaveY[ivalue - 2] +
                 WaveX[ivalue - 1] * WaveY[ivalue - 1] +
                 WaveX[ivalue] * WaveY[ivalue] +
                 WaveX[ivalue + 1] * WaveY[ivalue + 1] +
                 WaveX[ivalue + 2] * WaveY[ivalue + 2])
        mval = (5 * sumxy - sumx * sumy) / (5 * sumx2 - sumx**2)
        cval = (sumy - mval * sumx) / 5
        return mval * WaveX[ivalue] + cval
    return 0

def calculateErrors(SmErrors, SmIntensity, FitIntensity, Qvector):
    """
    Calculate errors using Pete's formulas.

    Parameters:
    SmErrors : array-like
        The smearing errors.
    SmIntensity : array-like
        The smearing intensities.
    FitIntensity : array-like
        The fitted intensities.
    Qvector : array-like
        The Q vector.

    Returns:
    array-like
        The calculated desmearing errors.
    """
    DsmErrors = FitIntensity * (SmErrors / SmIntensity)
    imax = len(FitIntensity)
    DsmErrors = np.resize(DsmErrors, imax)

    for i in range(2, imax - 2):
        if np.isfinite(FitIntensity[i - 2]) and np.isfinite(FitIntensity[i]) and np.isfinite(FitIntensity[i + 1]) and np.isfinite(FitIntensity[i + 2]):
            DsmErrors[i] += np.abs(calculateLineAverage(FitIntensity, Qvector, i) - FitIntensity[i])

    DsmErrors[0] = 3*DsmErrors[2]
    DsmErrors[1] = 2*DsmErrors[2]
    DsmErrors[imax - 1] = DsmErrors[imax - 3]
    DsmErrors[imax - 2] = DsmErrors[imax - 3]

    # Smooth the errors - do not smooth errors... 
    # DsmErrors = smooth_errors(DsmErrors)

    return DsmErrors

def smooth_errors(errors, window_len=3):
    """
    Smooth the errors using a simple moving average.

    Parameters:
    errors : array-like
        The errors to smooth.
    window_len : int
        The length of the smoothing window.

    Returns:
    array-like
        The smoothed errors.
    """
    smoothed = np.convolve(errors, np.ones(window_len)/window_len, mode='same')
    return smoothed


def oneDesmearIteration(SlitLength, QWave, DesmearIntWave, DesmearEWave, origSmearedInt, origSmearedErr, NormalizedError,ExtrapMethod,ExtrapQstart):
    """
    Perform one iteration of the desmearing process.

    Parameters:
    DesmearIntWave : array-like
        The intensity wave to desmear.
    DesmearQWave : array-like
        The Q wave corresponding to the intensity wave.
    DesmearEWave : array-like
        The error wave for the desmeared data.
    origSmearedInt : array-like
        The original smeared intensities.
    origSmearedErr : array-like
        The original smeared errors.
    NormalizedError : array-like
        The normalized error array.

    Returns:
    DesmearQWave, DesmearIntWave, DesmearEWave, NormalizedError
    """
    # Retrieve parameters from a hypothetical settings object
    # in Igor GUI system these are global poarameters. 
    BackgroundFunction = ExtrapMethod  # There are four choises of extension functions. 
    # 'flat', 'Power law', 'Porod', 'PowerLaw w flat'.
    if BackgroundFunction is None:
       BackgroundFunction = 'flat'
    #NumberOfIterations = 0  # Placeholder
    numOfPoints = len(DesmearIntWave)
    if ExtrapQstart is not None:
        BckgStartQ = ExtrapQstart
    else:
        BckgStartQ = QWave[-1] / 1.5

    if BckgStartQ > QWave[-1] / 1.5:
        BckgStartQ = QWave[-1] / 1.5

    SmFitIntensity = np.copy(DesmearIntWave)
    SmQWave = np.copy(QWave)
    SmErrors = np.copy(origSmearedErr)
    #here we extend the data to be at least slit length longer. 
    # it shoudl never fail as if it does with complex functions, inside will simply do flat extension.
    SmQWave, SmFitIntensity, SmErrors, ExtensionFailed = extendData(QWave,DesmearIntWave,SmErrors, SlitLength, BckgStartQ, BackgroundFunction)
    if ExtensionFailed:
        return 1
    #this is now smeared version of Intensity
    if SlitLength > 0:
        SmFitIntensity = smearIntensityArray(SmFitIntensity, SmQWave, SlitLength)

    # Resize arrays back to original length
    SmFitIntensity = SmFitIntensity[:numOfPoints]
    DesmearIntWave = DesmearIntWave[:numOfPoints]
    QWave = SmQWave[:numOfPoints]
    SmErrors = SmErrors[:numOfPoints]
    NormalizedError = NormalizedError[:numOfPoints]

    # Calculate normalized error
    # Handle division by zero for origSmearedErr
    with np.errstate(divide='ignore', invalid='ignore'):
        NormalizedError = np.where(origSmearedErr != 0, (origSmearedInt - SmFitIntensity) / origSmearedErr, 0)


    # Fast and slow convergence
    FastFitIntensity = DesmearIntWave * (origSmearedInt / SmFitIntensity)
    #SlowFitIntensity = DesmearIntWave + (origSmearedInt - SmFitIntensity)

    # Update DesmearIntWave based on normalized error
    for i in range(len(DesmearIntWave)):
        if abs(NormalizedError[i]) > 0.5:
            DesmearIntWave[i] = FastFitIntensity[i]
        else:
            DesmearIntWave[i] = DesmearIntWave[i]

    #NumberOfIterations += 1

    # Remove extremes from normalized error
    #min_loc = np.argmin(NormalizedError)
    #max_loc = np.argmax(NormalizedError)
    #NormalizedError[min_loc] = np.nan
    #NormalizedError[max_loc] = np.nan

    # Calculate errors
    DesmearEWave.fill(0)
    DesmearEWave = calculateErrors(origSmearedErr, origSmearedInt, DesmearIntWave, QWave)

    return QWave, DesmearIntWave, DesmearEWave, NormalizedError

# def removeNaNsFromArray(wave):
#     """
#     Remove NaNs from a wave by replacing them with zero.

#     Parameters:
#     wave : array-like
#         The wave to process.
#     """
#     wave[np.isnan(wave)] = 0

def desmearData(SMR_Qvec, SMR_Int, SMR_Error, SMR_dQ, slitLength=None, MaxNumIter = None, ExtrapMethod='flat',ExtrapQstart=None):
    """
    Perform the desmearing process on the provided data.

    Parameters:
    SlitLength : float
        The length of the slit. in Q units, [1/A] most likely
    MaxNumIter : int
        The number of iterations for the desmearing process.Typically 5-10, if not provided usign automatic method and max number set to 50. 
    SMR_Int : array-like
        The smeared intensities.
    SMR_Error : array-like
        The smeared errors.
    SMR_Qvec : array-like
        The Q vector for the smeared data.
    SMR_dQ : array-like
        The dQ values for the smeared data.
    ExtrapMethod: 
        "flat", "Power law", "Power law with flat", "Porod"
    ExtrapQstart
        Where to start extrapolation, if not provided set to Qmax/1.5

    Returns:
    tuple
        DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ
    """
    if SMR_Int is None or len(SMR_Int) == 0 or slitLength is None:
        return None, None, None, None, None

    tmpWork_Int = np.copy(SMR_Int)
    tmpWork_Error = np.copy(SMR_Error)
    tmpWork_Qvec = np.copy(SMR_Qvec)
    tmpWork_dQ = np.copy(SMR_dQ)
 
    DesmNormalizedError = np.zeros(len(SMR_Int))
    absNormalizedError = np.zeros(len(SMR_Int))
    #numOfPoints = len(SMR_Int)
    endme = 0
    oldendme = 0
    DesmearAutoTargChisq = 0.5
    NumIterations = 0
    if MaxNumIter is None:
        MaxNumIter = 20

    while True:
        tmpWork_Qvec, tmpWork_Int, tmpWork_Error, DesmNormalizedError = oneDesmearIteration(slitLength,tmpWork_Qvec, tmpWork_Int, tmpWork_Error, SMR_Int, SMR_Error, DesmNormalizedError, ExtrapMethod,ExtrapQstart)
        absNormalizedError = np.abs(DesmNormalizedError)
        endme = np.average(absNormalizedError)
        #this is difference in convergence between iterations
        difff = 1 - (oldendme / endme)
        oldendme = endme
        NumIterations += 1
        #Conditions under which we will end 
        if (endme < DesmearAutoTargChisq or abs(difff) < 0.02 or NumIterations > MaxNumIter):
            # print("Convergence reached")
            # print("Number of iterations (>", MaxNumIter,"): ", NumIterations)
            # print("Final average error (<0.5): ", endme)
            # print("Final convergence (<0.02): ", abs(difff))
            break

    DSM_Int = np.copy(tmpWork_Int)
    DSM_Qvec = np.copy(tmpWork_Qvec)
    DSM_Error = np.abs(tmpWork_Error)  # Remove negative values
    DSM_dQ = np.copy(tmpWork_dQ)
    #print(NumIterations)
    return DSM_Qvec, DSM_Int, DSM_Error, DSM_dQ
