import numpy  as np 

import numpy.polynomial.polynomial as poly

from scipy.optimize  import curve_fit
from scipy.signal    import savgol_filter
from scipy.stats     import t
from scipy.linalg    import svd

# R gas constant in kcal/(kelvin mol)
R_gas  = 1.987 / 1000 

def temperature_to_kelvin(T):

    # To kelvin, if required
    if np.max(T) < 270:
        T = T + 273.15

    return(T)

# not used in the app
def compute_bootstrap_regression_slope(x_data,y_data,num_bootstrap_samples=60):

    bootstrap_slopes = []

    for _ in range(num_bootstrap_samples):

        # Generate a bootstrap sample by resampling with replacement
        bootstrap_indices = np.random.choice(len(x_data), len(x_data), replace=True)
        bootstrap_x = x_data[bootstrap_indices]
        bootstrap_y = y_data[bootstrap_indices]

        # Fit the curve using poly.Polynomial.fit on the bootstrap sample
        poly_fit = poly.Polynomial.fit(bootstrap_x,bootstrap_y,1)
        slope    = poly_fit.coef[1]

        bootstrap_slopes.append(slope)

    # Specify the desired confidence level (e.g., 95%)
    # We must select this extreme value because for low n (i.e., 6 points) we will be underestimating the confidence interval
    # In other words, the 99% confidence interval will only contain the true slope in ~ 90 % of the cases
    confidence_level = 0.99

    # Calculate the lower and upper percentiles for each parameter
    lower_percentile = (1 - confidence_level) / 2
    upper_percentile = 1 - lower_percentile

    conf_interval = np.percentile(bootstrap_slopes, [lower_percentile * 100, upper_percentile * 100], axis=0)

    return conf_interval

def estimate_useful_wavelength_based_on_snr(wavelength,signal,wavelength_window_length=6,window_size=7,threshold=25):

    '''

    Subset the CD data based on the wavelength that show a good enough signal-to-noise (snr) ratio

    The formula for SNR is (a**2) / (b**2), where 
        a is the value of the smoothed curve that depends on the 'wavelength_window_length', and
        b is the value of the rolling standard deviation computed using a certain 'window_size'

    '''

    # Set the window size (odd number) and polynomial order
    deltaWL         = ( max(wavelength) - min(wavelength) ) / (len(wavelength) - 1) 

    odd_n_data_points_window_len         = np.ceil(wavelength_window_length / deltaWL) // 2 * 2 + 1

    # Apply Savitzky-Golay filter for smoothing
    smoothed_y = savgol_filter(signal,axis=1,
        window_length=odd_n_data_points_window_len,polyorder=3,mode="nearest")
    
    # Compute the residuals
    residuals = smoothed_y - signal
    
    # Rolling standard deviation
    rolling_std_by_column = np.apply_along_axis(lambda x: np.convolve(x, np.ones(window_size) / (window_size-1), mode='valid') ** 0.5, 
        axis=0, arr=residuals**2)
    
    length_difference = int((residuals.shape[0] - rolling_std_by_column.shape[0]) / 2)

    snr = (smoothed_y[length_difference:-length_difference,:]**2) / (rolling_std_by_column**2)

    useful_snr = snr > threshold

    useful_wl = np.all(useful_snr, axis=1)

    wavelength_useful = wavelength[length_difference:-length_difference][useful_wl]

    return wavelength_useful 

def estimate_useful_wavelength_based_on_amplitude(wavelength,signal,measurement_factor,relative_threshold = 0.5):

    """
    Find which wavelength data points have a change of at least X % (regarding the maximum amplitude)

    Requires:

        - wavelength           :   1D numpy array of length n
        - signal               :   2D numpy array with n rows and m columns
        - measurement_factor   :   1D numpy array of length m (an example would be temperature measurements)

    """

    # Step 1: Subtract the signal based on the measurement_factor
    lowest_factor_signal      = signal[:, np.argmin(measurement_factor)]
    highest_factor_signal     = signal[:, np.argmax(measurement_factor)]

    subtracted_signal              = highest_factor_signal - lowest_factor_signal

    # Step 2: Calculate the absolute values of the elements in the resulting vector
    absolute_values   = np.abs(subtracted_signal)

    # Step 3: Find the largest absolute value
    largest_absolute_value = np.max(absolute_values)

    threshold              = relative_threshold * largest_absolute_value

    indices_meeting_condition = (absolute_values > threshold)

    wavelength_useful = wavelength[indices_meeting_condition]

    return wavelength_useful

def chemical_unfolding_with_linear_dependence_one_curve(D,T,D50,M,folded_ellipticity,m1,unfolded_ellipticity,m2):

    T = temperature_to_kelvin(T)

    dG = M * (D50 - D)

    Keq                 = np.exp(-dG / (R_gas * T))  
    unfolded_fraction   =   Keq / (Keq + 1)
    folded_fraction     =  1 - unfolded_fraction

    # Here m1 and m2 are the slopes of the linear pretransition and post-transition changes, respectively. 
    Y = folded_fraction * ( folded_ellipticity + m1*D) + unfolded_fraction * ( unfolded_ellipticity + m2*D ) 

    return Y

def chemical_unfolding_one_curve(D,T,D50,M,folded_ellipticity,unfolded_ellipticity):
    # Remove the linear dependence on the denaturant agent concentration

    return chemical_unfolding_with_linear_dependence_one_curve(D,T,D50,M,folded_ellipticity,0,unfolded_ellipticity,0)

def chemical_unfolding_one_curve_slope_native(D,T,D50,M,folded_ellipticity,m1,unfolded_ellipticity):
    # Remove the linear dependence of the unfolded state on the denaturant agent concentration

    return chemical_unfolding_with_linear_dependence_one_curve(D,T,D50,M,folded_ellipticity,m1,unfolded_ellipticity,0)

def chemical_unfolding_one_curve_slope_unfolded(D,T,D50,M,folded_ellipticity,unfolded_ellipticity,m2):
    # Remove the linear dependence of the folded state on the denaturant agent concentration

    return chemical_unfolding_with_linear_dependence_one_curve(D,T,D50,M,folded_ellipticity,0,unfolded_ellipticity,m2)

def fit_chemical_unfolding(listOfChemConcentration,listOfSignals,temperature,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded):

    '''
    Fit the chemical unfolding profile of many curves at the same time
    Useful function to do global fitting of local and global parameters

    Requires:

        - The 'listOfChemConcentration' (one per dataset)
        - The 'listOfSignals'           (one per dataset)

    Returns:

        The fitting of the chemical unfolding curves based on the parameters D50, m, folded_ellipticity and unfolded_ellipticity
    '''

    # We need 1D arrays for numpy

    try:
        allSignal       = np.concatenate(listOfSignals)
    except:
        allSignal       = listOfSignals
    
    def chemical_unfolding(dummyVariable,*args):

        '''
        Calculate the chemical unfolding profile of many curves at the same time

        Requires:

            - The 'listOfChemConcentration' containing each of them a single dataset

        The other arguments have to be in the following order:

            - Global m 
            - Global D50 
            - Single intercepts, folded
            - Single slopes, folded
            - Single intercepts, unfolded
            - Single slopes, unfolded

        Returns:

            The melting curves based on the parameters m, D50,
                slopes and intercept of the folded and unfolded states
        '''
 
        totalDataSets = len(listOfChemConcentration)
        M             = args[0] # M
        D50           = args[1] # Concentration at which the sample is 50 % unfolded

        interceptsFolded    = args[2:(2+totalDataSets)]
        interceptsUnfolded  = args[(2+totalDataSets):(2+totalDataSets*2)]

        if fitSlopeNative:

            slopesFolded   = args[(2+totalDataSets*2):(2+totalDataSets*3)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,D in enumerate(listOfChemConcentration):

            bN, bU = interceptsFolded[i],   interceptsUnfolded[i]

            if fitSlopeNative and fitSlopeUnfolded:

                Y = chemical_unfolding_with_linear_dependence_one_curve(D,temperature,D50,M,bN,slopesFolded[i],bU,slopesUnfolded[i])
            
            if fitSlopeNative and not fitSlopeUnfolded:

                Y = chemical_unfolding_one_curve_slope_native(D,temperature,D50,M,bN,slopesFolded[i],bU)

            if not fitSlopeNative and  fitSlopeUnfolded:

                Y = chemical_unfolding_one_curve_slope_unfolded(D,temperature,D50,M,bN,bU,slopesUnfolded[i])

            if not fitSlopeNative and not fitSlopeUnfolded:

                Y = chemical_unfolding_one_curve(D,temperature,D50,M,bN,bU)

            signal.append(Y)

        return np.array(signal).flatten()

    global_fit_params, cov = curve_fit(chemical_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds))

    return global_fit_params, cov

    return None

def thermal_unfolding_one_curve(T,Tm,dh,bN,kN,bU,kU,Cp):

    T  = temperature_to_kelvin(T)
    Tm = temperature_to_kelvin(Tm)
    
    dG = dh * (1 - T / Tm) - Cp * (Tm - T + T * np.log(T / Tm))
    
    Ku = np.exp(-dG / (R_gas * T))
    
    Y = (Ku / (1 + Ku)) * (kU * T + bU) + (1 / (1 + Ku)) * (kN * T + bN)
    
    return Y

def fit_thermal_unfolding(listOfTemperatures,listOfSignals,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded):

    '''
    Fit the thermal unfolding profile of many curves at the same time
    Useful function to do global fitting of local and global parameters

    Requires:

        - The 'listOfTemperatures' (one per dataset)
        - The 'listOfSignals'      (one per dataset)

    Returns:

        The fitting of the melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
            slopes and intercept of the folded and unfolded states

    '''

    # We need 1D arrays for numpy
    try:
        allSignal       = np.concatenate(listOfSignals)
    except:
        allSignal       = listOfSignals
    
    def thermal_unfolding(dummyVariable,*args):

        '''
        Calculate the thermal unfolding profile of many curves at the same time

        Requires:

            - The 'listOfTemperatures' containing each of them a single dataset

        The other arguments have to be in the following order:

            - Global melting temperature    
            - Global enthalpy of unfolding 
            - Single intercepts, folded
            - Single slopes, folded
            - Single intercepts, unfolded
            - Single slopes, unfolded

        Returns:

            The melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
                slopes and intercept of the folded and unfolded states

        '''

        # Heat capacity of unfolding
        Cp = 0 

        totalDataSets = len(listOfTemperatures)
        Tm            = args[0] # Temperature of melting
        dh            = args[1] # Enthalpy of unfolding

        interceptsFolded    = args[2:(2+totalDataSets)]
        interceptsUnfolded  = args[(2+totalDataSets):(2+totalDataSets*2)]

        kN, kU = 0, 0
        if fitSlopeNative:

            slopesFolded   = args[(2+totalDataSets*2):(2+totalDataSets*3)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,T in enumerate(listOfTemperatures):

            bN, bU = interceptsFolded[i],   interceptsUnfolded[i]
            
            if fitSlopeNative:
                kN = slopesFolded[i]
            
            if fitSlopeUnfolded:
                kU = slopesUnfolded[i]

            Y = thermal_unfolding_one_curve(T,Tm,dh,bN,kN,bU,kU,Cp)
            signal.append(Y)

        return np.array(signal).flatten()

    global_fit_params, cov = curve_fit(thermal_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds))

    return global_fit_params, cov

def generate_initial_params_combinations(n_params,paramsSign,
    logSpaceLimitLeft=-3,logSpaceLimitRight=3):

    all_combis = []

    # Set the max number of possible initial combinations

    combisLimitLst = 420 * (1.55 ** np.arange(1,6))

    n_params_temp = n_params
    if n_params > 5:
        n_params_temp = 5

    combisLimit = combisLimitLst[n_params_temp-1]

    combisLimit = combisLimit * ( 1.7** (len([x for x in paramsSign if x != 'unkown']) ) )

    NlogNumbers = int(10**(np.log10(combisLimit) /  n_params) )

    logSeq = np.logspace(logSpaceLimitLeft, logSpaceLimitRight, NlogNumbers)

    logSeq  = np.concatenate((logSeq, -logSeq))

    def recursive_generation_of_combinations(bs=[]):

        if len(bs) == n_params:

            all_combis.append(bs)

        else:

            for logNumber in logSeq:

                recursive_generation_of_combinations(bs + [logNumber])

    recursive_generation_of_combinations()

    for i,sign in enumerate(paramsSign):

        if sign == 'positive':

            all_combis = [combi for combi in all_combis if combi[i] >= 0]

        if sign == 'negative':

            all_combis = [combi for combi in all_combis if combi[i] < 0]

    return all_combis

def apply_svd(X):

    U, S, VT = np.linalg.svd(X)

    # Calculate the total variance or correlation
    total_variance      = np.sum(S ** 2) 
    cumulative_variance = np.cumsum(S ** 2) 

    # The matrix V contains the variation of each component against the temperature / measurement dimension

    a_is = []

    for i in range(U.shape[1]):

        def coefficients_bi(column):
        # Your custom logic here
            return U[:,i].dot(column)

        a_i = np.apply_along_axis(coefficients_bi, axis=0, arr=X)

        a_is.append(a_i)

    coefficients = np.array(a_is)        

    # Basis spectra
    basis_spectra       = U
    
    # Cumulated explained variance of the components
    explained_variance  = cumulative_variance / total_variance * 100

    return explained_variance, basis_spectra, coefficients

def align_basis_spectra_and_coefficients(X,basis_spectra,coefficients):

    # Align basis spectra peaks to the original data
    # In other words, we want that if the original spectra has a peak with positive values of the CD signal,
    # so does our basis spectra

    # Fix the n_cutoff to remove the first n and last n rows of X before finding the peak

    n_cutoff = 5

    maxV_abs = np.abs( np.max(X[n_cutoff:-n_cutoff,:]) )
    minV_abs = np.abs( np.min(X[n_cutoff:-n_cutoff,:]) )

    positive_peak = maxV_abs > minV_abs

    k = basis_spectra.shape[1]

    for i in range(k):

        prcomp_i = basis_spectra[:,i]
        
        maxV_abs_prcomp_i = np.abs( np.max(prcomp_i[n_cutoff:-n_cutoff]) )
        minV_abs_prcomp_i = np.abs( np.min(prcomp_i[n_cutoff:-n_cutoff]) )

        positive_peak_prcomp_i = maxV_abs_prcomp_i > minV_abs_prcomp_i 

        if positive_peak_prcomp_i != positive_peak:

            coeff_i  = coefficients[i,:]
            
            basis_spectra[:,i]  = - prcomp_i
            coefficients[i,:]   = - coeff_i

    return basis_spectra, coefficients

def angleFromCathets(adjacentLeg,oppositeLeg):

    '''
    Input:  length of each leg

    Output: angle between the hypotenuse and the adjacent leg
    '''
 
    hypotenuse = np.sqrt(adjacentLeg**2 + oppositeLeg**2)

    return np.arccos(adjacentLeg / hypotenuse)

def get_2d_counterclockwise_rot_matrix(angle_in_radians):

    '''
    Obtain the rotation matrix for a 2d coordinates system using a counterclockwise direction

    '''

    rotM = np.array([[np.cos(angle_in_radians), np.sin(angle_in_radians)],
        [-np.sin(angle_in_radians), np.cos(angle_in_radians)]])

    return rotM

def get_3d_counterclockwise_rot_matrix_around_z_axis(angle_in_radians):

    rotM =  np.array([[np.cos(angle_in_radians), np.sin(angle_in_radians), 0],
        [-np.sin(angle_in_radians), np.cos(angle_in_radians), 0],
        [0, 0, 1]])

    return rotM

def get_3d_clockwise_rot_matrix_around_y_axis(angle_in_radians):

    rotM =  np.array([[np.cos(angle_in_radians), 0, np.sin(angle_in_radians)],
        [0, 1, 0],
        [-np.sin(angle_in_radians), 0, np.cos(angle_in_radians)]])

    return rotM

def rotate_two_basis_spectra(X,basis_spectra,pca_based=False):

    """
    Create a new basis spectra using a linear combination of the first and second basis spectra

    Requires:

        X             : the raw data matrix of size n*m, where 'n' is the number of measured wavelengths 
                        and 'm' is the number of acquired spectra

        basis_spectra : the matrix containing the set of basis spectra

        pca_based     : boolean to decide if we need to center the matrix X

    Returns:

        basis_spectraNew       : the new set of basis spectra
        coefficients           : the new set of associated coefficients

    """

    if pca_based:

        X_mean    = np.mean(X, axis=1,keepdims=True)
        X         = X - X_mean

    first_spectrum = X[:,0]

    c1        = first_spectrum.dot(basis_spectra[:,0])
    c2        = first_spectrum.dot(basis_spectra[:,1])

    rotAngle      = angleFromCathets(c1,c2)

    rotM          = get_2d_counterclockwise_rot_matrix(rotAngle)

    basis_spectraNew = np.dot(basis_spectra[:,:2], rotM)
    coefficients     = np.dot(basis_spectraNew.T, X)

    return basis_spectraNew, coefficients

def rotate_three_basis_spectra(X,basis_spectra,pca_based=False):

    """
    Create a new basis spectra using a linear combination from the first, second and third basis spectra

    Requires:

        X             : the raw data matrix of size n*m, where 'n' is the number of measured wavelengths 
                        and 'm' is the number of acquired spectra

        basis_spectra : the matrix containing the set of basis spectra

        pca_based     : boolean to decide if we need to center the matrix X

    Returns:

        basis_spectra       : the new set of basis spectra
        coefficients_subset : the new set of associated coefficients

    """

    if pca_based:

        X_mean    = np.mean(X, axis=1,keepdims=True)
        X         = X - X_mean

    first_spectrum = X[:,0]

    c1        = first_spectrum.dot(basis_spectra[:,0])
    c2        = first_spectrum.dot(basis_spectra[:,1])
    c3        = first_spectrum.dot(basis_spectra[:,2])

    zAngle = angleFromCathets(c1,c2)
    yAngle = angleFromCathets(np.sqrt(c1**2+c2**2),c3)

    rotZaxis   = get_3d_counterclockwise_rot_matrix_around_z_axis(zAngle)
    rotYaxis   = get_3d_clockwise_rot_matrix_around_y_axis(yAngle)

    basisZrot         = np.dot(basis_spectra[:,:3], rotZaxis)
    basis_spectraNew  = np.dot(basisZrot, rotYaxis)
    coefficients      = np.dot(basis_spectraNew.T, X)

    return basis_spectraNew, coefficients

def reconstruct_spectra(basis_spectra,coefficients,X=None,pca_based=False):

    """
    Reconstruct the original spectra based on the set of basis spectra and the associated coefficients 

    Requires:


        basis_spectra           : the matrix containing the set of basis spectra
        coefficients_subset     : the associated coefficients of each basis spectrum

        pca_based     : boolean to decide if we need to extract the mean from the the X raw data matrix
        X             : only used pca_based equals TRUE! 
                        X is the raw data matrix of size n*m, where 
                        'n' is the number of measured wavelengths and 
                        'm' is the number of acquired spectra

    Returns:

        fitted       : the reconstructed matrix which should be close the original raw data

    """

    fitted =   (basis_spectra @ coefficients)

    # Add the mean, if needed
    if pca_based:

        X_mean    = np.mean(X, axis=1,keepdims=True)
        fitted    = fitted + X_mean

    return fitted

def explained_variance_from_orthogonal_vectors(vectors,coefficients,total_variance):

    """
    Useful to get the percentage of variance, not in the coordinate space provided by PCA/SVD, 
    but against a different set of (rotated) vectors.

    Input:
            - vectors        :   numpy matrix of size n*m, the columns contain a set of orthogonal vectors
            - coefficients   :   numpy matrix of size m*z, the rows    contain a set of associated coefficients
            - total_variance :   float, total variance of the original data (mean subtracted if we performed PCA...)

    Output:

            - a list containing the amout of explained variance by each orthogonal vector

    """

    explained_variance = []

    for i in range(vectors.shape[1]):

        a = np.linalg.norm(coefficients[i,:])**2
        b = np.linalg.norm(vectors[:, i])**2  

        explained_variance.append(a / b)

    return 100 * np.cumsum(explained_variance) / total_variance

def apply_pca(X):

    X_mean    = np.mean(X, axis=0)
    X         = X - X_mean

    # Decide if the spectra have a maximum or a minimum as a peak 
    maxV_abs = np.abs( np.max(X[:,4:-4]) )
    minV_abs = np.abs( np.min(X[:,4:-4]) )

    positive_peak = maxV_abs > minV_abs

    # compute the covariance matrix
    cov_mat   = np.cov(X , rowvar = False)

    # find the eigen vectors and associated eigen values
    eigen_values , eigen_vectors = np.linalg.eigh(cov_mat)

    #sort the eigenvalues in descending order
    sorted_index = np.argsort(eigen_values)[::-1]
     
    sorted_eigenvalue = eigen_values[sorted_index]

    #similarly sort the eigenvectors 
    sorted_eigenvectors = eigen_vectors[:,sorted_index]

    # compute the total variance
    total_eigenvalues = np.sum(sorted_eigenvalue)

    # compute the explained variance
    exp_var_pca = (sorted_eigenvalue / total_eigenvalues * 100)

    # compute the cumulative explained variance
    cum_sum_eigenvalues = np.cumsum(exp_var_pca)

    principal_components = sorted_eigenvectors

    a_is = []

    for i in range(principal_components.shape[1]):

        def coefficients_bi(column):
        # Your custom logic here
            return principal_components[:,i].dot(column)

        a_i = np.apply_along_axis(coefficients_bi, axis=1, arr=X)

        a_is.append(a_i)

    coefficients = np.array(a_is)

    return cum_sum_eigenvalues, principal_components, coefficients