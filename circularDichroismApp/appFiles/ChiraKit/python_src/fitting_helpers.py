import numpy  as np 
import pandas as pd

from scipy.signal    import savgol_filter
from scipy.linalg    import pinv
from scipy.stats     import linregress

from helpers import *

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

    # Remove consecutive wavelengths
    if len(wavelength_useful) > 8:

        subset = []

        for wavelength in wavelength_useful:

            if not subset or wavelength != subset[-1] + 1:

                subset.append(wavelength)

        wavelength_useful = subset

    return wavelength_useful

def fit_baselines(signal_matrix,measurement_factor,threshold):

    '''
    Estimate the temperature / chemical dependence of the native/unfolded state signal.
    If you use the temperature as measurement_factor variable, you should use Kelvin units!

    Input:
        - signal_matrix      of dimensions n*m
        - measurement_factor of length m                             (e.g., temperature data)
        - threshold (float), to select certain columns of the matrix (e.g., 8°C)
    Output:
        - intercept and slope of the native and unfolded states
    '''

    max_native   = np.min(measurement_factor) + threshold
    min_unfolded = np.max(measurement_factor) - threshold

    native_signal       = filter_matrix_by_vector((signal_matrix).T,measurement_factor,
        np.min(measurement_factor),max_native)

    native_measurement_factor      = filter_vector_by_values(measurement_factor,np.min(measurement_factor),max_native)

    unfolded_signal     = filter_matrix_by_vector((signal_matrix).T,measurement_factor,
        min_unfolded,np.max(measurement_factor))

    unfolded_measurement_factor    = filter_vector_by_values(measurement_factor,min_unfolded,np.max(measurement_factor))

    fitNative = []

    for i in range(native_signal.shape[1]):

        mask = ~np.isnan(native_signal[:,i])
        x    = native_measurement_factor[mask]
        y    = native_signal[:,i][mask]

        fitNative.append(linregress(x,y))

    fitUnfolded = []

    for i in range(unfolded_signal.shape[1]):

        mask = ~np.isnan(unfolded_signal[:,i])
        x    = unfolded_measurement_factor[mask]
        y    = unfolded_signal[:,i][mask]

        fitUnfolded.append(linregress(x,y))

    kN           = np.array( [c.slope     for c in fitNative] )
    bN           = np.array( [c.intercept for c in fitNative] )
    
    kU           = np.array( [c.slope     for c in fitUnfolded] )
    bU           = np.array( [c.intercept for c in fitUnfolded] )

    return bN, kN, bU, kU

def signal_2d_matrices_to_full_matrix(signal_lst2d,rows,cols):

    matrix = np.full((rows, cols), np.nan) 

    row, col = 0,0

    for mat in signal_lst2d:

        mat_rows, mat_cols = mat.shape

        matrix[row:(row+mat_rows),col:(col+mat_cols)] = mat

        row += mat_rows
        col += mat_cols

    return matrix

def update_params_and_bounds(p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded,n1,n2,n3):

    '''
    To delete the pre and post transition slopes parameters
    n1 is the number of shared parameters
    n2 is the number of local  parameters (without the slopes)
    n3 is the total number of datasets
    '''

    if not fitSlopeUnfolded:

        p0, low_bound, high_bound = p0[:-n3], low_bound[:-n3], high_bound[:-n3]

    if not fitSlopeNative:

        l1, l2 = (n1+n3*n2), (n1+n3*(n2+1))

        indices_to_delete = [x for x in range(l1,l2)]

        p0           = np.delete(p0, indices_to_delete) 
        low_bound    = np.delete(low_bound, indices_to_delete) 
        high_bound   = np.delete(high_bound, indices_to_delete)

    return p0, low_bound, high_bound

def extend_bounds(fit_params,n_params_to_discard,low_bound,high_bound,threshold,extension):

    diffMin = (fit_params - low_bound)[n_params_to_discard:]
    diffMax = (high_bound - fit_params)[n_params_to_discard:]

    re_fit = False

    for i in range(len(diffMin)):

        if diffMin[i] < threshold:
            low_bound[i+n_params_to_discard] -= extension
            re_fit                        = True

        if diffMax[i] < threshold:
            high_bound[i+n_params_to_discard] += extension
            re_fit                         = True

    return re_fit, low_bound, high_bound

def generate_initial_params_combinations(n_params,paramsSign,
    logSpaceLimitLeft=-3,logSpaceLimitRight=3):

    all_combis = []

    # Set the max number of possible initial combinations

    combisLimitLst = 420 * (1.55 ** np.arange(1,6))

    n_params_temp = n_params
    if n_params > 5:
        n_params_temp = 5

    combisLimit = combisLimitLst[n_params_temp-1]
    combisLimit = combisLimit * ( 1.7** (len([x for x in paramsSign if x != 'unknown']) ) )

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

def model_n(model):

    model_map = {'dimer': 2,'trimer': 3,'tetramer': 4}

    n = 1  # default value

    for key in model_map:
        if key in model.lower():
            n = model_map[key]
            break  # exit loop once a match is found

    return n

def fn_two_state_monomer(K):
    '''
    Given the equilibrium constant K, of N <-> U, return the fraction of folded protein
    '''
    return (1/(1 + K))

def fu_two_state_dimer(K,C):
    '''
    Given the equilibrium constant K, of N2 <-> 2U, 
    and the concentration of dimer equivalent C, return the fraction of unfolded protein
    '''
    return solve_one_root_quadratic(4*C, K, -K)

def fu_two_state_trimer(K,C):
    '''
    Given the equilibrium constant K, of N3 <-> 3U, 
    and the concentration of trimer equivalent C, return the fraction of unfolded protein
    '''
    p = K/27/np.square(C)
    return solve_one_root_depressed_cubic(p,-p)

def fu_two_state_tetramer(K,C):
    '''
    Given the equilibrium constant K, of N4 <-> 4U, 
    and the concentration of tetramer equivalent C, return the fraction of folded protein
    '''

    A = 1
    D = K/256/np.power(C,3)
    E = -D

    b = D/A
    c = E/A

    P = -c
    Q = -np.square(b)/8

    R = -Q/2 + np.sqrt(np.square(Q)/4+P**3/27)

    U = np.cbrt(R)
    y = U-P/(3*U)
    W = np.sqrt(2*y)

    x4 = 0.5*(-W+np.sqrt(-(2*y-2*b/W)))  

    x4_sel = np.logical_and(np.greater(x4,0),np.less(x4,1.01))

    fu = x4_sel*np.nan_to_num(x4,nan=0.0)

    return fu

def fi_three_state_tetramer_monomeric_intermediate(K1,K2,Ct):
    '''
    Given the equilibrium constant K1, of N4 <-> 4I, K2, of I <-> U,
    and the concentration of tetramer equivalent Ct, return the fraction of intermediate
    '''
    Pt = Ct*4

    A = 4*(Pt**3)/K1
    D = 1+K2
    E = -1

    b = D/A
    c = E/A

    P = -c
    Q = -np.square(b)/8

    R = -Q/2 + np.sqrt(np.square(Q)/4+P**3/27)

    U = np.cbrt(R)
    y = U-P/(3*U)
    W = np.sqrt(2*y)

    x4 = 0.5*(-W+np.sqrt(-(2*y-2*b/W)))  

    x4_sel = np.logical_and(np.greater(x4,0),np.less(x4,1.01))

    fi = x4_sel*np.nan_to_num(x4,nan=0.0)

    return fi

def fi_three_state_dimer_monomeric_intermediate(K1,K2,C):
    '''
    Given the equilibrium constant K1, of N2 <-> 2I, K2, of 2I <-> 2U
    and the concentration of dimer equivalent C, return the fraction of intermediate
    '''
    return solve_one_root_quadratic(4*C,K1*(1+K2),-K1)

def fu_three_state_dimer_dimeric_intermediate(K1,K2,C):
    '''
    Given the equilibrium constant K1, of N2 <-> I2, K2, of I2 <-> 2U
    and the concentration of dimer equivalent C, return the fraction of unfolded protein
    '''
    return solve_one_root_quadratic(4*C*(1+K1), K1*K2, -K1*K2)

def fi_three_state_dimer_dimeric_intermediate(fu,K2,C):
    '''
    Given the fraction of unfolded protein fu, the equilibrium constant K2, of I2 <-> 2U,
    and the concentration of dimer equivalent C, return the fraction of intermediate
    '''
    return 4*np.square(fu)*C/K2

def fi_three_state_trimer_monomeric_intermediate(K1,K2,C):
    '''
    Given the equilibrium constant K1, of N3 <-> 3I, K2, of 3I <-> 3U
    and the concentration of trimer equivalent C, return the fraction of unfolded protein
    '''
    p = K1*(1+K2)/27/np.square(C)
    q = -K1/27/np.square(C)

    return solve_one_root_depressed_cubic(p,q)

def fu_three_state_trimer_trimeric_intermediate(K1,K2,C):
    '''
    Given the equilibrium constant K1, of N3 <-> I3, K2, of I3 <-> 3U
    and the concentration of trimer equivalent C, return the fraction of unfolded protein
    '''
    p = K1*K2 / (27*np.square(C)*(1+K1))
    q = -p

    return solve_one_root_depressed_cubic(p,q)

def fi_three_state_trimer_trimeric_intermediate(fu,K2,C):
    '''
    Given the fraction of unfolded protein fu, the equilibrium constant K2, of I3 <-> 3U,
    and the concentration of trimer equivalent C, return the fraction of intermediate
    '''
    return 27*np.square(C)*(fu**3) / K2

def linear_signal(T,intercept,slope):

    '''
    For the linear dependence of the pre and post transition states. 
    '''

    return intercept + temperature_to_kelvin(T)*slope

def generate_params_dfs(column_names,fit_params,fit_errors,fitSlopeNative,fitSlopeUnfolded):

    fit_params = np.array(fit_params)
    fit_errors = np.array(fit_errors)

    # Convert the NumPy array to a Pandas DataFrame with column names
    df_fit_params = pd.DataFrame(fit_params, columns=column_names)
    df_fit_errors = pd.DataFrame(fit_errors, columns=column_names)

    # Remove non necessary columns
    if not fitSlopeNative:

        df_fit_params = df_fit_params.drop('kN', axis=1)
        df_fit_errors = df_fit_errors.drop('kN', axis=1)
        column_names.remove('kN')

    if not fitSlopeUnfolded:

        df_fit_params = df_fit_params.drop('kU', axis=1)
        df_fit_errors = df_fit_errors.drop('kU', axis=1)
        column_names.remove('kU')

    # Convert to numeric all columns except the last two
    df_fit_params[ column_names[:-2] ] = df_fit_params[ column_names[:-2] ].apply(pd.to_numeric)
    df_fit_errors[ column_names[:-2] ] = df_fit_errors[ column_names[:-2] ].apply(pd.to_numeric)

    return df_fit_params, df_fit_errors

def split_all_params_two_state(fit_params,totalDataSets,fitSlopeNative,fitSlopeUnfolded):

    param1            = fit_params[0] # Temperature of melting /  M value
    param2            = fit_params[1] # Enthalpy of unfolding  /  D50 - concentration at which the sample is 50 % unfolded

    interceptsFolded    = fit_params[2:(2+totalDataSets)]
    interceptsUnfolded  = fit_params[(2+totalDataSets):(2+totalDataSets*2)]

    slopesFolded   = fit_params[(2+totalDataSets*2):(2+totalDataSets*3)] if fitSlopeNative   else np.full(len(interceptsFolded), 0)
    slopesUnfolded = fit_params[(len(fit_params)-totalDataSets):]        if fitSlopeUnfolded else np.full(len(interceptsFolded), 0)

    return [param1, param2, interceptsFolded, slopesFolded, interceptsUnfolded, slopesUnfolded]

def split_all_params_three_state(fit_params,totalDataSets,fitSlopeNative,fitSlopeUnfolded,withCp=False):

    p1, p2, p3, p4   = fit_params[:4] 
    start            = 4

    if withCp:

        p5     = fit_params[4] 
        start += 1

    interceptsFolded       = fit_params[(start+totalDataSets*0):(start+totalDataSets*1)]
    interceptsUnfolded     = fit_params[(start+totalDataSets*1):(start+totalDataSets*2)]
    interceptsIntermediate = fit_params[(start+totalDataSets*2):(start+totalDataSets*3)]

    slopesFolded   = fit_params[(start+totalDataSets*3):(start+totalDataSets*4)] if fitSlopeNative   else np.full(len(interceptsFolded), 0)
    slopesUnfolded = fit_params[(len(fit_params)-totalDataSets):]                if fitSlopeUnfolded else np.full(len(interceptsFolded), 0)

    if withCp:

        return [p1, p2, p3, p4, p5, interceptsFolded, slopesFolded, interceptsUnfolded, slopesUnfolded, interceptsIntermediate]
    
    else:

        return [p1, p2, p3, p4, interceptsFolded, slopesFolded, interceptsUnfolded, slopesUnfolded, interceptsIntermediate]

def generate_bounds_df(
    shared_params_names,local_params_names,wavelengths,fitSlopeNative,fitSlopeUnfolded,
    low_bound,high_bound,params_splt,errors_splt,
    last_col_name,dataset_name,df_fit_params,df_fit_errors):

    bounds_df2 = generate_bounds_df_indiv(
        shared_params_names,local_params_names,wavelengths,fitSlopeNative,fitSlopeUnfolded,
        low_bound,high_bound,errors_splt,
        last_col_name,dataset_name,df_fit_errors,True)

    bounds_df1 = generate_bounds_df_indiv(
            shared_params_names,local_params_names,wavelengths,fitSlopeNative,fitSlopeUnfolded,
        low_bound,high_bound,params_splt,
        last_col_name,dataset_name,df_fit_params,False)

    return pd.merge(bounds_df1, bounds_df2, how='left')

def generate_bounds_df_indiv(
    shared_params_names,local_params_names,wavelengths,fitSlopeNative,fitSlopeUnfolded,
    low_bound,high_bound,params_splt,
    last_col_name,dataset_name,df_fit_params,errors=False):

    nGlob = len(shared_params_names)
    nLoc  = len(local_params_names)

    # Create a copy to avoid in place modification
    names  = shared_params_names[:] 

    for lpn in local_params_names:
        names += [lpn  for _ in wavelengths]
    
    if fitSlopeNative:
        names += ['kN' for _ in wavelengths]
    
    if fitSlopeUnfolded:
        names += ['kU' for _ in wavelengths]

    data_ori  = np.concatenate((["All" for _ in range(nGlob)],np.tile(wavelengths, nLoc + fitSlopeNative + fitSlopeUnfolded)))
    bounds_df = pd.DataFrame({'Param': names,'Lower bound':low_bound,'Upper bound':high_bound,last_col_name:data_ori,'Dataset':dataset_name})
    
    # Create a mask for the rows that should not be sorted

    # Split the DataFrame into unchanged and to be sorted parts
    unchanged_df = bounds_df[ bounds_df.index < nGlob  ].copy()
    to_sort_df   = bounds_df[ bounds_df.index >= nGlob ].copy()

    if errors:
        colname = 'Relative error (%)'
    else:
        colname = 'Fitted value'

    unchanged_df[colname] = params_splt[:nGlob]

    # Sort the part to be sorted by the 'WL / Basis spectrum' column
    sorted_df = to_sort_df.sort_values(by=last_col_name)

    # Melt the DataFrame to reshape it
    df_long = pd.melt(df_fit_params, 
                       id_vars   = [last_col_name, 'Dataset'], 
                       var_name  = 'Param',
                       value_name= colname)

    df_merged = pd.merge(sorted_df, df_long, how='left')

    # Concatenate the unchanged part with the sorted part
    bounds_df = pd.concat([unchanged_df, df_merged], ignore_index=True)

    # Get the list of columns
    columns = bounds_df.columns.tolist()

    # Move the last column to the third place
    # Remove the last column and insert it at the third position
    last_column = columns.pop()     # Remove the last column
    columns.insert(2, last_column)  # Insert it at the third position (index 2)

    # Reorder the DataFrame
    return bounds_df[columns]

def merge_fractions_dictionaries(fractions_lst,unique_dimer_concs):

    merged_dict = {}

    for index, d in enumerate(fractions_lst):

        conc = np.round(unique_dimer_concs[index]*1e6,2)

        # Modify keys by adding the index
        updated_dict = {f'{key}_{conc}': value for key, value in d.items()}
        # Merge the updated dictionary into the final dictionary
        merged_dict.update(updated_dict)

    return merged_dict

def generate_predicted_matrix(predictedLst,experimental_data):

    '''
    Requires:
                - 'predictedLst'      : List of y-values, row-wise
                - 'experimental_data' : numpy matrix, as many rows as length(predictedLst). It may contain NAs
    '''

    signal_predicted      = np.array(predictedLst) # List of y-values, row-wise

    # Create a mask for positions where experimental_data has NaN values
    nan_mask = np.isnan(experimental_data)

    # Assign NaN to the same positions in signal_predicted
    signal_predicted[nan_mask] = np.nan

    return signal_predicted
