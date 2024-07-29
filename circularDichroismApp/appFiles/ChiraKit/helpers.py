import pandas as pd 
import numpy  as np 
import math
import re

def round_to_significant_digits(arr, digits):
    rounded_arr = np.zeros_like(arr)
    for i, value in np.ndenumerate(arr):
        if np.abs(value) < np.finfo(float).eps:
            rounded_arr[i] = value  # Keep small numbers close to zero as is
        else:
            magnitude = np.power(10, digits - np.floor(np.log10(np.abs(value))) - 1)
            rounded_arr[i] = np.around(value * magnitude) / magnitude
    return rounded_arr

def string_to_units(string):

    '''
    Convert the string into one of the followings:

        'milliabsorbance', 'molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity',
        'meanUnitMolarEllipticity', 'meanUnitMolarExtinction'

    Useful function to reverse the output from the function 'workingUnits2ProperLabel' from 'helpers_plotting.R'

    E.g. 
    Input - 'Millidegrees (m°)'      /   Output - millidegrees
    Input - 'Milliabsorbance (mΔA)'  /   Output - milliabsorbance
    Input - 'Δε ...'                /   Output - meanUnitMolarExtinction
    Input - '... MRE ...'           /   Output - meanUnitMolarEllipticity

    '''

     # Set milidegrees as default output
    units = 'millidegrees'

    string = string.lower()

    unit_mappings = {
    'absorbance'            : 'absorbance',
    'degrees'               : 'degrees',
    'mue'                   : 'meanUnitMolarEllipticity',
    'δε'                    : 'meanUnitMolarExtinction',
    'delta epsilon'         : 'meanUnitMolarExtinction',
    'molar extinction'      : 'molarExtinction',
    'molar ellipticity'     : 'molarEllipticity'
    }

    for key, value in unit_mappings.items():
        if key in string:
            units = value
            break  # Exit the loop once a match is found

    if 'milli' in string or 'mili' in string:

        units = 'milli' + units

    return units

def guess_input_units_from_metadata_dictionary(metadata_dic):

    for key, value in metadata_dic.items():

        key_lower = key.lower()

        # Known case, exported files from ChiraKit
        if 'Units of the CD signal' in key:

            return string_to_units(value)

    # Set 'Millidegrees' as default output
    return 'millidegrees'

def guess_parameter_from_metadata_dictionary(metadata_dic,parameter_names):

    '''
    parameter_names should be a list of strings 
    E.g., ['concentration','conc'], 
    '''

    for key, value in metadata_dic.items():

        key_lower = key.lower()

        for parameter in parameter_names:

            if parameter in key_lower:

                try:

                    return float(value)

                except:

                    pass

    return 0

def filter_matrix_by_vector(np_matrix,np_vector,min_value,max_value):

    """
    Filter the matrix using a certain range

    Requires: 
        
        1) The matrix 'np_matrix' of dimensions
    n*m where n matches the length of np_vector 
        2) The vector 'np_vector' of length n
        3) The lower bound 'min_value'
        4) The upper bound 'max_value'

    Returns the filtered matrix 
    """

    np_tog = np.column_stack((np_matrix, np_vector))
    tot    = np_matrix.shape[1]
    np_tog = np_tog[np_tog[:,tot]  >= min_value]
    np_tog = np_tog[np_tog[:,tot]  <= max_value]
    np_tog = np.delete(np_tog, tot, 1)

    return(np_tog)

def filter_vector_by_values(np_vector,min_value,max_value):

    """
    Filter the vector using a certain range

    Requires: 
        
        1) The vector 'np_vector' of length n.
        2) The lower bound 'min_value'
        3) The upper bound 'max_value'


    Returns the filtered vector 
    """

    np_temp = np_vector[np_vector   >= min_value]
    np_temp = np_temp[np_temp       <= max_value]

    return(np_temp)

def get_temp_at_maximum_of_derivative(temps,signal_derivative):

    tms_derivative = np.take(temps, np.argmax(signal_derivative,axis=0)) 
    return tms_derivative

def get_temp_at_minimum_of_derivative(temps,signal_derivative):

    tms_derivative = np.take(temps, np.argmin(signal_derivative,axis=0)) 
    return tms_derivative

def extract_words(input_string):

    # Pattern to match words with only alphabet characters or underscores (without hyphens)
    pattern = r'\b[A-Za-z_]+\b'

    words = re.findall(pattern, input_string)

    return words

def clean_function_text(text):

    # Use numpy notation
    cleaned_text = text.replace('e^', 'exp')
    cleaned_text = cleaned_text.replace('exp(','np.exp(') 
    cleaned_text = cleaned_text.replace('^', '**') 
    cleaned_text = cleaned_text.replace('log(','np.log(')
    cleaned_text = cleaned_text.replace('sqrt(','np.sqrt(')

    return cleaned_text

## check or remove code below!!!
def check_good_parameters_estimation(params,low_bound,high_bound,params_name):

    """
        
    Check that the estimated parameters are far from the bounds, 
    if not there is a problem with the fitting.

    For the first  4 params 'kN', 'bN', 'kU', 'bU' we will normalize low_bound - high_bound range to 0-1
    and then verify that the parameters lie in the interval 0.02-0.98.
        
    """

    low_bound  = np.array(low_bound)
    high_bound = np.array(high_bound)
    params     = np.array(params)

    params_normalized = (params[:4] - low_bound[:4] ) / (high_bound[:4] - low_bound[:4])

    lie_in_correct_interval   = np.logical_and(params_normalized < 0.98,
        params_normalized > 0.02)

    """

    For Tm or T1 and T2 and T_onset or T_onset1 and T_onset2  we will check that they are 1 degree from the boundaries

    """

    for temp_param_name in ["Tm","T1","T2","T_onset","T_onset1","T_onset2"]:

        if temp_param_name in params_name:

            position_of_tm = params_name.index(temp_param_name) 
            tm             = params[position_of_tm]
            tm_bound_low   = low_bound[position_of_tm]
            tm_bound_up    = high_bound[position_of_tm]

            lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([tm > (tm_bound_low+1) and tm < (tm_bound_up-1)]))

    """

    For Tm or dHm, dHm1 and dHm2 we will check that they are between 2.5 and 750 kcal/mol. 10466 and 3098230 in Joules

    """

    for dh_param_name in ["dHm","dHm1","dHm2"]:

        if dh_param_name in params_name:

            position_of_dh = params_name.index(dh_param_name) 
            dh             = params[position_of_dh]
            lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([dh > 10466 and dh < 3098230]))

    """

    For Ki we will check that it lies between 1e-3 and 1e3

    """

    if "Ki" in params_name:
        position_of_ki = params_name.index("Ki") 
        ki             = params[position_of_ki]
        lie_in_correct_interval = np.append(lie_in_correct_interval, np.array([ki >= 0.001 and ki < 1000]))


    return all(lie_in_correct_interval)
