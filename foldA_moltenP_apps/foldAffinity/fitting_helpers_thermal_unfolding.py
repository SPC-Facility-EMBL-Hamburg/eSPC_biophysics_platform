import numpy  as np
import pandas as pd
import scipy

from scipy.optimize     import curve_fit

def fit_single_thermal_curve(temps, fluo, cp):
    '''
    Local fit of thermal curve

    Args:
        temps: Temperatures
        fluo:  Fluorescence intensities
        cp:    DeltaCp

    Returns: None
    '''

    init_dH = 150
    window_size = int(len(temps) / 10)
    (init_bottom_slope, init_bottom_intercept, junk1, junk2, junk3) = scipy.stats.linregress(
        temps[0:window_size] + 273.15, fluo[0:window_size])
    tse = len(temps) - 1
    (init_top_slope, init_top_intercept, junk1, junk2, junk3) = scipy.stats.linregress(
        temps[tse - window_size:tse] + 273.15, fluo[tse - window_size:tse])
    fluo_midpoint = (np.max(fluo) + np.min(fluo)) / 2

    # Determine initial melting temperatures
    init_Tm = temps[np.argmin(np.abs(fluo - fluo_midpoint))].squeeze()

    # Starting values for fit
    p0 = (init_Tm, init_dH, init_top_intercept, init_bottom_intercept, init_top_slope, init_bottom_slope)

    # Define boundaries for fit
    low_bound, high_bound = [-np.inf] * 6, [np.inf] * 6
    # Define temps
    low_bound[0], high_bound[0]  = np.min(temps) + 0, np.max(temps) - 0
    # Define DeltaH
    low_bound[1] = 0        

    # Function with fixed Cp
    model = single_thermal_curve_woCp(cp)

    # Fit
    par, cov = curve_fit(model, temps, fluo, p0=p0, bounds=(low_bound, high_bound), max_nfev=1E4, method='trf')
    errors = np.sqrt(np.diag(cov))

    return par, errors

def single_thermal_curve(temperatures, Tm, dH, unfolded_intercept, folded_intercept, unfolded_slope, folded_slope, Cp):
    '''
    Uses parameters from local/global fit and outputs thermal unfolding curve
    Args:
        temperatures:
        Tm: melting temperatures
        dH: unfolding enthalpy in kcal/mol
        unfolded_intercept:
        folded_intercept:
        unfolded_slope:
        folded_slope:
        Cp: deltaCp of thermal unfolding in kcal/(mol*K)
    Returns: None
    '''
    # Tm, dH, unfolded_intercept, folded_intercept, unfolded_slope, folded_slope, Cp = parameters
    T = temperatures + 273.15
    # Gas constant in kcal/(mol*K)
    R = 1.987E-3
    # dG and dH in kcal/mol; Cp in kcal/(mol*K)
    dG = dH * (1 - T / (Tm + 273.15)) - Cp * (Tm + 273.15 - T + T * np.log(T / (Tm + 273.15)))
    Ku = np.exp(-dG / (R * T))
    Y = (Ku / (1 + Ku)) * (unfolded_slope * T + unfolded_intercept) + (1 / (1 + Ku)) * (
                folded_slope * T + folded_intercept)
    return Y

def single_thermal_curve_woCp(Cp):
    '''
    Wrapper script that returns single_thermal_curve function with a fixed Cp
    Args:
        Cp: deltaCp of thermal unfolding
    Returns: function
    '''

    # This returns of single_thermal_curve function with a fixed Cp
    def return_func(temperatures, Tm, dH, unfolded_intercept, folded_intercept, unfolded_slope, folded_slope):
        return single_thermal_curve(temperatures, Tm, dH, unfolded_intercept, folded_intercept, unfolded_slope,
                                         folded_slope, Cp)

    return return_func

def global_thermal_curves(temperatures, *args):  # Tms, dHs, unfolded_intercepts, folded_intercepts, unfolded_slopes, folded_slopes, Cp):
    '''
    This is a function to calculate all curves at once
    It is used for the global fit
    The arguments have to be in the following order:
    Single Tms
    Single dHs
    Single unfolded_intercepts
    Single folder_intercepts
    Global unfolded_slope
    Global folded_slope
    '''
    # Number of "repeats"
    num_datasets = np.sum(temperatures == temperatures[0])
    Tms = list(args[:num_datasets])
    dHs = args[num_datasets:2 * num_datasets]
    unfolded_intercepts = args[2 * num_datasets:3 * num_datasets]
    folded_intercepts = args[3 * num_datasets:4 * num_datasets]
    unfolded_slope = args[4 * num_datasets]
    folded_slope = args[4 * num_datasets + 1]
    Cp = args[4 * num_datasets + 2]
    # New temperature
    single_temp = temperatures[:len(temperatures) // num_datasets]
    # Data is appended to a list
    all_data = []
    for i in range(num_datasets):
        T = single_temp + 273.15
        R = 1.987 / 1000
        if Tms[i] < 0:
            Tms[i] = 0
        dG = dHs[i] * (1 - T / (Tms[i] + 273.15)) - Cp * (Tms[i] + 273.15 - T + T * np.log(T / (Tms[i] + 273.15)))
        Ku = np.exp(-dG / (R * T))
        Y = (Ku / (1 + Ku)) * (unfolded_slope * T + unfolded_intercepts[i]) + (1 / (1 + Ku)) * (
                    folded_slope * T + folded_intercepts[i])
        all_data.append(Y)
    return np.hstack(all_data)

def global_thermal_curves_woCp(Cp):
    '''
    Wrapper script that returns global_thermal_curves function with a fixed Cp
    Args:
        Cp: deltaCp of thermal unfolding
    Returns: function
    '''

    def return_func(temperatures, *args):
        args = (*args, Cp)
        return global_thermal_curves(temperatures, *args)

    return return_func

def single_thermal_curve_Kd(temperatures, Tm, dH0, unfolded_intercept, folded_intercept, unfolded_slope, folded_slope, Cp, Kd, conc):
    '''
    Uses parameters from local/global fit and outputs thermal unfolding curve
    for a particular Kd and Ku0
    Args:
        temperatures:
        Tm: melting temperatures
        dH: unfolding enthalpy in kcal/mol
        unfolded_intercept:
        folded_intercept:
        unfolded_slope:
        folded_slope:
        Cp: deltaCp of thermal unfolding in kcal/(mol*K)
        Kd: dissociation constant in M^-1
        dH0: unfolding enthalpy for pure protein (no ligand)
        
    Returns: None
    '''
    # Tm, dH, unfolded_intercept, folded_intercept, unfolded_slope, folded_slope, Cp = parameters
    T = temperatures + 273.15
    # Gas constant in kcal/(mol*K)
    R = 1.987E-3
    # dG and dH in kcal/mol; Cp in kcal/(mol*K)
    dG0 = dH0 * (1 - T / (Tm + 273.15)) - Cp * (Tm + 273.15 - T + T * np.log(T / (Tm + 273.15)))
    Ku0 = np.exp(-dG0 / (R * T))
    Y = (Ku0 / (Ku0 + conc/Kd  + 1)) * (unfolded_slope * T + unfolded_intercept) + \
        (1 - Ku0 / (Ku0 + conc/Kd  + 1)) * (folded_slope * T + folded_intercept)
    return Y