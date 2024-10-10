import numpy  as np 

from scipy.optimize  import curve_fit
from scipy.integrate import solve_ivp

from constants       import *
from helpers         import *
from fitting_helpers import *

'''
Useful references for unfolding models:
    - Rumfeldt, Jessica AO, et al. "Conformational stability and folding mechanisms of dimeric proteins." Progress in biophysics and molecular biology 98.1 (2008): 61-84.
    - Bedouelle, Hugues. "Principles and equations for measuring and interpreting protein stability: From monomer to tetramer." Biochimie 121 (2016): 29-37.
    - Mazurenko, Stanislav, et al. "Exploration of protein unfolding by modelling calorimetry data from reheating." Scientific reports 7.1 (2017): 16321.

All thermodynamic parameters are used in kcal mol units

Unfolding functions for monomers have an argument called 'extra_arg' that is not used. 
This is because unfolding functions for oligomers require the protein concentration in that position
'''

def eq_constant(T,DH1,T1,Cp):

    '''
    T1 is the temperature at which ΔG(T) = 0
    ΔH1, the variation of enthalpy between the two considered states at T1
    Cp the variation of calorific capacity between the two states
    '''

    T  = temperature_to_kelvin(T)
    T1 = temperature_to_kelvin(T1)

    DG = DH1*(1 - T/T1) - Cp*(T1 - T + T*np.log(T/T1)) 
    K  = np.exp(-DG / (R_gas * T))

    return K

def arrhenius(T, Tf, Ea):
    """
    Arrhenius equation: defines dependence of reaction rate constant k on temperature
    In this version of the equation we use Tf (a temperature of k=1)
    to get rid of instead of pre-exponential constant A

    T, Tf, must be given in Kelvin, Ea in kcal units
    """

    T  = temperature_to_kelvin(T)
    Tf = temperature_to_kelvin(Tf)

    return np.exp(-Ea / R_gas * (1 / T - 1 / Tf))

def two_state_thermal_unfold_curve(T,Tm,dh,bN,kN,bU,kU,extra_arg,Cp=0):

    '''
    N ⇔ U
    '''

    K   = eq_constant(T,dh,Tm,Cp)
    fn  = fn_two_state_monomer(K)

    return fn*linear_signal(T,bN,kN) + (1-fn)*linear_signal(T,bU,kU)

def two_state_thermal_unfold_curve_dimer(T,Tm,dh,bN,kN,bU,kU,C,Cp=0):

    """
    N2 ⇔ 2U   C is the total concentration (M) of the protein in dimer equivalent.
    """

    K  = eq_constant(T,dh,Tm,Cp)
    fu = fu_two_state_dimer(K,C)

    return (1-fu)*linear_signal(T,bN,kN) + fu*linear_signal(T,bU,kU)*2

def two_state_thermal_unfold_curve_trimer(T,Tm,dh,bN,kN,bU,kU,C,Cp=0):

    """
    N3 ⇔ 3U   C is the total concentration (M) of the protein in trimer equivalent.
    """

    K  = eq_constant(T,dh,Tm,Cp)
    fu = fu_two_state_trimer(K,C)

    return (1-fu)*linear_signal(T,bN,kN) + fu*linear_signal(T,bU,kU)*3

def two_state_thermal_unfold_curve_tetramer(T,Tm,dh,bN,kN,bU,kU,C,Cp=0):

    """
    N4 ⇔ 4U   C is the total concentration (M) of the protein in tetramer equivalent.
    """

    K  = eq_constant(T,dh,Tm,Cp)
    fu = fu_two_state_tetramer(K,C)

    return (1-fu)*linear_signal(T,bN,kN) + fu*linear_signal(T,bU,kU)*4

def map_two_state_model_to_signal_fx(model):

    signal_fx_map = {
    'Monomer':  two_state_thermal_unfold_curve,
    'Dimer':    two_state_thermal_unfold_curve_dimer,
    'Trimer':   two_state_thermal_unfold_curve_trimer,
    'Tetramer': two_state_thermal_unfold_curve_tetramer
    }

    return signal_fx_map.get(model)

def fit_thermal_unfolding(listOfTemperatures,listOfSignals,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,signal_fx,Cp=0,listOfOligomerConc=None):

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
    allSignal       = concat_signal_lst(listOfSignals)
    
    def thermal_unfolding(dummyVariable,*args):

        '''
        Calculate the thermal unfolding profile of many curves at the same time

        Requires:

            - The 'listOfTemperatures' containing each of them a single dataset

        The other arguments have to be in the following order:

            - Global melting temperature    
            - Global enthalpy of unfolding 
            - Single intercepts, folded
            - Single intercepts, unfolded
            - Single slopes, folded
            - Single slopes, unfolded

        Returns:

            The melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
                slopes and intercept of the folded and unfolded states

        '''

        totalDataSets = len(listOfTemperatures)
        Tm,dh         = args[:2] # Temperature of melting, Enthalpy of unfolding

        interceptsFolded    = args[2:(2+totalDataSets)]
        interceptsUnfolded  = args[(2+totalDataSets):(2+totalDataSets*2)]

        if fitSlopeNative:

            slopesFolded   = args[(2+totalDataSets*2):(2+totalDataSets*3)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,T in enumerate(listOfTemperatures):

            bN, bU = interceptsFolded[i],   interceptsUnfolded[i]
            
            kN = 0 if not fitSlopeNative          else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded        else slopesUnfolded[i]
            C  = 0 if listOfOligomerConc is None  else listOfOligomerConc[i]

            Y = signal_fx(T,Tm,dh,bN,kN,bU,kU,C,Cp)

            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(thermal_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds))

    return global_fit_params, cov

def unfolding_curve_monomer_monomeric_intermediate(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,extra_arg,Cp1=0,CpTh=0):

    """
    Three states reversible unfolding N <-> I <-> U
    """

    A = eq_constant(T,DH1,T1,Cp1)
    B = eq_constant(T,DH2,T2,CpTh-Cp1)

    den = (1+A+A*B)

    xN, xI, xU = 1/den, A/den, A*B/den

    return xN*linear_signal(T,bN,kN) + xI*bI + xU*linear_signal(T,bU,kU)

def unfolding_curve_dimer_monomeric_intermediate(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1=0,CpTh=0):

    """
    N2 ⇔ 2Ι ⇔ 2U Three-state unfolding with a monomeric intermediate
    C = concentration in dimer equivalent
    CpTotal = Cp1 + 2*Cp2
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,(CpTh-Cp1)/2)
    
    fi = fi_three_state_dimer_monomeric_intermediate(K1,K2,C)
    fu = fi*K2

    return (1-fu-fi)*linear_signal(T,bN,kN) + fi*bI*2 + fu*linear_signal(T,bU,kU)*2

def unfolding_curve_trimer_monomeric_intermediate(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1=0,CpTh=0):

    """
    N3 ⇔ 3Ι ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration of the trimer equivalent
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,(CpTh-Cp1)/3)   # We should actually find how Cp2 depends on CpTh
    
    fi    = fi_three_state_trimer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return (1-fu-fi)*linear_signal(T,bN,kN) + fi*bI*3 + fu*linear_signal(T,bU,kU)*3

def unfolding_curve_tetramer_monomeric_intermediate(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1=0,CpTh=0):

    """
    N4 ⇔ 4Ι ⇔ 4U Three-state unfolding with a monomeric intermediate
    C = concentration of the trimer equivalent
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,(CpTh-Cp1)/4)  
    
    fi    = fi_three_state_tetramer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return (1-fu-fi)*linear_signal(T,bN,kN) + fi*bI*4 + fu*linear_signal(T,bU,kU)*4

def unfolding_curve_trimer_trimeric_intermediate(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1=0,CpTh=0):

    """
    N3 ⇔ Ι3 ⇔ 3U Three-state unfolding with a trimeric intermediate
    C = concentration of the trimer equivalent
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,CpTh-Cp1)   
    
    fu = fu_three_state_trimer_trimeric_intermediate(K1,K2,C)
    fi = fi_three_state_trimer_trimeric_intermediate(fu,K2,C)

    return (1-fu-fi)*linear_signal(T,bN,kN) + fi*bI + fu*linear_signal(T,bU,kU)*3

def unfolding_curve_dimer_dimeric_intermediate(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1=0,CpTh=0):

    """
    N2 ⇔ Ι2 ⇔ 2U Three-state unfolding with a monomeric intermediate
    C       = molar concentration in dimer equivalent
    CpTotal = Cp1 + Cp2
    """

    K1  = eq_constant(T,DH1,T1,Cp1)
    K2  = eq_constant(T,DH2,T2,CpTh-Cp1)
    
    fu  = fu_three_state_dimer_dimeric_intermediate(K1,K2,C)
    fi  = fi_three_state_dimer_dimeric_intermediate(fu,K2,C)

    return (1-fu-fi)*linear_signal(T,bN,kN) + fi*bI + fu*linear_signal(T,bU,kU)*2

def map_three_state_model_to_signal_fx(model):

    signal_fx_map = {
    'Monomer':                         unfolding_curve_monomer_monomeric_intermediate,
    'Dimer_monomeric_intermediate':    unfolding_curve_dimer_monomeric_intermediate,
    'Dimer_dimeric_intermediate':      unfolding_curve_dimer_dimeric_intermediate,
    'Trimer_monomeric_intermediate':   unfolding_curve_trimer_monomeric_intermediate,
    'Trimer_trimeric_intermediate':    unfolding_curve_trimer_trimeric_intermediate,
    'Tetramer_monomeric_intermediate': unfolding_curve_tetramer_monomeric_intermediate
    }

    return signal_fx_map.get(model)

def fit_thermal_unfolding_three_species(listOfTemperatures,listOfSignals,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,signal_fx,listOfOligomerConc=None,fixed_t=False,t1=0,t2=0):

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

    allSignal       = concat_signal_lst(listOfSignals)

    if fixed_t:

        paramsInit  = np.delete(initialParameters, [0, 2])
        lowBoundsI  = np.delete(lowBounds,         [0, 2])
        highBoundsI = np.delete(highBounds,        [0, 2])
    
    else:
        
        paramsInit  = initialParameters
        lowBoundsI  = lowBounds
        highBoundsI = highBounds

    def thermal_unfolding(dummyVariable,*args):

        '''
        Calculate the thermal unfolding profile of many curves at the same time

        Requires:

            - The 'listOfTemperatures' containing each of them a single dataset

        The other arguments have to be in the following order:

            - Global melting temperature   for the first transition   
            - Global enthalpy of unfolding for the first transition
            - Global melting temperature   for the second transition   
            - Global enthalpy of unfolding for the second transition
            - Single intercepts, folded
            - Single intercepts, unfolded
            - Single intercepts, intermediate

            - Single slopes,     folded
            - Single slopes,     unfolded

        Returns:

            The melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
                slopes and intercept of the folded, intermediate and unfolded states

        '''

        if not fixed_t:
            start             = 4
            T1, DH1, T2, DH2  = args[:4]   
            
        else:
            start     = 2
            DH1, DH2  = args[:2]  
            T1,  T2   = t1, t2

        totalDataSets = len(listOfTemperatures)

        interceptsFolded       = args[(start+totalDataSets*0):(start+totalDataSets*1)]
        interceptsUnfolded     = args[(start+totalDataSets*1):(start+totalDataSets*2)]
        interceptsIntermediate = args[(start+totalDataSets*2):(start+totalDataSets*3)]

        if fitSlopeNative:

            slopesFolded   = args[(start+totalDataSets*3):(start+totalDataSets*4)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,T in enumerate(listOfTemperatures):

            bN, bU, bI = interceptsFolded[i], interceptsUnfolded[i], interceptsIntermediate[i]
            
            kN = 0 if not fitSlopeNative          else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded        else slopesUnfolded[i]
            C  = 0 if listOfOligomerConc is None  else listOfOligomerConc[i]

            Y = signal_fx(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,0,0)

            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(thermal_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds))

    return global_fit_params, cov

def fit_thermal_unfolding_three_species_with_cp(listOfTemperatures,listOfSignals,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,CpTh,signal_fx,listOfOligomerConc=None,fixed_t=False,t1=0,t2=0):

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
    allSignal       = concat_signal_lst(listOfSignals)

    if fixed_t:

        paramsInit  = np.delete(initialParameters, [0, 2])
        lowBoundsI  = np.delete(lowBounds,         [0, 2])
        highBoundsI = np.delete(highBounds,        [0, 2])
    
    else:
        paramsInit  = initialParameters
        lowBoundsI  = lowBounds
        highBoundsI = highBounds

    def thermal_unfolding(dummyVariable,*args):

        '''
        Calculate the thermal unfolding profile of many curves at the same time

        Requires:

            - The 'listOfTemperatures' containing each of them a single dataset

        The other arguments have to be in the following order:

            - Global melting temperature   for the first transition   
            - Global enthalpy of unfolding for the first transition
            - Global melting temperature   for the second transition   
            - Global enthalpy of unfolding for the second transition
            - Single intercepts, folded
            - Single intercepts, unfolded
            - Single intercepts, intermediate

            - Single slopes,     folded
            - Single slopes,     unfolded

        Returns:

            The melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
                slopes and intercept of the folded, intermediate and unfolded states

        '''

        totalDataSets = len(listOfTemperatures)

        if not fixed_t:
            start                  = 5
            T1, DH1, T2, DH2, Cp1  = args[:5]   
            
        else:
            start          = 3
            DH1, DH2, Cp1  = args[:3]  
            T1,  T2        = t1, t2

        interceptsFolded       = args[(start+totalDataSets*0):(start+totalDataSets*1)]
        interceptsUnfolded     = args[(start+totalDataSets*1):(start+totalDataSets*2)]
        interceptsIntermediate = args[(start+totalDataSets*2):(start+totalDataSets*3)]

        if fitSlopeNative:

            slopesFolded   = args[(start+totalDataSets*3):(start+totalDataSets*4)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,T in enumerate(listOfTemperatures):

            bN, bU, bI = interceptsFolded[i], interceptsUnfolded[i], interceptsIntermediate[i]
            
            kN = 0 if not fitSlopeNative   else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded else slopesUnfolded[i]
            C  = 0 if listOfOligomerConc is None  else listOfOligomerConc[i]

            Y = signal_fx(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1,CpTh)

            signal.append(Y)

        return np.concatenate(signal).flatten()

    global_fit_params, cov = curve_fit(thermal_unfolding,1,allSignal,
        p0=paramsInit,bounds=(lowBoundsI, highBoundsI))

    return global_fit_params, cov

def ode_irrev(T, xn,Tf, Ea,scan_rate_v):
    
    """
    ordinary differential equation for the native fraction native versus temperature
    dxn/dT = -1/v*k(T)*xn

    Tf is the temperature at which k = 1/min
    """

    return -1 / scan_rate_v * arrhenius(T, Tf, Ea) * xn

def get_ode_function_irrev(scan_rate_v):

    # Return ode function with a fixed scan_rate_v
    def return_func(T, xn,Tf, Ea):

        return ode_irrev(T, xn, Tf, Ea,scan_rate_v)

    return return_func 

def ode_rev_irrev(T, xf, dh,t1,Tf, Ea,scan_rate_v,Cp1=0):
    
    """
    ordinary differential equation for irreversibly denaturated fraction versus temperature
    dxf/dT = -1/v*k(T)*K(T)/(K(T) + 1)*(1 - xf)

    Here we call 'xf' to the final (irreversibly denaturated) state. Nomeclature from Jose M. Sanchez-Ruiz 1992

    Tf is the temperature at which k = 1/min

    """

    k = arrhenius(T, Tf, Ea)
    K = eq_constant(T,dh,t1,Cp1)
    
    return k*K*(1 - xf)/(scan_rate_v*(K + 1))

def get_ode_function_rev_irrev(scan_rate_v):

    # Return ode function with a fixed scan_rate_v
    def return_func(T, xf, dh,t1,Tf, Ea):

        return ode_rev_irrev(T, xf, dh,t1,Tf, Ea,scan_rate_v)

    return return_func

def irrev_signal(T, Tf,Ea,bF,bN,kN,kF,scan_rate_v):

    """
    Nomeclature from Jose M. Sanchez-Ruiz 1992
    T,t1, Tf should be in Kelvin!
    The variable
        'xf' refers to the final (irreversibly denaturated) state
        'xn' refers to the folded  state
        'xu' refers to the unfolded state
    """

    ivp_result = solve_ivp(ode_irrev,t_span=[np.min(T), np.max(T)],t_eval=T,y0=[1],args=(Tf,Ea,scan_rate_v),method="LSODA")

    xn = ivp_result.y[0, :]    

    return xn*linear_signal(T,bN,kN) + (1 - xn)*linear_signal(T,bF,kF)  

def rev_irrev_signal(T, dh,t1,Tf,Ea,bF,bN,bU,kN,kF,scan_rate_v,Cp1=0):

    """
    Nomeclature from Jose M. Sanchez-Ruiz 1992
    T,t1, Tf should be in Kelvin!
    The variable
        'xf' refers to the final (irreversibly denaturated) state
        'xn' refers to the folded  state
        'xu' refers to the unfolded state
    """

    ivp_result = solve_ivp(ode_rev_irrev,t_span=[np.min(T), np.max(T)],t_eval=T,y0=[0],args=(dh,t1,Tf,Ea,scan_rate_v),method="LSODA")

    xf = ivp_result.y[0, :]
    Ks = eq_constant(T,dh,t1,Cp1)
    
    xn = (1 - xf) / (1 + Ks)
    xu = 1 - xf - xn

    return xf*linear_signal(T,bF,kF) + xn*linear_signal(T,bN,kN) + xu*bU

def fit_irrev_unfolding(listOfTemperatures,listOfSignals,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,scan_rate_v,fixed_tf_ea=False,Tf_fix=0,Ea_fix=0):

    """
    k
    N  -> D.  One irreversible step

    where N, and D are the native, and final (irreversibly denatured) states of the protein

    """
    allSignal = concat_signal_lst(listOfSignals)
    ode       = get_ode_function_irrev(scan_rate_v)

    if fixed_tf_ea:

        paramsInit  = np.delete(initialParameters, [0, 1])
        lowBoundsI  = lowBounds[2:]
        highBoundsI = highBounds[2:]
    
    else:
        paramsInit  = initialParameters
        lowBoundsI  = lowBounds
        highBoundsI = highBounds

    def get_signal(T,Tf,Ea,bF,bN,kN,kF):

        """
        Nomeclature from Jose M. Sanchez-Ruiz 1992
        The variable
            'xf' refers to the final (irreversibly denaturated) state
            'xn' refers to the folded  state
            'xu' refers to the unfolded state
        """

        ivp_result = solve_ivp(ode,t_span=[np.min(T), np.max(T)],t_eval=T,y0=[1],args=(Tf,Ea),method="LSODA")

        xn = ivp_result.y[0, :]        
        xf = (1 - xn) 

        signal = xf*linear_signal(T,bF,kF) + xn*linear_signal(T,bN,kN)

        return signal

    def thermal_unfolding(dummyVariable,*args):

        '''
        Calculate the thermal unfolding profile of many curves at the same time

        Requires:

            - The 'listOfTemperatures' containing each of them a single dataset

        Returns:

            The melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
                slopes and intercept of the folded and unfolded states

        '''

        totalDataSets = len(listOfTemperatures)
  
        if not fixed_tf_ea:
            start         = 2
            Tf            = args[0]  
            Ea            = args[1] 
            
        else:
            start         = 0
            Tf            = Tf_fix 
            Ea            = Ea_fix
    
        interceptsFolded       = args[(start+totalDataSets*0):(start+totalDataSets*1)]
        interceptsUnfolded     = args[(start+totalDataSets*1):(start+totalDataSets*2)]

        if fitSlopeNative:

            slopesFolded   = args[(start+totalDataSets*2):(start+totalDataSets*3)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,T in enumerate(listOfTemperatures):

            bN, bU = interceptsFolded[i],  interceptsUnfolded[i]
            
            kN = 0 if not fitSlopeNative   else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded else slopesUnfolded[i]

            Y = get_signal(T, Tf, Ea,bU,bN,kN,kU) # Arguments are in the right order. Check the change in notation.

            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(thermal_unfolding,1,allSignal,
        p0=paramsInit,bounds=(lowBoundsI, highBoundsI))

    return global_fit_params, cov

def fit_1rev_2irrev_unfolding(listOfTemperatures,listOfSignals,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,scan_rate_v,fixed_tf_ea=False,Tf_fix=0,Ea_fix=0):

    """
        K     k
    N  <-> I -> D.  Reversible step followed by irreversible step

    where N, I, and D are the native, unfolded intermediate, and final (irreversibly denatured) states of the protein

    """

    allSignal = concat_signal_lst(listOfSignals)
    ode       = get_ode_function_rev_irrev(scan_rate_v)

    if fixed_tf_ea:

        paramsInit  = np.delete(initialParameters, [2, 3])
        lowBoundsI  = np.concatenate((lowBounds[:2],lowBounds[4:]))
        highBoundsI = np.concatenate((highBounds[:2],highBounds[4:]))
    
    else:
        paramsInit  = initialParameters
        lowBoundsI  = lowBounds
        highBoundsI = highBounds

    def get_signal(T, dh,t1,Tf,Ea,bF,bN,bU,kN,kF,Cp1=0):

        """
        Nomeclature from Jose M. Sanchez-Ruiz 1992
        The variable
            'xf' refers to the final (irreversibly denaturated) state
            'xn' refers to the folded  state
            'xu' refers to the unfolded state
        """

        ivp_result = solve_ivp(ode,t_span=[np.min(T), np.max(T)],t_eval=T,y0=[0],args=(dh,t1,Tf,Ea),method="LSODA")

        xf = ivp_result.y[0, :]
        Ks = eq_constant(T,dh,t1,Cp1)
        
        xn = (1 - xf) / (1 + Ks)
        xu = 1 - xf - xn

        signal = xf*linear_signal(T,bF,kF) + xn*linear_signal(T,bN,kN) + xu*bU

        return signal

    def thermal_unfolding(dummyVariable,*args):

        '''
        Calculate the thermal unfolding profile of many curves at the same time

        Requires:

            - The 'listOfTemperatures' containing each of them a single dataset

        Returns:

            The melting curves based on the parameters Temperature of melting, enthalpy of unfolding,
                slopes and intercept of the folded and unfolded states

        '''

        totalDataSets = len(listOfTemperatures)
        T1            = args[0] # 
        DH1           = args[1] #
  
        if not fixed_tf_ea:
            start         = 4
            Tf            = args[2]  
            Ea            = args[3] 
            
        else:
            start         = 2
            Tf            = Tf_fix 
            Ea            = Ea_fix
    
        interceptsFolded       = args[(start+totalDataSets*0):(start+totalDataSets*1)]
        interceptsUnfolded     = args[(start+totalDataSets*1):(start+totalDataSets*2)]
        interceptsIntermediate = args[(start+totalDataSets*2):(start+totalDataSets*3)]

        if fitSlopeNative:

            slopesFolded   = args[(start+totalDataSets*3):(start+totalDataSets*4)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,T in enumerate(listOfTemperatures):

            bN, bU, bI = interceptsFolded[i],   interceptsUnfolded[i], interceptsIntermediate[i]
            
            kN = 0 if not fitSlopeNative   else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded else slopesUnfolded[i]

            Y = get_signal(T, DH1,T1,Tf, Ea,bU,bN,bI,kN,kU) # Arguments are in the right order. Check the change in notation.

            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(thermal_unfolding,1,allSignal,
        p0=paramsInit,bounds=(lowBoundsI, highBoundsI))

    return global_fit_params, cov

def two_state_rev_unfolding_fractions(T,DHm,Tm,extra_arg,Cp=0):

    K  = eq_constant(T,DHm,Tm,Cp) 
    fn = fn_two_state_monomer(K)

    return {'Native': fn, 'Unfolded': (1-fn)}

def two_state_dimer_unfolding_fractions(T,DH1,T1,C,Cp=0):

    """
    N2 ⇔ 2U   where C is the total concentration (M) of the protein in dimer equivalent.
    """

    K  = eq_constant(T,DH1,T1,Cp)
    fu = fu_two_state_dimer(K,C)

    return {'Native dimer':(1-fu), 'Unfolded monomer':fu}

def two_state_trimer_unfolding_fractions(T,DH1,T1,C,Cp=0):

    """
    N3 ⇔ 3U   
    C is the total concentration (M) of the protein in trimer equivalent.
    """

    K  = eq_constant(T,DH1,T1,Cp)
    fu = fu_two_state_trimer(K,C)

    return {'Native trimer':(1-fu), 'Unfolded monomer':fu}

def two_state_tetramer_unfolding_fractions(T,DH1,T1,C,Cp=0):

    """
    N4 ⇔ 4U   C is the total concentration (M) of the protein in tetramer equivalent.
    """

    K  = eq_constant(T,DH1,T1,Cp)
    fu = fu_two_state_tetramer(K,C)

    return {'Native tetramer':(1-fu), 'Unfolded monomer':fu}

def three_state_dimer_monomeric_intermediate_unfolding_fractions(T,DH1,T1,DH2,T2,C,Cp1=0,CpTh=0):

    """
    N2 ⇔ 2Ι ⇔ 2U Three-state unfolding with a monomeric intermediate
    C = concentration of the dimer equivalent

    Note that the N2 state is considered to have twice the number of residues of the protein monomer
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,(CpTh-Cp1)/2)
    
    fi = fi_three_state_dimer_monomeric_intermediate(K1,K2,C)
    fu = fi*K2
    fn = 1-fu-fi

    return {'Native dimer':fn, 'Monomeric intermediate':fi, 'Unfolded monomer':fu}

def three_state_trimer_monomeric_intermediate_unfolding_fractions(T,DH1,T1,DH2,T2,C,Cp1=0,Cp2=0):

    """
    N3 ⇔ 3Ι ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration of the trimer equivalent
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,Cp2)
    
    fi    = fi_three_state_trimer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return {'Native trimer':(1-fi-fu), 'Monomeric intermediate':fi, 'Unfolded monomer':fu}

def three_state_tetramer_monomeric_intermediate_unfolding_fractions(T,DH1,T1,DH2,T2,C,Cp1=0,CpTh=0):

    """
    N4 ⇔ 4Ι ⇔ 4U Three-state unfolding with a monomeric intermediate
    C = concentration of the trimer equivalent
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,(CpTh-Cp1)/4)
    
    fi    = fi_three_state_tetramer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return {'Native tetramer':(1-fi-fu), 'Monomeric intermediate':fi, 'Unfolded monomer':fu}

def three_state_trimer_trimeric_intermediate_unfolding_fractions(T,DH1,T1,DH2,T2,C,Cp1=0,CpTh=0):

    """
    N3 ⇔ Ι3 ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration of the trimer equivalent
    """

    K1 = eq_constant(T,DH1,T1,Cp1)
    K2 = eq_constant(T,DH2,T2,CpTh-Cp1)
    
    fu = fu_three_state_trimer_trimeric_intermediate(K1,K2,C)
    fi = fi_three_state_trimer_trimeric_intermediate(fu,K2,C)

    return {'Native trimer':(1-fi-fu), 'Trimeric intermediate':fi, 'Unfolded monomer':fu}

def three_state_dimer_dimeric_intermediate_unfolding_fractions(T,DH1,T1,DH2,T2,C,Cp1=0,CpTh=0):

    """
    N2 ⇔ Ι2 ⇔ 2U Three-state unfolding with a monomeric intermediate
    C       = molar concentration in dimer equivalent
    """

    K1  = eq_constant(T,DH1,T1,Cp1)
    K2  = eq_constant(T,DH2,T2,CpTh-Cp1)
    
    fu  = fu_three_state_dimer_dimeric_intermediate(K1,K2,C)
    fi  = fi_three_state_dimer_dimeric_intermediate(fu,K2,C)

    return {'Native dimer':(1-fu-fi), 'Dimeric intermediate':fi, 'Unfolded monomer':fu}

def three_state_monomer_monomeric_intermediate_unfolding_fractions(T,DH1,DH2,T1,T2,extra_arg=None,Cp1=0,CpTh=0):

    '''
    Compute the folded, intermediate and unfolded fraction for a system with
    two reversible equilibria. F <-> I <-> U 
    '''

    A = eq_constant(T,DH1,T1,Cp1)
    B = eq_constant(T,DH2,T2,CpTh-Cp1)

    den = (1+A+A*B)

    return {'Native':1/den, 'Intermediate':A/den, 'Unfolded':A*B/den}

def rev_irrev_fractions(T, dh,t1,Tf,Ea,scan_rate_v,Cp=0):

    """
    Nomeclature from Jose M. Sanchez-Ruiz 1992
    The variable
        'xf' refers to the final (irreversibly denaturated) state
        'xn' refers to the folded  state
        'xu' refers to the unfolded state
    """

    ivp_result = solve_ivp(ode_rev_irrev,t_span=[np.min(T), np.max(T)],t_eval=T,y0=[0],args=(dh,t1,Tf,Ea,scan_rate_v),method="LSODA")

    xf = ivp_result.y[0, :]
    Ks = eq_constant(T,dh,t1,Cp)
    
    xn = (1 - xf) / (1 + Ks)
    xu = 1 - xf - xn

    return {'Native':xn, 'Unfolded (reversible)':xu, 'Unfolded (irreversible)':xf}

def irrev_fractions(T, Tf,Ea,scan_rate_v):

    # N -> U

    ivp_result = solve_ivp(ode_irrev,t_span=[np.min(T), np.max(T)],t_eval=T,y0=[1],args=(Tf,Ea,scan_rate_v),method="LSODA")

    xn = ivp_result.y[0, :]    

    return {'Native':xn, 'Unfolded (irreversible)':(1 - xn)}

def predict_all_signal_irrev_three_state(
    T, T1,DH1,Tf,Ea,interceptsFolded,slopesFolded,
    interceptsUnfolded,slopesUnfolded,interceptsIntermediate,scan_rate):

    predicted = []

    for bN, kN, bU, kU, bI in zip(interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded,interceptsIntermediate):

        Y = rev_irrev_signal(T, DH1,T1,Tf,Ea,bU,bN,bI,kN,kU,scan_rate)
        predicted.append(Y)

    return np.array(predicted)

def predict_all_signal_irrev_two_state(
    T,Tf,Ea,interceptsFolded,slopesFolded,
    interceptsUnfolded,slopesUnfolded,scan_rate):

    predicted = []

    for bN, kN, bU, kU in zip(interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded):

        Y = irrev_signal(T,Tf,Ea,bU,bN,kN,kU,scan_rate)
        predicted.append(Y)

    return np.array(predicted)

def predict_all_signal_nmer_with_intermediate(
    T, T1,DH1,T2,DH2,interceptsFolded,slopesFolded,
    interceptsUnfolded,slopesUnfolded,interceptsIntermediate,oligoConcLst,signal_fx,Cp1=0,CpTh=0):

    predicted = []

    for C,bN, kN, bU, kU, bI in zip(oligoConcLst,interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded,interceptsIntermediate):

        Y = signal_fx(T,T1,DH1,T2,DH2,bN,kN,bU,kU,bI,C,Cp1,CpTh)

        predicted.append(Y)

    return np.array(predicted)

def map_two_state_model_to_fractions_fx(model):

    fractions_fx_map = {
    'Monomer':  two_state_rev_unfolding_fractions,
    'Dimer':    two_state_dimer_unfolding_fractions,
    'Trimer':   two_state_trimer_unfolding_fractions,
    'Tetramer': two_state_tetramer_unfolding_fractions
    }

    return fractions_fx_map.get(model)

def map_three_state_model_to_fractions_fx(model):

    fractions_fx_map = {
    'Monomer':                         three_state_monomer_monomeric_intermediate_unfolding_fractions,
    'Dimer_monomeric_intermediate':    three_state_dimer_monomeric_intermediate_unfolding_fractions,
    'Dimer_dimeric_intermediate':      three_state_dimer_dimeric_intermediate_unfolding_fractions,
    'Trimer_monomeric_intermediate':   three_state_trimer_monomeric_intermediate_unfolding_fractions,
    'Trimer_trimeric_intermediate':    three_state_trimer_trimeric_intermediate_unfolding_fractions,
    'Tetramer_monomeric_intermediate': three_state_tetramer_monomeric_intermediate_unfolding_fractions
    }

    return fractions_fx_map.get(model)
