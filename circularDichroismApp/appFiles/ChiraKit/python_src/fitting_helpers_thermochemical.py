import numpy  as np 

from scipy.optimize  import curve_fit

from constants       import *
from helpers         import *
from fitting_helpers import *

'''
The code here is partially based on  
Louise Hamborg et al., 2020. Global analysis of protein stability by temperature and chemical
denaturation

https://github.com/KULL-Centre/ProteinUnfolding2D (MIT License, accesed on September 2024)

'''
def eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1):

    '''
    Input:
            - T:    temperature 					                                  (1D numpy array)
            - D:    denaturant agent concentration                                   (1D numpy array)
            - Tm:   temperature at which the equilibrium constant equals one         (float)
            - DHm:  variation of enthalpy between the two considered states at Tm    (float)
            - Cp0:  variation of calorific capacity between the two states           (float)
            - m0:   m-value at the reference temperature (Tref)                      (float)
            - m1:   variation of calorific capacity between the two states           (float)
    Output:

            - K: the equilibrium constant at a certain temperature and denaturant agent concentration  
    '''
    T   = temperature_to_kelvin(T)

    Tm  = temperature_to_kelvin(Tm)

    DT  = T - Tref_cst

    DG   = DHm*(1 - T/Tm) + Cp0*(T - Tm - T*np.log(T/Tm)) - D*(m0 + m1*DT)

    DG_RT = -DG / (R_gas * T)

    K     = np.exp(DG_RT)

    return K

def signal_two_state_tc_unfolding_monomer(T,D,DHm,Tm,Cp0,m0,m1,bN,bU,a0,c0,a1,c1,extra_arg=None):

    '''
    T,D,DHm,Tm,Cp0,m0,m1 are defined in the function 'eq_constant_termochem'
    bN and bU are the baseline of the native and folded states (at the temperature Tref_cst defined in constants.py)
    a0 and a1 are the slopes for the temperature 
    c0 and c1 are the slopes for the denaturant agent concentration
    extra_arg is not used but required
    '''

    K   = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fn  = fn_two_state_monomer(K)

    return  fn*(linear_signal(T,bN,a0) + c0*D) + (1- fn)*(linear_signal(T,bU,a1) + c1*D)

def signal_two_state_tc_unfolding_dimer(T,D,DHm,Tm,Cp0,m0,m1,bN,bU,a0,c0,a1,c1,C):

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fu = fu_two_state_dimer(K,C)

    return  (1-fu)*(linear_signal(T,bN,a0) + c0*D) + 2*fu*(linear_signal(T,bU,a1) + c1*D)

def signal_two_state_tc_unfolding_trimer(T,D,DHm,Tm,Cp0,m0,m1,bN,bU,a0,c0,a1,c1,C):

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fu = fu_two_state_trimer(K,C)

    return  (1-fu)*(linear_signal(T,bN,a0) + c0*D) + 3*fu*(linear_signal(T,bU,a1) + c1*D)

def signal_two_state_tc_unfolding_tetramer(T,D,DHm,Tm,Cp0,m0,m1,bN,bU,a0,c0,a1,c1,C):

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fu = fu_two_state_tetramer(K,C)

    return  (1-fu)*(linear_signal(T,bN,a0) + c0*D) + 4*fu*(linear_signal(T,bU,a1) + c1*D)

def fit_tc_unfolding(listOfChemConcentration,listOfSignals,listOfTemperatures,initialParameters,
    lowBounds,highBounds,fitSlopeNativeT,fitSlopeUnfoldedT,
    fitSlopeNativeD,fitSlopeUnfoldedD,signal_fx,listOfOligomerConc=None):

    '''
    Fit the thermal and chemical unfolding profile of many curves at the same time
    Useful function to do global fitting of local and global parameters
    '''

    allSignal       = concat_signal_lst(listOfSignals)
        
    def tc_unfolding(dummyVariable,*args):

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
 
        totalDataSets     = len(listOfChemConcentration)
        DHm,Tm,Cp0,m0,m1  = args[:5] 
        start             = 5

        interceptsFolded    = args[start:(start+totalDataSets)]
        interceptsUnfolded  = args[(start+totalDataSets):(start+totalDataSets*2)]

        if fitSlopeNativeT:
            slopesFoldedT = args[start + totalDataSets * 2 : start + totalDataSets * 3]

        if fitSlopeUnfoldedT:
            offset = 3 if fitSlopeNativeT else 2
            slopesUnfoldedT = args[start + totalDataSets * offset : start + totalDataSets * (offset + 1)]

        if fitSlopeNativeD:
            offset = 4 if fitSlopeNativeT and fitSlopeUnfoldedT else 3 if fitSlopeUnfoldedT else 2
            slopesFoldedD = args[start + totalDataSets * offset : start + totalDataSets * (offset + 1)]

        if fitSlopeUnfoldedD:

            slopesUnfoldedD = args[(len(args)-totalDataSets):]

        signal = []

        for i,D in enumerate(listOfChemConcentration):

            bN, bU = interceptsFolded[i],   interceptsUnfolded[i]

            kNt = 0 if not fitSlopeNativeT         else slopesFoldedT[i]
            kUt = 0 if not fitSlopeUnfoldedT       else slopesUnfoldedT[i]
            kNd = 0 if not fitSlopeNativeD         else slopesFoldedD[i]
            kUd = 0 if not fitSlopeUnfoldedD       else slopesUnfoldedD[i]
            C   = 0 if listOfOligomerConc is None  else listOfOligomerConc[i]

            T   = listOfTemperatures[i]

            Y   = signal_fx(T,D,DHm,Tm,Cp0,m0,m1,bN,bU,kNt,kNd,kUt,kUd,C)
                    
            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(tc_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds))

    return global_fit_params, cov

def tc_monomer_two_state_rev_unfolding_fractions(T,D,DHm,Tm,Cp0,m0,m1,extra_arg=None):

    '''
    Compute the folded and unfolded fraction for a system with one reversible equilibrium (LEM model)

    F <-> U, K1 = [U] / [F] 
    '''

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fn =  fn_two_state_monomer(K)

    return {'Native':fn, 'Unfolded':(1-fn)}

def tc_dimer_two_state_rev_unfolding_fractions(T,D,DHm,Tm,Cp0,m0,m1,C):

    '''
    Compute the folded and unfolded fraction for a system with one reversible equilibrium (LEM model)
    F2 <-> 2U
    '''

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fu = fu_two_state_dimer(K,C)

    return {'Native dimer':(1-fu), 'Unfolded':fu}

def tc_trimer_two_state_rev_unfolding_fractions(T,D,DHm,Tm,Cp0,m0,m1,C):

    # F3 <-> 3U

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fu = fu_two_state_trimer(K,C)

    return {'Native trimer':(1-fu), 'Unfolded monomer':fu}

def tc_tetramer_two_state_rev_unfolding_fractions(T,D,DHm,Tm,Cp0,m0,m1,C):

    # F4 <-> 4U

    K  = eq_constant_termochem(T,D,DHm,Tm,Cp0,m0,m1)
    fn = fn_two_state_tetramer(K,C)

    return {'Native trimer':fn, 'Unfolded monomer':(1-fn)}

def map_two_state_tc_model_to_signal_fx(model):

    signal_fx_map = {
    'Monomer':  signal_two_state_tc_unfolding_monomer,
    'Dimer':    signal_two_state_tc_unfolding_dimer,
    'Trimer':   signal_two_state_tc_unfolding_trimer,
    'Tetramer': signal_two_state_tc_unfolding_tetramer
    }

    return signal_fx_map.get(model)

def map_two_state_tc_model_to_fractions_fx(model):

    fractions_fx_map = {
    'Monomer':  tc_monomer_two_state_rev_unfolding_fractions,
    'Dimer':    tc_dimer_two_state_rev_unfolding_fractions,
    'Trimer':   tc_trimer_two_state_rev_unfolding_fractions,
    'Tetramer': tc_tetramer_two_state_rev_unfolding_fractions
    }

    return fractions_fx_map.get(model)