
from scipy.optimize  import curve_fit

from constants       import *
from fitting_helpers import *

'''
All thermodynamic parameters are used in kcal mol units

Unfolding functions for monomers have an argument called 'extra_arg' that is not used. 
This is because unfolding functions for oligomers require the protein concentration in that position
'''

def eq_constant_chem(D,D50,M,T):

    '''
    D50 is the denaturant concentration at which ΔG(D) = 0
    M represents the slope for the linear dependence of ΔG(D) on D
    T is the temperature in Kelvin units
    '''

    T  = temperature_to_kelvin(T)
    dG = M * (D50 - D)
    K  = np.exp(-dG / (R_gas * T))

    return K

def two_state_chemical_unfolding_curve_monomer(D,T,D50,M,bN,kN,bU,kU,extra_arg=None):
    #N2 ⇔ 2U

    K  = eq_constant_chem(D,D50,M,T)
    fn = fn_two_state_monomer(K)

    return fn*(D*kN+bN) + (1-fn)*(D*kU+bU)

def two_state_chemical_unfolding_curve_dimer(D,T,D50,M,bN,kN,bU,kU,C):
    #N2 ⇔ 2U

    K  = eq_constant_chem(D,D50,M,T)
    fu = fu_two_state_dimer(K,C)

    return (1-fu)*(D*kN+bN) + fu*(D*kU+bU)*2

def two_state_chemical_unfolding_curve_trimer(D,T,D50,M,bN,kN,bU,kU,C):

    """
    N3 ⇔ 3U   C is the total concentration (M) of the protein in trimer equivalent.
    """

    K  = eq_constant_chem(D,D50,M,T)
    fu = fu_two_state_trimer(K,C)

    return (1-fu)*(D*kN+bN) + fu*(D*kU+bU)*3

def two_state_chemical_unfolding_curve_tetramer(D,T,D50,M,bN,kN,bU,kU,C):

    """
    N4 ⇔ 4U   C is the total concentration (M) of the protein in tetramer equivalent.
    """

    K  = eq_constant_chem(D,D50,M,T)
    fu = fu_two_state_tetramer(K,C)

    return (1-fu)*(D*kN+bN) + fu*(D*kU+bU)*4

def fit_chemical_unfolding(listOfChemConcentration,listOfSignals,temperature,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,signal_fx,listOfOligomerConc=None):

    '''
    Fit the chemical unfolding profile of many curves at the same time
    Useful function to do global fitting of local and global parameters

    Requires:

        - The 'listOfChemConcentration' (one per dataset)
        - The 'listOfSignals'           (one per dataset)

    Returns:

        The fitting of the chemical unfolding curves 
    '''

    allSignal       = concat_signal_lst(listOfSignals)
        
    temperature = temperature_to_kelvin(temperature)

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

            kN = 0 if not fitSlopeNative          else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded        else slopesUnfolded[i]
            C  = 0 if listOfOligomerConc is None  else listOfOligomerConc[i]

            Y = signal_fx(D,temperature,D50,M,bN,kN,bU,kU,C)

            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(chemical_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds))

    return global_fit_params, cov

def chemical_unfolding_curve_monomer_monomeric_intermediate(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,extra_arg=None):

    """
    Three states reversible unfolding: N <-> I <-> U
    """

    A = eq_constant_chem(D,D50v1,M1,T)
    B = eq_constant_chem(D,D50v2,M2,T)

    den = (1+A+A*B)

    fN, fI, fU = 1/den, A/den, A*B/den
    
    return fN*(kN*D + bN) + fI*bI + fU*(kU*D + bU) 

def chemical_unfolding_curve_dimer_monomeric_intermediate(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C):

    """
    N2 ⇔ 2Ι ⇔ 2U Three-state unfolding with a monomeric intermediate
    C = concentration in dimer equivalent
    mTotal = m1 + m2*2
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fi = fi_three_state_dimer_monomeric_intermediate(K1,K2,C)
    fu = fi*K2

    return (1-fu-fi)*(kN*D + bN) + fi*bI*2 + fu*(kU*D + bU)*2

def chemical_unfolding_curve_trimer_monomeric_intermediate(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C):

    """
    N3 ⇔ 3Ι ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration in dimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fi    = fi_three_state_trimer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return (1-fi-fu)*(kN*D + bN) + fi*bI*3 + fu*(kU*D + bU)*3

def chemical_unfolding_curve_tetramer_monomeric_intermediate(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C):

    """
    N3 ⇔ 3Ι ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration in dimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fi    = fi_three_state_tetramer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return (1-fi-fu)*(kN*D + bN) + fi*bI*4 + fu*(kU*D + bU)*4

def chemical_unfolding_curve_trimer_trimeric_intermediate(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C):

    """
    N3 ⇔ Ι3 ⇔ 3U Three-state unfolding with a trimeric intermediate
    C = concentration in dimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fu = fu_three_state_trimer_trimeric_intermediate(K1,K2,C)
    fi = fi_three_state_trimer_trimeric_intermediate(fu,K2,C)

    return (1-fi-fu)*(kN*D + bN) + fi*bI + fu*(kU*D + bU)*3

def chemical_unfolding_curve_dimer_dimeric_intermediate(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C):

    """
    N2 ⇔ Ι2 ⇔ 2U Three-state unfolding with a monomeric intermediate
    C = concentration in dimer equivalent
    mTotal = m1 + m2
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fu  = fu_three_state_dimer_dimeric_intermediate(K1,K2,C)
    fi  = fi_three_state_dimer_dimeric_intermediate(fu,K2,C)

    return (1-fu-fi)*(kN*D + bN) + fi*bI + fu*(kU*D + bU)*2

def fit_chemical_unfolding_three_species(listOfChemConcentration,listOfSignals,temperature,initialParameters,
    lowBounds,highBounds,fitSlopeNative,fitSlopeUnfolded,signal_fx,listOfOligomerConc=None,
    fixed_d=False,d1=None,d2=None):

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

    allSignal       = concat_signal_lst(listOfSignals)
    
    T = temperature_to_kelvin(temperature)

    if fixed_d:

        paramsInit  = np.delete(initialParameters, [1, 3])
        lowBoundsI  = np.delete(lowBounds, [1, 3])
        highBoundsI = np.delete(highBounds, [1, 3])
    
    else:
        
        paramsInit  = initialParameters
        lowBoundsI  = lowBounds
        highBoundsI = highBounds

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

        if not fixed_d:
            start                 = 4
            M1, D50v1, M2, D50v2  = args[:4]   
            
        else:
            start           = 2
            M1, M2          = args[:2]  
            D50v1,  D50v2   = d1, d2

        interceptsFolded       = args[(start+totalDataSets*0):(start+totalDataSets*1)]
        interceptsUnfolded     = args[(start+totalDataSets*1):(start+totalDataSets*2)]
        interceptsIntermediate = args[(start+totalDataSets*2):(start+totalDataSets*3)]

        if fitSlopeNative:

            slopesFolded   = args[(start+totalDataSets*3):(start+totalDataSets*4)]

        if fitSlopeUnfolded:

            slopesUnfolded = args[(len(args)-totalDataSets):]

        signal = []

        for i,D in enumerate(listOfChemConcentration):

            bN, bU, bI = interceptsFolded[i], interceptsUnfolded[i], interceptsIntermediate[i]

            kN = 0 if not fitSlopeNative   else slopesFolded[i]
            kU = 0 if not fitSlopeUnfolded else slopesUnfolded[i]
            C  = 0 if listOfOligomerConc is None  else listOfOligomerConc[i]

            Y = signal_fx(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C)

            signal.append(Y)

        return np.concatenate(signal,axis=0)

    global_fit_params, cov = curve_fit(chemical_unfolding,1,allSignal,
        p0=initialParameters,bounds=(lowBounds, highBounds),method='trf')

    return global_fit_params, cov

def chem_monomer_two_state_rev_unfolding_fractions(T,D,M,D50,extra_arg=None):

    '''
    Compute the folded and unfolded fraction for a system with one reversible equilibrium (LEM model)

    F <-> U, K1 = [U] / [F] 
    '''

    K  = eq_constant_chem(D,D50,M,T)  
    fn =  fn_two_state_monomer(K)

    return {'Native':fn, 'Unfolded':(1-fn)}

def chem_dimer_two_state_rev_unfolding_fractions(T,D,M,D50,C):

    '''
    Compute the folded and unfolded fraction for a system with
    one reversible equilibrium (LEM model)
    F2 <-> 2U
    '''

    K  = eq_constant_chem(D,D50,M,T)
    fu = fu_two_state_dimer(K,C)

    return {'Native dimer':(1-fu), 'Unfolded':fu}

def chem_trimer_two_state_rev_unfolding_fractions(T,D,M,D50,C):

    # F3 <-> 3U

    K  = eq_constant_chem(D,D50,M,T)
    fu = fu_two_state_trimer(K,C)

    return {'Native trimer':(1-fu), 'Unfolded monomer':fu}

def chem_tetramer_two_state_rev_unfolding_fractions(T,D,M,D50,C):

    # F4 <-> 4U

    K  = eq_constant_chem(D,D50,M,T)
    fu = fu_two_state_tetramer(K,C)

    return {'Native trimer':(1-fu), 'Unfolded monomer':fu}

def chemical_unfolding_monomer_monomeric_intermediate_fractions(D,T,D50v1,M1,D50v2,M2,extra_arg=None):

    '''
    Compute the folded, intermediate and unfolded fraction for a system with
    two reversible equilibria, chemical unfolding (LEM model)

    F <-> I <-> U ; K1 = [I] / [F] ; K2 = [U] / [I] 

    Input:

        Parameter name  Detail                                      Units
        - 'T'           temperature                                 in celsius or kelvin
        - 'D'           concentration of the denaturant agent       in molar    
        - 'M1'          dependence of the DG1 on D                  in kcal/mol/M
        - 'D50v1'       concentration of D where DG1 equals zero    in molar
        - 'M2'          dependence of the DG2 on D                  in kcal/mol/M
        - 'D50v2'       concentration of D where DG2 equals zero    in molar
    '''

    A = eq_constant_chem(D,D50v1,M1,T)
    B = eq_constant_chem(D,D50v2,M2,T)

    den = (1+A+A*B)

    xN, xI, xU = 1/den, A/den, A*B/den

    return {'Native':xN, 'Intermediate':xI, 'Unfolded':xU}

def chemical_unfolding_dimer_monomeric_intermediate_fractions(D,T,D50v1,M1,D50v2,M2,C):

    """
    N2 ⇔ 2Ι ⇔ 2U Three-state unfolding with a monomeric intermediate
    C    = concentration in dimer equivalent
    mTot = m1 + 2*m2
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fi = fi_three_state_dimer_monomeric_intermediate(K1,K2,C)
    fu = fi*K2
    fn = 1-fu-fi

    return {'Native dimer':fn, 'Monomeric intermediate':fi, 'Unfolded':fu}

def chemical_unfolding_trimer_monomeric_intermediate_fractions(D,T,D50v1,M1,D50v2,M2,C):

    """
    N3 ⇔ 3Ι ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration in trimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fi    = fi_three_state_trimer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return {'Native trimer':(1-fi-fu), 'Monomeric intermediate':fi, 'Unfolded monomer':fu}

def chemical_unfolding_tetramer_monomeric_intermediate_fractions(D,T,D50v1,M1,D50v2,M2,C):

    """
    N3 ⇔ 3Ι ⇔ 3U Three-state unfolding with a monomeric intermediate
    C = concentration in trimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fi    = fi_three_state_tetramer_monomeric_intermediate(K1,K2,C)
    fu    = fi*K2

    return {'Native tetramer':(1-fi-fu), 'Monomeric intermediate':fi, 'Unfolded monomer':fu}

def chemical_unfolding_trimer_trimeric_intermediate_fractions(D,T,D50v1,M1,D50v2,M2,C):

    """
    N3 ⇔ Ι3 ⇔ 3U Three-state unfolding with a trimeric intermediate
    C = concentration in trimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fu = fu_three_state_trimer_trimeric_intermediate(K1,K2,C)
    fi = fi_three_state_trimer_trimeric_intermediate(fu,K2,C)

    return {'Native trimer':(1-fi-fu), 'Trimeric intermediate':fi, 'Unfolded monomer':fu}

def chemical_unfolding_dimer_dimeric_intermediate_fractions(D,T,D50v1,M1,D50v2,M2,C):

    """
    N2 ⇔ 2Ι ⇔ 2U Three-state unfolding with a monomeric intermediate
    C = concentration in dimer equivalent
    """

    K1 = eq_constant_chem(D,D50v1,M1,T)
    K2 = eq_constant_chem(D,D50v2,M2,T)
    
    fu  = fu_three_state_dimer_dimeric_intermediate(K1,K2,C)
    fi  = fi_three_state_dimer_dimeric_intermediate(fu,K2,C) 

    return {'Native dimer':(1-fu-fi), 'Dimeric intermediate':fi, 'Unfolded':fu}

def predict_all_signal_nmer_with_intermediate_chemical(
    D,T,D50v1,M1,D50v2,M2,interceptsFolded,slopesFolded,
    interceptsUnfolded,slopesUnfolded,interceptsIntermediate,oligoConcLst,signal_fx):

    predicted = []

    for C,bN, kN, bU, kU, bI in zip(oligoConcLst,interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded,interceptsIntermediate):

        Y = signal_fx(D,T,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI,C)

        predicted.append(Y)

    return np.array(predicted)

def map_two_state_chem_model_to_signal_fx(model):

    signal_fx_map = {
    'Monomer':  two_state_chemical_unfolding_curve_monomer,
    'Dimer':    two_state_chemical_unfolding_curve_dimer,
    'Trimer':   two_state_chemical_unfolding_curve_trimer,
    'Tetramer': two_state_chemical_unfolding_curve_tetramer
    }

    return signal_fx_map.get(model)

def map_three_state_chem_model_to_signal_fx(model):

    signal_fx_map = {
    'Monomer':                         chemical_unfolding_curve_monomer_monomeric_intermediate,
    'Dimer_monomeric_intermediate':    chemical_unfolding_curve_dimer_monomeric_intermediate,
    'Dimer_dimeric_intermediate':      chemical_unfolding_curve_dimer_dimeric_intermediate,
    'Trimer_monomeric_intermediate':   chemical_unfolding_curve_trimer_monomeric_intermediate,
    'Trimer_trimeric_intermediate':    chemical_unfolding_curve_trimer_trimeric_intermediate,
    'Tetramer_monomeric_intermediate': chemical_unfolding_curve_tetramer_monomeric_intermediate
    }

    return signal_fx_map.get(model)

def map_two_state_chem_model_to_fractions_fx(model):

    fractions_fx_map = {
    'Monomer':  chem_monomer_two_state_rev_unfolding_fractions,
    'Dimer':    chem_dimer_two_state_rev_unfolding_fractions,
    'Trimer':   chem_trimer_two_state_rev_unfolding_fractions,
    'Tetramer': chem_tetramer_two_state_rev_unfolding_fractions
    }

    return fractions_fx_map.get(model)

def map_three_state_chem_model_to_fractions_fx(model):

    fractions_fx_map = {
    'Monomer':                         chemical_unfolding_monomer_monomeric_intermediate_fractions,
    'Dimer_monomeric_intermediate':    chemical_unfolding_dimer_monomeric_intermediate_fractions,
    'Dimer_dimeric_intermediate':      chemical_unfolding_dimer_dimeric_intermediate_fractions,
    'Trimer_monomeric_intermediate':   chemical_unfolding_trimer_monomeric_intermediate_fractions,
    'Trimer_trimeric_intermediate':    chemical_unfolding_trimer_trimeric_intermediate_fractions,
    'Tetramer_monomeric_intermediate': chemical_unfolding_tetramer_monomeric_intermediate_fractions
    }

    return fractions_fx_map.get(model)