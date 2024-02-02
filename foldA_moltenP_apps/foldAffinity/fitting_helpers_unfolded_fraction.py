import numpy  as np
import pandas as pd

from scipy.optimize     import curve_fit
from scipy.optimize     import minimize_scalar
from scipy              import stats


def calculate_L_free(ligand_conc, Ku, Kd, pconc):
    '''
    Calculate the equilibrium concentration of free ligand
    in the equation: U + F <=>[Ku] L + F <=>[Kd] LF
    Args: 
       ligand_conc: initial ligand concentrations
       pconc: total protein concentration
       Ku: eqilibrum constant of unfolding
       Kd: dissociation constant in M
    '''
    # First calculate complex concentration by solving quadratic expression
    p = pconc/(Ku+1) + ligand_conc + Kd
    q = pconc*ligand_conc/(Ku+1)
    # if pconc > Kd:
    #     complex_conc = .5*p + np.sqrt(0.25*p**2 - q)
    # else:
    complex_conc = .5*p - np.sqrt(0.25*p**2 - q)
    # Calculate Lfree
    L_free = ligand_conc - complex_conc
    return L_free

def calculate_F_free(protein_conc, Ku, Kdim):
    '''
    Calculate the equilibrium concentration of free folded protein
    in the equation: U + F <=> [Ku] F + F <=>[Kdim] FF
    Args: 
       protein_conc: total protein concentration vector
       Ku: eqilibrum constant of unfolding
       Kdim: dimerization constant in M
    '''
    # First calculate complex concentration by solving quadratic expression
    a =  Kdim
    b =  1 + Ku
    c =  -protein_conc

    root1          = (-b-np.sqrt(b**2-4*a*c)) / (2*a)
    root2          = (-b+np.sqrt(b**2-4*a*c)) / (2*a)

    sel1  = np.greater(root1,np.zeros(len(root1)))
    sel2  = np.greater(root2,np.zeros(len(root2)))

    F              = root1 * sel1 + root2 * sel2

    return F

def get_ku(T, Tm, dH, Cp):
    '''
    Unfolding constant KU vs. T
    This approximation is only valid around Tm

    Input:
    T: Temperatures in °C
    Tm: Melting temperature of protein in °C
    dH: Enthalpy change in kcal/mol unfolding
    Cp: Change of heat capacity in kcal/molK for unfolding

    Returns: equilibrium constant KU
    '''
    # Convert to K
    T_K = T +  273.15
    Tm_K = Tm + 273.15
    # Gas constant in kcal/molK
    R = 1.987E-3 
    # Free enthalpy
    dG = dH*(1-T_K/Tm_K) - Cp*(Tm_K -T_K + T_K * np.log(T_K/Tm_K))
    # Equilibrium constant
    Ku = np.exp(-dG/(R*T_K))
    return Ku

def calculate_fraction_unfolded(binding_temperature, concs, fitting_params, Cp):
    '''
    Build the (isothermal) fraction folded curve at a specified temperature
    Args:
        binding_temperature: temperature to calculate fu
        concs: ligand concentrations
        fitting_params: results from local or global fit
        Cp: deltaCp of unfolding
    Returns:
    '''
    num_datasets = len(concs)
    fraction_unfolded = np.empty((concs).shape)
    for c_index in range(num_datasets):
        Tm = fitting_params[0, c_index]  # Tm
        dH = fitting_params[1, c_index]  # dH
        Ku = get_ku(binding_temperature, Tm, dH, Cp)
        fraction_unfolded[c_index] = (Ku / (1 + Ku))
    return fraction_unfolded

def calculate_fitted_isothermal(ligand_conc, Ku, Kd, pconc):
    '''
    Calculate binding curve vs. ligand concentration
    Args:
        ligand_conc: ligand concentrations
        Ku: equilibrium constant of thermal unfolding
        Kd: dissociation constant
        protein_conc:
    Returns:
    '''
    # Calculate L_free
    L_free = calculate_L_free(ligand_conc, Ku, Kd, pconc)
    # Calculate fraction unfolded
    fu = 1 / (1 + 1/Ku *(1+L_free/Kd))
    return fu

def calculate_fitted_isothermal_simple(protein_conc):
    '''
    This is a wrapper function that simplifies the isothermal curve
    to two parameters (Ku and Kd)
    This function can be used to fit Kd and Ku
    '''
    def return_func(ligand_conc, Ku, Kd):
        return calculate_fitted_isothermal(ligand_conc, Ku, Kd, protein_conc)
    return return_func

def calculate_fitted_isothermal_simple_onlyKu(protein_conc, Kd):
    '''
    This is a wrapper function that simplifies the isothermal curve
    to one parameter (Ku)
    This function can be used to fit Ku
    '''
    def return_func(ligand_conc, Ku):
        return calculate_fitted_isothermal(ligand_conc, Ku, Kd, protein_conc)
    return return_func

def calculate_fitted_isothermal_dim(ligand_conc, Ku, Kdim):
    '''
    Calculate binding curve vs. protein concentration
    Args:
        ligand_conc: ¡¡¡¡¡¡ PROTEIN CONCENTRATION !!!!!!
        Ku: equilibrium constant of thermal unfolding
        Kdim: dissociation constant of the protein dimer
    Returns:
    '''
    # Calculate folded amount
    F = calculate_F_free(ligand_conc, Ku, Kdim)
    # Calculate unfolded amount
    U  = F_free*Ku
    # Calculate dimer amount
    FF = Kdim * F**2
     # Calculate fraction unfolded
    fu = U / (U+F+FF)
    return fu

def calculate_fitted_isothermal_simple_dim():
    '''
    This is a wrapper function that simplifies the isothermal curve
    to two parameters (Ku and Kdim)
    This function can be used to fit Kdim and Ku
    '''
    def return_func(ligand_conc, Ku, Kdim):
        return calculate_fitted_isothermal_dim(ligand_conc, Ku, Kdim) #ligand_conc vector has protein concentrations!!! 
    return return_func

def calculate_fitted_isothermal_onlyKd(protein_conc, Ku):
    '''
    This is a wrapper function that simplifies the isothermal curve
    to one parameter (Kd)
    This function can be used to fit Kd
    '''
    def return_func(ligand_conc, Kd):
        return calculate_fitted_isothermal(ligand_conc, Ku, Kd, protein_conc)
    return return_func 

def calculate_fitted_isothermal_2kds(conc_lig, Ku, Kd1,Kd2,  conc_prot):
    #Kd1 = Kd2 = Kd Uncomment if you want to fit only 1 kd and two sites
    # Cooperativity (experimental)
    c = 1 
    B = Kd1 + Kd2 + (2*conc_prot - conc_lig) / c
    C = (conc_prot - conc_lig) * (Kd1+Kd2) + Kd1*Kd2 * (1+Ku)
    D = -conc_lig*Kd1*Kd2 * (1+Ku)
    # More variables
    p = B * c
    q = C * c
    r = D * c
    # Solutions
    m = (3*q - p**2)/3
    n = (2*p**3 - 9*p*q + 27*r)/27

    sel1  = np.greater(4*m**3 + 27*n**2,np.zeros(len(n)))
    sel2  = np.greater(np.zeros(len(n)),m)
    cond1 = np.logical_and(sel1,sel2)
    
    sel1  = np.greater(4*m**3 + 27*n**2,np.zeros(len(n)))
    sel2  = np.greater(m,np.zeros(len(n)))
    cond2 = np.logical_and(sel1,sel2) 
    
    cond3 = np.less(4*m**3 + 27*n**2,np.zeros(len(n))) 

    A_free1 =  oneRoot1(n,m,p)  *   cond1
    A_free2 =  oneRoot2(n,m,p)  *   cond2
    A_free3 =  three_real_aux(n,m,p,conc_lig) *   cond3

    A_free = A_free3 + A_free1 + A_free2 
        
    B_free = Kd1*Kd2*(conc_lig - A_free) / ( Kd1+Kd2+2*A_free/c) / A_free
    
    # Set B_free to protein concentration if ligand concentration is 0
    cond4  = np.equal(conc_lig,0)

    B_free = np.nan_to_num(B_free,nan=0)
    B_free = np.add(conc_prot*cond4,B_free)

    BA = A_free*B_free / Kd2
    AB = A_free*B_free / Kd1
    ABA = AB * A_free / ( Kd2*c)
    # then use L_free to get fraction unfolded
    fit_fraction_unfolded = Ku*B_free / (Ku*B_free + B_free + BA + AB + ABA)
    return fit_fraction_unfolded

def calculate_fitted_isothermal_2kds_simple2(conc_prot):
    '''
    This is a wrapper function that simplifies the isothermal curve
    to only two parameters (Kd and Ku)
    This function can be used to fit Kd and Ku
    '''
    def return_func(conc_lig, Ku, Kd):
        return calculate_fitted_isothermal_2kds(conc_lig, Ku, Kd, Kd, conc_prot)
    return return_func

def calculate_fitted_isothermal_2kds_simple2_onlyKu(protein_conc, Kd):
    '''
    This is a wrapper function that simplifies the isothermal curve
    to one parameter (Ku)
    This function can be used to fit Ku
    '''

    def return_func(ligand_conc, Ku):
        return calculate_fitted_isothermal_2kds(ligand_conc, Ku, Kd,Kd, protein_conc)
    return return_func

def calculate_fitted_isothermal_2different_kds_simple2(conc_prot):

    def return_func(conc_lig, Ku, Kd1, Kd2):
        return calculate_fitted_isothermal_2kds(conc_lig, Ku, Kd1, Kd2, conc_prot)
    return return_func


def three_real_aux(n,m,p,l):

    # Auxilary function to solve cubic equation for 3 real roots

    x0 = (2 * np.sqrt(-m/3) * np.cos(1/3*np.arccos(3*n/(2*m)*np.sqrt(-3/m)) ) ) - p/3
    x1 = (2 * np.sqrt(-m/3) * np.cos(1/3*np.arccos(3*n/(2*m)*np.sqrt(-3/m)) - 2*np.pi/3) ) - p/3
    x2 = (2 * np.sqrt(-m/3) * np.cos(1/3*np.arccos(3*n/(2*m)*np.sqrt(-3/m)) - 4*np.pi/3) ) - p/3

    x0 = np.nan_to_num(x0,nan=0.0)
    x1 = np.nan_to_num(x1,nan=0.0)
    x2 = np.nan_to_num(x2,nan=0.0)

        
    c1   = np.less(x0,l*1.02)
    c2   = np.greater(x0,np.zeros(len(x0)))
    cx0  = np.logical_and(c1,c2)
    
    c1   = np.less(x1,l*1.02)
    c2   = np.greater(x1,np.zeros(len(x1)))
    cx1  = np.logical_and(c1,c2)   
 
    c1   = np.less(x2,l*1.02)
    c2   = np.greater(x2,np.zeros(len(x2)))
    cx2  = np.logical_and(c1,c2)

    A_free1 = x0 * cx0
    A_free2 = x1 * cx1
    A_free3 = x2 * cx2

    A_free_f = A_free1 + A_free2 + A_free3
        
    return(A_free_f)

def oneRoot1(n,m,p):
    
    # Auxilary function to solve cubic equation for 1 real root, case 1

    a_t = (-2 * abs(n)/n * np.sqrt(-m/3)* np.cosh(1/3 * np.arccosh(-3*abs(n)/2/m * np.sqrt(-3/m)))) - p/3
    
    a_t = np.nan_to_num(a_t,nan=0.0)
    
    return(a_t)

def oneRoot2(n,m,p):

    # Auxilary function to solve cubic equation for 1 real root, case 2

    a_t = (-2 * np.sqrt(m/3)*np.sinh(1/3*np.arcsinh(3*n/2/m * np.sqrt(3/m))))  - p/3
    
    a_t = np.nan_to_num(a_t,nan=0.0)

    return(a_t)
    
# Helper functions to obtain the asymmetric confidence interval
# Vaida Paketuryte et al., 2021. "Uncertainty in protein–ligand binding constants: asymmetric
# confidence intervals versus standard errors"

#n is the number of data points, p the number of parameters, alfa is the desired confidence interval
def rss_p(rrs0,n,p,alfa):

    critical_value = stats.f.ppf(q=1-alfa, dfn=1, dfd=n-p) 

    return rrs0 * ( 1 + critical_value / (n-p) )
  
def get_rss(yVals,yPredicted):

    residuals = yVals - yPredicted

    return np.sum(residuals**2)

def get_desired_rss(yVals,yPredicted,n,p):

    rss = get_rss(yVals,yPredicted)

    return rss_p(rss,n,p,0.05)

def kd_to_residual_oneKd(fitting_function,fu,lig_conc,start_Ku):

    """

    fitting_function has a fixed value of Kd.
    We use this function to fit Ku and obtain the residual sum of squares

    fu is the measured fraction unfolded

    """

    params, pcov = curve_fit(fitting_function, lig_conc, fu, max_nfev=1E3, method='trf',
        p0=[start_Ku], bounds=((0), (np.inf)))

    fuPred      = fitting_function(lig_conc, *params)

    return get_rss(fu,fuPred)

def get_asymmetric_ci95(kd_estimated,fu,lig_conc,start_Ku,rss_desired,pconc,model):

    def f_to_optimize(kd):


        kd = kd / 1e9 # back to molar 

        if model == "One_Site":

            fit_function_for_asymmetric_ci = calculate_fitted_isothermal_simple_onlyKu(pconc,kd)

        if model == "Two_Sites_One_Kd":

            fit_function_for_asymmetric_ci = calculate_fitted_isothermal_2kds_simple2_onlyKu(pconc,kd)

        rss = kd_to_residual_oneKd(fit_function_for_asymmetric_ci,fu,lig_conc,start_Ku)

        return np.abs(rss - rss_desired)

    boundsMin = np.array([kd_estimated/500*1e9,kd_estimated*1e9]) # *1e9 to nanomolar to help minimizer   
    boundsMax = np.array([kd_estimated*1e9,kd_estimated*1e9*500]) # *1e9 to nanomolar to help minimizer   

    kd_min95 = minimize_scalar(f_to_optimize, bounds=boundsMin,method='bounded')
    kd_max95 = minimize_scalar(f_to_optimize, bounds=boundsMax,method='bounded')

    kd_min95, kd_max95 = kd_min95.x / 1e9, kd_max95.x / 1e9

    return np.array([kd_min95, kd_max95])
