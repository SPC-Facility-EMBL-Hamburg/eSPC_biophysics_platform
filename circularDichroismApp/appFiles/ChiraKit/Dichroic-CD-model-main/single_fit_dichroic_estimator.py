# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 13:23:05 2023

@author: urosz

This code was later edited by Osvaldo Burastero to allow integration in the ChiraKit online tool, European Molecular Biology Laboratory, Hamburg, 2024

"""

import numpy as np
import re
import os
from scipy import optimize


# Constant variables defined below!
#  v**n x w**m MATRIX
v_pow_max,w_pow_max = 5,30 #max powers of v, w for v**n x w**m matrix -SAME for all the polynomials


#Extract powers of v and w from a given partition function polynomial term.
def extract_powers(term): 
    v_power, w_power = 0, 0
    v_match = re.search(r'v\*\*(\d+)', term)
    w_match = re.search(r'w\*\*(\d+)', term)
    if v_match:
        v_power = int(v_match.group(1))
    if w_match:
        w_power = int(w_match.group(1))
    # If there's a v or w without power, assign power as 1
    if 'v' in term and v_power == 0:
        v_power = 1
    if 'w' in term and w_power == 0:
        w_power = 1
    return v_power, w_power

#Extract coefficient (integer) from a term.
def extract_coefficient(term):
    # Checking for terms that start with 'v' or 'w' without any numeral in front
    if term.startswith('v') or term.startswith('w'):
        return 1
    match = re.search(r'(\d+)\*', term)
    if match:
        return int(match.group(1))
    return 1

#Counting of single and double H-bonded states in single-helix states (up to v**5)
def vw_powers(v_st,w_st,n):
    if (v_st==2 or v_st==3 or v_st==4 or v_st==5) and w_st > 0:
        if w_st >3:
            h2= (w_st + 3) - 6
        else:
            h2 =0
        h1 = (w_st + 3) - h2
    else:
        h1=0
        h2=0
                
    c = (n+1) - (h1+h2)
        
    return (h1,h2,c)

#Counting of single and double H-bonded states in double-helix states
def vw_powers_double(v_st,w_st,n):
    h1_d,h2_d,c_d=0,0,0
    
    if w_st<5 and w_st>1:
        h2_d=0
        h1_d=w_st+6
        
    elif w_st==5:
        h2_d=0.5
        h1_d=10.5
        
    elif w_st>5:
        nw_2=w_st//2
        
        h2_d=((w_st-6)*(nw_2-2)+(w_st-5)+(w_st-4))/nw_2
        h1_d=((nw_2-2)*12 + 21)/nw_2
                
    c_d = (n+1) - (h1_d+h2_d)
        
    return (h1_d,h2_d,c_d)

"""
SPECTROSCOPIC WEIGHTS matrix generator
"""

def spectro_matrices(v_Pow_max,w_Pow_max): #RETURNS matrices to calculate contribution of each P term

    matrix_h1 = np.zeros((v_pow_max+1 ,w_pow_max+1), dtype=int) # matrix of single-bonded helix 
    matrix_h2 = np.zeros((v_pow_max+1 ,w_pow_max+1), dtype=int)
    matrix_double_h1 = np.zeros((v_pow_max+1 ,w_pow_max+1), dtype=float)
    matrix_double_h2 = np.zeros((v_pow_max+1 ,w_pow_max+1), dtype=float)
    
    for vv in range(v_Pow_max+1):
        for ww in range(w_Pow_max+1):
            h1,h2,c=vw_powers(vv,ww,100)
            if vv==4 or vv==5:
                h1_d,h2_d,c_d=vw_powers_double(vv,ww,100)
            else:
                h1_d,h2_d=0,0
                
            matrix_h1[vv,ww]=h1
            matrix_h2[vv,ww]=h2
            matrix_double_h1[vv,ww]=h1_d
            matrix_double_h2[vv,ww]=h2_d
            
    #CONTRIBUTIONS of EACH TERM in POLYNOMIAL
    M_H2=matrix_h2.transpose().flatten()
    M_H1=matrix_h1.transpose().flatten()
    M_double_H1=matrix_double_h1.transpose().flatten()
    M_double_H2=matrix_double_h2.transpose().flatten() 
    
    return (M_H2,M_H1,M_double_H2,M_double_H1)

# Partition function polynomial to MATRIX (array)
def polynomial_to_matrix(poly,v_Pow_max,w_Pow_max):
    terms = poly.split('+')
    max_v_power, max_w_power = v_Pow_max,w_Pow_max

    matrix = np.zeros((max_v_power + 1,max_w_power + 1), dtype=int)

    for term in terms:
        term=term.replace(" ","")
        v_power, w_power = extract_powers(term)
        coefficient = extract_coefficient(term)

        if v_power==0 and w_power==0:
            matrix[v_power,w_power] = int(term)
        else:
            matrix[v_power,w_power] = coefficient
    return matrix

def matrix_vw(v,w,v_pow_max,w_pow_max):
    """        w_j**k                 """
    w_powers = np.arange(w_pow_max + 1)
    w_pow = np.power(w[:, None], w_powers)
    
    """         v**i                  """
    v_powers = np.arange(v_pow_max + 1)
    v_pow = np.power(v, v_powers)   
       
    """     (v**i)*(w_j)**k           """
    wv = v_pow[:,None,None]*w_pow
    
    return wv.transpose(1, 2, 0).reshape(len(w), -1)

def run_helicity_estimation(poly_double_all,poly_total_all,n_input,t_input,Y):

    '''

    n_input = int   # define number of peptide bonds -> N_res-1 for unblocked and N_res+1 for blocked (N- and C-) peptides
    t_input = float # experimental temperature in degrees Celsius / Kelvin
    Y       = float # measured MRE (molar residue ellipticity)  [deg cm2 dmol-1 per peptide bond] 

    Careful, MRE is here normalised by the number of peptide bonds, not by the number of residues! 

    '''

    n_input = int(n_input)
    
    # convert to Celsius if input was in Kelvin
    if t_input > 273.15:

        t_input = t_input - 273.15

    poly_double = poly_double_all[n_input-4] # Because index '0' corresponds to four peptide bonds
    poly_total  = poly_total_all[n_input-4]  #

    t_C = np.array([t_input])
    Ns  = [n_input]

    coef_matrix_polys_total  = [polynomial_to_matrix(poly_total,v_pow_max,w_pow_max).transpose().flatten()] #overall coeficients
    coef_matrix_polys_double = [polynomial_to_matrix(poly_double,v_pow_max,w_pow_max).transpose().flatten()]#double helices
    coef_matrix_polys_single = [coef_matrix_polys_total[0]-coef_matrix_polys_double[0] ] #single helices

    M_H2,M_H1,M_double_H2,M_double_H1 = spectro_matrices(v_pow_max, w_pow_max)


    MH_matrices = M_H2,M_H1,M_double_H2,M_double_H1
    temps       = t_C
    data_cd     = Y

    ################
    #
    # Default model parameters obtained from global fit of CD data (see Article)
    #
    ################

    DICHROIC_MODEL_FIT_PARAMS=[0.07,-41000,100,3.4,2100,-45] # v, H2_inf, dH2/dT, k, Coil_T0, dCoil/dT

    def MODEL_CD(ww): # fitting function -> returns difference between input ellipticity (Y) and model-calculated value
        
        w= np.array([ww])
        v,H2,dH2,k,C,dC=DICHROIC_MODEL_FIT_PARAMS

        M_H2,M_H1,M_double_H2,M_double_H1 = MH_matrices

        VW = matrix_vw(v,w,v_pow_max,w_pow_max)

        #VSI
        matrix_polys_tot = [coef_matrix_polys_total[i]*VW for i in range(len(Ns))]
        P_matrixs_tot = [matrix_poly_tot/np.sum(matrix_poly_tot,axis=1)[:,None] for matrix_poly_tot in matrix_polys_tot] # each term divided by total sum of polynomial- normalization -> PROBABILITIES

        matrix_polys_2 = [coef_matrix_polys_double[i]*VW for i in range(len(Ns))]
        P_matrixs_2 = [matrix_poly_2/np.sum(matrix_poly_tot,axis=1)[:,None] for matrix_poly_2,matrix_poly_tot in zip(matrix_polys_2,matrix_polys_tot)] # each term divided by total sum of polynomial- normalization -> PROBABILITIES

        cd = [((P_matrix_tot-P_matrix_2)*( M_H2[None,:]*(H2 + dH2*temps)[:,None] +M_H1[None,:]*((H2 + dH2*temps)*(1-k/6))[:,None] +(np.full(M_H2.shape, ii+1, dtype=int)-(M_H2+M_H1))[None,:] *(C  + dC *temps)[:,None])+
               
               (P_matrix_2)*( M_double_H2[None,:]*(H2 + dH2*temps)[:,None] + M_double_H1[None,:]*((H2 + dH2*temps)*(1-k/6))[:,None] +(np.full(M_double_H2.shape, ii+1, dtype=int)-(M_double_H2+M_double_H1))[None,:] *(C  + dC *temps)[:,None])).sum(axis=1)/(ii+1) for P_matrix_tot,P_matrix_2,ii in zip(P_matrixs_tot,P_matrixs_2,Ns)]
        
        fh = [((P_matrix_tot-P_matrix_2)*(M_H2[None,:] + M_H1[None,:]) +(P_matrix_2)*(M_double_H2[None,:] + M_double_H1[None,:])).sum(axis=1)/(ii+1) for P_matrix_tot,P_matrix_2,ii in zip(P_matrixs_tot,P_matrixs_2,Ns)]

        return (cd[0][0]-data_cd, fh[0][0])

    def CD_opt(w_opt):
        return MODEL_CD(w_opt)[0]

    def FH(w_opt):
        return MODEL_CD(w_opt)[1]

    try:
        
        root = optimize.brentq(CD_opt,0.05,4.1,xtol=2e-6, rtol=8.8e-12) # minimize the function CD_opt to fit the input ellipticity Y
        
        return {"Helicity_fraction":round(FH(root),3),"Propagation_parameter_w":round(root,3)}

    except:
        
        #print("Input values of ellipticity out of reasonable range.")
        
        return {"Helicity_fraction":np.NaN,"Propagation_parameter_w":np.NaN}