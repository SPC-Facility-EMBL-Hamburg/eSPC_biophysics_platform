# -*- coding: utf-8 -*-
"""
Created on Wed Sep 27 13:23:05 2023

@author: urosz

This code was later minimally edited by Osvaldo Burastero to allow integration in the ChiraKit online tool, 
European Molecular Biology Laboratory, Hamburg, 2024

"""

from scipy         import optimize
from CD_model_Uros import *

def run_helicity_estimation(poly_double_all,poly_total_all,n_input,t_input,Y):

    n_input = int(n_input)
    
    # Convert to Celsius if input was in Kelvin
    if t_input > 273.15:

        t_input = t_input - 273.15

    poly_double = poly_double_all[n_input-4] # Because index '0' corresponds to four peptide bonds
    poly_total  = poly_total_all[n_input-4]  #

    t_C = np.array([t_input])
    Ns  = [n_input]

    coef_matrix_polys_total =[polynomial_to_matrix(poly_total,v_pow_max,w_pow_max).transpose().flatten()]   # Overall coeficients
    coef_matrix_polys_double=[polynomial_to_matrix(poly_double,v_pow_max,w_pow_max).transpose().flatten()] # Double  helices
    coef_matrix_polys_single=[coef_matrix_polys_total[0]-coef_matrix_polys_double[0] ]                     # Single  helices

    M_H2,M_H1,M_double_H2,M_double_H1 = spectro_matrices(v_pow_max, w_pow_max)

    MH_matrices = M_H2,M_H1,M_double_H2,M_double_H1

    temps   = t_C
    data_cd = Y

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
        return {"Helicity_fraction":np.nan,"Propagation_parameter_w":np.nan}
