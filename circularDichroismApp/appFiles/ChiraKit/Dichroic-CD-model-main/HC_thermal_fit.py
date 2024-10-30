# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 08:01:51 2023

@author: urosz

This code was later minimally edited by Osvaldo Burastero to allow integration in the ChiraKit online tool,
European Molecular Biology Laboratory, Hamburg, 2024
"""

from CD_model_Uros  import *

# FIXED CD PARAMETERS AND NUCLEATION CONSTANT
# These are baseline CD parameters and the nucleation constant from Zavrtanik et al. (2024)
cd_params_v_fix = [-41.0, 0.1, 3.4, 2.1, -0.045, 0.07]

# Given thermodynamic parameters (dG,dH,dCp) and T-range, this function returns molar residual ellipticity (MRE) as a function of T.
def calculate_cd_signal(dg, dh, dcp, t, n,poly_double_all,poly_total_all):
    """
    Calculate the circular dichroism (CD) model using thermodynamic parameters.

    Parameters:
    -----------
    dg : float
        Gibbs free energy at reference temperature (273.15 K) in units kcal/(mol res).
    dh : float
        Enthalpy at reference temperature (273.15 K)  in units kcal/(mol res).
    dcp : float
        Heat capacity at reference temperature (273.15 K) in units kcal/(mol K res).
    t : array
        Array of temperature values in °C or Kelvin.
    n : int
        Number of residues in the peptide.
    poly_total_all : array
        list of  coefficient matrices for total helices polynomial.
    poly_double_all : array
         list of coefficient matrices for double helices polynomial.

    Returns:
    --------
    cd_model : array
        Modeled CD values as a function of temperature in units 10^3 deg cm^2/(dmol res).
    """

    poly_double = poly_double_all[int(n)-4] # Because index '0' corresponds to four peptide bonds
    poly_total  = poly_total_all[int(n)-4]  # Because index '0' corresponds to four peptide bonds

    coef_matrix_polys_total = polynomial_to_matrix(poly_total,v_pow_max,w_pow_max).transpose().flatten()   # Overall coeficients
    coef_matrix_polys_double= polynomial_to_matrix(poly_double,v_pow_max,w_pow_max).transpose().flatten()  # Double  helices

    T0 = 273.15             # Reference temperature in Kelvin
    R = 1.987 * 10 ** (-3)  # Gas constant [kcal K-1 mol-1]

    h2, dh2, k, c, dc, v                 = cd_params_v_fix  # Unpack fixed parameters
    M_H2, M_H1, M_double_H2, M_double_H1 = spectro_matrices(v_pow_max, w_pow_max)  # Unpack spectroscopic matrices

    # Convert temperatures to Kelvin, if needed
    t_K = t + 273.15 if np.max(t) < 270 else t

    # Calculate temperature-dependent enthalpy
    H_t_w = dh + dcp * (t_K - T0)

    # Temperature dependence of the w parameter
    w = np.exp((-1 / R) * (dg / T0 + H_t_w * ((1 / t_K) - (1 / T0)) + dcp * (1 - (T0 / t_K) - np.log(t_K / T0))))

    # Construct the V x W matrix
    VW = matrix_vw(v, w, v_pow_max, w_pow_max)

    # Apply polynomials to the VW matrix
    matrix_polys_tot = coef_matrix_polys_total * VW
    matrix_polys_2   = coef_matrix_polys_double * VW

    # Calculate partition function (q) and normalize for probabilities -> P_matrixs_tot,P_matrixs_2
    q_sums_t = np.sum(matrix_polys_tot, axis=1)

    P_matrixs_tot = matrix_polys_tot / q_sums_t[:,
                                       None]  # each term divided by total sum of polynomial- normalization -> PROBABILITIES
    P_matrixs_2 = matrix_polys_2 / q_sums_t[:,
                                   None]  # each term divided by total sum of polynomial- normalization -> PROBABILITIES

    # Compute the CD signal
    cd_model = np.sum((P_matrixs_tot - P_matrixs_2) * (
                M_H2[None, :] * (h2 + dh2 * t)[:, None] + M_H1[None, :] * ((h2 + dh2 * t) * (1 - k / 6))[:, None] + (
                                                                                                                                np.full(
                                                                                                                                    M_H2.shape,
                                                                                                                                    n + 1,
                                                                                                                                    dtype=int) - (
                                                                                                                                            M_H2 + M_H1))[
                                                                                                                    None,
                                                                                                                    :] * (
                                                                                                                                     c + dc * t)[
                                                                                                                         :,
                                                                                                                         None]) +

                      (P_matrixs_2) * (M_double_H2[None, :] * (h2 + dh2 * t)[:, None] + M_double_H1[None, :] * ((
                                                                                                                            h2 + dh2 * t) * (
                                                                                                                            1 - k / 6))[
                                                                                                               :,
                                                                                                               None] + (
                                                                                                                                   np.full(
                                                                                                                                       M_double_H2.shape,
                                                                                                                                       n + 1,
                                                                                                                                       dtype=int) - (
                                                                                                                                               M_double_H2 + M_double_H1))[
                                                                                                                       None,
                                                                                                                       :] * (
                                                                                                                                        c + dc * t)[
                                                                                                                            :,
                                                                                                                            None]),
                      axis=1) / (n + 1)

    fh = ((P_matrixs_tot - P_matrixs_2) * (M_H2[None, :] + M_H1[None, :]) + (P_matrixs_2) * (
                M_double_H2[None, :] + M_double_H1[None, :])).sum(axis=1) / (n + 1)

    return (cd_model,fh)


# OBTAIN OBJECTIVE FUNCTION FOR OPTIMIZATION
def get_objective_fx_peptide(n,poly_double_all,poly_total_all):
    """
    Obtain the objective function for optimization.
    """

    def objective_function(tcd,dg, dh, dcp):

        return calculate_cd_signal(dg, dh, dcp, tcd, n,poly_double_all,poly_total_all)[0]

    return objective_function


