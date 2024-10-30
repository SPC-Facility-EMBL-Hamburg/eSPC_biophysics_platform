# -*- coding: utf-8 -*-
"""
Created on Thu Dec  7 08:01:51 2023

@author: urosz
"""

import numpy as np
import re

# Constant variables defined below!
#  v**n x w**m MATRIX
v_pow_max,w_pow_max = 5,38 #max powers of v, w for v**n x w**m matrix -SAME for all the polynomials (max=38 for peptide with 40 units)


# Extract powers of v and w from a given partition function polynomial term.
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


# Extract coefficient (integer) from a term.
def extract_coefficient(term):
    # Checking for terms that start with 'v' or 'w' without any numeral in front
    if term.startswith('v') or term.startswith('w'):
        return 1
    match = re.search(r'(\d+)\*', term)
    if match:
        return int(match.group(1))
    return 1


# Counting of single and double H-bonded microstates in single-helix states (up to v**5)
def vw_powers(v_st, w_st, n):
    if (v_st == 2 or v_st == 3 or v_st == 4 or v_st == 5) and w_st > 0:
        if w_st > 3:
            h2 = (w_st + 3) - 6
        else:
            h2 = 0
        h1 = (w_st + 3) - h2
    else:
        h1 = 0
        h2 = 0

    c = (n + 1) - (h1 + h2)

    return (h1, h2, c)


# Counting of single and double H-bonded states in double-helix microstates.
def vw_powers_double(v_st, w_st, n):
    h1_d, h2_d, c_d = 0, 0, 0

    if w_st < 5 and w_st > 1:
        h2_d = 0
        h1_d = w_st + 6

    elif w_st == 5:
        h2_d = 0.5
        h1_d = 10.5

    elif w_st > 5:
        nw_2 = w_st // 2

        h2_d = ((w_st - 6) * (nw_2 - 2) + (w_st - 5) + (w_st - 4)) / nw_2
        h1_d = ((nw_2 - 2) * 12 + 21) / nw_2

    c_d = (n + 1) - (h1_d + h2_d)

    return (h1_d, h2_d, c_d)


# RETURNS matrices to calculate CD contribution of each partition function polynomial term
def spectro_matrices(v_Pow_max, w_Pow_max):
    matrix_h1 = np.zeros((v_Pow_max + 1, w_Pow_max + 1), dtype=int)  # matrix of single helices
    matrix_h2 = np.zeros((v_Pow_max + 1, w_Pow_max + 1), dtype=int)  # matrix of double helices
    matrix_double_h1 = np.zeros((v_Pow_max + 1, w_Pow_max + 1), dtype=float)
    matrix_double_h2 = np.zeros((v_Pow_max + 1, w_Pow_max + 1), dtype=float)

    for vv in range(v_Pow_max + 1):
        for ww in range(w_Pow_max + 1):
            h1, h2, c = vw_powers(vv, ww, 100)
            if vv == 4 or vv == 5:
                h1_d, h2_d, c_d = vw_powers_double(vv, ww, 100)
            else:
                h1_d, h2_d = 0, 0

            matrix_h1[vv, ww] = h1
            matrix_h2[vv, ww] = h2
            matrix_double_h1[vv, ww] = h1_d
            matrix_double_h2[vv, ww] = h2_d

    # CONTRIBUTIONS of EACH TERM in POLYNOMIAL
    M_H2 = matrix_h2.transpose().flatten()
    M_H1 = matrix_h1.transpose().flatten()
    M_double_H1 = matrix_double_h1.transpose().flatten()
    M_double_H2 = matrix_double_h2.transpose().flatten()

    return (M_H2, M_H1, M_double_H2, M_double_H1)


# Partition function polynomial to MATRIX (array)
def polynomial_to_matrix(poly, v_Pow_max, w_Pow_max):
    terms = poly.split('+')
    max_v_power, max_w_power = v_Pow_max, w_Pow_max

    matrix = np.zeros((max_v_power + 1, max_w_power + 1), dtype=int)

    for term in terms:
        term = term.replace(" ", "")
        v_power, w_power = extract_powers(term)
        coefficient = extract_coefficient(term)

        if v_power == 0 and w_power == 0:
            matrix[v_power, w_power] = int(term)
        else:
            matrix[v_power, w_power] = coefficient
    return matrix


# returns VxW matrix to calculate all possible values of v**n x w**m states for n<v_pow_max and m<w_pow_max (maximum-lenght helix)
def matrix_vw(v, w, v_pow_max, w_pow_max):
    """        w_j**k                 """
    w_powers = np.arange(w_pow_max + 1)
    w_pow = np.power(np.array(w)[:, None], w_powers)

    """         v**i                  """
    v_powers = np.arange(v_pow_max + 1)
    v_pow = np.power(v, v_powers)

    """     (v**i)*(w_j)**k           """
    wv = v_pow[:, None, None] * w_pow

    return wv.transpose(1, 2, 0).reshape(len(w), -1)