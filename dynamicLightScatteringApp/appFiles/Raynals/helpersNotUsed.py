import numpy as np

from scipy import spatial

from scipy.optimize      import nnls
from scipy.optimize      import minimize
from scipy.optimize      import LinearConstraint

from helpers import *

"""
Cemetery of functions that were tested/explored but are not used in the finals Raynals tool
You can safely delete this file :)

All functions regarding the tikhonov_regularized_inversion are initially based on a 
Quora answer (https://scicomp.stackexchange.com/questions/10671/tikhonov-regularization-in-the-non-negative-least-square-nnls-pythonscipy) 
given by Dr. Brian Borchers (https://scicomp.stackexchange.com/users/2150/brian-borchers)
"""

def g1_from_g2_approximation(g2,beta):
        
    '''
    Approximate sqrt(x) with -x when x <0
    '''

    c1 = np.greater(g2,1)
    s1 = np.sqrt( (g2-1) / beta) # unitless

    s1 = np.nan_to_num(s1)

    c2 = np.logical_not(c1)
    s2 = - (g2-1) / beta # unitless, approximation!!!!

    return c1 * s1 + c2 * s2

def tikhonov_regularized_inversion(
    kernel: np.ndarray, alpha: float, data: np.ndarray) -> np.ndarray:

    """
    Solve x / minimize ||Ax - b|| + alpha||Ix|| 

    A is the kernel
    b is the vector with the measurements
    x is the unknown vector we want to estimate
    I is the identity matrix
    """

    data = data.reshape(-1, 1)
    I = alpha * np.eye(*kernel.shape)
    C = np.concatenate([kernel, I], axis=0)
    d = np.concatenate([data, np.zeros_like(data)])
    x, _ = nnls(C, d.flatten())
    
    return x

def tikhonov_regularized_inversion_withWeights(kernel,alpha,data):

    """    
    Solve x / minimize ||w*(Ax - b)|| + alpha||Ix|| 

    A is the kernel
    b is the vector with the measurements
    x is the unknown vector we want to estimate
    I is the identity matrix
    w is a vector of weights, in this case we will assign weights as

        sqrt(1), sqrt(1/2), sqrt(1/3), ... , sqrt(1/n)  

    """

    W      = np.arange(1,len(data)+1)**(1/2) 
    W      = W / np.max(W)    # Set max weight to one 
    W      = np.append(W,100) # Extrem weight to force sum of contributions equal to 1
    data   = np.sqrt(W) * np.append(data,1) 
    data   = data.reshape(-1, 1)
    kernel = np.vstack([kernel,np.ones(kernel.shape[1])])
    kernel = np.sqrt(W)[:, None] * kernel
    
    I      = alpha * np.eye(*kernel.shape)
    C      = np.concatenate([kernel, I], axis=0)
    d      = np.concatenate([data, np.zeros_like(data)])
    x, _   = nnls(C, d.flatten())

    return x

def tikhonov_regularized_inversion_withWeightsAndSmooth(kernel,alpha,data):

    """
    
    Solve x / minimize ||w*(Ax - b)|| + alpha||Mx|| 

    A is the kernel
    b is the vector with the measurements
    x is the unknown vector we want to estimate
    M is the second order derivative matrix
    w is a vector of weights, in this case we will assign weights as

        sqrt(1), sqrt(1/2), sqrt(1/3), ... , sqrt(1/n) 

    """

    W      = np.arange(1,len(data)+1)**(1/2) 
    W      = W / np.max(W)
    W      = np.append(W,100) # Extrem weight to force sum of contributions equal to 1
    data   = np.sqrt(W) * np.append(data,1) 
    data   = data.reshape(-1, 1)
    kernel = np.vstack([kernel,np.ones(kernel.shape[1])])
    kernel = np.sqrt(W)[:, None] * kernel
    
    cols   = kernel.shape[1]
    
    M = np.zeros((cols,cols))
    for i in range(1,M.shape[1]-1):
        M[i,i-1] = -1
        M[i,i]   =  2
        M[i,i+1] = -1
    
    L      = alpha * M
    C      = np.concatenate([kernel, L], axis=0)
        
    d      = np.concatenate([data, np.zeros(cols).reshape(-1, 1)])
    x, _   = nnls(C, d.flatten())

    return x

def tikhonov_regularized_inversion_withWeightsAndSmooth2(kernel,alpha,beta,data):

    """
    
    Solve x / minimize ||w*(Ax - b)|| + alpha||Mx|| + beta||Ix||

    A is the kernel
    b is the vector with the measurements
    x is the unknown vector we want to estimate
    M is the second order derivative matrix
    I is the identity matrix
    w is a vector of weights, in this case we will assign weights as

        sqrt(1), sqrt(1/2), sqrt(1/3), ... , sqrt(1/n)  

    Moreover, we add a last equation to the system such that sum(x) equals 1!

    """

    W      = np.arange(1,len(data)+1)**(1/2) * 0 + 1 # all weights are equal 
    W      = W / np.max(W)
    W      = np.append(W,100) # Extrem weight to force sum of contributions equal to 1

    data   = np.sqrt(W) * np.append(data,1) 
    data   = data.reshape(-1, 1)
    kernel = np.vstack([kernel,np.ones(kernel.shape[1])])
    kernel = np.sqrt(W)[:, None] * kernel
    
    cols   = kernel.shape[1]
    
    M = np.zeros((cols,cols))
    for i in range(1,M.shape[1]-1):
        M[i,i-1] = -1
        M[i,i]   =  2
        M[i,i+1] = -1
    
    L      = alpha * M
    C      = np.concatenate([kernel, L], axis=0)    
    d      = np.concatenate([data, np.zeros(cols).reshape(-1, 1)])
    
    I      = beta * np.eye(*kernel.shape)
    C      = np.concatenate([C, I], axis=0)
    d      = np.concatenate([d, np.zeros_like(data)])
    x, _   = nnls(C, d.flatten())

    return x


def fitAutocorrelation5(times,autocorrelation,contributionsGuess,s_space,betaPrior,smoothFactor = 0):

    decay_rates = 1 / s_space

    bounds     = [(0, 1) for _ in range(len(contributionsGuess))] # Relative contribution can't be higher than 1
    bounds[0]  = (0,0) # 0 at the ends
    bounds[-1] = (0,0) # 0 at the ends
    bounds    += [(betaPrior*0.9,betaPrior*1.1)] # Plus minus 10 %

    def g2residuals(decay_rates,times,alpha,g2measured,beta,contributions):
        
        error   = g2measured - g2_finite_aproximation(decay_rates,times,beta,contributions)
        error   = 0.5 * np.sum(error**2) / len(g2measured)
        
        delta   = 2*contributions[1:-1] - contributions[:-2] - contributions[2:] 
        
        regularizor = alpha**2 * np.sum((delta**2)) + 1e2*((np.sum(contributions)-1)**2)
                
        return error + regularizor 

    def toMinimize(contributionsAndBeta):

        contributions = contributionsAndBeta[:-1]
        beta          = contributionsAndBeta[-1]  # Last value is beta
    
        return g2residuals(decay_rates,times,smoothFactor,autocorrelation,beta,np.array(contributions))

    contributionsAndBeta = np.append(contributionsGuess,betaPrior)

    res = minimize(toMinimize,contributionsAndBeta,bounds=bounds)

    betaEst, contEst = res.x[-1], res.x[:-1]

    return betaEst, contEst 

def g2_residuals(decay_rates,times,alpha,g2measured,beta,contributions,ridgePenalty):
        
    """

    Useful to fit the experimental second order autocorrelation using a NNLS approach

    alpha is a regularization parameter for smoothing the obtained decay rates
    it affects the second derivative of the contributions

    Returns the error

    """

    error   = g2measured - g2_finite_aproximation(decay_rates,times,beta,contributions)

    # Assign weights for the fitting (based on the fact that the ith point in the autocorrelation function 
    # was estimated from ith measurements !)
    weights = np.arange(len(g2measured))[::-1] / len(g2measured)

    error   = np.sum(weights*(error**2))

    # Get second derivative
    x     = decay_rates      

    y     = contributions
    dydx  = np.gradient(y,1/x)    
    dy2dx = np.gradient(dydx,1/x)
            
    regularizor = alpha**2 * np.sum((dy2dx**2))
    ridge       = 10 * np.sum(contributions**2)   
    
    return error + regularizor + ridge

def fitAutocorrelation(times,autocorrelation,contributionsGuess,s_space,betaPrior,smoothFactor = 0.08):

    decay_rates = 1 / s_space

    bounds     = [(0, 1) for _ in range(len(contributionsGuess))] # Relative contribution can't be higher than 1
    bounds[0]  = (0,0) # 0 at the ends
    bounds[-1] = (0,0) # 0 at the ends
    bounds    += [(betaPrior*0.9,betaPrior*1.1)] # Plus minus 10 %

    # Total contributions should add 1 !!!
    linear_constraint = LinearConstraint([[1 for _ in range(len(contributionsGuess))]+[0]], [1], [1])

    def toMinimize(contributionsAndBeta):

        contributions = contributionsAndBeta[:-1]
        beta          = contributionsAndBeta[-1]  # Last value is beta
    
        return g2_residuals(decay_rates,times,smoothFactor,autocorrelation,beta,contributions)

    contributionsAndBeta = np.append(contributionsGuess,betaPrior)

    res = minimize(toMinimize,contributionsAndBeta, 
                    bounds=bounds,
                    constraints=linear_constraint,
                    options={'maxiter':1e4})

    betaEst = res.x[-1]
    contEst = res.x[:-1]

    return betaEst, contEst