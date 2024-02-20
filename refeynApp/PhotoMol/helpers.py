import numpy as np
from scipy.signal import find_peaks 
 
def findPeaks(x,histogram_centers,height=14, distance=4,prominence=8,masses=True):
    """
    x is 1D vector where each observation is equally spaced in a 2nd (not given) dimension
    x could be a time series, i.e. [23,21,16,22,19] where each value is the mean temperature of a day
    In our case, x are histogram counts of the observed masses (or constrasts) from a mass photometry experiment

    """

    if masses:

        sel1     = histogram_centers < 650
        sel2     = np.logical_and(histogram_centers >= 650, histogram_centers  < 1500)
        sel3     = np.logical_and(histogram_centers >= 1500, histogram_centers < 5000)
        
        pks1     = find_peaks(x[sel1], height=height,   distance=distance,   prominence=prominence)[0]
        pks2     = find_peaks(x[sel2], height=height*2, distance=distance*3, prominence=prominence)[0]
        pks3     = find_peaks(x[sel3], height=height*5, distance=distance*8, prominence=prominence*2)[0]

        pks = ( pks1, pks2+len(x[sel1]), pks3 + len(x[sel1]) + len(x[sel2]) )
        pks = np.concatenate(pks)

    # For contrasts
    else:

        sel1     = histogram_centers <   -0.1
        sel2     = histogram_centers >=  -0.1
        
        pks1     = find_peaks(x[sel1], height=height*3,   distance=distance*20,   prominence=prominence)[0]
        pks2     = find_peaks(x[sel2], height=height, distance=distance, prominence=prominence)[0]

        pks = ( pks1, pks2+len(x[sel1]))
        pks = np.concatenate(pks)

    return histogram_centers[pks]

def multi_gauss(x,*params):
    '''
    Multiple gaussian function
    Inputs x values and gaussian parameters
    Outputs y values
    '''

    y = np.zeros_like(x)
    for i in range(0, len(params), 3):
        ctr = params[i]
        amp = params[i+1]
        std = params[i+2]

        gaussian = amp*np.exp(-(x-ctr)**2 / (2*(std**2)))

        y = y + gaussian

    return y

def truncated_multi_gauss_with_baseline(x,lower_limit_of_resolution,baseline,*params):
    
    #lower_limit_of_resolution' is the lower limit of a mass that can be observed, i.e. 30 for refeyn model 1
    values_have_sense = np.greater(np.abs(x),lower_limit_of_resolution)
    y = multi_gauss(x,*params)
    y = y * values_have_sense

    return y + baseline

def r_squared(data, fit):
    '''
    R squared
    '''
    mean_data = np.mean(data)
    ss_tot = np.sum(np.power(fit - mean_data, 2))
    ss_res = np.sum(np.power(data - fit, 2))
    return 1 - (ss_res/ss_tot)

def compute_contrasts_to_mass(contrasts,slope,intercept):
        '''
        Function to convert masses from
        contrasts using known calibration parameters 

        Caution! slope and intercept are based on f(mass) = contrast !!!! 
        In other words, contrast = slope*mass + intercept

        '''

        interceptInverse = -intercept / slope 
        slopeInverse     = 1 / slope          

        masses_kDa   = np.polyval(np.array([slopeInverse,interceptInverse]), contrasts)

        return masses_kDa