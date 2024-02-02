import numpy as np
import miepython

def generateParticlesPopulation(means,sds,samplesSize):

    """
    Create a distribution of particles of different sizes using a gaussian function
    
    Input: 
        - means        : in nanometer, one for each population of particles
        - sds          : in nanometer, one for each population of particles
        - samplesSize  : one for each population of particles

    Output:

        - The generated sample 

    """

    if isinstance(means,float):
        
        return np.random.normal(means, sds, int(samplesSize)) / 1e9 # From nm to m

    hrDistributions = []

    for i in range(len(means)):

        hrDistribution     = np.random.normal(means[i], sds[i], int(samplesSize[i])) / 1e9 # From nm to m
        hrDistributions.append(hrDistribution)

    return np.concatenate(hrDistributions)

def mieIntensity(hr=np.arange(10,1000,1),angle=173,lambda0= 668,refractiveIndex=1.33-0.01j):

    """
    Obtain the scattered intensity using Mie Theory at a certain angle & wavelength
    For sizes < wavelength / 10, the scattered intensity follows the Rayleigh theory and 
    is proportional to the 6th power of the radius. 

    Input: 
        - hr        : hydrodynamic radius of the particles, in nanometers, numpy array
        - angle     : angle of the DLS detector, in degrees, float
        - lambda0   : the used wavelength, in nanometers, float
        - refractiveIndex : unitless, contains an imaginary part denoted with 'j', float

    Output:

        - The scattered intensity 
    """

    try:
        refractiveIndex = complex(refractiveIndex)
    except:
        return None

    mu      = np.cos(angle*np.pi/180) # from degrees    to radians
    hr      = hr/1e-9                 # from nanometers to meters
    lambda0 = lambda0/1e-9            # from nanometers to meters

    sigma_sca = []

    for x in hr:

        x *= 2       # Convert radius to diameter
        size_factor = np.pi/lambda0 * x
        geometric_cross_section = np.pi * x**2/4 * 1e4                          # cm**2
        sigma_sca.append(geometric_cross_section * miepython.i_unpolarized(refractiveIndex,size_factor,mu,'qsca'))

    sigma_sca = np.array(sigma_sca) 
    
    return sigma_sca