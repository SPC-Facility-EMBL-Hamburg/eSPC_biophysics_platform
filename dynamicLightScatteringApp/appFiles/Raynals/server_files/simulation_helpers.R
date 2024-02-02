source_python("./simulation_helpers.py")
# Take a vector of hydrodynamic radius and generate
# the corresponding number, volume and intensity distributions
# This code assumes sphere like particles and uses the Mie Theory  
# implemented in the miepython package
generate_distributions <- function(
    particles,angle=173,lambda0=668,refractiveIndex='1.33 - 0.01j') {
  
  # particles is a dataframe with a column called 'hr' (units in meters)
  # angle (of the DLS detector) must be given in degrees
  # lambda0 (wavelength) in nanometers
  # refractiveIndex is unitless and contains an imaginary part
  
  # Uncomment only if we change to a light absorbing medium
  # Remove all spaces to avoid error in the python code when converting to complex
  # refractiveIndex <- gsub("\\s", "", refractiveIndex)
  
  # Formula of a sphere
  particles$Volume <- 4/3 * 3.14159265 * (particles$hr**3)
  totalVolume      <- sum(particles$Volume)
  
  # Discretize the hydrodynamic radius space
  intervals <- 10**(seq(log10(0.09), log10(1e6), length.out = 200))*1e-9
  
  contributionsVolume     <- c()
  contributionsNumber     <- c()
  contributionsIntensity  <- c()
  intervalsMiddlePoint    <- c()
  
  for (i in 2:length(intervals)) {
    tempDF <- particles[particles$hr > intervals[i-1] & 
                          particles$hr <= intervals[i],    ]
    
    midHr <- 10**(seq(log10(intervals[i-1]), log10(intervals[i]), length.out = 3))[2]*1e9
    
    # Intensity using the miepython package
    tempDF$intensity       <- mieIntensity(np_array(tempDF$hr*1e9),angle,lambda0,refractiveIndex) # Input units are nanometers
    
    intervalsMiddlePoint   <- c(intervalsMiddlePoint,midHr)
    contributionsNumber    <- c(contributionsNumber,sum(nrow(tempDF)))
    contributionsVolume    <- c(contributionsVolume,sum(tempDF$Volume))
    contributionsIntensity <- c(contributionsIntensity,sum(tempDF$intensity))
    
  }
  
  contributionsNumber    <- contributionsNumber / nrow(particles)
  contributionsVolume    <- contributionsVolume / totalVolume
  contributionsIntensity <- contributionsIntensity / sum(contributionsIntensity)
  
  return(list('hr'=intervalsMiddlePoint,
              'contributionsNumber'=contributionsNumber,
              'contributionsVolume'=contributionsVolume,
              'contributionsIntensity'=contributionsIntensity))
  
}


