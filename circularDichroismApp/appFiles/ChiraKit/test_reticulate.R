rm(list=ls())
gc()
library(reticulate)
library(reshape2)
library(tidyverse)
library(plotly)
use_condaenv("r-reticulate",required = TRUE)

setwd('/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/')

source("server_files/helpers.R")
source_python("./cdAnalyzer.py")
source_python("./cdUnitsConverter.py")
source_python("./helpers.py")
source_python('./loadCDfilesHelpers.py')

t1 = cd_experiment_general()
t2 = cd_experiment_general()

t1$load_data('/home/os/Downloads/1ITa_peptide.csv','t')
t2$load_data('/home/os/Downloads/buffer_background.csv','t')

# List of signal dataframes 
signalSel <- list(t1$signalInput)
htSel     <- list(t1$signalHT)

temperatureSel <- list(t1$temperature)

# Assign the mean of the sample scans temperature to the rest of the curves
temps          <- unlist(temperatureSel)

if (all(is.null(temps))) {
  newTemp <- NA
} else {
  newTemp        <- mean(temps,na.rm = T)
}

signalSelVectorized <- list()
htSelVectorized     <- list()

spectra_counter    <- 0
experiment_counter <- 0

# Transform list of signal dataframes into a list of signal vectors
for (signalDF in signalSel) {
  
  experiment_counter <- experiment_counter + 1
  htDF               <- htSel[[experiment_counter]]
  
  print('here1')
  
  for (columnID in 1:ncol(signalDF)) {
    spectra_counter                        <- spectra_counter + 1
    signalSelVectorized[[spectra_counter]] <- signalDF[,columnID]
    htSelVectorized[[spectra_counter]]     <- htDF[,columnID]
  }
  
}

wlSelVectorized <- list(seq(280,180,length.out = 200),
                        seq(280,180,length.out = 200),
                        seq(280,180,length.out = 200)
)

interpolatedAverage <- averageListOfVectors(wlSelVectorized,signalSelVectorized)
signalNewSample     <- matrix(interpolatedAverage$average,ncol = 1)


