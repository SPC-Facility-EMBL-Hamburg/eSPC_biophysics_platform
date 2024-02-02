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
t3 = cd_experiment_general()
t4 = cd_experiment_general()
t5 = cd_experiment_general()
t6 = cd_experiment_general()

t1$load_data('/home/os/Downloads/R78907.d01','t')
t2$load_data('/home/os/Downloads/R78907.d02','t')
t3$load_data('/home/os/Downloads/R78907.d03','t')
t4$load_data('/home/os/Downloads/R78906.d01','t')
t5$load_data('/home/os/Downloads/R78906.d02','t')
t6$load_data('/home/os/Downloads/R78906.d03','t')

t = cd_experiment_comparison()

signal = do.call(
  cbind,
  list(t1$signalInput, t2$signalInput,t3$signalInput,
       t4$signalInput,t5$signalInput,t6$signalInput)
  )

t$signalDesiredUnit     = np_array(signal)

t$labels     = np_array(c('c','c','c','s','s','s'))
t$wavelength = t1$wavelength

t$summarise_signal_per_label()
t$generate_comparison_labels()
t$generate_difference_spectra()
t$find_distances()
t$workingUnits <- 'meanResidueMolarExtinction'

source("server_files/helpers_plotting.R")
source("server_files/plotFunctionsSpectraComparison.R")

plot_distances(t$comparison_labels,t$distances,'s')




