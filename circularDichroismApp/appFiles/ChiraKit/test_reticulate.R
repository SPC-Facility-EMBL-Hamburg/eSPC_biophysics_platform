rm(list=ls())
gc()
library(reticulate)
library(reshape2)
library(tidyverse)
library(plotly)
library(ggdendro)

#use_condaenv("r-reticulate",required = TRUE)

setwd('/home/osvaldo/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/')

source("server_files/helpers.R")
source_python("./cdUnitsConverter.py")
source_python("./helpers.py")
source_python('./loadCDfilesHelpers.py')

source_python("./cdAnalyzer.py")

cd = cdAnalyzer()
cd$loadExperiment('./www/Pkng_WT_0.7mgperml_12.11.2019.csv','t')

cd$loadExperiment('./www/R78907.d01','t2')
cd$initializeExperimentModif()



