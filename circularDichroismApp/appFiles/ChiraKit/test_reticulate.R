rm(list=ls())
gc()
library(reticulate)
library(reshape2)
library(tidyverse)
library(plotly)
library(factoextra)

appName     <- "ChiraKit"
user        <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

setwd('/home/osvaldo/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/')

source("server_files/helpers.R")


source_python("python_src/cdUnitsConverter.py")
source_python("python_src/helpers.py")
source_python('python_src/loadCDfilesHelpers.py')

source_python("python_src/cdAnalyzer.py")

cdAnalyzer <- cdAnalyzer()
group <- 'A'
cdAnalyzer$clean_experiments('chemical')
cdAnalyzer$experimentNamesChemical <- c(group)

cdAnalyzer$experimentsChemical[[group]]   <- cd_experiment_chemical_unfolding()

cdAnalyzer$experimentsChemical[[group]]$load_unfolding_data('/home/osvaldo/Downloads/test.csv')
cdAnalyzer$experimentsChemical[[group]]$reshape_signal_oligomer('Chemical')

source("server_files/helpers_unfolding.R")
df <- generate_chemical_unfolding_df(cdAnalyzer)
nrow(df)




