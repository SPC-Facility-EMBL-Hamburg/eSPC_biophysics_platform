# load libraries
library(reshape2)
library(tidyverse)
library(reticulate)
library(pracma)
library(data.table)
require(scales)
library(viridis)

setwd("/home/os/spc_shiny_servers/thermoA_app/thermoAffinity/")
source("server_files/global_variables.R")

source_python("helpers.py")
source_python("mst.py")

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
#source("server_files/simulation_helpers.R")
source("server_files/plot_functions.R")
#source("server_files/server_download_helpers.R")

scripts <- list.files("server_files/fitting_helpers/")
lapply(paste0("server_files/fitting_helpers/",scripts), source)

if (TRUE) {
  mst <- MST_fit()
  mst$load_MST_xlsx("./www/demo.xlsx")
  mst$set_signal("Raw Fluorescence")
  
  print(mst$concs)
  mst$get_cold_fluo(-1,0)
  mst$get_fnorm(0,1)
  
  mst$F_norm  
  
  fit <- fit_fluo_1_site(mst$F_norm,mst$concs*1e6,rep(0.025,length(mst$concs)),
                         39,31,47,F,
                         39,31,47,F,
                         0.8,0.01,1e2,F)
  augment(fit$fit_obj)
  fit$tidy_fit
  

}








  
  