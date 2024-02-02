library(reshape2)
library(tidyverse)
library(reticulate)
library(pracma)
library(data.table)
library(grid)
library(Cairo)

user      <- Sys.info()['user']
users_dir <- paste0("/home/",user,"/data_users/")

notebook_app  <- (Sys.info()["nodename"] == "osvaldo")

if (notebook_app) {
  use_python("/home/osvaldo/miniconda3/bin/python")
} else  {
  reticulate::use_condaenv("r-reticulate",required = TRUE)
} 

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/foldA_moltenP_apps/moltenprot/")

# set the corrrect path for the docker user
if (user == 'shiny') {
  base_dir <- "/home/shiny/moltenprot/"
}

global_chunck_n     <- 16 # should be global_plot_columns * global_plot_rows
global_plot_columns <- 4
global_plot_rows    <- 4

