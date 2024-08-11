library(reshape2)
library(tidyverse)
library(reticulate)
library(pracma)
library(data.table)
library(grid)
library(Cairo)

user      <- Sys.info()['user']

reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/foldA_moltenP_apps/moltenprot/")

#  path for the docker user
if (user == 'shiny') {
  base_dir <- "/home/shiny/moltenprot/"
}

global_chunck_n     <- 16 # should be global_plot_columns * global_plot_rows
global_plot_columns <- 4
global_plot_rows    <- 4

