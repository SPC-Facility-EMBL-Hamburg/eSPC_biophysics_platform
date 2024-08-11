# load libraries
library(reshape2)
library(Cairo)
library(tidyverse)
library(reticulate)
library(pracma)
library(data.table)
require(scales)
library(viridis)
library(minpack.lm)
library(broom)
library(grid)

user      <- Sys.info()['user']

notebook_app  <- (Sys.info()["nodename"] == "osvaldo")

if (notebook_app) {
  use_python("/home/osvaldo/miniconda3/bin/python")
} else  {
  reticulate::use_python(paste0("/home/",user,"/myenv/bin/python"), required = TRUE)
}

# developer path
base_dir <- paste0("/home/",user,"/spc_shiny_servers/foldA_moltenP_apps/foldAffinity/")

# set the corrrect path for the docker user
if (user == 'shiny') {
  base_dir <- "/home/shiny/foldAffinity/"
}

# Number of elements per plot in the fluorescence fit
# (the plots are plot_columns X plot_rows grids

global_n_rows_conditions_table <- 24
global_chunck_n <- 16    # should be global_plot_columns * global_plot_rows

global_plot_columns <- 4
global_plot_rows    <- 4

# viridis color palette
global_colors_palette_signal <- c("#440154","#471063","#481d6f","#472a7a",
                                  "#414487","#3c4f8a","#375a8c",
                                  "#32648e","#2a788e","#26828e",
                                  "#228b8d","#1f958b","#22a884",
                                  "#2cb17e","#3bbb75","#4ec36b",
                                  "#7ad151","#95d840","#b0dd2f","#cae11f",
                                  "#fde725")
