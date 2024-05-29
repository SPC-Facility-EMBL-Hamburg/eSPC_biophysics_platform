# load libraries
library(reshape2)
library(tidyverse)
library(reticulate)
library(minpack.lm)
library(nlstools)
library(broom)
library(data.table)
library(colourpicker)
library(RColorBrewer)
## Set the working directory according to the user name

user      <- Sys.info()['user']
users_dir <- paste0("/home/",user,"/data_users/")

notebook_app  <- (Sys.info()["nodename"] == "osvaldo")

if (notebook_app) {
  use_python("/home/osvaldo/miniconda3/bin/python")
} else  {
  reticulate::use_condaenv("r-reticulate",required = TRUE)
} 

# Developer path
base_dir <- paste0('/home/',Sys.info()["user"],'/spc_shiny_servers/thermoA_app', '/thermoAffinity/')

# Set the correct path for the docker user
if (user == 'shiny') {
  base_dir <- "/home/shiny/thermoAffinity/"
}


# To convert between units 
factorList <- list("molar"=1,"milimolar"=1e3,"micromolar"=1e6,'nanomolar'=1e9)

# viridis color palette for numeric gradient
global_colors_palette_signal <- c("#440154","#471063","#481d6f","#472a7a",
                                  "#414487","#3c4f8a","#375a8c",
                                  "#32648e","#2a788e","#26828e",
                                  "#228b8d","#1f958b","#22a884",
                                  "#2cb17e","#3bbb75","#4ec36b",
                                  "#7ad151","#95d840","#b0dd2f","#cae11f",
                                  "#fde725")

#Custom render function - used to colour the cell based on the cell value (hex code of a color)
myrenderer <- "function(instance, td, row, col, prop, value, cellProperties) {
                Handsontable.renderers.TextRenderer.apply(this, arguments);
                if (instance.params) {
                    hcols = instance.params.col_highlight
                    hcols = hcols instanceof Array ? hcols : [hcols]
                    hrows = instance.params.row_highlight
                    hrows = hrows instanceof Array ? hrows : [hrows]
                    
                    for (i = 0; i < hcols.length; i++) { 
                        if (hcols[i] == col && hrows[i] == row) {
                            td.style.background = instance.getDataAtCell(row, col);
                        }
                    }
                }
  }"     

#Custom render function - used to apply a checkbox style to a logical column
myrendererBoolean <- "function(instance, td, row, col, prop, value, cellProperties) {
            Handsontable.renderers.CheckboxRenderer.apply(this, arguments);}"
