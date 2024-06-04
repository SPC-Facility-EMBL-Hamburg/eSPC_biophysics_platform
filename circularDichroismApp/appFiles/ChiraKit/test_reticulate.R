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

t1 = cd_experiment_comparison()

t1$load_data('/home/os/Downloads/zero_new.csv','t')

t1$signalDesiredUnit <- t1$signalInput[,1:6]
t1$labels <- np_array(as.character(c(1,1,1,2,2,2)))
  
t1$labels 

t1$summarise_signal_per_label()
t1$generate_comparison_labels()

library(tidyverse)

means <- t1$means

plot_list   <- list()

fig <- plot_ly()

for (i in 1:(t1$labels_unique_N-1)) {
  
  for (ii in (i+1):t1$labels_unique_N) {
    
    print(i)
    print(ii)
    tempDf <- data.frame(x=means[,i],y=means[,ii])
    
    print(tempDf)
    
    fig    <- fig %>% add_trace(data=tempDf,x=~x,y=~y,
                                type = 'scatter', 
                                marker = list(size = 9,
                                              line = list(width = 0)))
    
  }
    
}

fig

yRange <- find_y_plotting_range(means,sds)

x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
          tickfont = list(size = axis_size),range = c(minWL <- min(wl) - 5, max(wl) + 5),
          showgrid = F)

y <- list(title = workingUnits2ProperLabel(selectedUnits),
          titlefont = list(size = axis_size), tickfont = list(size = axis_size),
          showgrid = F, range = c(yRange$yMin,yRange$yMax))

fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                      legend = list(font = list(size = axis_size-1)))







