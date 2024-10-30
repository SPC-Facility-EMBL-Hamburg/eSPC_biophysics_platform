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

setwd('/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/')

source("server_files/helpers.R")


source_python("python_src/cdUnitsConverter.py")
source_python("python_src/helpers.py")
source_python('python_src/loadCDfilesHelpers.py')

source_python("python_src/cdAnalyzer.py")

f <- read.csv('secondary_structure_estimation_files/AU-A128_PCDDB-Nov22.txt',sep='',header=FALSE)
head(f)
dim(f)

lbl <- read_file('secondary_structure_estimation_files/Labels-SMP180_PCDDBOrder.txt')

lbl <- strsplit(lbl,split = '\t')[[1]]

length(lbl)

colnames(f) <- lbl

f$wl <- 240:175

f2 <- reshape2::melt(f,id.vars = 'wl')


fig <- ggplot(f2[f2$variable == "FepA",],aes(wl,value,color=variable)) +
  geom_point() 

ggplotly(fig)

f3 <- (f2[f2$variable == "FepA",c(1,3)]) 
  
colnames(f3) <- c('wavelength','signal')
f3 <- f3[f3$signal != 0,]


write.csv(f3,'/home/os/Downloads/FepA.csv',row.names = FALSE)





