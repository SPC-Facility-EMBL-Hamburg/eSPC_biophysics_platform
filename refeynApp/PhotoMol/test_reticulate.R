rm(list=ls())
gc()
setwd("/home/osvaldo/spc_shiny_servers/refeynApp/PhotoMol/")
library(reticulate)
library(tidyverse)
library(plotly)
library(reshape2)
source("server_files/helpers.R")
source_python("refeyn.py")
source_python("helpers.py")
source("server_files/plot_functions.R")

refeyn <- Refeyn()
refeyn$load_data_h5("./www/demo.h5")
refeyn$create_histo()
refeyn$hist_mass
refeyn$hist_counts
refeyn$fit_histo()

refeynFit           <- as.data.frame(refeyn$fit)

baseline <- 0

for (i in 1:4) {
  
  tempV       <- refeynFit[,i+1]
  sel         <- tempV > (baseline+0.05)
  tempV       <- tempV / sum(sel)
  tempV[!sel] <- "NA"
  
  refeynFit[,i+1] <- as.numeric(tempV)
  
}

ggplot(refeynFit,aes(V1,V5))+
  geom_line()

m           <- read.csv("~/Mass_histogramNormalised_PhotoMol_2023-02-20.csv")
colnames(m) <- c("y","x")

ggplot(m,aes(x,y))+
  geom_bar(stat = "identity")

g <- read.csv("~/Fitted_gaussiansNormalised_PhotoMol_2023-02-20.csv")
head(g)

ggplot(m,aes(x,y))+
  geom_bar(stat = "identity")+
  geom_line(data = g,aes(Mass.kDa,Molecule..1 ))


