rm(list=ls())
gc()
library(reticulate)
library(reshape2)
library(tidyverse)
library(plotly)
library(factoextra)

use_condaenv("r-reticulate",required = TRUE)

setwd('/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/')

source("server_files/helpers.R")

source_python("./cdUnitsConverter.py")
source_python("./helpers.py")
source_python('./loadCDfilesHelpers.py')

source_python("./cdAnalyzer.py")

cd = cd_experiment_general()

cd$load_data('/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/www/example_files/default_spectra_GQuadruplex.csv','')
cd$spectraNames <- trimws(cd$spectraNames)

# Use new params
useNewParams   <- T
useNewParams2  <- F

if (useNewParams) {
  parametersSec <- read.csv('./www/example_files/new_fr.csv',
                            check.names = FALSE)
  
  parametersSec   <- parametersSec[rowSums(parametersSec[,2:6]) > 0,]
  
  cd$signalInput  <- cd$signalInput[,cd$spectraNames %in% parametersSec[,1]]
  cd$spectraNames <- cd$spectraNames[cd$spectraNames %in% parametersSec[,1]]
  
  parametersSec <- parametersSec[,2:7] 
  parametersSec <- sweep(parametersSec, 1, parametersSec$tbs, "/")
  parameters <- parametersSec %>% select(-tbs) 
  
} else if (useNewParams2) {
  
  parametersSec <- read.csv('./www/example_files/new_fr_loop.csv',
                            check.names = FALSE)[,1:9]
  
  colnames(parametersSec) <- c('id','aa','sa','as','dl','ss','o','tbs','dl2')
  
  parametersSec$o <- parametersSec$o - parametersSec$dl2 + parametersSec$dl
  parametersSec$dl <- parametersSec$dl2 
  
  parametersSec <- parametersSec %>% select(-dl2)
  
  parametersSec   <- parametersSec[rowSums(parametersSec[,2:7]) > 0,]
  
  cd$signalInput  <- cd$signalInput[,cd$spectraNames %in% parametersSec[,1]]
  cd$spectraNames <- cd$spectraNames[cd$spectraNames %in% parametersSec[,1]]
    
  parametersSec   <- parametersSec[,2:8] 
  parametersSec$o <- parametersSec$ss + parametersSec$o
  parametersSec   <- sweep(parametersSec, 1, parametersSec$tbs, "/")
  parameters      <- parametersSec %>% select(-tbs,-ss) 
  
} else {
  # Replace with original params
  parametersSec <- read.csv('./www/example_files/default_secondary_parameters.csv',
                            check.names = FALSE)
  
  colnames(parametersSec) <- c('id','aa','sa','as','dl','o')
  
  parameters        <- parametersSec[,-1]
}

head(parameters)

cd_reference           <- base::t(cd$signalInput)
colnames(cd_reference) <- cd$wavelength

source("server_files/helpers_G-Quadruplex.R")

all_dfs <- list()
rmsds   <- c()

for (i in 1:nrow(parameters)) {
  test_spectra            <- matrix(cd_reference[i,],nrow=1)
  
  # Leave-one-out method
  cd_reference2           <- cd_reference[-i,]
  parameters2             <- parameters[-i,]
  rownames(cd_reference2) <- cd$spectraNames[-i]

  res1          <- svd_spectra_analysis(cd_reference2,parameters2,test_spectra,cd$spectraNames[i])
  pred          <- (res1$predicted/100)
  
  real          <- parameters[i,]
  rownames(real)<- cd$spectraNames[i]
  
  colnames(real) <- paste0(colnames(real),'_real')
  
  rmsds <- c(rmsds,sqrt(sum((real - pred)**2)/ncol(real)))
  
  # Convert row names to a column
  df    <- cbind(real,pred)
  all_dfs[[i]] <- df

}

all_dfs <- do.call(rbind,all_dfs)

p1 <- ggplot(all_dfs,aes(aa,aa_real))+
  geom_point()+
  theme_classic(base_size = 14)+
  geom_abline(slope = 1,intercept = 0)+
  geom_smooth(method = 'lm')+
  xlim(c(0,0.7))+ylim(c(0,0.7))

p2 <- ggplot(all_dfs,aes(sa,sa_real))+
  geom_point()+
  theme_classic(base_size = 14)+
  geom_abline(slope = 1,intercept = 0)+
  xlim(c(0,0.7))+ylim(c(0,0.7))+
  geom_smooth(method = 'lm')

p3 <- ggplot(all_dfs,aes(as,as_real))+
  geom_point()+
  theme_classic(base_size = 14)+
  geom_abline(slope = 1,intercept = 0)+
  xlim(c(0,0.7))+ylim(c(0,0.7))+
  geom_smooth(method = 'lm')

p4 <- ggplot(all_dfs,aes(dl,dl_real))+
  geom_point()+
  theme_classic(base_size = 14)+
  geom_abline(slope = 1,intercept = 0)+
  xlim(c(0,0.7))+ylim(c(0,0.7))+
  geom_smooth(method = 'lm')

p5 <- ggplot(all_dfs,aes(o,o_real))+
  geom_point()+
  theme_classic(base_size = 14)+
  geom_abline(slope = 1,intercept = 0)+
  xlim(c(0,0.7))+ylim(c(0,0.7))+
  geom_smooth(method = 'lm')


p1
p2
p3
p4
p5

library(ggpubr)

ggarrange(p1, p2, p3, p4, p5, 
          ncol = 3, nrow = 2)


df       <- data.frame(rmsds,'name'=cd$spectraNames)
df$rmsds <- signif(df$rmsds,3)
df

(lm(aa_real ~ aa, data = all_dfs))
(lm(sa_real ~ sa, data = all_dfs))
(lm(as_real ~ as, data = all_dfs))
(lm(dl_real ~ dl, data = all_dfs))
(lm(o_real  ~ o,  data = all_dfs))

cor(all_dfs$aa_real,all_dfs$aa)
cor(all_dfs$sa_real,all_dfs$sa)
cor(all_dfs$as_real,all_dfs$as)
cor(all_dfs$dl_real,all_dfs$dl)
cor(all_dfs$o_real,all_dfs$o)









