rm(list=ls())
gc()

packages <- c("shinydashboard","shinycssloaders","rhandsontable","plotly","shinyalert","reticulate",
              "DT","reshape2","ggridges","scales",'tidyverse','pracma','ggpubr')

invisible(lapply(packages, library, character.only = TRUE))

#reticulate::use_condaenv("r-reticulate",required = TRUE)
reticulate::use_python("/home/osvaldo/miniconda3/bin/python")

setwd("/home/osvaldo/spc_shiny_servers/dynamicLightScatteringApp/appFiles/Raynals/")
#setwd("/home/osvaldo/spc_shiny_servers/dynamicLightScatteringApp/appFiles/Raynals/")

source("server_files/helpers.R")
source_python("./dlsAnalyzer.py")
source_python("./helpers.py")

expName <- 'CA_5_replicates'

dlsAnalyzer  <- dlsAnalyzer()
dlsAnalyzer$loadExperiment("www/220315CA_5_replicates.csv",expName)

dlsAnalyzer$experimentsOri[[expName]]$nMeasurements
dlsAnalyzer$experimentsOri[[expName]]$getQ()
dlsAnalyzer$experimentsOri[[expName]]$createFittingS_space(0.009,1e6,200) # in nm min and max Hr
dlsAnalyzer$experimentsOri[[expName]]$setAutocorrelationData()
dlsAnalyzer$experimentsOri[[expName]]$getBetaEstimate()
dlsAnalyzer$experimentsOri[[expName]]$getG1correlation()
dlsAnalyzer$experimentsOri[[expName]]$getInitialEstimates(0.01)

dlsAnalyzer$experimentsOri[[expName]]$predictAutocorrelationCurves()

estimatedContributionsFirst <- dlsAnalyzer$getExperimentProperties('contributionsGuess')[[1]]

estimatedContributions2 <- do.call(cbind,estimatedContributionsFirst)
data           <- data.frame(dlsAnalyzer$experimentsOri[[1]]$ds,estimatedContributions2)

head(data)

colnames(data) <- c('diff',paste0('C',1:25))

data <- reshape2::melt(data,id.vars='diff')

data <- data %>% filter(variable %in% c('C1','C2','C3','C4'))
str(data)

overlapFactor          <- 1
estimatedContributions <- data

head(estimatedContributions)

## Mapping
mapping <- data.frame('Condition'=unique(estimatedContributions$variable),'Group'=c('A','A','B','B'))

source('server_files/plot_functions.R')

fig <- plot_type_hr_distribution('histogram',
                                 estimatedContributions,
                                 16,
                                 0.5,
                                 'default',
                                 T,
                                 F,
                                 mapping
)

fig
#### 
expName <- 'testt'

dlsAnalyzer  <- dlsAnalyzer()
dlsAnalyzer$loadExperiment("www/test.csv",expName)

#file = "/home/os/arise/DLS_manuscript/analysisRaynals/SimulatedDataDLS_case2/autocorrelation_2023-02-08.csv"
#dlsAnalyzer$loadExperiment(file,expName)
dlsAnalyzer$experimentsOri[[expName]]$setAutocorrelationData()

dlsAnalyzer$experimentsOri[[expName]]$nMeasurements
dlsAnalyzer$experimentsOri[[expName]]$getQ()
dlsAnalyzer$experimentsOri[[expName]]$createFittingS_space(0.009,1e6,200) # in nm min and max Hr
dlsAnalyzer$experimentsOri[[expName]]$getBetaEstimate()
dlsAnalyzer$experimentsOri[[expName]]$getG1correlation()
alphaVec <- (5**seq(-6,2,0.1))**2
dlsAnalyzer$experimentsOri[[expName]]$getInitialEstimatesManyAlpha(alphaVec,1e8)
dlsAnalyzer$experimentsOri[[expName]]$getOptimalAlphaLcurve()
dlsAnalyzer$experimentsOri[[expName]]$getInitialEstimatesOptimalAlphaLcurve(1e8)

fidelity <- dlsAnalyzer$experimentsOri[[expName]]$curvesResidualNorm # One column per curve
penalty  <- dlsAnalyzer$experimentsOri[[expName]]$curvesPenaltyNorm  # One column per curve

idx <- T

expN   <- dlsAnalyzer$experimentNames[idx]
sNames <- dlsAnalyzer$getExperimentProperties("sampleInfoRelevant")[idx]
rn     <- dlsAnalyzer$getExperimentProperties("curvesResidualNorm")[idx]
pn     <- dlsAnalyzer$getExperimentProperties("curvesPenaltyNorm")[idx]
idSel  <- dlsAnalyzer$getExperimentProperties("alphaOptIdx")[idx]

lCurveData  <- formatNorms(expN,sNames,rn,pn,idSel)

ggplot(df,aes(fidelity,penalty))+
  geom_point()+
  geom_point(data=dfSel,
             shape=4,color='blue',size=8) +
  scale_x_log10()+
  scale_y_log10()+
  xlab('Fidelity term')+
  ylab('Penalty term')+
  facet_wrap(~variable)+
  theme_bw(base_size = 12)







