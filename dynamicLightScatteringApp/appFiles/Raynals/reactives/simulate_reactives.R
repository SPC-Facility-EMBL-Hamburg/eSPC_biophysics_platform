output$simulation_params_table <- renderRHandsontable({
  
  popMeans    <- c("5","10","20","5 50","10 200")
  popSDs      <- c("0.1","0.1","0.1","0.1 1","0.1 5")
  particleNum <- c("1e2","1e2","1e2","1e4 10","1e4 5") 
    
  nReps <- ceiling(input$SimNumber / 5)
  
  table <- data.frame(
    "Simulation name"           = paste0("Sim. ",1:(5*nReps)),
    "Population mean(s) (nm)"   = rep(popMeans,nReps),
    "Population sd(s) (nm)"     = rep(popSDs,nReps),
    "Number of particles"       = rep(particleNum,nReps),
    "Wavelength (nm)"           = 817,
    "Angle (º)"                 = 150,
    "Refractive index"          = 1.33,
    "Temperature (°C)"          = 20,
    "Viscosity (pascal-second * 1e4)" = 8.9,
    "Show"                      = as.logical(T),
    check.names = F)
  
  table <- table[1:input$SimNumber,]
  
  return(rhandsontable(table,rowHeaders = F,stretchH = "all") %>% 
           hot_col(c(9),format = "0[.]00") %>% 
           hot_col(c("Wavelength (nm)","Angle (º)","Temperature (°C)"),
                   format = 1))
})

simulatedDistributions <- eventReactive(input$launchSimulation,{
  
  simParams        <- hot_to_r(input$simulation_params_table) 
  allDistributions <- list()
  
  #To store the info of the selected distributions
  simNames  <- simParams[,"Simulation name"][simParams$Show] 
  simAngles <- simParams[,"Angle (º)"][simParams$Show]         
  simWL     <- simParams[,"Wavelength (nm)"][simParams$Show]
  simTemp   <- simParams[,"Temperature (°C)"][simParams$Show]
  simRefInd <- simParams[,"Refractive index"][simParams$Show]
  simVisco  <- simParams[,"Viscosity (pascal-second * 1e4)"][simParams$Show]
  
  # Create table to export (with the used values)!
  reactives$simulatedMetaData               <- simParams[simParams$Show,]
  reactives$simulatedMetaData$GaussianError <- input$autocorrelationError
  reactives$simulatedMetaData$Intercept     <- input$autocorrelationIntercept
  
  counter <- 0
  # Iterate the different distributions
  for (i in 1:nrow(simParams)) {
    
    # Simulate if selected
    if (simParams[i,"Show"]) {
      
      counter <- counter + 1
      
      # These user inputs are strings and values are separated by spaces
      means       <- get_positions(simParams[i,"Population mean(s) (nm)"])
      sdeviations <- get_positions(simParams[i,"Population sd(s) (nm)"])
      sampleSizes <- get_positions(simParams[i,"Number of particles"])
      
      if (sum(sampleSizes) > 1e7) {
        shinyalert(text = "<b>Total sample size can't be greater than 1e7</b>",
                   type = 'warning',html = T)
        Sys.sleep(1)
        return(NULL)
      }
      
      vectorSizes <- sapply(list(means,sdeviations,sampleSizes), length)
      
      if (length(unique(vectorSizes)) != 1) {
        shinyalert(text = paste0("<b>The number of means, standard deviations and 
                             sample sizes are not equal = ",
                                 paste0(vectorSizes,collapse=" "),'</b>'),
                   type = 'warning',html = T)
        Sys.sleep(1)
        return(NULL)
      }
      
      # Generate a sample of hydrodynamic radius using gaussian distributions
      particles <- data.frame(
        'hr'=generateParticlesPopulation(means,sdeviations,sampleSizes))
      
      distributions <- generate_distributions(
        particles,simParams[i,"Angle (º)"],
        simParams[i,"Wavelength (nm)"],simParams[i,"Refractive index"])
      
      distributions$variable <- simParams[i,"Simulation name"]
      
      allDistributions[[counter]] <- data.frame(distributions)
      
    }
    
  }

  allDistributions <- do.call(rbind,allDistributions)
  
  return(list("allDistributions"=allDistributions,"simNames"=simNames,
              "simAngles"=simAngles,"simWL"=simWL,
              "simTemp"=simTemp,"simRefInd"=simRefInd,
              "simViscosity"=simVisco))
  
})

