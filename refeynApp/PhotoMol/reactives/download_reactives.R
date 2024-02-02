output$download_fitting_params_table <- downloadHandler(filename = function() {
  paste0("Fitting_Parameters_PhotoMol_",Sys.Date(),".csv")
  },content = function(file) {
    
    parameter <- c("Bin width","Min. observed mass","Upper limit for the standard deviation",
                   "Tolerance to the initial guesses","Starting values","Baseline",
                   "Window range")
    
    value    <-  c(input$bin_width,input$min_observed_mass,input$upper_limit_std,
                   input$position_tolerance,input$starting_values,input$baseline,
                   paste0(input$window_range[1],"-",input$window_range[2]))
    
    df <- data.frame(parameter,value)
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_params_table <- downloadHandler(filename = function() {
  paste0("Fitted_Parameters_PhotoMol_",Sys.Date(),".csv")
  },content = function(file) {
    
    req(refeyn$massesLoaded)
    table <- refeyn$fit_table
    df    <- table[,1:(ncol(table)-1)]
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_fitting_params_table_calibration <- downloadHandler(filename = function() {
  paste0("Fitting_Parameters_Calibration_PhotoMol_",Sys.Date(),".csv")
  },content = function(file) {
    
    parameter <- c("Bin width","Initial guesses * 1e3","Known masses",
                   "Window range * 1e3")
    
    value    <-  c(input$bin_widthContrast,input$starting_valuesContrast,input$knownMasses,
                   paste0(input$window_rangeContrast[1],"-",input$window_rangeContrast[2]))
    
    df <- data.frame(parameter,value)
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_params_table_calibration <- downloadHandler(filename = function() {
  paste0("Fitted_Parameters_Calibration_PhotoMol_",Sys.Date(),".csv")
  },content = function(file) {
    
    if (is.null(modify_refeynCalibration_data())) {return(NULL)}
    
    table <- refeynCalib$fit_table
    table[,1] <- table[,1]*factorForContrast
    table[,2] <- table[,2]*factorForContrast
    table[,1] <- paste0(signif(table[,1],2), " / 1e3")
    table[,2] <- paste0(signif(table[,2],2), " / 1e3")
    
    df <- table[,1:(ncol(table)-1)]
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_mass_histogram <- downloadHandler(filename = function() {
  paste0("Mass_histogram_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {
  
  req(refeyn$massesLoaded)

  counts <- refeyn$hist_counts
  mass   <- refeyn$hist_mass 
  
  df <- data.frame("Histogram counts"=counts,"Mass (kDa)"=mass)
  
  write.csv(df,file,row.names = F,quote = F)
})

output$download_fitted_gaussians <- downloadHandler(filename = function() {
  paste0("Fitted_gaussians_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {
  
  req(refeyn$massesLoaded)
  refeynFit           <- as.data.frame(refeyn$fit)
  legends             <- get_legend_from_rhandTable(input$legendInfo)
  
  if (length(legends) == 1) {
    colnames(refeynFit) <- c("Mass kDa",legends)
  } else {
    colnames(refeynFit) <- c("Mass kDa",legends[2:length(legends)],legends[1])
  }
  
  write.csv(refeynFit,file,row.names = F,quote = F)
})

output$downloadH5file <- downloadHandler(
  
  filename <- function() {
    paste0('withMass-',input$massPhotometryFile$name)
  },
  
  content <- function(file) {
    refeyn$export_h5_dataset('output.h5')
    file.copy("output.h5", file)
  }
)

output$download_mass_histogramNormalised <- downloadHandler(filename = function() {
  paste0("Mass_histogramNormalised_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {
  
  req(refeyn$massesLoaded)
  
  counts <- refeyn$hist_counts
  mass   <- refeyn$hist_mass 
  counts <- counts / sum(counts)
  
  df <- data.frame("Normalised histogram counts"=counts,"Mass (kDa)"=mass)
  
  write.csv(df,file,row.names = F,quote = F)
})

output$download_fitted_gaussiansNormalised <- downloadHandler(filename = function() {
  paste0("Fitted_gaussiansNormalised_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {
  
  req(refeyn$massesLoaded)
  
  # Retrieve the mass data that was used for the fitting
  dfMass <- data.frame("mass"=refeyn$masses_kDa)
  dfMass <- dfMass %>% 
    filter(mass >= refeyn$hist_window[1]) %>% 
    filter(mass <= refeyn$hist_window[2])
  
  # Retrieve the fitted gaussians 
  refeynFit           <- as.data.frame(refeyn$fit)
  legends             <- get_legend_from_rhandTable(input$legendInfo)
  
  baseline <- input$baseline / nrow(dfMass) # Normalise the baseline
  
  if (length(legends) == 1) {
    colnames(refeynFit) <- c("Mass kDa",legends)
  } else {
    colnames(refeynFit) <- c("Mass kDa",legends[2:length(legends)],legends[1])
  }
  
  # Normalise data
  for (i in 1:length(legends)) {
    
    tempV       <- refeynFit[,i+1]
    sel         <- tempV > (baseline+0.05)
    # Same division as applied to the histogram dataset!
    tempV       <- tempV / nrow(dfMass)
    # Set NA because we don't want to plot the gaussian when it's almost flat
    tempV[!sel] <- "NA"
    
    refeynFit[,i+1] <- tempV
    
  }
  
  write.csv(refeynFit,file,row.names = F,quote = F)
})

