output$download_fitting_params_table <- downloadHandler(filename = function() {
  paste0("Fitting_Parameters_PhotoMol_",Sys.Date(),".csv")
  },content = function(file) {
    
    parameter <- c("Bin width","Min. observed mass","Upper limit for the standard deviation",
                   "Tolerance to the initial guesses","Baseline",
                   "Window range","Starting values1","Starting values2",
                   "Starting values3","Starting values4","Starting values5",
                   "Starting values6","Starting values7","Starting values8")
    
    value    <-  c(
      input$bin_width,input$min_observed_mass,input$upper_limit_std,
      input$position_tolerance,input$baseline,
      paste0(input$window_range[1],"-",input$window_range[2]),
      input$starting_values1,input$starting_values2,input$starting_values3,input$starting_values4,
      input$starting_values5,input$starting_values6,input$starting_values7,input$starting_values8)
    
    df <- data.frame(parameter,value)
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_params_table <- downloadHandler(filename = function() {
  paste0("Fitted_Parameters_PhotoMol_",Sys.Date(),".csv")
  },content = function(file) {

  req(photoMolModels$allMassesLoaded)

  fit_tables <- list()

  i     <- 0
  cnt1  <- 0

  for (model in photoMolModels$models) {

    i <- i + 1

    table  <- as.data.frame(model$fit_table)

    if (length(photoMolModels$models) >1) table$File <- model$name

    if (nrow(table) > 0) {
      cnt1 <- cnt1 + 1
      fit_tables[[cnt1]] <- table
    }

  }

  if (length(fit_tables) > 0) {
    fit_table <- do.call(rbind,fit_tables)[,-5]
  } else {
    fit_table <- NULL
  }
    
  write.csv(fit_table,file,row.names = F,quote = F)
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
    table[,1] <- table[,1]*csfFactorForContrast
    table[,2] <- table[,2]*csfFactorForContrast
    table[,1] <- paste0(signif(table[,1],2), " / 1e3")
    table[,2] <- paste0(signif(table[,2],2), " / 1e3")
    
    df <- table[,1:(ncol(table)-1)]
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_mass_histogram <- downloadHandler(filename = function() {
  paste0("Mass_histogram_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {
  
  req(photoMolModels$allMassesLoaded)

  dfs <- list()

  i <- 0
  for (refeyn in photoMolModels$models) {

    i <- i + 1

    counts <- refeyn$hist_counts
    mass   <- refeyn$hist_mass

    df <- data.frame("Histogram counts"=counts,"Mass (kDa)"=mass,"File"=refeyn$name)

    dfs[[i]] <- df
  }

  df <- do.call(rbind,dfs)

  write.csv(df,file,row.names = F,quote = F)

})


output$download_fitted_gaussians <- downloadHandler(filename = function() {
  paste0("Fitted_gaussians_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {

  req(photoMolModels$allMassesLoaded)

  legendsAll             <- get_legend_from_rhandTable(input$legendInfo)

  id_start  <- 1

  refeynFits <- list()

  for (refeyn in photoMolModels$models) {

    if (!is.null(refeyn$fit)) {

      refeynFit  <- as.data.frame(refeyn$fit)

      id_end <- ifelse(ncol(refeynFit) <= 3,id_start, id_start + ncol(refeynFit) - 2)
      legends      <- legendsAll[id_start:id_end]
      id_start     <- id_end + 1

      if (length(legends) == 1) {
        refeynFit           <- refeynFit[,-ncol(refeynFit)]
        colnames(refeynFit) <- c("Mass kDa",legends)
      } else {
        colnames(refeynFit) <- c("Mass kDa",legends[2:length(legends)],legends[1])
      }

      colnames(refeynFit) <- paste0(colnames(refeynFit),' ',refeyn$name)

      refeynFits[[length(refeynFits)+1]] <- refeynFit

    }

  }

  max_rows <- max(sapply(refeynFits,nrow))

  for (i in 1:length(refeynFits)) {
    refeynFit <- refeynFits[[i]]
    if (nrow(refeynFit) < max_rows) {

      # Number of rows to add
      N <- max_rows - nrow(refeynFit)

      # Create a dataframe with N rows of NA values
      na_rows <- data.frame(matrix(NA, nrow = N, ncol = ncol(refeynFit)))

      # Set column names to match the original dataframe
      colnames(na_rows) <- colnames(refeynFit)

      # Add the NA rows to the original dataframe
      refeynFit <- rbind(refeynFit, na_rows)

      refeynFits[[i]] <- refeynFit
    }
  }

  refeynFit   <- do.call(cbind, refeynFits)

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
  
  req(photoMolModels$allMassesLoaded)

  dfs <- list()

  i <- 0
  for (refeyn in photoMolModels$models) {

    i <- i + 1

    counts <- refeyn$hist_counts
    mass   <- refeyn$hist_mass
    counts <- counts / sum(counts)

    df <- data.frame("Normalized histogram counts"=counts,
                     "Mass (kDa)"=mass,"File"=refeyn$name)

    dfs[[i]] <- df
  }

  df <- do.call(rbind,dfs)

  write.csv(df,file,row.names = F,quote = F)

})

output$download_fitted_gaussiansNormalised <- downloadHandler(filename = function() {
  paste0("Fitted_gaussiansNormalised_PhotoMol_",Sys.Date(),".csv")
},content = function(file) {

  req(photoMolModels$allMassesLoaded)

  legendsAll             <- get_legend_from_rhandTable(input$legendInfo)

  id_start  <- 1

  refeynFits <- list()

  for (refeyn in photoMolModels$models) {

    if (!is.null(refeyn$fit)) {

        # Retrieve the mass data that was used for the fitting
      dfMass <- data.frame("mass"=refeyn$masses_kDa)
      dfMass <- dfMass %>%
        filter(mass >= refeyn$hist_window[1]) %>%
        filter(mass <= refeyn$hist_window[2])

      refeynFit  <- as.data.frame(refeyn$fit)

      id_end <- ifelse(ncol(refeynFit) <= 3,id_start, id_start + ncol(refeynFit) - 2)
      legends      <- legendsAll[id_start:id_end]
      id_start     <- id_end + 1

      if (length(legends) == 1) {
        refeynFit           <- refeynFit[,-ncol(refeynFit)]
        colnames(refeynFit) <- c("Mass kDa",legends)
      } else {
        colnames(refeynFit) <- c("Mass kDa",legends[2:length(legends)],legends[1])
      }

      colnames(refeynFit) <- paste0(colnames(refeynFit),' ',refeyn$name)

      refeynFit[,2:ncol(refeynFit)] <- refeynFit[,2:ncol(refeynFit)] / nrow(dfMass)

      refeynFits[[length(refeynFits)+1]] <- refeynFit

    }

  }

  max_rows <- max(sapply(refeynFits,nrow))

  for (i in 1:length(refeynFits)) {
    refeynFit <- refeynFits[[i]]
    if (nrow(refeynFit) < max_rows) {

      # Number of rows to add
      N <- max_rows - nrow(refeynFit)

      # Create a dataframe with N rows of NA values
      na_rows <- data.frame(matrix(NA, nrow = N, ncol = ncol(refeynFit)))

      # Set column names to match the original dataframe
      colnames(na_rows) <- colnames(refeynFit)

      # Add the NA rows to the original dataframe
      refeynFit <- rbind(refeynFit, na_rows)

      refeynFits[[i]] <- refeynFit
    }
  }

  refeynFit   <- do.call(cbind, refeynFits)

  write.csv(refeynFit,file,row.names = F,quote = F)

})

