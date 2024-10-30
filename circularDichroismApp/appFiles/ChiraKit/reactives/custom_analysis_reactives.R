# Create the Table to fill with the experimental parameters data
observeEvent(list(input$legendInfo,input$workingUnits,input$numberOfExpParams),{
  
  req(reactives$data_loaded)
  # Retrieve the dataframe from the legends table (1. Load Input Tab)
  legendDf <- getLegendDF(input$legendInfo)
  id       <- legendDf$Internal.ID
  
  # Initialize the dataframe with the CD curves 
  df           <- data.frame(id)
  colnames(df) <- c('CD_curve')
  
  # accept no more than two experimental parameters
  nParams <- min(input$numberOfExpParams,2)
  
  # abort if we don't have one or experimental parameters
  if (!nParams %in% c(1,2)) return(NULL)
  
  paramsNames <- c(input$firstParamName,input$secondParamName)
  
  # add one column per experimental parameter. Only one is accepted for now
  for (i in 1:nParams) {
    df[,paramsNames[i]] <- 0
  }
  
  df$Dataset_name <- 'D1'
  
  # Remove non-selected CD curves
  df <- df[legendDf$Show,]
  
  # Remove experiments with non-matching units
  id_to_keep         <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
  internalID_all     <- cdAnalyzer$get_experiment_properties('internalID')
  internalID_to_keep <- unlist(internalID_all[id_to_keep])
  
  df <- df[df$CD_curve %in% internalID_to_keep,]
  
  # Assign the created dataframe to the Table custom_parameters_data (available at the 2d. Custom analysis Tab)
  output$custom_parameters_data <- renderRHandsontable({
    rhandsontable(df,rowHeaders=NULL)    %>% 
      hot_col(col = c(1), readOnly=TRUE) %>% 
      hot_table(stretchH='all')
  })
  
})

# Update the experimental parameters name
observeEvent(list(input$firstParamName,input$secondParamName),{
  
  req(reactives$data_loaded)
  req(input$custom_parameters_data)
  
  # Retrieve the custom dataframe from the experimental parameters Table
  customDF <- hot_to_r(input$custom_parameters_data)
  
  current_colnames <- colnames(customDF)
  
  if (is.na(input$firstParamName) || input$firstParamName == '') {
    return(NULL)
  }
  
  colnames(customDF)[2] <- substr(input$firstParamName, 1, 20) # Use the first 20 characters
  
  # not used for now...
  if (ncol(customDF) == 4) colnames(customDF)[3] <- substr(input$secondParamName, 1, 20) # Use the first 20 characters
  
  # Assign the created dataframe to the Table custom_parameters_data (available at the 2d. Custom analysis Tab)
  output$custom_parameters_data <- renderRHandsontable({
    rhandsontable(customDF,rowHeaders=NULL)    %>% 
      hot_col(col = c(1), readOnly=TRUE) %>% 
      hot_table(stretchH='all')
  })
  
})

observeEvent(input$btn_create_custom_dataset,{
  
  req(reactives$data_loaded)
  
  reactives$customDatasetCreated          <- NULL
  reactives$spectra_was_decomposed_custom <- FALSE
  
  cdAnalyzer$clean_experiments('custom')
  
  # Retrieve which CD experiments should be used to build
  # a new dataset for the custom model analysis
  df_ids2find        <- hot_to_r(input$custom_parameters_data)
  groups             <- unique(df_ids2find$Dataset_name)
  
  append_record_to_logbook(c('Creating a custom analysis dataset with the following data',df_to_lines(df_ids2find)))
  
  # Create one thermal dataset per group
  for (group in groups) {
    
    df_temp <- df_ids2find[df_ids2find$Dataset_name == group,]
    
    relevantSpectra  <- df_temp$CD_curve 
    relevantF1       <- df_temp[,2]
    
    merged <- get_signal_dfs_from_selected_spectra(relevantSpectra,cdAnalyzer)
    
    sorted_indexes  <- order(relevantF1)
    relevantF1      <- relevantF1[sorted_indexes]
    
    sorted_signal   <- as.matrix(merged[,-1][,sorted_indexes],drop = FALSE)
    
    # Assign the signal and temperature data to the new thermal unfolding experiment
    cdAnalyzer$experimentsCustom[[group]]                    <- CdExperimentCustomAnalysis()
    
    cdObj <- cdAnalyzer$experimentsCustom[[group]]
    
    cdObj$wavelength         <- np_array(merged[,1])
    cdObj$signalDesiredUnit  <- np_array(sorted_signal)
    cdObj$name               <- group
    
    cdObj$first_measurement_dimension  <- np_array(relevantF1)
    
    cdObj$first_exp_param_name         <- colnames(df_temp)[2] 
    
    # cdObj$second_exp_param_name is not used, will be set as 'Dataset_name'
    # we assume that the user will not use 'Dataset_name' as one parameter to be fitted  ... 
    # If you want this code to be used, 
    # change the hidden UI elements in 2d_ui_load_custom_data.R
    cdObj$second_exp_param_name        <- colnames(df_temp)[3] 
        
    # Convert the string of selected wavelengths to a numeric vector
    selected_wl <- parse_selected_wavelengths(input$selected_wavelength_thermal_unfolding)
    
    cdObj$assign_useful_signal(selected_wl)
    
    # add the second measurement dimension, if needed
    if (ncol(df_temp) == 4) {
      
      relevantF2  <- df_temp[,3]
      relevantF2  <- relevantF2[sorted_indexes]
      cdObj$second_measurement_dimension <- np_array(relevantF2)
    } 
    
  }
  
  cdAnalyzer$experimentNamesCustom <- groups
  
  reactives$customWorkingUnits   <- input$workingUnits
  reactives$customDatasetCreated <- TRUE
  
})

observeEvent(input$btn_find_wl_custom,{
  
  req(reactives$customDatasetCreated)
  
  custom_exps            <- cdAnalyzer$experimentNamesCustom
  
  wavelength_filtered_all <- list()
  counter                 <- 0
  
  for (exp in custom_exps) {
    
    counter    <- counter + 1
    relevantF1 <- cdAnalyzer$experimentsCustom[[exp]]$first_measurement_dimension
    cdAnalyzer$experimentsCustom[[exp]]$estimate_useful_signal_based_on_snr_and_amplitude(relevantF1)
    wavelength_filtered_all[[counter]] <- cdAnalyzer$experimentsCustom[[exp]]$wavelength_filtered
    
  }
  
  wavelength_filtered_common <- Reduce(intersect, wavelength_filtered_all)
  if (length(wavelength_filtered_common) == 0) {
    shinyalert(text = "The automatic selection algorithm didn't work. 
               Please select the wavelength manually.",
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    wavelength_filtered_common <- 220
  }
  
  append_record_to_logbook('Running the automatic wavelength selection')
  
  updateTextInput(session, "selected_wavelength_custom", 
                  value = paste(wavelength_filtered_common,collapse=" "))
  
})

observeEvent(list(input$selected_wavelength_custom,input$analysis_model_custom),{
  
  req(reactives$customDatasetCreated)
  
  if (input$analysis_model_custom == 'fixedWL') {
    append_record_to_logbook(paste0("Setting the 'Selected wavelength(s)' to: ",
                                    input$selected_wavelength_custom))
    
    reactives$spectra_decomposition_method_custom <- 'None'
    
  }
  
  reactives$data_loaded <- FALSE
  # Convert the string of selected wavelengths to a numeric vector
  selected_wl <- parse_selected_wavelengths(input$selected_wavelength_custom)
  
  custom_exps <- cdAnalyzer$experimentNamesCustom
  
  for (exp in custom_exps) {
    
    cdAnalyzer$experimentsCustom[[exp]]$assign_useful_signal(selected_wl)
    
  }
  
  reactives$data_loaded <- TRUE
  
})

output$customCurves <- renderPlotly({
  
  req(reactives$data_loaded)
  req(reactives$customDatasetCreated)
  
  df <- generate_custom_df(cdAnalyzer)
  
  plot_unfolding_exp(df,
                     reactives$customWorkingUnits,
                     input$plot_width_custom, input$plot_height_custom, 
                     input$plot_type_custom, input$plot_axis_size_custom,
                     'fixedWL',input$use_log_axis_custom)
  
})

# Get initial parameter estimates
observeEvent(input$btn_get_initial_params,{
  
  req(reactives$customDatasetCreated)
  
  showModal(modalDialog(
    
    tags$h3('Please select the limits of the log spaced grid search:'),
    
    numericInput('leftLimitLogSearch','Lower limit (10**x)',-3,-6,0),
    
    numericInput('rightLimitLogSearch','Higher limit (10**x)',3,0,6),
    
    footer=tagList(
      actionButton('submitLogSearchLimits', 'Submit'),
      modalButton('Cancel')
    )
  ))
  
})
  
renderInitialParamsTable <- function() {
  
  output$initialParamsValuesSVD <- NULL
  output$initialParamsValues    <- NULL
  
  custom_exps <- cdAnalyzer$experimentNamesCustom
  nExps       <- length(custom_exps)
  
  p0s        <- c()
  p0s_names  <- c()
  lowBounds  <- c()
  highBounds <- c()
  
  c1 <- input$analysis_model_custom != 'fixedWL'

  if (c1) req(reactives$spectra_was_decomposed_custom)
  
  for (exp in custom_exps) {
    
    p0s        <- c(p0s,unlist(cdAnalyzer$experimentsCustom[[exp]]$p0))
    
    lowBounds  <- c(lowBounds,unlist(cdAnalyzer$experimentsCustom[[exp]]$low_bound))
    highBounds <- c(highBounds,unlist(cdAnalyzer$experimentsCustom[[exp]]$high_bound))
    
    p0s_names_to_add <- c(
      unlist(cdAnalyzer$experimentsCustom[[exp]]$global_params_names),
      unlist(cdAnalyzer$experimentsCustom[[exp]]$local_params_names_extended)
    )
    
    if (nExps > 1) {
      p0s_names_to_add <- paste0(exp,' ',p0s_names_to_add)
    }
    
    p0s_names <- c(p0s_names,p0s_names_to_add)
    
  }
  
  df_p0  <- data.frame('parameter'   = p0s_names, 'initial_value' = p0s,
                       'lower_limit' = lowBounds,  'high_limit'    = highBounds)
  
  rtable <- rhandsontable(df_p0,rowHeaders=NULL) %>% 
    hot_col('parameter',     readOnly=TRUE) %>% 
    hot_col(c('initial_value','lower_limit','high_limit')) %>% 
    hot_table(stretchH='all')
  
  if (c1) {
    output$initialParamsValuesSVD <- renderRHandsontable({rtable})
  } else {
    output$initialParamsValues <- renderRHandsontable({rtable})
  }
  
  return(NULL)
}

observeEvent(input$submitLogSearchLimits,{
  
  removeModal()
  
  withBusyIndicatorServer("fitCustomHidden",{
    
    custom_exps <- cdAnalyzer$experimentNamesCustom
    nExps       <- length(custom_exps)
    
    if (nExps == 0) return(NULL)
    
    c1 <- input$analysis_model_custom != 'fixedWL'

    if (c1) req(reactives$spectra_was_decomposed_custom)
    
    for (exp in custom_exps) {
      
      # Assign the signal to the change in the PCA/SVD coefficients
      if (c1) {
        cdAnalyzer$experimentsCustom[[exp]]$assign_useful_signal_svd(input$selectedK_custom)
      }
      
      cdAnalyzer$experimentsCustom[[exp]]$generate_model_function(input$custom_model)
      cdAnalyzer$experimentsCustom[[exp]]$get_initial_values(input$leftLimitLogSearch,
                                                             input$rightLimitLogSearch)
      
    }

    renderInitialParamsTable()
    
    Sys.sleep(0.5)

  })
  
  shinyalert(text = paste("<b>The search for the initial values has finished.</b>"),
             type = "success",closeOnEsc = T,closeOnClickOutside = T,
             html=T)
  
})

updateInitialEstimatesOrFittingBoundaries <- function(change) {
  
  # Change has four elements whose  index start at 0: 
  # selected row, selected column, previous value, new value
  
  changesRow <- change[[1]] + 1
  changesCol <- change[[2]] + 1
  
  custom_exps <- cdAnalyzer$experimentNamesCustom
  
  params_counter <- 0
  # find which values was selected by the user and change it in the python class
  for (exp in custom_exps) {
    
    pyObj <- cdAnalyzer$experimentsCustom[[exp]]
    
    n_params <- length(pyObj$global_params_names) + length(pyObj$local_params_names_extended)
    
    params_counter <- params_counter + n_params
    
    if (changesRow <= params_counter) {
      
      selN   <- changesRow - params_counter + n_params
      selExp <- exp
      break
    } 
    
  }
  
  pyObj <- cdAnalyzer$experimentsCustom[[selExp]]
  
  updatedVal <- change[[4]]
  
  #change initial value
  if (changesCol == 2) {
    pyObj$p0[selN] <- updatedVal
    # change low limit (fitting threshold)
  } else if (changesCol == 3) {
    pyObj$low_bound[selN] <- updatedVal
    # change high limit (fitting threshold)
  } else {
    pyObj$high_bound[selN] <- updatedVal
  }
  
  return(NULL)
}

# update the initial values or fitting  boundaries based on the user input
observeEvent(input$initialParamsValues$changes$changes, {
  
  # Get the changes made to the table
  changes <- input$initialParamsValues$changes$changes
  
  for (change in changes) {
    updateInitialEstimatesOrFittingBoundaries(change)
  }
  
})

# update the initial values or fitting  boundaries based on the user input
observeEvent(input$initialParamsValuesSVD$changes$changes, {
  
  # Get the changes made to the table
  changes <- input$initialParamsValuesSVD$changes$changes
  
  for (change in changes) {
    updateInitialEstimatesOrFittingBoundaries(change)
  }
  
})

observeEvent(input$btn_fit_custom_data,{
  
  req(reactives$customDatasetCreated)
  custom_exps <- cdAnalyzer$experimentNamesCustom
  
  withBusyIndicatorServer("fitCustomHidden",{
    
    reactives$custom_data_was_fitted <- FALSE

    if (length(custom_exps) == 0) return(NULL)
    
    for (exp in custom_exps) {
      
      cdAnalyzer$experimentsCustom[[exp]]$fit_signal()
      
      if (!is.null(cdAnalyzer$experimentsCustom[[exp]]$fit_params)) {
        
        shinyalert(text = paste("<b>Fitting done for the experiment ",exp,"</b>"),
                   type = "success",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        
      } else {
        
        shinyalert(text = 
                     paste("<b>Fitting failed for the experiment ",exp,".
              You may try a) changing the initial estimates or fitting boundaries,
              b) reducing the number of parameters to be fitted,
              c) forcing certain parameters to be positive or negative by including the 
              pattern 'pos' or 'neg' in the parameter name,
              d) scaling all the function parameters to lie between -100 and 100 and
              setting the limits of the log spaced grid search to -2 and 2.</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        return(NULL)
        
      }
      
    }
    
    Sys.sleep(0.5)
    reactives$custom_data_was_fitted <- TRUE
    
    reactives$fitted_coefficients_method_custom <- 'fixedWL'
    
    append_record_to_logbook("Fitting the CD signal versus custom factor curve")
    
    # Set SVD / PCA fit to FALSE
    reactives$custom_data_was_fitted_svd_or_pca <- FALSE
    
    boundaries_updated <- sapply(custom_exps,function(exp) cdAnalyzer$experimentsCustom[[exp]]$boundaries_updated)
    boundaries_updated <- any(boundaries_updated)
    
    if (boundaries_updated) {
      renderInitialParamsTable()
      shinyalert(text = paste("<b>One or more of the fitted parameters was close to the fitting boundaries. 
                              Consequently, the 'initial parameter estimates' Table was automatically updated, 
                              and the data was fitted again.</b>"),
                 type = "info",closeOnEsc = T,closeOnClickOutside = T,
                 html=T)
    }
      
  })
  
  
})

output$fittedParams_custom <- renderTable({
  
  req(reactives$custom_data_was_fitted)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,type = 'Custom')

  return(df)
}, digits = 4)

output$fittedErrors_custom <- renderTable({
  
  req(reactives$custom_data_was_fitted)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,'Custom',errors=TRUE)
  
  return(df)
})

output$fittedCustomCurves <- renderPlotly({
  
  req(reactives$custom_data_was_fitted)
  
  df    <- generate_custom_df(cdAnalyzer)
  dfFit <- generate_custom_df(cdAnalyzer,'signal_predicted')

  fig   <- plot_unfolding_fitting(
    df,dfFit,
    reactives$customWorkingUnits,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom, input$plot_axis_size_custom,
    'fixedWL',input$use_log_axis_custom)
  
  return(fig)
})

output$residualsCustomCurves <- renderPlot({
  
  req(reactives$custom_data_was_fitted)
  
  df      <- generate_custom_df(cdAnalyzer)
  dfFit   <- generate_custom_df(cdAnalyzer,signal_type='signal_predicted')
  
  join_cols <- colnames(df)
  join_cols <- join_cols[join_cols != 'value']
  
  tog           <- inner_join(df,dfFit,by=join_cols) 
  tog$residuals <- tog$value.y - tog$value.x
  
  fig <- plot_residuals(
    tog,
    input$plot_axis_size_custom,
    svd_or_pca_based=FALSE, xlab=cdAnalyzer$experimentsCustom[[1]]$first_exp_param_name,
    input$use_log_axis_custom)
  
  return(fig)
  
})

## Start of SVD/PCA decomposition reactives

output$customSpectra <- renderPlotly({
  
  req(reactives$customDatasetCreated)
  
  df  <- generate_custom_df(cdAnalyzer,signal_type='signalDesiredUnit')
  
  fig <- plot_unfolding_exp_spectra(
    df,
    reactives$customWorkingUnits,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom,  input$plot_axis_size_custom,
    plot_mode=input$plot_style_custom)
  
  return(fig)
  
})

observeEvent(input$btn_decompose_spectra_custom,{
  
  reactives$spectra_was_decomposed_custom <- FALSE
  
  exps <- cdAnalyzer$experimentNamesCustom
  
  if (length(exps) == 0) return(NULL)
  
  explained_variance_threshold <- input$explained_variance_threshold_custom
  
  # vector to store the number of useful components
  ks <- c()
  for (exp in exps) {
    
    if (input$analysis_model_custom == 'spectraDecompositionSVD') {
      cdAnalyzer$experimentsCustom[[exp]]$decompose_spectra_svd()
    }
    
    if (input$analysis_model_custom == 'spectraDecompositionPCA') {
      cdAnalyzer$experimentsCustom[[exp]]$decompose_spectra_pca()
    }
    
    if (!cdAnalyzer$experimentsCustom[[exp]]$decompositionDone) {
      
      shinyalert(text = 
                   paste("<b>The spectra decomposition algorithm did
                   not converge. Please remove noisy data."),
                 type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                 html=T)
      return(NULL)
      
    }
    
    cdAnalyzer$experimentsCustom[[exp]]$filter_basis_spectra(explained_variance_threshold)
    cdAnalyzer$experimentsCustom[[exp]]$align_basis_spectra_and_coefficients()
    cdAnalyzer$experimentsCustom[[exp]]$reconstruct_spectra()
    
    ks <- c(ks,cdAnalyzer$experimentsCustom[[exp]]$k)
  }
  
  reactives$show_basis_change_option <- all(ks %in% c(2,3) )
  
  reactives$spectra_decomposition_method_custom <- gsub(
    'spectraDecomposition','',input$analysis_model_custom)
  
  append_record_to_logbook(paste0("Decomposing the CD spectra",
                                  '. Method: ',reactives$spectra_decomposition_method_custom,
                                  '. Explained variance threshold: ',explained_variance_threshold))
  
  updateNumericInput(session,'selectedK_custom',NULL,min(ks),1,min(ks))
  
  reactives$spectra_was_decomposed_custom <- TRUE
  
})

observeEvent(input$explained_variance_threshold_custom,{
  
  req(reactives$spectra_was_decomposed_custom)
  
  reactives$spectra_was_decomposed_custom <- FALSE
  
  exps <- cdAnalyzer$experimentNamesCustom
  
  explained_variance_threshold <- input$explained_variance_threshold_custom
  
  ks <- c()
  for (exp in exps) {
    
    cdAnalyzer$experimentsCustom[[exp]]$filter_basis_spectra(explained_variance_threshold)
    cdAnalyzer$experimentsCustom[[exp]]$reconstruct_spectra()
    ks <- c(ks,cdAnalyzer$experimentsCustom[[exp]]$k)
  }
  
  reactives$show_basis_change_option <- all(ks %in% c(2,3) )
  
  append_record_to_logbook(paste0("Setting the variance threshold to: ",
                                  explained_variance_threshold))
  
  updateNumericInput(session,'selectedK_custom',NULL,min(ks),1,min(ks))
  reactives$spectra_was_decomposed_custom <- TRUE
  
})

observeEvent(input$btn_change_basis_custom,{
  
  req(reactives$spectra_was_decomposed_custom)
  reactives$spectra_was_decomposed_custom <- FALSE
  
  exps <- cdAnalyzer$experimentNamesCustom
  
  for (exp in exps) {
    
    cdAnalyzer$experimentsCustom[[exp]]$rotate_basis_spectra()
  }
  
  append_record_to_logbook("Applying change of basis.")
  
  reactives$spectra_was_decomposed_custom <- TRUE
})

observeEvent(input$btn_flip_spectrum_custom,{
  
  req(reactives$spectra_was_decomposed_custom)
  
  showModal(modalDialog(
    
    tags$h3('Please select the basis spectrum to invert:'),
    
    selectInput('experiment_to_flip_spectrum_custom','Dataset name',
                choices = c(cdAnalyzer$experimentNamesCustom)),
    
    selectInput('selected_k_spectrum_custom','Selected basis spectrum',
                choices = 1:input$selectedK_custom),
    
    footer=tagList(
      actionButton('submitInversionCustom', 'Submit'),
      modalButton('Cancel')
    )
  ))
  
})

observeEvent(input$submitInversionCustom,{
  
  removeModal()
  
  reactives$spectra_was_decomposed_custom <- FALSE
  
  exp <- input$experiment_to_flip_spectrum_custom
  k   <- as.numeric(input$selected_k_spectrum_custom) 
  
  append_record_to_logbook(paste0("Inverting the  ",
                                  k,"th basis spectrum of the dataset ",
                                  exp))
  
  cdAnalyzer$experimentsCustom[[exp]]$invert_selected_spectrum(k-1) #python index starts at 0
  
  reactives$spectra_was_decomposed_custom <- TRUE
  
})

output$customBasisSpectra <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_custom)
  
  df  <- get_basis_spectra_df(cdAnalyzer,'Custom')
  fig <- plot_basis_spectra(
    df,
    reactives$customWorkingUnits,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom, input$plot_axis_size_custom)
  
  return(fig)
  
})

output$customFittedSpectra <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_custom)
  
  df     <- generate_custom_df(cdAnalyzer,signal_type='signalDesiredUnit')
  dfFit  <- generate_custom_df(cdAnalyzer,signal_type='fitted_spectra')
  fig    <- plot_unfolding_exp_spectra(
    df,
    reactives$customWorkingUnits,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom, input$plot_axis_size_custom,
    dfFit)
  
  return(fig)
})

output$customExplainedVariance <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_custom)
  df  <- get_explained_variance_df(cdAnalyzer,'Custom')
  fig <- plot_explained_variance(
    df,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom, input$plot_axis_size_custom)
  
  return(fig)
  
})

output$customSVDCoefficients <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_custom)
  df  <- get_coefficients_df(cdAnalyzer,'Custom')
  
  fig <- plot_unfolding_exp(
    df,
    reactives$customWorkingUnits,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom, input$plot_axis_size_custom,
    reactives$spectra_decomposition_method_custom,
    input$use_log_axis_custom)
  
  return(fig)
  
})

observeEvent(input$btn_fit_custom_data_svd,{
  
  req(reactives$spectra_was_decomposed_custom)
  exps <- cdAnalyzer$experimentNamesCustom
  
  withBusyIndicatorServer("fitCustomHidden",{
    
    reactives$custom_data_was_fitted_svd_or_pca <- FALSE

    if (length(exps) == 0) return(NULL)
    
    for (exp in exps) {
      cdAnalyzer$experimentsCustom[[exp]]$fit_signal()
      
      if (!is.null(cdAnalyzer$experimentsCustom[[exp]]$fit_params)) {
        
        shinyalert(text = paste("<b>Fitting done for the experiment ",exp,"</b>"),
                   type = "success",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        
      } else {
        
        shinyalert(text = 
                     paste("<b>Fitting failed for the experiment ",exp,".
              You may try a) changing the initial estimates or fitting boundaries,
              b) reducing the number of parameters to be fitted,
              c) forcing certain parameters to be positive or negative by including the 
              pattern 'pos' or 'neg' in the parameter name,
              d) scaling all the function parameters to lie between -100 and 100 and
              setting the limits of the log spaced grid search to -2 and 2.</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        return(NULL)
        
      }
      
    }
    
    Sys.sleep(0.5)
    reactives$custom_data_was_fitted_svd_or_pca <- TRUE
    
    reactives$fitted_coefficients_method_custom <- reactives$spectra_decomposition_method_custom
    
    append_record_to_logbook(paste0("Fitting the CD coefficients versus custom experimental parameter curve",
                                    '. Coefficients of interest: ',input$selectedK_custom))
    
    # Set fixed wavelength fit to FALSE
    reactives$custom_data_was_fitted    <- FALSE
    
    boundaries_updated <- sapply(custom_exps,function(exp) cdAnalyzer$experimentsCustom[[exp]]$boundaries_updated)
    boundaries_updated <- any(boundaries_updated)
    
    if (boundaries_updated) {
      renderInitialParamsTable()
      shinyalert(text = paste("<b>One or more of the fitted parameters was close to the fitting boundaries. 
                              Consequently, the 'initial parameter estimates' Table was automatically updated, 
                              and the data was fitted again.</b>"),
                 type = "info",closeOnEsc = T,closeOnClickOutside = T,
                 html=T)
    }
    
    Sys.sleep(0.5)
    
  })
  
})

output$customFittedSVDCoefficients <- renderPlotly({
  
  req(reactives$custom_data_was_fitted_svd_or_pca)
  
  df    <- generate_custom_df(cdAnalyzer)
  dfFit <- generate_custom_df(cdAnalyzer,'signal_predicted')
  fig   <- plot_unfolding_fitting(
    df,dfFit,
    reactives$customWorkingUnits,
    input$plot_width_custom, input$plot_height_custom, 
    input$plot_type_custom, input$plot_axis_size_custom,
    reactives$fitted_coefficients_method_custom,
    input$use_log_axis_custom)
  
  return(fig)
  
})

output$customResidualsSVDCoefficients <- renderPlot({
  
  req(reactives$custom_data_was_fitted_svd_or_pca)
  
  df      <- generate_custom_df(cdAnalyzer)
  dfFit   <- generate_custom_df(cdAnalyzer,signal_type='signal_predicted')
  
  join_cols <- colnames(df)
  join_cols <- join_cols[join_cols != 'value']
  
  tog           <- inner_join(df,dfFit,by=join_cols) 
  tog$residuals <- tog$value.y - tog$value.x
  
  fig <- plot_residuals(
    tog,
    input$plot_axis_size_custom,
    svd_or_pca_based=TRUE, xlab=cdAnalyzer$experimentsCustom[[1]]$first_exp_param_name,
    input$use_log_axis_custom)
  
  return(fig)
  
})


output$fittedParams_customSVD <- renderTable({
  
  req(reactives$custom_data_was_fitted_svd_or_pca)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,'Custom')
  
  return(df)
  
}, digits = 4)

output$fittedErrors_customSVD <- renderTable({
  
  req(reactives$custom_data_was_fitted_svd_or_pca)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,'Custom',errors=TRUE)
  
  return(df)
  
})
