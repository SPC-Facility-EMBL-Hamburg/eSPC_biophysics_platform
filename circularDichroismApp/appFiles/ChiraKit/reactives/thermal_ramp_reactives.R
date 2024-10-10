# Fit a certain thermal experiment using the unfolding model selected by the user
# Returns TRUE (boolean) if the fitting could be run, otherwise it returns FALSE
fitThermalExperiment <- function(exp) {
  
  pyExp <- cdAnalyzer$experimentsThermal[[exp]]

  user_sel_model <- input$thermal_unfolding_model
  fitSlopeN      <- input$fitSlopeNative
  fitSlopeU      <- input$fitSlopeUnfolded
  
  append_record_to_logbook(paste0("Selected model: ", user_sel_model))
  
  # Models that require the user to input some value for the parameter Cp
  models_req_cp <- c('twoStateDimer','twoState','twoStateTrimer','twoStateTetramer',
                     'threeStateCp','threeStateDimerMICp','threeStateDimerDICp')
  
  if (user_sel_model %in% models_req_cp) {
    
    Cp             <- suppressWarnings(as.numeric(input$Cp))
    is_not_numeric <- is.na(Cp)
    
    if (is_not_numeric) {
      popUpWarning("Please set a valid Cp.")
      return(F)
    }
  }  
  
  if (user_sel_model %in% c('threeStateCp','threeStateDimerMICp','threeStateDimerDICp')) {
    
    if (Cp == 0) {
      popUpWarning("If you want to fit the Cp of the first equilibrium constant, please set a global Cp greater than 0.")
      return(F)
    }
  }  
  
  if (user_sel_model %in% c('twoStateDimer','twoState','twoStateTrimer','twoStateTetramer')) {
    
    selModel <- switch(user_sel_model,
                       'twoStateDimer'    = 'Dimer',
                       'twoStateTrimer'   = 'Trimer',
                       'twoStateTetramer' = 'Tetramer',
                       'Monomer'  # Default value if none of the cases match
    )
      
    pyExp$fit_signal(fitSlopeN,fitSlopeU,Cp,selModel)
    append_record_to_logbook(paste0("Cp (set by the user): ", Cp))
  } 
  
  if (user_sel_model %in% c('threeState','threeStateDimerMI','threeStateDimerDI',
                            'threeStateTrimerMI','threeStateTrimerTI','threeStateTetramerMI')) {
    
    selModel <- switch(user_sel_model,
                       'threeStateDimerMI' = 'Dimer_monomeric_intermediate',
                       'threeStateDimerDI' = 'Dimer_dimeric_intermediate',
                       'threeStateTrimerMI' = 'Trimer_monomeric_intermediate',
                       'threeStateTrimerTI' = 'Trimer_trimeric_intermediate',
                       'threeStateTetramerMI' = 'Tetramer_monomeric_intermediate',
                       'Monomer'  # Default value if none of the cases match
    )
    
    pyExp$fit_signal_three_state(fitSlopeN,fitSlopeU,input$T1_init,input$T2_init,selModel)
    append_record_to_logbook(paste0("Cp1 and Cp2: ", 0))
    append_record_to_logbook(paste0("T1 init: ", input$T1_init, ". T2 init: ", input$T2_init))
  }  
  
  if (user_sel_model %in% c('threeStateCp','threeStateDimerMICp','threeStateDimerDICp')) {
    
    selModel <- 'Monomer'
    if (user_sel_model == 'threeStateDimerMICp') selModel <- 'Dimer_monomeric_intermediate'
    if (user_sel_model == 'threeStateDimerDICp') selModel <- 'Dimer_dimeric_intermediate'
    
    pyExp$fit_signal_three_state_with_cp(fitSlopeN,fitSlopeU,input$T1_init,input$T2_init,Cp,selModel)
    append_record_to_logbook(paste0("Global Cp: ", Cp))
    append_record_to_logbook(paste0("T1 init: ", input$T1_init, ". T2 init: ", input$T2_init))
  }
  
  if (user_sel_model %in% c('twoState_irrev','threeState_rev_irrev')) {
    
    v              <- input$scan_rate
    is_not_numeric <- is.na(suppressWarnings(as.numeric(v)))
    
    if (is_not_numeric) {
      popUpWarning("Please set a valid scan rate. For example, '0.25'. The units are degrees per minute.")
      return(F)
    }
    
    v <- as.numeric(v)
    
    if (v == 0) {
      popUpWarning("The scan rate can not be zero.")
      return(F)
    }
    
    all_temperatures_lst <- pyExp$xAxis_lst
    for (vec in all_temperatures_lst) {
      has_duplicates <- any(duplicated(vec))
      if (has_duplicates) {
        popUpWarning("Please remove temperature duplicates")
        return(F)
      }
    }
    
    append_record_to_logbook(paste0("Scan rate: ", v))
  }
  
  if (user_sel_model == 'twoState_irrev') pyExp$fit_signal_two_state_irrev(v,fitSlopeN,fitSlopeU)

  if (user_sel_model == 'threeState_rev_irrev') {
    
    if (input$T1_init_irrev <= 0) {
      popUpWarning("The initial estimate for T1 must be greater than zero. Fitting process canceled.")
      return(F)
    }
    
    popUpInfo("Fitting started, it might take a few minutes.")
    pyExp$fit_signal_three_state_irrev(v,input$fitSlopeNative,input$fitSlopeUnfolded,input$T1_init_irrev)
    append_record_to_logbook(paste0("T1 init: ", input$T1_init_irrev))
  } 
  
  return(T)
}

# Create the Table to fill with the temperature data
observeEvent(list(input$legendInfo,input$workingUnits,input$oligomeric_state_term),{
  
  req(reactives$data_loaded)
  
  # Retrieve the dataframe from the legends table (1. Load Input Tab)
  legendDf <- getLegendDF(input$legendInfo)
  id       <- legendDf$Internal.ID
  
  temperature <- cdAnalyzer$get_experiment_properties('temperature')
  temperature <- unlist(temperature)

  oligo_state <- input$oligomeric_state_term
  
  if (oligo_state == 'monomer') {
    df           <- data.frame(id,temperature,'A')
    colnames(df) <- c('CD_curve','Temperature (°C/K)','Dataset')
  } else {
    df           <- data.frame(id,0,temperature,'A')
    
    oligoColumnName <- '[Dimer] (μM)'
    if (oligo_state == 'trimer')   oligoColumnName <- '[Trimer] (μM)'
    if (oligo_state == 'tetramer') oligoColumnName <- '[Tetramer] (μM)'
    
    colnames(df) <- c('CD_curve',oligoColumnName,'Temperature (°C/K)','Dataset')
  }  
  
  # Remove non-selected CD curves
  df <- df[legendDf$Show,]
  
  # Remove experiments with non-matching units
  id_to_keep         <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)

  internalID_all     <- cdAnalyzer$get_experiment_properties('internalID')
  internalID_to_keep <- unlist(internalID_all[id_to_keep])
  
  df <- df[df$CD_curve %in% internalID_to_keep,]

  # Assign the created dataframe to the Table thermal_denaturation_data (available at the 2a. Thermal analysis Tab)
  output$thermal_denaturation_data <- renderRHandsontable({
    rhandsontable(df,rowHeaders=NULL)    %>% 
      hot_col(col = c(1), readOnly=TRUE) %>% 
      hot_table(stretchH='all')
  })
  
})

observeEvent(input$btn_create_thermal_dataset,{
  
  req(reactives$data_loaded)
  
  reactives$thermalDatasetCreated  <- NULL
  reactives$spectra_was_decomposed <- NULL
  
  cdAnalyzer$clean_experiments('thermal')
  
  # Retrieve which CD experiments should be used to build
  # a new dataset for thermal denaturation
  df_ids2find        <- hot_to_r(input$thermal_denaturation_data)
  
  # Reassign the column name
  idx <- which(grepl('Temperature',colnames(df_ids2find)))
  colnames(df_ids2find)[idx] <- 'Temperature'
  
  groups             <- unique(df_ids2find$Dataset)
  
  append_record_to_logbook(c('Creating a thermal dataset with the following data',df_to_lines(df_ids2find)))
  
  if (nrow(df_ids2find) < 2) {
    popUpWarning("Please load more CD data")
    return(NULL)
  }
  
  we_have_nas <- sum(is.na(as.numeric(df_ids2find$Temperature))) > 0
  
  if (we_have_nas) {
    popUpWarning("Please check the temperature column, there is missing data or non-numeric values")
    return(NULL)
  }
  
  # Create one thermal dataset per group
  for (group in groups) {
    
    df_temp <- df_ids2find[df_ids2find$Dataset == group,]
    
    relevantSpectra      <- df_temp$CD_curve 
    relevantTemperature  <- df_temp$Temperature 
    
    # To degree Celsius, required latter for plotting... 
    if (max(relevantTemperature) > 250) {
      relevantTemperature <- relevantTemperature - 273.15 
    }
    
    merged <- get_signal_dfs_from_selected_spectra(relevantSpectra,cdAnalyzer)
    
    sorted_indexes      <- order(relevantTemperature)
    relevantTemperature <- relevantTemperature[sorted_indexes]
    
    sorted_signal   <- as.matrix(merged[,-1][,sorted_indexes],drop = FALSE)
    
    # Assign the signal and temperature data to the new thermal unfolding experiment
    cdAnalyzer$experimentsThermal[[group]]                    <- CdExperimentThermalRamp()
    cdAnalyzer$experimentsThermal[[group]]$wavelength         <- np_array(merged[,1])
    cdAnalyzer$experimentsThermal[[group]]$signalDesiredUnit  <- np_array(sorted_signal)
    cdAnalyzer$experimentsThermal[[group]]$temperature        <- np_array(relevantTemperature)
    cdAnalyzer$experimentsThermal[[group]]$temperature_ori    <- np_array(relevantTemperature)
    cdAnalyzer$experimentsThermal[[group]]$name               <- group
    
    # Convert the string of selected wavelengths to a numeric vector
    selected_wl <- parse_selected_wavelengths(input$selected_wavelength_thermal_unfolding)
    
    cdAnalyzer$experimentsThermal[[group]]$assign_useful_signal(selected_wl)
    
    if (input$oligomeric_state_term != 'monomer') {
      relevantOligoConc  <- df_temp[,2]
      relevantOligoConc  <- relevantOligoConc[sorted_indexes]/1e6 # Transform micromolar (user input) to molar
      cdAnalyzer$experimentsThermal[[group]]$oligo_conc_molar <- np_array(relevantOligoConc)

    } else {
      
      updateSelectInput(session,'thermal_unfolding_model',NULL,
                        choices = c(
                          'N ⇌ U'                  = 'twoState',
                          'N ⇌ I ⇌ U'              = 'threeState',
                          'N ⇌ I ⇌ U (fit Cp1)'    = 'threeStateCp',
                          'N → U '                 = 'twoState_irrev',
                          'N ⇌ I → U'              = 'threeState_rev_irrev'))
      
    }
    
    if (input$oligomeric_state_term == 'dimer') {
      updateSelectInput(session,'thermal_unfolding_model',NULL,
                        choices = c(
                          'N2 ⇌ 2U'                     = 'twoStateDimer',
                          'N2 ⇌ 2I ⇌ 2U'                = 'threeStateDimerMI',
                          'N2 ⇌ I2 ⇌ 2U'                = 'threeStateDimerDI',
                          'N2 ⇌ 2M ⇌ 2U (fit Cp of K1)' = 'threeStateDimerMICp',
                          'N2 ⇌ I2 ⇌ 2U (fit Cp of K1)' = 'threeStateDimerDICp'))
    }
    
    if (input$oligomeric_state_term == 'trimer') {
      updateSelectInput(session,'thermal_unfolding_model',NULL,
                        choices = c(
                          'N3 ⇌ 3U'       = 'twoStateTrimer',
                          'N3 ⇌ 3I ⇌ 3U'  = 'threeStateTrimerMI',
                          'N3 ⇌ I3 ⇌ 3U'  = 'threeStateTrimerTI'))
    }
    
    if (input$oligomeric_state_term == 'tetramer') {
      updateSelectInput(session,'thermal_unfolding_model',NULL,
                        choices = c(
                          'N4 ⇌ 4U'       = 'twoStateTetramer',
                          'N4 ⇌ 4I ⇌ 4U'  = 'threeStateTetramerMI'))
    }
    
    cdAnalyzer$experimentsThermal[[group]]$reshape_signal_oligomer('Thermal')    
    
  }

  cdAnalyzer$experimentNamesThermal <- groups
  
  reactives$thermalWorkingUnits   <- input$workingUnits
  reactives$thermalDatasetCreated <- TRUE
  
})

observeEvent(input$thermalUnfoldingFile,{
  
  req(input$thermalUnfoldingFile)
  
  reactives$thermalDatasetCreated  <- NULL
  reactives$spectra_was_decomposed <- NULL
  
  cdAnalyzer$clean_experiments('thermal')
  group <- 'A'
    # Assign the signal and temperature data to the new thermal unfolding experiment
  cdAnalyzer$experimentsThermal[[group]]                    <- CdExperimentThermalRamp()
  cdAnalyzer$experimentsThermal[[group]]$name               <- group
    

  cdAnalyzer$experimentNamesThermal <- c(group)

  is_monomer <- !are_headers_numeric(input$thermalUnfoldingFile$datapath)

  if (is_monomer) {

    loadState <- cdAnalyzer$experimentsThermal[[group]]$load_unfolding_data_monomer(input$thermalUnfoldingFile$datapath)

    updateSelectInput(session,'thermal_unfolding_model',NULL,
                      choices = c(
                        'N ⇌ U'                  = 'twoState',
                        'N ⇌ I ⇌ U'              = 'threeState',
                        'N ⇌ I ⇌ U (fit Cp1)'    = 'threeStateCp'))

    reactives$thermalWorkingUnits   <- 'yMonomer'

  } else {

    loadState <- cdAnalyzer$experimentsThermal[[group]]$load_unfolding_data_oligomer(input$thermalUnfoldingFile$datapath)

    updateSelectInput(session,'thermal_unfolding_model',NULL,
                      choices = c(
                        'N2 ⇌ 2U'                     = 'twoStateDimer',
                        'N2 ⇌ 2I ⇌ 2U'                = 'threeStateDimerMI',
                        'N2 ⇌ I2 ⇌ 2U'                = 'threeStateDimerDI',
                        'N2 ⇌ 2M ⇌ 2U (fit Cp of K1)' = 'threeStateDimerMICp',
                        'N2 ⇌ I2 ⇌ 2U (fit Cp of K1)' = 'threeStateDimerDICp',
                        'N3 ⇌ 3U'       = 'twoStateTrimer',
                        'N3 ⇌ 3I ⇌ 3U'  = 'threeStateTrimerMI',
                        'N3 ⇌ I3 ⇌ 3U'  = 'threeStateTrimerTI',
                        'N4 ⇌ 4U'       = 'twoStateTetramer',
                        'N4 ⇌ 4I ⇌ 4U'  = 'threeStateTetramerMI'))

    reactives$thermalWorkingUnits   <- 'yOligomer'

  }

  cdAnalyzer$experimentsThermal[[group]]$reshape_signal_oligomer('Thermal')

    if (!loadState) {
      popUpWarning("The file could not be loaded. Please check the format.")
      return(NULL)
    } else {
      popUpSuccess("The file was loaded successfully. Navigate to the Thermal unfolding panel")
    }

  reactives$thermalDatasetCreated <- TRUE
  
})

observeEvent(input$btn_find_wl,{
  
  req(reactives$thermalDatasetCreated)
  
  thermal_exps            <- cdAnalyzer$experimentNamesThermal
  
  wavelength_filtered_all <- list()
  counter                 <- 0
  
  for (exp in thermal_exps) {
    
    counter <- counter + 1
    temperature <- cdAnalyzer$experimentsThermal[[exp]]$temperature
    cdAnalyzer$experimentsThermal[[exp]]$estimate_useful_signal_based_on_snr_and_amplitude(temperature)
    wavelength_filtered_all[[counter]] <- cdAnalyzer$experimentsThermal[[exp]]$wavelength_filtered
    
  }
  
  wavelength_filtered_common <- Reduce(intersect, wavelength_filtered_all)
  if (length(wavelength_filtered_common) == 0) {
    
    popUpWarning("The automatic selection algorithm didn't work, select the wavelength manually")
    wavelength_filtered_common <- 220
  }
    
  append_record_to_logbook('Running the automatic wavelength selection')
  
  updateTextInput(session, "selected_wavelength_thermal_unfolding", 
                  value = paste(wavelength_filtered_common,collapse=" "))
    
})

observeEvent(list(input$selected_wavelength_thermal_unfolding,input$analysis_model_thermal),{
  
  req(reactives$thermalDatasetCreated)
  
  if (input$analysis_model_thermal == 'fixedWL') {
    append_record_to_logbook(paste0("Setting the 'Selected wavelength(s)' to: ",
                                    input$selected_wavelength_thermal_unfolding))
    
    reactives$spectra_decomposition_method_thermal <- 'None'
    
  }

  reactives$thermalDatasetCreated <- FALSE
  Sys.sleep(0.1)
  # Convert the string of selected wavelengths to a numeric vector
  selected_wl <- parse_selected_wavelengths(input$selected_wavelength_thermal_unfolding)

  thermal_exps <- cdAnalyzer$experimentNamesThermal
  
  for (exp in thermal_exps) {
    
    cdAnalyzer$experimentsThermal[[exp]]$assign_useful_signal(selected_wl)
    cdAnalyzer$experimentsThermal[[exp]]$reshape_signal_oligomer('Thermal')
    
  }
  
  reactives$thermalDatasetCreated <- TRUE
  
})

output$meltingCurves <- renderPlotly({
  
  req(reactives$thermalDatasetCreated)
  
  thermal_ramp_df <- generate_thermal_ramp_df(cdAnalyzer)
  
  plot_unfolding_exp(thermal_ramp_df,
                     reactives$thermalWorkingUnits,
                     input$plot_width_melt, input$plot_height_melt, 
                     input$plot_type_melt, input$plot_axis_size_melt,
                     reactives$spectra_decomposition_method_thermal,
                     F,input$x_legend_pos,input$y_legend_pos,input$show_title
                     )  
  
})

observeEvent(input$btn_fit_melting_data,{
  
  req(reactives$thermalDatasetCreated)
  
  append_record_to_logbook("Fitting procedure started")
  
  withBusyIndicatorServer("fitThermalHidden",{
    
    reactives$melting_data_was_fitted <- FALSE
    thermal_exps <- cdAnalyzer$experimentNamesThermal
    
    if (length(thermal_exps) == 0) return(NULL)
    
    for (exp in thermal_exps) {
      
      fit <- fitThermalExperiment(exp)
      if (!fit) return(NULL)
    }
    
    Sys.sleep(0.5)
    reactives$melting_data_was_fitted <- TRUE
    
    reactives$fitted_coefficients_method_thermal <- 'fixedWL'
    
    append_record_to_logbook(
      paste0('. Fit native slope mode: ',input$fitSlopeNative,
             '. Fit unfolded slope mode: ',input$fitSlopeUnfolded))
    
    # Set SVD / PCA fit to FALSE
    reactives$melting_data_was_fitted_svd_or_pca <- FALSE
    
  })
  
  popUpSuccess("<b>Fitting done!</b>")
    
})

output$fittedParams_melting <- renderTable({
  
  req(reactives$melting_data_was_fitted)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,type='Thermal')
  
  return(df)
})

output$fittingBounds_melting <- renderTable({
  
  req(reactives$melting_data_was_fitted)
  df <- get_fitting_bounds_unfolding(cdAnalyzer,type='Thermal')

  return(df)
})


output$fittedErrors_melting <- renderTable({
  
  req(reactives$melting_data_was_fitted)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,errors=TRUE)
  
  return(df)
  
})

output$fittedMeltingCurves <- renderPlotly({
 
  req(reactives$melting_data_was_fitted)
  
  df      <- generate_thermal_ramp_df(cdAnalyzer)
  dfFit   <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signal_predicted')
  fig     <- plot_unfolding_fitting(
    df,dfFit,
    reactives$thermalWorkingUnits,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    reactives$fitted_coefficients_method_thermal,
    F,input$x_legend_pos,input$y_legend_pos,input$show_title
    )
  
  return(fig)
  
})

output$fractions_melting <- renderPlotly({
  
  req(reactives$melting_data_was_fitted)
  
  fractions_df <- generate_fractions_df(cdAnalyzer)
  
  fig <- plot_unfolding_fractions(
    fractions_df,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    'Temperature (°C)',
    xLegend   = input$x_legend_pos,
    yLegend   = input$y_legend_pos,
    showTitle = input$show_title)
  
  return(fig)
  
})

output$residualsMeltingCurves <- renderPlot({
  
  req(reactives$melting_data_was_fitted)
  
  df      <- generate_thermal_ramp_df(cdAnalyzer)
  dfFit   <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signal_predicted')
  
  tog           <- inner_join(df,dfFit,by=c('wavelength','temperature','legend'),relationship = "many-to-many") 
  tog$residuals <- tog$value.y - tog$value.x
    
  fig <- plot_residuals(
    tog,
    input$plot_axis_size_melt)
                        
  return(fig)
  
})

output$meltingSpectra <- renderPlotly({
  
  req(reactives$thermalDatasetCreated)

  df  <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signalDesiredUnit')
  fig <- plot_unfolding_exp_spectra(
    df,
    reactives$thermalWorkingUnits,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    plot_mode = input$plot_style_melt,
    xLegend   = input$x_legend_pos,
    yLegend   = input$y_legend_pos,
    showTitle = input$show_title)
  
  return(fig)
})

observeEvent(input$btn_decompose_spectra,{
  
  reactives$spectra_was_decomposed <- NULL
  
  thermal_exps <- cdAnalyzer$experimentNamesThermal
  
  if (length(thermal_exps) == 0) return(NULL)
  
  explained_variance_threshold <- input$explained_variance_threshold

  # vector to store the number of useful components
  ks <- c()
  for (exp in thermal_exps) {

    if (input$analysis_model_thermal == 'spectraDecompositionSVD') {
      cdAnalyzer$experimentsThermal[[exp]]$decompose_spectra_svd()
    }
    
    if (input$analysis_model_thermal == 'spectraDecompositionPCA') {
      cdAnalyzer$experimentsThermal[[exp]]$decompose_spectra_pca()
    }

    if (!cdAnalyzer$experimentsThermal[[exp]]$decompositionDone) {
      
      popUpWarning("The spectra decomposition algorithm did not converge. Please remove noisy data.")
      return(NULL)
      
    }
    
    cdAnalyzer$experimentsThermal[[exp]]$filter_basis_spectra(explained_variance_threshold)
    cdAnalyzer$experimentsThermal[[exp]]$align_basis_spectra_and_coefficients()
    cdAnalyzer$experimentsThermal[[exp]]$reconstruct_spectra()
    
    ks <- c(ks,cdAnalyzer$experimentsThermal[[exp]]$k)
  }
  
  reactives$show_basis_change_option <- all(ks %in% c(2,3) )
  
  updateNumericInput(session,'selectedK',NULL,min(ks),1,min(ks))
  
  reactives$spectra_decomposition_method_thermal <- gsub(
    'spectraDecomposition','',input$analysis_model_thermal)
  
  append_record_to_logbook(paste0("Decomposing the CD spectra",
                                  '. Method: ',reactives$spectra_decomposition_method_thermal,
                                  '. Explained variance threshold: ',explained_variance_threshold))
  
  reactives$spectra_was_decomposed <- TRUE
  
})

observeEvent(input$explained_variance_threshold,{
  
  req(reactives$spectra_was_decomposed)
  
  reactives$spectra_was_decomposed <- NULL
  
  thermal_exps <- cdAnalyzer$experimentNamesThermal
  
  explained_variance_threshold <- input$explained_variance_threshold
  
  ks <- c()
  for (exp in thermal_exps) {
    
    cdAnalyzer$experimentsThermal[[exp]]$filter_basis_spectra(explained_variance_threshold)
    cdAnalyzer$experimentsThermal[[exp]]$reconstruct_spectra()
    ks <- c(ks,cdAnalyzer$experimentsThermal[[exp]]$k)
  }
  
  reactives$show_basis_change_option <- all(ks %in% c(2,3) )
  
  append_record_to_logbook(paste0("Setting the variance threshold to: ",
                                  explained_variance_threshold))
  
  updateNumericInput(session,'selectedK',NULL,min(ks),1,min(ks))
  reactives$spectra_was_decomposed <- TRUE
  
})

observeEvent(input$btn_flip_spectrum,{
  
  req(reactives$spectra_was_decomposed)

  showModal(modalDialog(
    
    tags$h3('Please select the basis spectrum to invert:'),

    selectInput('experiment_to_flip_spectrum','Dataset name',
                choices = c(cdAnalyzer$experimentNamesThermal)),
    
    selectInput('selected_k_spectrum','Selected basis spectrum',
                choices = 1:input$selectedK),
    
    footer=tagList(
      actionButton('submitInversion', 'Submit'),
      modalButton('Cancel')
    )
  ))
  
})

observeEvent(input$submitInversion,{
  
  removeModal()
  
  reactives$spectra_was_decomposed <- NULL
  
  exp <- input$experiment_to_flip_spectrum
  k   <- as.numeric(input$selected_k_spectrum) - 1
  
  append_record_to_logbook(paste0("Inverting the  ",
                                  k+1,"th basis spectrum of the dataset ",
                                  exp))
  
  cdAnalyzer$experimentsThermal[[exp]]$invert_selected_spectrum(k)
  
  reactives$spectra_was_decomposed <- TRUE
  
})

observeEvent(input$btn_change_basis,{
  req(reactives$spectra_was_decomposed)
  reactives$spectra_was_decomposed <- NULL
  
  thermal_exps <- cdAnalyzer$experimentNamesThermal
  
  for (exp in thermal_exps) {
    
    cdAnalyzer$experimentsThermal[[exp]]$rotate_basis_spectra()
  }
  
  append_record_to_logbook("Applying change of basis.")
  
  reactives$spectra_was_decomposed <- TRUE
})

output$basisSpectra <- renderPlotly({
  
  req(reactives$spectra_was_decomposed)
  
  df  <- get_basis_spectra_df(cdAnalyzer)
  fig <- plot_basis_spectra(
    df,
    reactives$thermalWorkingUnits,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    xLegend   = input$x_legend_pos,
    yLegend   = input$y_legend_pos,
    showTitle = input$show_title)
  
  return(fig)
})

output$fittedSpectra <- renderPlotly({
  
  req(reactives$spectra_was_decomposed)
  
  df     <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signalDesiredUnit')
  dfFit  <- generate_thermal_ramp_df(cdAnalyzer,signal_type='fitted_spectra')
  fig    <- plot_unfolding_exp_spectra(
    df,
    reactives$thermalWorkingUnits,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    dfFit,
    xLegend   = input$x_legend_pos,
    yLegend   = input$y_legend_pos,
    showTitle = input$show_title
    )
  
  return(fig)
})

# To debug
output$explainedVariance <- renderPlotly({
  
  req(reactives$spectra_was_decomposed)
  df  <- get_explained_variance_df(cdAnalyzer,'Thermal')
  fig <- plot_explained_variance(
    df,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    xLegend   = input$x_legend_pos,
    yLegend   = input$y_legend_pos,
    showTitle = input$show_title)
  
  return(fig)
})

output$svdCoefficients <- renderPlotly({
  
  req(reactives$spectra_was_decomposed)
  df  <- get_coefficients_df(cdAnalyzer)
  
  fig <- plot_unfolding_exp(
    df, reactives$thermalWorkingUnits,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    reactives$spectra_decomposition_method_thermal,
    F,input$x_legend_pos,input$y_legend_pos,input$show_title)
  
  return(fig)
})

observeEvent(input$btn_fit_melting_data_svd,{
  
  req(reactives$spectra_was_decomposed)
  
  withBusyIndicatorServer("fitThermalHidden",{
    
    reactives$melting_data_was_fitted_svd_or_pca <- FALSE
    thermal_exps <- cdAnalyzer$experimentNamesThermal
    
    if (length(thermal_exps) == 0) return(NULL)
    
    for (exp in thermal_exps) {
      
      cdAnalyzer$experimentsThermal[[exp]]$assign_useful_signal_svd(input$selectedK)
      cdAnalyzer$experimentsThermal[[exp]]$reshape_signal_oligomer('Thermal')
      
      Sys.sleep(0.1)
      fitThermalExperiment(exp)
      
    }
    
    Sys.sleep(0.5) 
    
    reactives$melting_data_was_fitted_svd_or_pca <- TRUE
    
    append_record_to_logbook(paste0("Fitting the CD coefficients versus Temperature curve",
                                    '. Coefficients of interest: ',input$selectedK,
                                    '. Fit native slope mode: ',input$fitSlopeNative,
                                    '. Fit unfolded slope mode: ',input$fitSlopeUnfolded))
    
    reactives$fitted_coefficients_method_thermal <- reactives$spectra_decomposition_method_thermal
    
    # Set fixed wavelength fit to FALSE
    reactives$melting_data_was_fitted <- FALSE
    
    Sys.sleep(0.5)    
  })
  
  popUpSuccess("<b>Fitting done!</b>")
  
})

output$fittedSVDCoefficients <- renderPlotly({
  
  req(reactives$melting_data_was_fitted_svd_or_pca)
  
  df      <- generate_thermal_ramp_df(cdAnalyzer)
  dfFit   <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signal_predicted')
  fig     <- plot_unfolding_fitting(
    df,dfFit,
    reactives$thermalWorkingUnits,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    reactives$fitted_coefficients_method_thermal,
    F,input$x_legend_pos,input$y_legend_pos,input$show_title)
  
  return(fig)
  
})

output$residualsSVDCoefficients <- renderPlot({
  
  req(reactives$melting_data_was_fitted_svd_or_pca)
  
  df      <- generate_thermal_ramp_df(cdAnalyzer)
  dfFit   <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signal_predicted')

  tog           <- inner_join(df,dfFit,by=c('wavelength','temperature','legend'),relationship = "many-to-many") 
  tog$residuals <- tog$value.y - tog$value.x
  
  fig <- plot_residuals(
    tog,
    input$plot_axis_size_melt,
    svd_or_pca_based=TRUE)
  
  return(fig)
  
})

output$fittedParams_meltingSVD <- renderTable({
  
  req(reactives$melting_data_was_fitted_svd_or_pca)
  
  df <- get_fitted_params_unfolding(cdAnalyzer)
  
  return(df)
})

output$fittingBoundsSVD_melting <- renderTable({
  
  req(reactives$melting_data_was_fitted_svd_or_pca)
  df <- get_fitting_bounds_unfolding(cdAnalyzer)
  
  return(df)
})

output$fittedErrors_meltingSVD <- renderTable({
  
  req(reactives$melting_data_was_fitted_svd_or_pca)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,errors=TRUE)
  
  return(df)
  
})

output$fractions_melting_svd <- renderPlotly({
  
  req(reactives$melting_data_was_fitted_svd_or_pca)
  
  fractions_df <- generate_fractions_df(cdAnalyzer)
  
  fig <- plot_unfolding_fractions(
    fractions_df,
    input$plot_width_melt, input$plot_height_melt, 
    input$plot_type_melt, input$plot_axis_size_melt,
    'Temperature (°C)',
    xLegend   = input$x_legend_pos,
    yLegend   = input$y_legend_pos,
    showTitle = input$show_title)
  
  return(fig)
  
})
