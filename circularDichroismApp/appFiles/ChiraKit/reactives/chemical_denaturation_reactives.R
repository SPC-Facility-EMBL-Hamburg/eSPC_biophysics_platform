# Fit a certain chemical experiment using the unfolding model selected by the user
fitChemicalExperiment <- function(exp) {
  
  chem_conc <- cdAnalyzer$experimentsChemical[[exp]]$chem_concentration
  cdAnalyzer$experimentsChemical[[exp]]$estimate_baselines_parameters(chem_conc,2)
  
  if (input$chemical_unfolding_model == 'twoState') {
    cdAnalyzer$experimentsChemical[[exp]]$fit_signal(
      input$fitSlopeNative_chemical,input$fitSlopeUnfolded_chemical)
  } 
  
  if (input$chemical_unfolding_model == 'threeState') {
    cdAnalyzer$experimentsChemical[[exp]]$fit_signal_three_state(
      input$fitSlopeNative_chemical,input$fitSlopeUnfolded_chemical,
      input$D50v1_init,input$D50v2_init)
  }
  
  return(NULL)
}

# Create the Table to fill with the concentration of the chemical denaturation agent
observeEvent(list(input$legendInfo,input$workingUnits),{
  
  req(reactives$data_loaded)
  
  # Retrieve the dataframe from the legends table (1. Load Input Tab)
  legendDf <- getLegendDF(input$legendInfo)
  id  <- legendDf$Internal.ID
  
  # Initialize the dataframe with the CD curves and the denaturant agent concentration
  df           <- data.frame(id,0,'A',0)
  colnames(df) <- c('CD_curve','[Denaturant agent] (M)','Dataset_name','Temperature (°C or K)')
  
  # Remove non-selected CD curves
  df <- df[legendDf$Show,]
  
  # Remove experiments with non-matching units
  id_to_keep         <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
  internalID_all     <- cdAnalyzer$getExperimentProperties('internalID')
  internalID_to_keep <- unlist(internalID_all[id_to_keep])
  
  df <- df[df$CD_curve %in% internalID_to_keep,]
  
  # Assign the created dataframe to the Table chemical_denaturation_conc (available at the 2b. Chemical unfolding Tab)
  output$chemical_denaturation_conc <- renderRHandsontable({
    rhandsontable(df,rowHeaders=NULL)    %>% 
      hot_col(col = c(1), readOnly=TRUE) %>% 
      hot_table(stretchH='all')
  })
  
}) 

observeEvent(input$btn_create_chemical_dataset,{
  
  req(reactives$data_loaded)
  
  reactives$chemicalDatasetCreated          <- NULL
  reactives$spectra_was_decomposed_chemical <- NULL
  
  cdAnalyzer$clean_experiments('chemical')
  
  # Retrieve which CD experiments should be used to build a new dataset for chemical denaturation
  df_ids2find        <- hot_to_r(input$chemical_denaturation_conc)
 
  # Replace column name
  colnames(df_ids2find)[4] <- 'Temperature'
  
  # Remove CD curves with denaturant concentration equal to zero
  df_ids2find        <- df_ids2find[df_ids2find[,2] > 0,]
  groups             <- unique(df_ids2find$Dataset_name)
  
  # Stop if there are no CD curves of interest  
  if (nrow(df_ids2find) == 0 ) return(NULL)
  
  append_record_to_logbook(c('Creating a chemical dataset with the following data',df_to_lines(df_ids2find)))
  
  # Check that there is only one temperature per group
  nested_df <- df_ids2find %>% 
    group_by(Dataset_name) %>% 
    summarize(unique_count = n_distinct(Temperature))
  
  if (any(nested_df$unique_count > 1)) {
    
    shinyalert(text = "Please verify that all the CD curves sharing the same dataset name
               have also the same temperature.",
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    
    return(NULL)
  }
  
  # Create one thermal dataset per group
  for (group in groups) {
    
    df_temp <- df_ids2find[df_ids2find$Dataset_name == group,]
    
    relevantSpectra    <- df_temp$CD_curve 
    relevantDenatConc  <- df_temp[,2]
    
    merged <- get_signal_dfs_from_selected_spectra(relevantSpectra,cdAnalyzer)
    
    sorted_indexes      <- order(relevantDenatConc)
    relevantDenatConc   <- relevantDenatConc[sorted_indexes]
    
    sorted_signal   <- as.matrix(merged[,-1][,sorted_indexes],drop = FALSE)
    
    # Assign the signal and temperature data to the new chemical unfolding experiment
    cdAnalyzer$experimentsChemical[[group]]                    <- cd_experiment_chemical_unfolding()
    cdAnalyzer$experimentsChemical[[group]]$wavelength         <- np_array(merged[,1])
    cdAnalyzer$experimentsChemical[[group]]$signalDesiredUnit  <- np_array(sorted_signal)
    cdAnalyzer$experimentsChemical[[group]]$temperature        <- mean(df_temp$Temperature)
    cdAnalyzer$experimentsChemical[[group]]$name               <- group
    cdAnalyzer$experimentsChemical[[group]]$chem_concentration <- np_array(relevantDenatConc)
    
    # Convert the string of selected wavelengths to a numeric vector
    selected_wl <- parse_selected_wavelengths(input$selected_wavelength_thermal_unfolding)
    
    cdAnalyzer$experimentsChemical[[group]]$assign_useful_signal(selected_wl)
    
  }
  
  cdAnalyzer$experimentNamesChemical <- groups
  
  reactives$chemicalWorkingUnits   <- input$workingUnits
  reactives$chemicalDatasetCreated <- TRUE
  
})

observeEvent(input$btn_find_wl_chemical,{
  
  req(reactives$chemicalDatasetCreated)
  
  chemical_exps           <- cdAnalyzer$experimentNamesChemical
  wavelength_filtered_all <- list()
  counter                 <- 0
  
  for (exp in chemical_exps) {
    
    counter  <- counter + 1
    chem_con <- cdAnalyzer$experimentsChemical[[exp]]$chem_concentration
    cdAnalyzer$experimentsChemical[[exp]]$estimate_useful_signal_based_on_snr_and_amplitude(chem_con)
    wavelength_filtered_all[[counter]] <- cdAnalyzer$experimentsChemical[[exp]]$wavelength_filtered
    
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
  
  updateTextInput(session, "selected_wavelength_chemical_unfolding", 
                  value = paste(wavelength_filtered_common,collapse=" "))
  
})

observeEvent(list(input$selected_wavelength_chemical_unfolding,input$analysis_model_chemical),{
  
  req(reactives$chemicalDatasetCreated)
  
  if (input$analysis_model_chemical == 'fixedWL') {
    append_record_to_logbook(paste0("Setting the 'Selected wavelength(s)' to: ",
                                    input$selected_wavelength_chemical_unfolding))
    
    reactives$spectra_decomposition_method_chemical <- 'None'
    
  }
  
  reactives$data_loaded <- FALSE
  # Convert the string of selected wavelengths to a numeric vector
  selected_wl   <- parse_selected_wavelengths(input$selected_wavelength_chemical_unfolding)
  
  chemical_exps <- cdAnalyzer$experimentNamesChemical
  
  for (exp in chemical_exps) {
    
    cdAnalyzer$experimentsChemical[[exp]]$assign_useful_signal(selected_wl)
    
  }
  
  reactives$data_loaded <- TRUE
  
})
 
output$chemicalCurves <- renderPlotly({
  
  req(reactives$data_loaded)
  req(reactives$chemicalDatasetCreated)

  df <- generate_chemical_unfolding_df(cdAnalyzer)
  
  fig <- plot_unfolding_exp(
    df,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem)
  
  return(fig)
  
})

observeEvent(input$btn_fit_chemical_data,{
  
  req(reactives$chemicalDatasetCreated)

  withBusyIndicatorServer("fitChemicalHidden",{
    
    reactives$chemical_data_was_fitted <- FALSE
    chemical_exps <- cdAnalyzer$experimentNamesChemical
    
    if (length(chemical_exps) == 0) return(NULL)
    
    for (exp in chemical_exps) {
      
      fitChemicalExperiment(exp)
                                              
    }
    
    Sys.sleep(0.5)
    reactives$chemical_data_was_fitted <- TRUE
    
    reactives$fitted_coefficients_method_chemical <- 'fixedWL'
    
    append_record_to_logbook(paste0("Fitting the CD signal versus chemical agent concentration curve",
                                    '. Fit native slope mode: ',input$fitSlopeNative_chemical,
                                    '. Fit unfolded slope mode: ',input$fitSlopeUnfolded_chemical))
    
    # Set SVD / PCA fit to FALSE
    reactives$chemical_data_was_fitted_svd_or_pca <- FALSE
    Sys.sleep(0.5)
    
  })
  
  shinyalert(text = paste("<b>Fitting done!</b>"),
             type = "success",closeOnEsc = T,closeOnClickOutside = T,
             html=T)
  
})


# To debug why it randomly collapses shiny when clicking it...
output$fittedParams_chemical <- renderTable({
  
  req(reactives$chemical_data_was_fitted)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,'Chemical')
  
  return(df)
})

output$fittedErrors_chemical <- renderTable({
  
  req(reactives$chemical_data_was_fitted)

  df <- get_fitted_params_unfolding(cdAnalyzer,'Chemical',errors=TRUE)
  
  return(df)
})


output$fittedChemicalCurves <- renderPlotly({
  
  req(reactives$chemical_data_was_fitted)
  
  df    <- generate_chemical_unfolding_df(cdAnalyzer)
  dfFit <- generate_chemical_unfolding_df(cdAnalyzer,'signal_predicted')
  fig   <- plot_unfolding_fitting(
    df,dfFit,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem)
  
  return(fig)
})

output$fractions_chemical <- renderPlotly({
  
  req(reactives$chemical_data_was_fitted)
  
  fractions_df <- generate_fractions_df(cdAnalyzer,'Chemical')
  
  fig <- plot_unfolding_fractions(
    fractions_df,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem,
    'Denaturant agent concentration (M)')
  
  return(fig)
  
})

output$residualsChemicalCurves <- renderPlot({
  
  req(reactives$chemical_data_was_fitted)
  
  df      <- generate_chemical_unfolding_df(cdAnalyzer)
  dfFit   <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='signal_predicted')
  
  tog           <- inner_join(df,dfFit,by=c('wavelength','chem_conc','legend'),relationship = "many-to-many") 
  tog$residuals <- tog$value.y - tog$value.x
  
  fig <- plot_residuals(
    tog,
    input$plot_axis_size_chem,
    svd_or_pca_based=FALSE,xlab="[Denaturant agent] (M)")
  
  return(fig)
  
})

## Start of SVD/PCA decomposition reactives

output$chemUnfoldingSpectra <- renderPlotly({
  
  req(reactives$chemicalDatasetCreated)
  
  df  <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='signalDesiredUnit')
  fig <- plot_unfolding_exp_spectra(
    df,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem,
    plot_mode=input$plot_style_chem)
  
  return(fig)
  
})

observeEvent(input$btn_decompose_spectra_chemical,{
  
  reactives$spectra_was_decomposed_chemical <- NULL
  
  chemical_exps <- cdAnalyzer$experimentNamesChemical
  
  if (length(chemical_exps) == 0) return(NULL)
  
  explained_variance_threshold <- input$explained_variance_threshold_chemical
  
  # vector to store the number of useful components
  ks <- c()
  for (exp in chemical_exps) {
    
    if (input$analysis_model_chemical == 'spectraDecompositionSVD') {
      cdAnalyzer$experimentsChemical[[exp]]$decompose_spectra_svd()
    }
    
    if (input$analysis_model_chemical == 'spectraDecompositionPCA') {
      cdAnalyzer$experimentsChemical[[exp]]$decompose_spectra_pca()
    }
    
    if (!cdAnalyzer$experimentsChemical[[exp]]$decompositionDone) {
      
      shinyalert(text = 
                   paste("<b>The spectra decomposition algorithm did
                   not converge. Please remove noisy data."),
                 type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                 html=T)
      return(NULL)
      
    }
    
    cdAnalyzer$experimentsChemical[[exp]]$filter_basis_spectra(explained_variance_threshold)
    cdAnalyzer$experimentsChemical[[exp]]$align_basis_spectra_and_coefficients()
    cdAnalyzer$experimentsChemical[[exp]]$reconstruct_spectra()
    
    ks <- c(ks,cdAnalyzer$experimentsChemical[[exp]]$k)
  }
  
  reactives$show_basis_change_option <- all(ks %in% c(2,3) )
  
  reactives$spectra_decomposition_method_chemical <- gsub(
    'spectraDecomposition','',input$analysis_model_chemical)
  
  append_record_to_logbook(paste0("Decomposing the CD spectra",
                                  '. Method: ',reactives$spectra_decomposition_method_chemical,
                                  '. Explained variance threshold: ',explained_variance_threshold))
  
  updateNumericInput(session,'selectedK_chemical',NULL,min(ks),1,min(ks))
  
  reactives$spectra_was_decomposed_chemical <- TRUE
  
})

observeEvent(input$explained_variance_threshold_chemical,{
  
  req(reactives$spectra_was_decomposed_chemical)
  
  reactives$spectra_was_decomposed_chemical <- NULL
  
  chemical_exps <- cdAnalyzer$experimentNamesChemical
  
  explained_variance_threshold <- input$explained_variance_threshold_chemical
  
  ks <- c()
  for (exp in chemical_exps) {
    
    cdAnalyzer$experimentsChemical[[exp]]$filter_basis_spectra(explained_variance_threshold)
    cdAnalyzer$experimentsChemical[[exp]]$reconstruct_spectra()
    ks <- c(ks,cdAnalyzer$experimentsChemical[[exp]]$k)
  }
  
  reactives$show_basis_change_option <- all(ks %in% c(2,3) )
  
  append_record_to_logbook(paste0("Setting the variance threshold to: ",
                                  explained_variance_threshold))
  
  updateNumericInput(session,'selectedK_chemical',NULL,min(ks),1,min(ks))
  reactives$spectra_was_decomposed_chemical <- TRUE
  
})

observeEvent(input$btn_change_basis_chemical,{
  
  req(reactives$spectra_was_decomposed_chemical)
  reactives$spectra_was_decomposed_chemical <- NULL
  
  chemical_exps <- cdAnalyzer$experimentNamesChemical
  
  for (exp in chemical_exps) {
    
    cdAnalyzer$experimentsChemical[[exp]]$rotate_basis_spectra()
  }
  
  append_record_to_logbook("Applying change of basis.")
  
  reactives$spectra_was_decomposed_chemical <- TRUE
})

observeEvent(input$btn_flip_spectrum_chemical,{
  
  req(reactives$spectra_was_decomposed_chemical)
  
  showModal(modalDialog(
    
    tags$h3('Please select the basis spectrum to invert:'),
    
    selectInput('experiment_to_flip_spectrum_chemical','Dataset name',
                choices = c(cdAnalyzer$experimentNamesChemical)),
    
    selectInput('selected_k_spectrum_chemical','Selected basis spectrum',
                choices = 1:input$selectedK_chemical),
    
    footer=tagList(
      actionButton('submitInversionChemical', 'Submit'),
      modalButton('Cancel')
    )
  ))
  
})

observeEvent(input$submitInversionChemical,{
  
  removeModal()
  
  reactives$spectra_was_decomposed_chemical <- NULL
  
  exp <- input$experiment_to_flip_spectrum_chemical
  k   <- as.numeric(input$selected_k_spectrum_chemical) 
  
  append_record_to_logbook(paste0("Inverting the  ",
                                  k,"th basis spectrum of the dataset ",
                                  exp))
  
  cdAnalyzer$experimentsChemical[[exp]]$invert_selected_spectrum(k-1) # python index starts at 0
  
  reactives$spectra_was_decomposed_chemical <- TRUE
  
})

output$chemBasisSpectra <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_chemical)
  
  df  <- get_basis_spectra_df(cdAnalyzer,'Chemical')
  fig <- plot_basis_spectra(
    df,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem)
  
  return(fig)
  
})

output$chemFittedSpectra <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_chemical)
  
  df     <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='signalDesiredUnit')
  dfFit  <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='fitted_spectra')
  fig    <- plot_unfolding_exp_spectra(
    df,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem,
    dfFit)
  
  return(fig)
})



output$chemExplainedVariance <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_chemical)
  df  <- get_explained_variance_df(cdAnalyzer,'Chemical')
  fig <- plot_explained_variance(
    df,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem)
  
  return(fig)
  
})

output$chemSVDCoefficients <- renderPlotly({
  
  req(reactives$spectra_was_decomposed_chemical)
  df  <- get_coefficients_df(cdAnalyzer,'Chemical')
  
  fig <- plot_unfolding_exp(
    df,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem,
    reactives$spectra_decomposition_method_chemical)
  
  return(fig)
  
})

observeEvent(input$btn_fit_chemical_data_svd,{
  
  req(reactives$spectra_was_decomposed_chemical)
  
  withBusyIndicatorServer("fitChemicalHidden",{
    
    reactives$chemical_data_was_fitted_svd_or_pca <- FALSE
    chemical_exps <- cdAnalyzer$experimentNamesChemical
    
    if (length(chemical_exps) == 0) return(NULL)
    
    for (exp in chemical_exps) {
      
      cdAnalyzer$experimentsChemical[[exp]]$assign_useful_signal_svd(input$selectedK_chemical)
      fitChemicalExperiment(exp)
      
    }
    
    Sys.sleep(0.5)
    reactives$chemical_data_was_fitted_svd_or_pca <- TRUE
    
    reactives$fitted_coefficients_method_chemical <- reactives$spectra_decomposition_method_chemical
    
    append_record_to_logbook(paste0("Fitting the CD coefficients versus chemical agent concentration curve",
                                    '. Coefficients of interest: ',input$selectedK_chemical,
                                    '. Fit native slope mode: ',input$fitSlopeNative_chemical,
                                    '. Fit unfolded slope mode: ',input$fitSlopeUnfolded_chemical))
    
    # Set fixed wavelength fit to FALSE
    reactives$chemical_data_was_fitted_    <- FALSE
    Sys.sleep(0.5)
    
  })
  
  shinyalert(text = paste("<b>Fitting done!</b>"),
             type = "success",closeOnEsc = T,closeOnClickOutside = T,
             html=T)
  
})

output$chemFittedSVDCoefficients <- renderPlotly({
  
  req(reactives$chemical_data_was_fitted_svd_or_pca)
  
  df    <- generate_chemical_unfolding_df(cdAnalyzer)
  dfFit <- generate_chemical_unfolding_df(cdAnalyzer,'signal_predicted')
  fig   <- plot_unfolding_fitting(
    df,dfFit,
    reactives$chemicalWorkingUnits,
    input$plot_width_chem, input$plot_height_chem, 
    input$plot_type_chem, input$plot_axis_size_chem,
    reactives$fitted_coefficients_method_chemical)
  
  return(fig)
  
})

output$chemResidualsSVDCoefficients <- renderPlot({
  
  req(reactives$chemical_data_was_fitted_svd_or_pca)
  
  df      <- generate_chemical_unfolding_df(cdAnalyzer)
  dfFit   <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='signal_predicted')
  
  tog           <- inner_join(df,dfFit,by=c('wavelength','chem_conc','legend'),relationship = "many-to-many") 
  tog$residuals <- tog$value.y - tog$value.x
  
  fig <- plot_residuals(
    tog,
    input$plot_axis_size_chem,
    svd_or_pca_based=TRUE,xlab="[Denaturant agent] (M)")
  
  return(fig)
  
})

output$fittedParams_chemicalSVD <- renderTable({
  
  req(reactives$chemical_data_was_fitted_svd_or_pca)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,'Chemical')

  return(df)
  
})

output$fittedErrors_chemicalSVD <- renderTable({
  
  req(reactives$chemical_data_was_fitted_svd_or_pca)
  
  df <- get_fitted_params_unfolding(cdAnalyzer,'Chemical',errors=TRUE)
  
  return(df)
  
})


