# Load example dataset
observeEvent(input$GoLoadExample,{
  dsf$load_example_data()
  updateSelectInput(session, "which",choices  = c("350nm"))
  updateSliderInput(session, "sg_range", value = c(30,50),
                    min = 30, max = 50, step = 1)
  
  tables <- get_renderRHandsontable_list_filled(
    dsf$conditions_original,global_n_rows_conditions_table,
    dsf$conditions_original,rep(T,length(dsf$conditions_original)))
  
  output$table4 <- tables[[4]]
  output$table3 <- tables[[3]]
  output$table2 <- tables[[2]]
  output$table1 <- tables[[1]]
  
  updateSelectInput(
    session, "conc_units",choices  = c(
                                       "Micromolar" = "Micromolar",
                                       "Molar"="Molar",
                                       "Milimolar" = "Milimolar",
                                       "Nanomolar"="Nanomolar"))
  
  updateNumericInput(session, "protein_conc", value = 12)
  
})

# load data into the DSF_binding() when the file is uploaded
#observeEvent(input$Go,{
observeEvent(input$FLf,{ 
  
  output$table1 <- output$table2 <- output$table3 <- output$table4 <- NULL
  
  withBusyIndicatorServer("Go",{
    
    #system(paste0("rm -f ",users_dir,total_folders,"/*"))
    
    #check that we have the input file
    if (!(is.null(input$FLf))) {
      
      fileExtension <- getFileNameExtension(input$FLf$datapath)
      newFileName   <- paste0("0.",fileExtension)
      file.copy(input$FLf$datapath,newFileName,overwrite=TRUE)
      
      if (fileExtension == "json") {dsf$load_JSON_file(newFileName)}
      
      if ((fileExtension == "zip")) {
        
        system(paste0("rm -f *xlsx"))
        unzip(newFileName)
        xlsx_files <- list.files(".",pattern = "xlsx",recursive=T)
        
        dsf_objects_from_xlsx_files <- dsf_objects_from_xlsx_files(xlsx_files)
        
        dsf_objects   <- dsf_objects_from_xlsx_files$dsf_objects
        signal_keys   <- dsf_objects_from_xlsx_files$signal_keys
        
        signal_values <- c()
        temp_values   <- c()
        
        for (signal in signal_keys) {
          mergedSignal  <- get_merged_signal_dsf(dsf_objects,signal)
          signal_values <- c(signal_values,np_array(mergedSignal$signal))
          temp_values   <- c(temp_values,np_array(mergedSignal$temp))
        }
        
        signal_data_dictionary <- py_dict(signal_keys, signal_values, convert = F)
        temp_data_dictionary   <- py_dict(signal_keys, temp_values, convert = F)
        
        conditions_original    <- mergedSignal$conditions_ori
        
        dsf$signals                <- signal_keys
        dsf$signal_data_dictionary <- signal_data_dictionary
        dsf$temp_data_dictionary   <- temp_data_dictionary
        dsf$conditions_original    <- conditions_original
        
        shinyalert(text=paste("The capillaries were loaded in the following order:",
                         paste0(conditions_original,collapse = ' ')),
                   imageWidth = 800,imageHeight = 400,closeOnClickOutside=TRUE,closeOnEsc=TRUE,
                   type="info")
        
      }
      
      if (fileExtension == "txt") {

        fileType <- detect_txt_file_type(newFileName)
        
        if (fileType == 'MX3005P')     dsf$load_Agilents_MX3005P_qPCR_txt(newFileName)
        if (fileType == 'QuantStudio') dsf$load_QuantStudio_txt(newFileName)
        
      }
      
      if (fileExtension == "csv") {dsf$load_csv_file(newFileName)}
      
      if (fileExtension == "xlsx" | fileExtension == "xls") {
        
        # Get file type: DSF or nDSF
        sheet_names <- get_sheet_names_of_xlsx(newFileName)
        # ... Load the data to the python class ...
        if ("RFU" %in% sheet_names)  {
          dsf$load_Thermofluor_xlsx(newFileName)
        } else if ("Data Export" %in% sheet_names || "melting-scan" %in% sheet_names) {
          dsf$load_panta_xlsx(newFileName)
        } else if ("Profiles_raw" %in% sheet_names) {
          dsf$load_tycho_xlsx(newFileName)
        } else {
          dsf$load_nanoDSF_xlsx(newFileName,sheet_names)
        }
        
      }
      
      conditions_original <- dsf$conditions_original
      
      if (fileExtension == "json") {
        maxTpre <- ceiling(max(dsf$temps))
        minTpre <- floor(min(dsf$temps))
        dsf$set_signal(dsf$signal_type)
        
      } else {
        dsf$set_signal(dsf$signals[1])
      }
      
      maxT <- ceiling(max(dsf$temps))
      minT <- floor(min(dsf$temps))
      
      updateSelectInput(session, "which",choices  = dsf$signals)
      
      if (fileExtension == "json") {
        updateSliderInput(session, "sg_range", value = c(minTpre,maxTpre),min = minT, max = maxT, step = 1)
        # Load concentrations in micromolar
        tables <- get_renderRHandsontable_list_filled(
          conditions_original,global_n_rows_conditions_table,
          dsf$concentration_vector*1e6,dsf$include_vector)
        
        updateSelectInput(session, "conc_units",selected='Micromolar')
        
        # Load previous data into the tables and plots -fluorescence fit
        if (!is.null(dsf$fit_fluo_params)) {
          
          updateSelectInput(session, "model_selected",selected=dsf$fitting_option)
          updateNumericInput(session,'delta_cp',value=dsf$cp)
          
          output$params_table <- renderTable({
            get_params_df(dsf$fit_fluo_params,dsf$concentrations)
            
            })
          
          output$params_table_errors <- renderTable({
            get_params_errors_df(dsf$fit_fluo_params,
                                 dsf$concentrations,
                                 dsf$fit_fluo_errs)
          })
          
          output$fluo_fit_plot <- renderPlot({
            
            req(input$select_fitting_plot)
            
            selected <- get_selected_from_choice_label(input$select_fitting_plot,global_chunck_n)
            
            fluo_m      <- make_list_df4plot(dsf$fluo,dsf$concentrations,dsf$temps,global_chunck_n)
            fluo_m_pred <- make_list_df4plot(dsf$fit_fluo_pred,dsf$concentrations,dsf$temps,global_chunck_n)
            
            return( plot_fluorescence_fit(fluo_m,fluo_m_pred,selected) )
          })
          
        }
        
        # Load previous data into the tables and plots - isothermal fit
        if (!is.null(dsf$isothermal_fits)) {
          
          updateSliderInput(session, "if_range",value = c(min(dsf$isothermal_ts),max(dsf$isothermal_ts)))
          updateNumericInput(session,'protein_conc',value=round(dsf$pconc*1e6,2))
          
          output$isf <- renderPlotly({
            
            ist_df_data <- format_ist_data_exp_and_pred(dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts,
                                                        dsf$kd_models,dsf$kd_models_lower,dsf$kd_models_upper,
                                                        dsf$kd_model_conc,dsf$bind_params,dsf$bind_errors,
                                                        dsf$bind_ci95_asymmetric_low, dsf$bind_ci95_asymmetric_up)
            
            return( plot_isothermal_fitting(ist_df_data$data_exp,ist_df_data$data_pred,
                                            input$plot_width_isf, input$plot_height_isf, 
                                            input$plot_type_isf,input$plot_axis_size_isf,
                                            input$logScaleType_isf) )
          }) 
          
          }
          
          if (!is.null(dsf$tms_fit_pred)) {
            
            updateSelectInput(session, "model_selected_tm_shift",selected=dsf$tms_fit_model)
            
            output$tm_shift_fitting <- renderPlotly({
              
              p <- plot_tm_shift_fit(dsf$concentrations,dsf$tmsFromDerivative,
                                     data.frame(dsf$tms_fit_pred), # list to df
                                     data.frame(dsf$tms_fit_info), # list to df
                                     dsf$tms_fit_asymmetricCI95,
                                     input$plot_width_tmShift, input$plot_height_tmShift, 
                                     input$plot_type_tmShift,  input$plot_axis_size_tmShift,
                                     input$logScaleType_tmShift)
              return( p )
            })
            
            output$params_table_tm_shift <- renderTable({
              
              df <- format_tm_fit_df(data.frame(dsf$tms_fit_info)[,1:3],
                                     dsf$tms_fit_model,
                                     dsf$tms_fit_asymmetricCI95)
              
              return(df)
            })
            
          }
          
      } else {
        updateSliderInput(session, "sg_range", value = c(minT,maxT),min = minT, max = maxT, step = 1)
        tables <- get_renderRHandsontable_list(conditions_original,global_n_rows_conditions_table)
        
      }
      
      output$table4 <- tables[[4]]
      output$table3 <- tables[[3]]
      output$table2 <- tables[[2]]
      output$table1 <- tables[[1]]
      
    }})
},priority = 10)

# Autocomplete the concentration versus position table
observe({
  req(input$table1)
  if (input$fill_table) {
    tables <- get_renderRHandsontable_list_autofill(
      dsf$conditions_original,global_n_rows_conditions_table,
      input$initial_ligand,input$n_replicates,input$dil_factor,input$rev_order)
    
    output$table4 <- tables[[4]]
    output$table3 <- tables[[3]]
    output$table2 <- tables[[2]]
    output$table1 <- tables[[1]]
  }
})

# reactive expression to select the fluorescence and concentration data
modify_fluo_temp_cond <- reactive({
  
  condition_include_list <- get_include_vector(input$table1,input$table2,input$table3,input$table4,
                                               dsf$conditions_original,global_n_rows_conditions_table,input$conc_units)
  
  concentration_vector   <- condition_include_list$concentration_vector
  include_vector         <- as.logical(condition_include_list$include_vector)
 
  dsf$set_signal(input$which)
  
  # Store original data - useful to save and load the data using JSON
  dsf$include_vector       <- include_vector
  dsf$concentration_vector <- concentration_vector
  
  # Get signal window range
  sg_range_min_celsius <- input$sg_range[1]
  sg_range_max_celsius <- input$sg_range[2]
  
  left_bound    <-   max(sg_range_min_celsius,min(dsf$temps))
  right_bound   <-   min(sg_range_max_celsius,max(dsf$temps))
  
  # ... Modify in place the python class fluorescence signal according to the selected signal window range ...

  dsf$fluo  <- filter_fluo_by_temp(dsf$fluo,dsf$temps,left_bound,right_bound)
  dsf$temps <- filter_temp_by_temp(dsf$temps,left_bound,right_bound)
  
  
  
  # ... use only the conditions selected by the user ...
  dsf$concentrations <- np_array(concentration_vector[include_vector])
  
  dsf$fluo                    <- dsf$fluo[,include_vector]
  
  median_filter <- get_median_filter(input$median_filter)
  if (median_filter > 0) {dsf$median_filter(median_filter)}
  
  dsf$estimate_fluo_derivates(10) 
  
  updateSelectInput(session, "select_fitting_plot",
                    choices  = get_choices_fluo_fits(ncol(dsf$fluo),global_chunck_n))
  
  return(condition_include_list$concentration_vector)
  
})
