reactives <- reactiveValues(data_loaded=FALSE,report_was_created=FALSE,
                            global_max_conditions          = 384,
                            global_n_rows_conditions_table = 96)

observeEvent(input$FLf,{

  reactives$data_loaded             <- FALSE
  output$table1 <- output$table2 <- output$table3 <- output$table4 <- NULL
  
  withBusyIndicatorServer("Go",{
    
    system(paste0("rm -f ",users_dir,total_folders,"/*"))
    
    # Check that we have the xlsx file
    if (!(is.null(input$FLf))) {
      
      fileExtension <- getFileNameExtension(input$FLf$datapath)
      if ((fileExtension == "zip")) {
        
        system(paste0("rm -f *xlsx"))
        file.copy(input$FLf$datapath,"0.zip",overwrite=TRUE)
        unzip("0.zip")
        xlsx_files <- list.files(".",pattern = "xlsx")
        
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
        
        # We assign default time such that the scan rate is 1 degree per minute
        conditions             <- mergedSignal$conditions
        conditions_original    <- mergedSignal$conditions_ori
        
        dsf$signal_data_dictionary <- signal_data_dictionary
        dsf$temp_data_dictionary   <- temp_data_dictionary
        dsf$conditions_original    <- conditions_original
        dsf$signals                <- signal_keys
        dsf$conditions             <- conditions

      }

      newFileName <- paste0("0.",fileExtension)
      if (fileExtension == "txt") {
        # Copy txt file to the current folder
        file.copy(input$FLf$datapath,newFileName,overwrite=TRUE)
        
        fileType <- detect_txt_file_type(newFileName)
        
        if (fileType == 'MX3005P')     dsf$load_Agilents_MX3005P_qPCR_txt(newFileName)
        if (fileType == 'QuantStudio') dsf$load_QuantStudio_txt(newFileName)
        
      }
      
      if (fileExtension == "csv") {
        # Copy txt file to the current folder
        file.copy(input$FLf$datapath,newFileName,overwrite=TRUE)
        dsf$load_csv_file(newFileName)
        
      }
      
      if (fileExtension == "xlsx" | fileExtension == "xls") {
        
        # Copy xlsx file to the current folder
        file.copy(input$FLf$datapath,newFileName,overwrite=TRUE)
        
        # Get file type: DSF or nDSF
        sheet_names <- get_sheet_names_of_xlsx(newFileName)
        # ... Load the data to the python class ...
        if ("RFU" %in% sheet_names)  {
          dsf$load_Thermofluor_xlsx(newFileName)
        } else if ("Data Export" %in% sheet_names || "melting-scan" %in% sheet_names) {
          dsf$load_panta_xlsx(newFileName)
        } else if ("Profiles_raw" %in% sheet_names) {
          dsf$load_tycho_xlsx(newFileName)
        }else {
          dsf$load_nanoDSF_xlsx(newFileName,sheet_names)
        }
      
      }
      
      conditions <- dsf$conditions
      dsf$set_signal(dsf$signals[1])
      updateSelectInput(session, "which",choices  = dsf$signals)
      
      nConditions <- length(conditions)

      reactives$global_max_conditions          <- findClosestHigherValue(nConditions)
      reactives$global_n_rows_conditions_table <- reactives$global_max_conditions / 4

      tables <- get_renderRHandsontable_list(conditions,
                                             reactives$global_n_rows_conditions_table)
      
      output$table4 <- tables[[4]]
      output$table3 <- tables[[3]]
      output$table2 <- tables[[2]]
      output$table1 <- tables[[1]]
      
      min_temp <- round(min(dsf$temps) - 273.15) # To degree celsius
      max_temp <- round(max(dsf$temps) - 273.15) # To degree celsius
      
      updateSliderInput(session,"sg_range",NULL,min = min_temp, max = max_temp,value = c(min_temp+3,max_temp-3))
      
      reactives$data_loaded             <- TRUE
      
    }})
},priority = 10)

output$data_loaded             <- reactive({
  return(reactives$data_loaded)})

output$report_was_created             <- reactive({
  return(reactives$report_was_created)})

outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
outputOptions(output, "report_was_created", suspendWhenHidden = FALSE)

observe({
  
  req(input$table1)
  req(reactives$data_loaded)
  
  uniqueSeries <- unique(
    get_include_vector(input$table1,input$table2,input$table3,input$table4,
                       length(dsf$conditions_original),
                       reactives$global_n_rows_conditions_table,
                       reactives$global_max_conditions)$series_vector)

  updateSelectInput(session, "selected_cond_series",choices  = c("All",uniqueSeries))
  
})

observeEvent(input$layout_file$datapath,{
  
  req(input$table1)
  if (!(is.null(input$layout_file))) {
    
    conditions <- load_layout(input$layout_file$datapath)
    tot_cond <- length(dsf$conditions)
    dsf$conditions          <- conditions[1:tot_cond]
    dsf$conditions_original <- conditions[1:tot_cond]
    tables <- get_renderRHandsontable_list(conditions,reactives$global_n_rows_conditions_table)
    
    output$table4 <- tables[[4]]
    output$table3 <- tables[[3]]
    output$table2 <- tables[[2]]
    output$table1 <- tables[[1]]
  }
  
})

observeEvent(input$sort_conditions,{
  req(input$table1)
  
  output$table1 <- output$table2 <- output$table3 <- output$table4 <- NULL
  dsf$sort_by_conditions_name(input$sort_conditions)
  tables <- get_renderRHandsontable_list(dsf$conditions_original,
                                         reactives$global_n_rows_conditions_table)
  
  output$table4 <- tables[[4]]
  output$table3 <- tables[[3]]
  output$table2 <- tables[[2]]
  output$table1 <- tables[[1]]
  
})

modify_fluo_temp_cond <- reactive({
  
  nconditions <- length(dsf$conditions_original)
  condition_include_list <- get_include_vector(input$table1,input$table2,input$table3,input$table4,
                                               nconditions,reactives$global_n_rows_conditions_table,
                                               reactives$global_max_conditions)
  
  conditions_vector      <- as.character(condition_include_list$conditions_vector)
  include_vector         <- as.logical(condition_include_list$include_vector)
  series_vector          <- condition_include_list$series_vector
  
  dsf$set_signal(input$which)
  
  # Get signal window range
  sg_range_min_celsius <- input$sg_range[1]+273.15
  sg_range_max_celsius <- input$sg_range[2]+273.15
  
  left_bound    <-   max(sg_range_min_celsius,min(dsf$temps))
  right_bound   <-   min(sg_range_max_celsius,max(dsf$temps))
  
  # ... Modify in place the python class fluorescence signal according to the selected signal window range ...
  dsf$fluo  <- filter_fluo_by_temp(dsf$fluo,dsf$temps,left_bound,right_bound)
  dsf$temps <- filter_temp_by_temp(dsf$temps,left_bound,right_bound)
  
  if (input$selected_cond_series != "All") {
    include_vector         <- include_vector & (series_vector == input$selected_cond_series)
  }
  
  # ... use only the conditions selected by the user ...
  dsf$conditions         <- conditions_vector[include_vector]
  
  dsf$select_signal_columns(c(include_vector))
  
  median_filter <- get_median_filter(input$median_filter)
  if (median_filter > 0) {dsf$median_filter(median_filter)}
  
  dsf$fluo <- normalize_fluo_matrix_by_option(input$normalization_type,dsf$fluo,dsf$temps)
  
  dsf$estimate_fluo_derivates(input$SG_window2)

  return(condition_include_list$conditions_vector)
  
})

# Render signal plot
output$signal <- renderPlotly({
  
  req(input$table1)
  req(reactives$data_loaded)
  
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  fluo_m <- make_df4plot(dsf$fluo,dsf$conditions,
                         dsf$temps)
  
  if (!(is.null(fluo_m))) {
    p <- plot_fluo_signal(fluo_m,dsf$signal_type, 
                          input$plot_width, input$plot_height, 
                          input$plot_type,input$plot_font_size,input$plot_axis_size)
    return(p)
  }
  return(NULL)
  
}
)

# Render 1st derivative plot
output$signal_der1 <- renderPlotly({
  
  req(input$table1)
  req(reactives$data_loaded)
  
  if (length(dsf$conditions) > reactives$global_max_conditions - reactives$global_n_rows_conditions_table*1) {req(input$table4)}
  
  modify_fluo_temp_cond()
  
  fluo_m <- make_df4plot(dsf$derivative,dsf$conditions,dsf$temps)
  if (!(is.null(fluo_m))) {
    p <- plot_fluo_signal(fluo_m,"First derivative",
                          input$plot_width, input$plot_height, 
                          input$plot_type,input$plot_font_size,input$plot_axis_size)
    return(p)
  }
  return(NULL)
  
}
)

# Render 2nd derivative plot
output$signal_der2 <- renderPlotly({
  
  #req(fluo_signal_loaded())
  req(input$table1)
  req(reactives$data_loaded)
  
  modify_fluo_temp_cond()
  
  fluo_m <- make_df4plot(dsf$derivative2,dsf$conditions,dsf$temps)
  
  if (!(is.null(fluo_m))) {
    p <- plot_fluo_signal(fluo_m,"Second derivative",
                          input$plot_width, input$plot_height, 
                          input$plot_type,input$plot_font_size,input$plot_axis_size)
    return(p)
  }
  return(NULL)
  
}
)

# Render maximum of derivative plot
output$tm_derivative <- renderPlotly({
  
  #req(fluo_signal_loaded())
  req(input$table1)
  req(reactives$data_loaded)
  
  modify_fluo_temp_cond()
  
  p <- generate_max_der_plot(dsf$tms_from_deriv,dsf$conditions,
                             input$plot_width, input$plot_height, 
                             input$plot_type,input$plot_font_size,input$plot_axis_size)
  
  return(p)
  
}
)

# Fit when the user presses the button
fluo_fit_data <- eventReactive(input$btn_cal, {
  
  #req(fluo_signal_loaded())
  req(input$table1)
  # check we have data to fit
  if (ncol(dsf$fluo)>0)   {
    
    dsf$estimate_baselines_parameters(input$temp_range_baseline_estimation)
    
    withBusyIndicatorServer("btn_cal",{
      # ... Fit according to selection ...
      if (input$model_selected == "EquilibriumTwoState"  )   {
        dsf$cp <- input$delta_cp * 4184 # convert units
        dsf$EquilibriumTwoState()}
      
      if (input$model_selected == "EquilibriumThreeState")   {
        dsf$EquilibriumThreeState(input$t1min+273.15,input$t1max+273.15,
                                  input$t2min+273.15,input$t2max+273.15) # To kelvin 
        }
      
      if (input$model_selected == "EmpiricalTwoState"    )   {dsf$EmpiricalTwoState()     }
      
      if (input$model_selected == "EmpiricalThreeState"  )   {
        dsf$EmpiricalThreeState(input$t1min+273.15,input$t1max+273.15,
                                input$t2min+273.15,input$t2max+273.15)   
        }
      
      if (input$model_selected == "IrreversibleTwoState" )   {
        dsf$scan_rate <- input$scan_rate
        dsf$IrreversibleTwoState()  }
    })
    
    fluo_m           <- make_list_df4plot(dsf$fitted_fluo,dsf$fitted_conditions,
                                          dsf$temps,global_chunck_n)
    
    if (length(dsf$fluo_predictions_all) == 0 ) return(NULL)
    
    fluo_pred_matrix <- t(do.call(rbind,dsf$fluo_predictions_all))
    fluo_m_pred      <- make_list_df4plot(fluo_pred_matrix,dsf$fitted_conditions,
                                          dsf$temps,global_chunck_n)
    
    return(list("fluo_fit_real"=fluo_m,"fluo_fit_pred"=fluo_m_pred))
    
  } else {
    return(NULL)
  }
  
})

output$three_state_model_selected             <- reactive({
  return( grepl("Three",input$model_selected) ) 
})

observeEvent(input$model_selected,{
  
  req(input$table1)
  
  if (grepl("Three",input$model_selected)) {
    
    temps   <- dsf$temps - 273.15
    tempMax <- round(max(temps)) 
    tempMin <- round(min(temps))
    tempMiddle <- round( (max(temps) + min(temps)) * 0.5 )
    
    updateNumericInput(session, "t1max", value = tempMiddle,  min = tempMin, max = tempMax)
    updateNumericInput(session, "t1min", value = tempMin+5,   min = tempMin, max = tempMax)
    updateNumericInput(session, "t2max", value = tempMax-5,   min = tempMin, max = tempMax)
    updateNumericInput(session, "t2min", value = tempMiddle,  min = tempMin, max = tempMax)
    
  }
  
  return(NULL)
})

output$data_was_fitted <- reactive({
  req(input$table1)
  # re evaluate expression if the user fitted data again
  req(input$btn_cal)
  return(length(dsf$fitted_conditions)>1)

})

outputOptions(output, "data_was_fitted", suspendWhenHidden = FALSE)
outputOptions(output, "three_state_model_selected", suspendWhenHidden = FALSE)

output$params_table <- renderTable({
  req(fluo_fit_data())
  get_sorted_params_table(dsf$params_all,dsf$fitted_conditions,
                          global_chunck_n,dsf$params_name,input$sort_table_parameter)
})

output$params_table_errors <- renderTable({
  req(fluo_fit_data())
  return(get_sorted_params_table_errors(dsf$errors_percentage_all,
                                        dsf$fitted_conditions,
                                        dsf$params_name,input$sort_table_parameter))
})

observe({
  
  req(fluo_fit_data())
  updateSelectInput(session, "select_fitting_plot",
                    choices  = get_choices_fluo_fits(length(dsf$fitted_conditions),global_chunck_n))
  
})

observe({
  
  req(fluo_fit_data())
  updateSelectInput(session, "sort_table_parameter",
                    choices  = dsf$params_name)
  
})

# Plot the fluorescence fits
output$fluo_fit_plot <- renderPlot({
  
  req(fluo_fit_data())
  req(input$select_fitting_plot)
  
  selected <- get_selected_from_choice_label(input$select_fitting_plot,global_chunck_n)
  
  real_data  <- fluo_fit_data()$fluo_fit_real
  model_data <- fluo_fit_data()$fluo_fit_pred
  
  p <- plot_fluorescence_fit(real_data,model_data,selected)
  return(p)
}
)

# Plot the fluorescence fits standarized residuals
output$fluo_residuals_plot <- renderPlot({
  
  req(fluo_fit_data())
  req(input$select_fitting_plot)
  
  selected <- get_selected_from_choice_label(input$select_fitting_plot,global_chunck_n)
  
  real_data  <- fluo_fit_data()$fluo_fit_real
  model_data <- fluo_fit_data()$fluo_fit_pred
  
  p <- plot_fluorescence_residuals(real_data,model_data,selected,dsf$std_error_estimate_all,
                                   dsf$fitted_conditions)
  return(p)
}
)

# End of Plot the fluorescence fits

output$fitted_conditions_table <- renderTable({
  req(fluo_fit_data())
  return(get_fitted_conditions_table(dsf$conditions,
                                     dsf$fitted_conditions))
})

filter_conditions <- reactive({
  
  req(fluo_fit_data())
  selected_indexes <- get_selected_conditions_index(dsf$fitted_conditions,
                                                    input$sd_factor_bool,input$bs_factor_bool,
                                                    input$far_from_bounds,
                                                    dsf$errors_percentage_all,
                                                    dsf$baseline_factor_all,
                                                    dsf$parameters_far_from_bounds)
  
  return(selected_indexes)
})

get_score_table <- reactive({
  
  # Check we have fitted the fluorescence data
  req(fluo_fit_data())
  
  selected_indexes <- filter_conditions()
  if (!(any(selected_indexes))){return(NULL)}
  
  score_table <- NULL
  if ((dsf$model_name == "EquilibriumTwoState")) {
    dG_std <- sapply(dsf$dG_std, function(x) x*0.000239006) # to kcal
    
    params_all  <- map2(dG_std,dsf$dCp_component,function(x,y) c(x,y))
    params_name <- c("dg_std","cp_comp")
    
    score_table <- get_all_params_df(params_all[selected_indexes],
                                     dsf$fitted_conditions[selected_indexes],
                                     global_chunck_n,params_name)
    
    score_table$dg_std <- as.numeric(score_table$dg_std) 
    score_table        <- score_table %>% dplyr::arrange(dg_std)
    colnames(score_table) <- c("ΔG_unfolding (25°C) (kcal/mol) ","ΔCp_component(K)","Condition")
  } 
  
  if ((dsf$model_name == "EquilibriumThreeState")) {
    
    params_all <- map(dsf$dG_comb_std,function(x) c(x*0.000239006,99)) # to kcal, 99 is a place holder
    params_name <- c("score","temp_value")
    
    score_table <- get_all_params_df(params_all[selected_indexes],
                                     dsf$fitted_conditions[selected_indexes],
                                     global_chunck_n,params_name)
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(score) %>% select(-temp_value)
    colnames(score_table) <- c("ΔG_comb (25°C) (kcal/mol)","Condition")
    
  }
  
  if ((dsf$model_name == "EmpiricalTwoState")) {
    
    params_all <- map(dsf$score,function(x) c(x,99))
    params_name <- c("score","temp_value")
    
    score_table <- get_all_params_df(params_all[selected_indexes],
                                     dsf$fitted_conditions[selected_indexes],
                                     global_chunck_n,params_name)
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(-score) %>% select(-temp_value)
    colnames(score_table) <- c("sqrt( T_onset**2 + Tm**2 )","Condition")
    
  }
  
  if ((dsf$model_name == "EmpiricalThreeState")) {
    
    params_all <- map(dsf$T_eucl_comb,function(x) c(x,99))
    params_name <- c("score","temp_value")
    
    score_table <- get_all_params_df(params_all[selected_indexes],
                                     dsf$fitted_conditions[selected_indexes],
                                     global_chunck_n,params_name)
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(-score) %>% select(-temp_value)
    colnames(score_table) <- c("sqrt( T_onset_1**2 + T1**2 ) + sqrt( T_onset_1**2 + T2**2 )","Condition")
    
  }
  
  if ((dsf$model_name == "IrreversibleTwoState")) {
    
    params_all <- map(dsf$pkd,function(x) c(x,99))
    params_name <- c("score","temp_value")
    
    score_table <- get_all_params_df(params_all[selected_indexes],
                                     dsf$fitted_conditions[selected_indexes],
                                     global_chunck_n,params_name)
    
    score_table$score    <-  as.numeric(score_table$score) 
    score_table          <-  score_table %>% dplyr::arrange(-score) %>% select(-temp_value)
    colnames(score_table) <- c("pkd","Condition")
    
  }
  
  if (is.null(score_table)) {return(NULL)} 
  return(score_table)
  
})

output$score_table <- renderTable({
  
  req(get_score_table())
  return(get_score_table())
  
})

observe({
  
  req(fluo_fit_data())
  updateSelectInput(session, "select_plot_type",
                    choices  = get_choices_result_plot(dsf$model_name))
  
})


get_results_plot <- reactive({
  
  selected_indexes <- filter_conditions()
  
  if (!(any(selected_indexes))){return(NULL)}

  # Get the value of the parameter Tm
  if (input$select_plot_type %in% c("Unfolded fraction","The 25 highest Tms",
                                    "The 25 highest Tms versus Tonset")) {
    
    tm_position <- which(dsf$params_name == "Tm")  
    tms <- sapply(dsf$params_all, function(x) x[tm_position])
    
  }
  
  #Get the value of the parameter T_onset 
  if (input$select_plot_type %in% c("Unfolded fraction","The 25 highest Tms versus Tonset")) {
    
    if (input$model_selected == "EquilibriumTwoState") {t_onset <- dsf$T_onset}
    if (input$model_selected == "EmpiricalTwoState")   {
      t_onset_position <- which(dsf$params_name == "T_onset") 
      t_onset          <- sapply(dsf$params_all, function(x) x[t_onset_position])
    }
  }
  
  if (input$select_plot_type == "Unfolded fraction") {
    
    dh_position <- which(dsf$params_name == "dHm")  
    dhs <- sapply(dsf$params_all, function(x) x[dh_position])
    
    if (input$model_selected == "EquilibriumTwoState") {
      fig <- generate_fractions_plot(dhs[selected_indexes],
                                     tms[selected_indexes],
                                     dsf$fitted_conditions[selected_indexes],
                                     t_onset[selected_indexes], # useful only to set colors
                                     input$plot_width_results, input$plot_height_results, 
                                     input$plot_type_results,input$plot_font_size_results,
                                     input$plot_axis_size_results,FALSE) 
                                                         #set to TRUE  to use a different color for each position  
                                                         #set to FALSE to use a different color for each condition                                      
    } 
    
    return(fig)
  }
  
  if (input$select_plot_type == "The 25 highest Tms") {
    
    fig <- generate_tm_plot(tms[selected_indexes],dsf$fitted_conditions[selected_indexes],
                            input$plot_width_results, input$plot_height_results, 
                            input$plot_type_results,input$plot_axis_size_results)
    return(fig)
  }
  
  if (input$select_plot_type == "The 25 highest Tms versus Tonset") {
    
    fig <- generate_tm_tonset_plot(tms[selected_indexes],
                                   dsf$fitted_conditions[selected_indexes],
                                   t_onset[selected_indexes],
                                   input$plot_width_results, input$plot_height_results, 
                                   input$plot_type_results,input$plot_font_size_results,input$plot_axis_size_results)
    return(fig)
  }
  
  if (input$select_plot_type == "The 25 highest combinations of T1 and T2") {
    
    t1_position <- which(dsf$params_name == "T1") 
    t1s <- sapply(dsf$params_all, function(x) x[t1_position])
    
    t2_position <- which(dsf$params_name == "T2") 
    t2s <- sapply(dsf$params_all, function(x) x[t2_position])
    
    if (input$model_selected == "EquilibriumThreeState") {score <- dsf$dG_comb_std}
    if (input$model_selected == "EmpiricalThreeState")   {score <- dsf$T_eucl_comb}
    
    fig <- generate_2t_plot(t1s[selected_indexes],
                            t2s[selected_indexes],
                            dsf$fitted_conditions[selected_indexes],
                            score[selected_indexes],
                            input$plot_width_results, input$plot_height_results, 
                            input$plot_type_results,input$plot_font_size_results,input$plot_axis_size_results)
    return(fig)
    
  }

  if (input$select_plot_type == "Score versus condition") {
    
    req(get_score_table())
    return(generate_score_plot(get_score_table(),
                               input$plot_width_results,input$plot_height_results, 
                               input$plot_type_results,input$plot_axis_size_results))
    
  }
  
})

output$results_plot <- renderPlotly({
  
  req(fluo_fit_data())
  req(input$select_plot_type)

  plot <- get_results_plot()
  return(plot)
  
})

