reactives <- reactiveValues(data_loaded=FALSE,using_simulated_data=FALSE,
                            logScaleType_fit_values=c("micromolar"),is_csv=FALSE,
                            dataWasFitted=FALSE)

# Make available the value of reactives$data_loaded to the code in ui.R
# useful to have conditional panels
output$data_loaded <- reactive({reactives$data_loaded})
outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

output$is_csv <- reactive({reactives$is_csv})
outputOptions(output, "is_csv", suspendWhenHidden = FALSE)

# Load example dataset
observeEvent(input$GoLoadExample,{
  
  reactives$is_csv        <- FALSE
  reactives$dataWasFitted <- FALSE
  
  mst$load_example_data()

  mst$conditions_original <- mst$concs
  
  mst$set_signal("Raw Fluorescence")
  updateSliderInput(session,"cold_range",NULL,min = -4, max = 0,value = c(-1,0))
  updateSliderInput(session,"hot_range" ,NULL,min = 0, max = 22,value = c(0,1))
  
  updateSelectInput(session, "conc_units",
                    choices  = c(
                      "Micromolar" = "micromolar","Nanomolar"="nanomolar","Molar"="molar","Milimolar" = "milimolar"))
  
  updateNumericInput(session, "protein_conc", value = 25)
  updateSelectInput(session, "prot_conc_units", choices  = c(
    "Nanomolar"="Nanomolar","Micromolar" = "Micromolar"))

  tables <- get_renderRHandsontable_list(mst$concs)
  
  output$table4 <- tables[[4]]
  output$table3 <- tables[[3]]
  output$table2 <- tables[[2]]
  output$table1 <- tables[[1]]
  
  Sys.sleep(1)
  reactives$using_simulated_data <- TRUE
  reactives$data_loaded          <- TRUE

})

# load data into the MST_fit() when the user uploads data
observeEvent(input$FLf,{

  req(input$FLf)
  
  reactives$data_loaded   <- FALSE
  reactives$dataWasFitted <- FALSE
  
  output$table1 <- output$table2 <- output$table3  <- output$table4             <- NULL
  output$experiments2fit1 <- output$experiments2fit2 <- output$experiments2fit3 <- NULL
  
  if (!(is.null(input$FLf))) {
    fileExt <- getFileNameExtension(input$FLf$datapath)
  } 
    
  withBusyIndicatorServer("Go",{
    
    if (fileExt == "csv") {
      
      reactives$is_csv               <- TRUE
      
      file.copy(input$FLf$datapath,"0.csv",overwrite=TRUE)
      
      sep <- getSepCharacter("./0.csv")
      if (is.null(sep)) return(NULL) # No delimiter found
      
      header <- csvHasHeader("./0.csv",sep)

      mst$load_MST_csv("./0.csv",sep,header)
      mst$set_signal("Raw Fluorescence")
      
    } else {
      
      reactives$is_csv               <- FALSE
      
      if (fileExt == "zip") {
        
        system(paste0("rm -f *xlsx"))
        file.copy(input$FLf$datapath,"0.zip",overwrite=TRUE)
        unzip("0.zip")
        
        xlsx_filesOri   <- list.files(".",pattern = "xlsx")
        # Avoid μ in file names breaking our python code :)
        xlsx_files      <- gsub("μ", "micro", xlsx_filesOri)
        file.rename(xlsx_filesOri,xlsx_files)
        
        mst_objects <- mst_objects_from_xlsx_files(xlsx_files)
        
        merged      <- get_merged_signal_mst(mst_objects,xlsx_files)
       
        signal_keys   <- time_keys    <- c("Raw Fluorescence")
        signal_values <- c(np_array(merged$signal))
        time_values   <- c(np_array(merged$times))
                         
        signal_data_dictionary <- py_dict(signal_keys, signal_values, convert = F)
        time_data_dictionary   <- py_dict(time_keys, time_values,     convert = F)
        concs                  <- as.numeric(merged$conditions_ori)

        mst$signal_data_dictionary <- signal_data_dictionary
        mst$time_data_dictionary   <- time_data_dictionary
        mst$concs                  <- concs
        mst$protConc               <- as.numeric(merged$protConcVec)
        
        # Assing one experimental ID per dataset
        n_ids                        <- lapply(mst_objects, function (x) length(x$experimentID))
        
        experimentID <- c()
        
        for (i in 1:length(n_ids)) {
          
          fname <- basename(merged$file_order[i]) # Remove the path of the file
          fname <- sub('\\..[^\\.]*$', '', fname) # Remove the file extension
            
          experimentID <- c(experimentID,rep(fname,n_ids[i]))
        }
        
        mst$experimentID           <- experimentID

      } else {
        
        file.copy(input$FLf$datapath,"0.xlsx",overwrite=TRUE)
        mst$load_MST_xlsx("./0.xlsx")
        
      }
      
      mst$set_signal("Raw Fluorescence")
      maxT <- ceiling(max(mst$times))
      minT <- floor(min(mst$times))
      updateSliderInput(session,"cold_range",NULL,min = minT, max = 0,value = c(-1,0))
      updateSliderInput(session,"hot_range" ,NULL,min = 0, max = maxT,value = c(0,1))
      
    }
    
    # Automatic detection & transformation to micromolar
    
    # Highest concentration check
    c1 <- max(mst$concs[mst$concs != 0]) < 1e-2 
    # Lowest concentration check
    c2 <- min(mst$concs[mst$concs != 0]) < 1e-7
    
    if (c1 | c2) {
      updateSelectInput(
        session, "conc_units",
        choices  = c("Micromolar" = "micromolar","Nanomolar"="nanomolar",
                     "Molar"="molar","Milimolar" = "milimolar"))
      
      mst$concs    <- mst$concs*1e6
      mst$protConc <- mst$protConc*1e6
    }
    
    mst$conditions_original <- mst$concs
    tables                  <- get_renderRHandsontable_list(mst$concs,mst$protConc,
                                                            mst$experimentID)
    
    output$table4 <- tables[[4]]
    output$table3 <- tables[[3]]
    output$table2 <- tables[[2]]
    output$table1 <- tables[[1]]
    
    Sys.sleep(0.1)
    
    # Check if protein concentration is a vector of ones
    protInfo  <- !( sum(mst$protConc)            == length(mst$protConc) )
    expIDinfo <- !( sum(mst$experimentID == "A") == length(mst$experimentID) )
    dataLoadedMessage(protInfo,expIDinfo)
      
    reactives$using_simulated_data <- FALSE
    reactives$data_loaded          <- TRUE
    
  })
},priority = 10)

# Autocomplete the concentration versus position table
observe({
  req(input$table1)
  
  if (!(reactives$data_loaded)) {return(NULL)}
  
  if (input$fill_table) {
    
    # Quick check to verify that all inputs are numeric
    wrongInput  <- sapply(c(input$initial_ligand,input$n_replicates,input$dil_factor), is.na)
    if (any(wrongInput)) {return(NULL)}
    
    tables <- get_renderRHandsontable_list_autofill(
      mst$conditions_original,
      input$initial_ligand,input$n_replicates,input$dil_factor,input$rev_order)
    
    output$table4 <- tables[[4]]
    output$table3 <- tables[[3]]
    output$table2 <- tables[[2]]
    output$table1 <- tables[[1]]
  }
  
})

# reactive expression to select the fluorescence and concentration data
modify_fluo_temp_cond <- reactive({
  
  req(reactives$data_loaded) 
  condition_include_list <- get_include_vector(input$table1,input$table2,input$table3,input$table4,
                                               mst$conditions_original,input$conc_units)
  
  concentration_vector   <- condition_include_list$concentration_vector
  protein_conc_vector    <- condition_include_list$protConc_vector
  include_vector         <- as.logical(condition_include_list$include_vector)
  expID_vector           <- condition_include_list$expID_vector
  
  # Reset selection of signal
  mst$set_signal("Raw Fluorescence")
  
  # Quick check - the number of conditions should be the same as the length of
  # the boolean vector that containts which should be selected
  if (length(include_vector) != dim(mst$signal)[2]) {return(NULL)}

  # Avoid using only one data point! - somehow it breaks the python code :/
  if (sum(include_vector) < 2) return(NULL)

  # ... use only the conditions selected by the user ...
  mst$concs                 <-  concentration_vector[include_vector]
  mst$signal                <-  mst$signal[,include_vector]
  mst$expID_vector          <-  expID_vector[include_vector]        # Split experiments according to ID 
  mst$protConc              <-  protein_conc_vector[include_vector]
  
  if (!reactives$is_csv) {
    median_filter <- get_median_filter(input$median_filter)
    if (median_filter > 0) mst$median_filter(median_filter)
    
    if (input$normalization_type == "Divide_by_init_value") {
      mst$signal <- normalize_with_cold_region_mean(mst$signal,mst$times) 
    }
    
    mst$get_cold_fluo(input$cold_range[1],input$cold_range[2])
    mst$get_fnorm(input$hot_range[1],input$hot_range[2])
    
  } else {
    # Trick to select the signal if a csv was loaded (no time data!)
    mst$get_cold_fluo(-1,0) 
    mst$F_norm <- mst$F_cold
  }
  
  # Load the tables that contain which conditions should be fitted
  tables <- get_renderRHandsontable_listFit(unique(mst$expID_vector))
  output$experiments2fit1 <- tables[[1]]
  output$experiments2fit2 <- tables[[2]]
  output$experiments2fit3 <- tables[[3]]
  
  return( condition_include_list$concentration_vector )
  
})

observe({
  
  req(input$table1)
  req(reactives$data_loaded)
  
  if (reactives$using_simulated_data) {
    
    updateSelectInput(session, "signal_selected",
                      choices  = c( "F_hot / F_cold  "       = "f_norm",
                                    "Initial Fluorescence"   = "init_fluo"))
    
  } else {
    
    is_csv <- getFileNameExtension(input$FLf$datapath) == "csv"
    
    if (is_csv) {
      updateSelectInput(session, "signal_selected",
                        choices  = c("Initial Fluorescence"   = "init_fluo"))
    } else {
      updateSelectInput(session, "signal_selected",
                        choices  = c( "F_hot / F_cold  "       = "f_norm",
                                      "Initial Fluorescence"   = "init_fluo"))
    }
    
  }
  
}) 

output$ui_signal_tab_box <- renderUI({

  req(input$table1)
  # re evaluate expression if the user loaded a new file

  if (!(reactives$data_loaded)) {return(NULL)}
  
  if (reactives$using_simulated_data) {
    
    tabBox(title = "", width = 10,id = "tabset2",
           
           tabPanel("Signal",withSpinner(plotlyOutput("signal"))),
           tabPanel("Ligand Fluorescence",withSpinner(plotlyOutput("lig_fluo"))),
           tabPanel("Initial Fluorescence vs [Ligand]",withSpinner(plotlyOutput("cold_signal"))),
           tabPanel("F_hot / F_cold vs [Ligand]",withSpinner(plotlyOutput("fnorm_signal")))
    )
    
  } else {
    
    is_csv <- getFileNameExtension(input$FLf$datapath) == "csv"  
    
    if (is_csv) {
      
      tabBox(title = "", width = 10,id = "tabset1",
             
             tabPanel("Initial Fluorescence vs [Ligand]",withSpinner(plotlyOutput("cold_signal"))),
             tabPanel("Ligand Fluorescence",withSpinner(plotlyOutput("lig_fluo")))
             
      ) 
    } else {
      
      tabBox(title = "", width = 10,id = "tabset2",
             
             tabPanel("Signal",withSpinner(plotlyOutput("signal"))),
             tabPanel("Ligand Fluorescence",withSpinner(plotlyOutput("lig_fluo"))),
             tabPanel("Initial Fluorescence vs [Ligand]",withSpinner(plotlyOutput("cold_signal"))),
             tabPanel("F_hot / F_cold vs [Ligand]",withSpinner(plotlyOutput("fnorm_signal")))
      )
    }
    
  }
  
})

# Render signal plot
output$signal <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  fluo_m <- make_df4plot(mst$signal,mst$concs,mst$times)

  if (!(is.null(fluo_m))) {
    return( plot_fluo_signal(fluo_m,input$cold_range[1],input$cold_range[2],
                             input$hot_range[1],input$hot_range[2],"Fluorescence",
                             input$plot_width, input$plot_height, 
                             input$plot_type,input$plot_axis_size))
  }
  
  return(NULL)
})

output$lig_fluo <- renderPlotly({
  
  req(input$table1)

  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  p <- plot_lig_fluo(mst$F_cold,mst$concs,mst$expID_vector,
                     input$plot_width, input$plot_height, 
                     input$plot_type,input$plot_axis_size,
                     input$plot_legend_text_size)
  
  return(p)
  
})

output$cold_signal <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  p <- plot_fluo_vs_lig(mst$F_cold,mst$concs,mst$expID_vector,"cold",
                        input$plot_width, input$plot_height, 
                        input$plot_type,input$plot_axis_size,
                        input$plot_legend_text_size,input$logScaleType,F)
  return(p)
  
})

output$fnorm_signal <- renderPlotly({
  
  req(input$table1)
  if (is.null(modify_fluo_temp_cond())) {return(NULL)}
  
  p <- plot_fluo_vs_lig(mst$F_norm,mst$concs,mst$expID_vector,"f_norm",
                        input$plot_width, input$plot_height, 
                        input$plot_type,input$plot_axis_size,
                        input$plot_legend_text_size,
                        input$logScaleType,F)
  
  return(p)
  
})