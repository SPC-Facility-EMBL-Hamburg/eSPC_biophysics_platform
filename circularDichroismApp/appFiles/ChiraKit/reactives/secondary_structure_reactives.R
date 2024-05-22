# Reactives to process user reference sets
observeEvent(input$matrixF,{
  
  req(input$matrixF)

  matrixF <- read_matrix_file(input$matrixF$datapath)
  
  nrows <- nrow(matrixF)
  
  reactives$secStrRefMatrixF <- np_array(matrixF)
    
  df <- data.frame('Sec_str_element'=paste0('Element',1:nrows))
  
  output$secondary_structure_elements_names <- renderRHandsontable(
    {rhandsontable(df,rowHeaders=NULL)})
  
  append_record_to_logbook(paste0("Loading matrix F to ChiraKit : ",input$matrixF$name))
  
})

# Reactives to process user reference sets
observeEvent(input$matrixC,{
  
  req(input$matrixC)
  
  matrixC <- read_matrix_file(input$matrixC$datapath)
  
  reactives$secStrRefMatrixC <- np_array(matrixC)
  
  append_record_to_logbook(paste0("Loading matrix C to ChiraKit : ",input$matrixC$name))
  
})


# End of Reactives to process user reference sets

observeEvent(input$runSecStrEstimation,{
  
  reactives$secStrFittingWasDone <- FALSE
  
  # Clear previous TabPanels
  for (tabPanelTargetName in reactives$secStrFittingTabsNames) {
    
    removeTab(inputId = "secondary_structure_tabBox", 
              target = tabPanelTargetName)
    
  }
  
  reactives$secStrFittingTabsNames <- c()
  
  withBusyIndicatorServer("hiddenBtnSecStr",{
    
    expNames <- cdAnalyzer$experimentNames
    
    # We need to use lapply instead of a for-loop to avoid R rendering the TabPanels only with the last value
    
    append_record_to_logbook(paste0("Starting the secondary structure estimation algorithm.",
                                    'Lower WL: ',input$lower_wl_secStr))
    
    # use values for default reference sets
    if (input$useDefaultReferenceSet) {
      high_wl <- 240
      step_wl <- 1
      
      append_record_to_logbook(paste0("Using the default reference sets."))
      
    } else {
      
      matC    <- reactives$secStrRefMatrixC
      matF    <- reactives$secStrRefMatrixF
      high_wl <- input$maxWL_refSet
      step_wl <- input$stepWL_refSet
      
      # start of - check validity of the analysis
      
      low_wl  <- high_wl - (nrow(as.matrix(matC)-1) * step_wl)
      
      if (low_wl > 190) {
        
        shinyalert(text = "Please provide an F matrix that goes down to 
                 at least 190 nm.",type = "warning",
                   closeOnEsc = T,closeOnClickOutside = T,html=T)
        return(NULL)
      }
      
      if ((low_wl > 180 | input$lower_wl_secStr > 180 ) & nrow(as.matrix(matF)) > 4) {
        
        shinyalert(text = "Please provide an F matrix with four or less
                 secondary structure elements.",type = "warning",
                   closeOnEsc = T,closeOnClickOutside = T,html=T)
        return(NULL)
      }
      
      if ((low_wl > 170 | input$lower_wl_secStr > 170 ) & nrow(as.matrix(matF)) > 6) {
        
        shinyalert(text = "Please provide an F matrix with six or less
                 secondary structure elements.",type = "warning",
                   closeOnEsc = T,closeOnClickOutside = T,html=T)
        return(NULL)
      }
      
      append_record_to_logbook(paste0("Using the custom reference sets.",
                                      ' High wl: ',high_wl,
                                      ' nm. Step wl: ',step_wl, ' nm.'))
      
      # end of - check validity of the analysis
      
    }
    
    lapply(expNames, function(exp){
      
      check <- cdAnalyzer$experimentsOri[[exp]]$init_and_check_secondary_str_method(input$lower_wl_secStr)
      
      if (!check) return(NULL)
      
      if (input$useDefaultReferenceSet) {
        cdAnalyzer$experimentsOri[[exp]]$set_secondary_structure_method_references_default()
      } else {
        
        ss_lbls <- hot_to_r(input$secondary_structure_elements_names)[,1]
        
        cdAnalyzer$experimentsOri[[exp]]$set_secondary_structure_method_references_user(
          matF,matC,high_wl,step_wl,ss_lbls
        )
        
      }
      
      Sys.sleep(0.1)
      cdAnalyzer$experimentsOri[[exp]]$estimate_secondary_structure()      
      
      secondary_structure_content_lst <- cdAnalyzer$experimentsOri[[exp]]$secondary_structure_content
      
      if (length(secondary_structure_content_lst) > 0) {
        
        spectraNames <- cdAnalyzer$experimentsOri[[exp]]$spectraNames
        
        for (i in 1:length(spectraNames)) {
          
          sec_str_df     <- secondary_structure_content_lst[i]
          
          saneCD_curveName   <- gsub(':','',spectraNames[i])
          tabPanelTargetName <- paste0(saneCD_curveName,'secStrTarget')
          
          tabP <- tabPanel(title=saneCD_curveName,
                           value=tabPanelTargetName,
                           fluidRow(column(12,tableOutput(paste0('secStr_',saneCD_curveName)))))
          
          reactives$secStrFittingTabsNames <- c(
            reactives$secStrFittingTabsNames,tabPanelTargetName)
          
          appendTab("secondary_structure_tabBox",tabP,select=TRUE)
          
          output[[paste0('secStr_',saneCD_curveName)]] <- renderTable({sec_str_df})
          
          Sys.sleep(0.2)
          
        }
      }
      
      return(NULL)
    })
    
    list_of_signals  <- list()
    list_of_fittings <- list()
    
    output$fitted_CD_spectra_Sec_Str <- NULL
    includePlot <- FALSE
    
    for (exp in expNames) {
      
      cdAnalyzer$experimentsOri[[exp]]$query_spectra_sec_str
      
      secondary_structure_content_lst <- cdAnalyzer$experimentsOri[[exp]]$secondary_structure_content
      
      if (length(secondary_structure_content_lst) > 0) {
        
        includePlot      <- TRUE
        list_of_signals  <- append(list_of_signals,list(cdAnalyzer$experimentsOri[[exp]]$query_spectra_sec_str))
        list_of_fittings <- append(list_of_fittings,list(cdAnalyzer$experimentsOri[[exp]]$fitted_spectra_sec_str))
        
      }
      
    }
    
    if (includePlot) {
      
      reactives$secStrFittingWasDone <- TRUE
      
      reactives$list_of_signals  <- list_of_signals
      reactives$list_of_fittings <- list_of_fittings
      
      output$fitted_CD_spectra_Sec_Str <- renderPlotly({
        
        req(reactives$data_loaded)
        
        fig <- plot_fitted_spectra_sec_str(list_of_signals,list_of_fittings,
                                           high_wl,step_wl)
        
        return(fig)
        
      })
      
    }

  })
  
})

observeEvent(input$pdbFiles,{

  cd_data_files   <- input$pdbFiles$datapath
  names           <- input$pdbFiles$name
  
  # Clear previous TabPanels
  for (tabPanelTargetName in reactives$secStrCalculationTabsNames) {
    
    removeTab(inputId = "secondary_structure_calc_tabBox", 
              target = tabPanelTargetName)
    
  }
  
  reactives$secStrCalculationTabsNames <- c()
  
  reactives$secStrCalcWasDone <- FALSE
  
  withBusyIndicatorServer("hiddenBtnPDBfile",{
    
    # We need to use lapply instead of a for-loop to avoid R rendering the TabPanels only with the last value
    
    append_record_to_logbook(paste0("Calculating the secondary structure components of: ",
                                    paste0(names,collapse = '')))
    
    lapply(1:length(cd_data_files), function(i){
      
      f    <- cd_data_files[i]
      name <- remove_file_extension(names[i])
      
      df <- run_dssp_workflow(f)
      
      if (!is.null(df)) {
        
        reactives$secStrCalcWasDone <- TRUE
        
        saneFileName <- paste0('PDB ',gsub(':','',name))
        
        tabPanelTargetName <- paste0(saneFileName,'secStrCalcTarget')
        
        tabP <- tabPanel(title=saneFileName,
                         value=tabPanelTargetName,
                         fluidRow(column(12,tableOutput(paste0('secStrCalc_',saneFileName)))))
        
        reactives$secStrCalculationTabsNames <- c(
          reactives$secStrCalculationTabsNames,tabPanelTargetName)
        
        appendTab("secondary_structure_calc_tabBox",tabP,select=TRUE)
        
        output[[paste0('secStrCalc_',saneFileName)]] <- renderTable({df})
        
        Sys.sleep(0.2)
        
      }
      
      return(NULL)
    })
    
  })
  
},ignoreNULL = TRUE,ignoreInit = TRUE)




