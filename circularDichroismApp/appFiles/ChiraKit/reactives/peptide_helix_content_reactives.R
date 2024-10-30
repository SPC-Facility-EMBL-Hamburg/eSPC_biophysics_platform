# Create the Table to fill with the temperature data
observeEvent(input$btn_create_peptide_dataset,{
  
  req(reactives$data_loaded)
  
  temperatures <- c()
  cd_curves    <- c()
  pep_bonds    <- c()
  mre222       <- c()

  for (exp in cdAnalyzer$experimentNames) {
    
    cdExp <- cdAnalyzer$experimentsOri[[exp]]
    check <- cdExp$init_and_check_helicity_method()
    
    if (!check) next

    cdExp$get_mre_222nm()

    temperatures <- c(temperatures,cdExp$temperature)
    cd_curves    <- c(cd_curves,cdExp$internalID)
    mre222       <- c(mre222,cdExp$MRE_222nm)

    x            <- cdExp$numberOfCroms
    pep_bonds    <- c(pep_bonds,rep(x,length(cdExp$internalID)))
    
  }
  
  df_peptide <- data.frame(cd_curves,temperatures,pep_bonds,mre222) 
  
  if (nrow(df_peptide) == 0 ) return(NULL) 
  
  # Remove the curves that were not selected in the 1. Import data Tab
  legendDf <- getLegendDF(input$legendInfo)
  sel_ids  <- legendDf$Internal.ID[legendDf$Show]
  
  df_peptide <- df_peptide[df_peptide$cd_curves %in% sel_ids,]
  
  # Remove experiments with non-matching units
  id_to_keep         <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
  internalID_all     <- cdAnalyzer$get_experiment_properties('internalID')
  internalID_to_keep <- unlist(internalID_all[id_to_keep])
  
  df_peptide <- df_peptide[df_peptide$cd_curves %in% internalID_to_keep,]
  
  colnames(df_peptide) <- c('CD_curve','Temperature (°C or K)','#Peptide_bonds','MRE_222nm')
  
  append_record_to_logbook(c('Creating a peptide dataset with the following data',df_to_lines(df_peptide)))
  
  # Assign the created dataframe to the Table thermal_denaturation_data (available at the 2a. Thermal analysis Tab)
  output$thermal_peptide_data <- renderRHandsontable({
    rhandsontable(df_peptide,rowHeaders=NULL)    %>% 
      hot_col(col = c(1,3,4), readOnly=TRUE) %>% 
      hot_table(stretchH='all')
  })
  
})

observeEvent(input$btn_fit_helicity,{
  
  df     <- hot_to_r(input$thermal_peptide_data)
  temps  <- df[,2]
  Ys     <- df[,4] 
  nPepBs <- df[,3]
    
  nan_temps <- sum(is.na(temps))
  
  if (nan_temps > 0) {
    
    shinyalert(text = "Please fill the missing temperature data.",
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    
    return(NULL)
  }  
  
  if (max(temps) > 273) temps <- temps - 273.15 # If required, convert from Kelvin to Celsius
  
  results <- list()
  
  for (row in 1:nrow(df)) {
    
    # The variables poly_double_all and poly_total_all are load when sourcing cdAnalyzer.py (in server.R)
    
    Y     <- Ys[row]
    temp  <- temps[row]
    nPepB <- nPepBs[row]
      
    helicity <- run_helicity_estimation(poly_double_all,poly_total_all,
                                        nPepB,temp,Y)
    
    results[[row]] <- data.frame(helicity)
    
  }
  
  results <- do.call(rbind,results)
  
  results$cd_curve <- df[,1]
  
  output[['helicityTable']] <- renderTable({results})
  output$helicityPlot       <- renderPlotly({plot_helicity(temps,results[,1],as.character(nPepBs))})
  
  append_record_to_logbook('Estimating the peptide helicity using the method from U. Zavrtanik et al., 2024')
  
  
})





