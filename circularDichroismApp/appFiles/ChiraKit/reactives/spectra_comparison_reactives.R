# Create the Table to fill with the categorical labels of the spectra
observeEvent(list(input$legendInfo,input$workingUnits),{
  
  req(reactives$data_loaded)
  
  # Retrieve the dataframe from the legends table (1. Load Input Tab)
  legendDf <- getLegendDF(input$legendInfo)
  id  <- legendDf$Internal.ID
  
  # Initialize the dataframe with the CD curves and the labels
  df           <- data.frame(id,'Control')
  colnames(df) <- c('CD_curve','Label')
  
  # Remove non-selected CD curves
  df <- df[legendDf$Show,]
  
  # Remove experiments with non-matching units
  id_to_keep         <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
  internalID_all     <- cdAnalyzer$getExperimentProperties('internalID')
  internalID_to_keep <- unlist(internalID_all[id_to_keep])
  
  df <- df[df$CD_curve %in% internalID_to_keep,]
  
  # Assign the created dataframe to the Table spectra_labels (available at the 2e. Spectra comparison Tab)
  output$spectra_labels <- renderRHandsontable({
    rhandsontable(df,rowHeaders=NULL)    %>% 
      hot_col(col = c(1), readOnly=TRUE) %>% 
      hot_table(stretchH='all')
  })
  
}) 

observeEvent(input$btn_create_compare_dataset,{
  
  req(reactives$data_loaded)
  
  reactives$compareDatasetCreated          <- NULL
  
  #reactives$spectra_was_decomposed_chemical <- NULL
  
  # Retrieve which CD experiments should be used to build a new dataset for the compare spectra workflow
  df_ids2find        <- hot_to_r(input$spectra_labels)
  
  # Stop if there are no CD curves of interest  
  if (nrow(df_ids2find) == 0 ) return(NULL)
  
  append_record_to_logbook(c('Creating a comparison dataset with the following data',df_to_lines(df_ids2find)))
  
  relevantSpectra    <- df_ids2find$CD_curve 
  relevantLabels     <- df_ids2find[,2]
  
  if (length(unique(relevantLabels)) > 6) {
    
    shinyalert(text = 
    paste("<b>Please use six or less categories."),
    type = "warning",closeOnEsc = T,closeOnClickOutside = T,
    html=T)
    
    return(NULL)
  }
  
  merged <- get_signal_dfs_from_selected_spectra(relevantSpectra,cdAnalyzer)
  
  sorted_indexes      <- order(relevantLabels)
  relevantLabels      <- as.character(relevantLabels[sorted_indexes])
  
  sorted_signal   <- as.matrix(merged[,-1][,sorted_indexes],drop = FALSE)
  
  compareSpectraPyClass$wavelength          <- np_array(merged[,1])
  compareSpectraPyClass$signalDesiredUnit   <- np_array(sorted_signal)
  compareSpectraPyClass$labels              <- np_array(relevantLabels)

  compareSpectraPyClass$workingUnits        <- input$workingUnits
  
  compareSpectraPyClass$summarise_signal_per_label()
  compareSpectraPyClass$generate_comparison_labels()
  compareSpectraPyClass$generate_difference_spectra()
  compareSpectraPyClass$find_distances()
  
  Sys.sleep(2)
  reactives$compareDatasetCreated <- TRUE
  
})

output$cdSpectraAvg <- renderPlotly({
  
  req(reactives$data_loaded)
  req(reactives$compareDatasetCreated)

  plot_average_spectra(compareSpectraPyClass$means,compareSpectraPyClass$sds,
                       compareSpectraPyClass$labels_unique,compareSpectraPyClass$wavelength,
                       compareSpectraPyClass$workingUnits,
                       input$plot_width_compare,input$plot_height_compare,
                       input$plot_type_compare,input$plot_axis_size_compare)
  
})

output$cdSpectraDiff <- renderPlotly({
  
  req(reactives$data_loaded)
  req(reactives$compareDatasetCreated)
  
  dS    <- compareSpectraPyClass$difference_spectra
  dSe   <- compareSpectraPyClass$difference_spectra_sd
  dsLbl <- compareSpectraPyClass$difference_spectra_lbl
  
  plot_difference_spectra(dS,dSe,dsLbl,
                          compareSpectraPyClass$wavelength,compareSpectraPyClass$workingUnits,
                          input$plot_width_compare,input$plot_height_compare,
                          input$plot_type_compare,input$plot_axis_size_compare)
  
})

output$cdSpectraDist <- renderPlotly({
  
  req(reactives$data_loaded)
  req(reactives$compareDatasetCreated)
  
  plot_distances(compareSpectraPyClass$comparison_labels,compareSpectraPyClass$distances,
                 'boxplot',
                 input$plot_width_compare,input$plot_height_compare,
                 input$plot_type_compare,input$plot_axis_size_compare)
  
})



