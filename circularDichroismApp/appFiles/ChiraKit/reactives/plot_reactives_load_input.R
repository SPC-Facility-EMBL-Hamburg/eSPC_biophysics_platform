observeEvent(input$wavelengthRange,{
  
  req(reactives$data_loaded)
  
  reactives$data_loaded <- FALSE
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])
  
  append_record_to_logbook(paste0('Setting wavelength range (nm) to ',wlRange[1],' - ',wlRange[2]))
  
  wlRangeLeft <- wlRange[1]
  
  if (wlRangeLeft < 175) wlRangeLeft <- 175
    
  if (wlRangeLeft <= 190) {
    updateNumericInput(session,'lower_wl_secStr',NULL,
                       wlRangeLeft, min = wlRangeLeft, max = 190) 
  }

  Sys.sleep(0.2)
  reactives$data_loaded <- TRUE
  
})

output$cdSpectra <- renderPlotly({
  
  req(reactives$data_loaded)
  req(input$legendInfo)
  
  legends      <- get_legend_from_rhandTable(input$legendInfo)
  colorPalette <- get_colors_from_rhandTable(input$legendInfo)
  sels         <- get_sel_from_rhandTable(input$legendInfo)
  
  df_list       <- cdAnalyzer$get_experiment_properties('spectraNames')
  total_columns <- sum(sapply(df_list, length))
  
  req(length(legends) == total_columns)
  
  fig <- plotCDexperiments(cdAnalyzer, input$workingUnits, 
                           input$plot_width,  input$plot_height,
                           input$plot_type,  input$plot_axis_size, 
                           legends, colorPalette, sels,
                           showGridX=input$showGridX,
                           showGridY=input$showGridY)
  
  return(fig)
  
})

output$htSpectra <- renderPlotly({
  
  req( reactives$data_loaded, input$legendInfo, length(cdAnalyzer$experimentNames)>0 )
  
  legends      <- get_legend_from_rhandTable(input$legendInfo)
  colorPalette <- get_colors_from_rhandTable(input$legendInfo)
  sels         <- get_sel_from_rhandTable(input$legendInfo)
  
  df_list       <- cdAnalyzer$get_experiment_properties('spectraNames')
  total_columns <- sum(sapply(df_list, length))
  
  req(length(legends) == total_columns)
  
  fig <- plotCDexperimentsHT(cdAnalyzer,
                             input$plot_width,  input$plot_height,
                             input$plot_type,  input$plot_axis_size,
                             legends,colorPalette,sels,
                             showGridX=input$showGridX,
                             showGridY=input$showGridY)
  
  return(fig)
})

output$cd_ht_spectra <- renderPlotly({
  
  req( reactives$data_loaded, input$legendInfo, length(cdAnalyzer$experimentNames)>0 )
  
  legends      <- get_legend_from_rhandTable(input$legendInfo)
  colorPalette <- get_colors_from_rhandTable(input$legendInfo)
  sels         <- get_sel_from_rhandTable(input$legendInfo)
  
  df_list       <- cdAnalyzer$get_experiment_properties('spectraNames')
  total_columns <- sum(sapply(df_list, length))
  
  req(length(legends) == total_columns)
  
  fig <- plot_cd_and_voltage(cdAnalyzer, input$workingUnits,
                             input$plot_width,  input$plot_height,
                             input$plot_type,  input$plot_axis_size,
                             legends, colorPalette, sels,
                             showGridX=input$showGridX,
                             showGridY=input$showGridY)
  
  return(fig)
})

output$cdSpectraMiliDeg <- renderPlotly({
  
  req(reactives$data_loaded)
  req(input$legendInfo)
  
  legends      <- get_legend_from_rhandTable(input$legendInfo)
  colorPalette <- get_colors_from_rhandTable(input$legendInfo)
  sels         <- get_sel_from_rhandTable(input$legendInfo)
  
  df_list       <- cdAnalyzer$get_experiment_properties('spectraNames')
  total_columns <- sum(sapply(df_list, length))
  
  req(length(legends) == total_columns)
  
  fig <- plotCDexperiments(cdAnalyzer, input$workingUnits, 
                           input$plot_width,  input$plot_height,
                           input$plot_type,  input$plot_axis_size, 
                           legends, colorPalette, sels,useMilliDeg = TRUE,
                           showGridX=input$showGridX,
                           showGridY=input$showGridY)
  return(fig)
  
})


