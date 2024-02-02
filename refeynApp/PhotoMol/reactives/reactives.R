reactives <- reactiveValues(data_loaded=FALSE,data_loadedCalibration=FALSE,
                            calibrationMethod="calibrationFile")

updateInputBox <- function(){
  limits <- get_mass_limits(refeyn$hist_counts,refeyn$hist_mass) 
  updateSliderInput(session,"window_range",NULL,min = limits$min, 
                    max = limits$max,value = c(0,limits$max),step = 1)
  # Let's improve the initial guess for the peaks
  updateTextInput(session, "starting_values", value = paste(refeyn$pks_initial,collapse=" "))
  
  updateNumericInput(session,"leftLimitWindowRange",
                     value = limits$min, min = -1e6, max = 0, step = 1)
  updateNumericInput(session,"rightLimitWindowRange",
                     value = limits$max, min = 0, max = 1e6, step = 1)
  
  if (max(refeyn$pks_initial) > 500) {
    
    updateSliderInput(session,"upper_limit_std",NULL,min = 5, max = 200,value = 200)
    updateSliderInput(session,"position_tolerance",NULL,min = 1, max = 200,value = 200)
    
  }
  
}

# Load example dataset
observeEvent(input$GoLoadExample,{
  
  reactives$data_loaded <- FALSE
  refeyn$load_data_h5("demoTestFile")
  
  updateInputBox()
  reactives$data_loaded <- TRUE
  Sys.sleep(1)
})

observeEvent(input$massPhotometryFile,{
  
  req(input$massPhotometryFile)
  
  reactives$data_loaded <- FALSE
  
  withBusyIndicatorServer("Go",{
    
    fileExtension <- getFileNameExtension(input$massPhotometryFile$datapath)
    
    if (fileExtension == "h5") {
      file.copy(input$massPhotometryFile$datapath,"0.h5",overwrite=TRUE)
      refeyn$load_data_h5("0.h5")
    }
    
    if (fileExtension == "csv") {
      file.copy(input$massPhotometryFile$datapath,"0.csv",overwrite=TRUE)
      refeyn$load_data_csv("0.csv")
    }

    if (refeyn$massesLoaded) {
      updateInputBox()
      reactives$data_loaded <- TRUE
    } else  {
      shinyalert(title="Masses were not found in the input file. 
                 Please do the calibration and use the fitted intercept and slope to convert 
                 from contrasts to masses.", type = "warning")
    }
    
    
    Sys.sleep(1)
    
  })
},priority = 10)

output$dataLoaded <- reactive({reactives$data_loaded})
outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)

observeEvent(list(input$leftLimitWindowRange,input$rightLimitWindowRange),{
  
  updateSliderInput(session,"window_range",NULL,min = input$leftLimitWindowRange, 
                    max = input$rightLimitWindowRange,value = c(input$leftLimitWindowRange,input$rightLimitWindowRange),
                    step = 1)
})

modify_refeyn_data <- reactive({
  
  if (!(reactives$data_loaded)) {return(NULL)}
  
  Sys.sleep(0.2)
  
  if (refeyn$massesLoaded) {
    
    lower_limit_histogram <- input$window_range[1]
    upper_limit_histogram <- input$window_range[2]
    window                <- c(lower_limit_histogram,upper_limit_histogram)
    
    req(lower_limit_histogram < upper_limit_histogram) 
    
    refeyn$create_histo(window=window,bin_width=input$bin_width)
    starting_values <- get_guess_positions(input$starting_values)
    
    req(all(abs(py_to_r(starting_values)) > input$min_observed_mass))

    refeyn$fit_histo(guess_pos=starting_values,tol=input$position_tolerance,
                     max_std=input$upper_limit_std,
                     min_observed_mass=input$min_observed_mass,
                     baseline=input$baseline)
  }
  
  # We want to evaluate this expression everytime the user changes sth in the UI 
  return(refeyn) 
})

observeEvent(length(get_guess_positions(input$starting_values)),{
  
  req(reactives$data_loaded)
  legends <- paste0("Molecule #",1:length(get_guess_positions(input$starting_values)))
  
  if (length(legends) > 1) {legends <- c("Gaussian sum",legends)}
  numberOfLegends <- length(legends)
  
  colorPalette <- colorPalette9
  if(numberOfLegends >= 10) colorPalette <- colorPalette12

  legendDf <- data.frame(legends = legends,color=colorPalette[1:numberOfLegends],
                         select  = as.logical(rep(TRUE,numberOfLegends)))
  
  color_cells <- data.frame(col=2,row=1:numberOfLegends)
  
  output$legendInfo <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,colHeaders=NULL,
                                                          col_highlight = color_cells$col - 1,
                                                          row_highlight = color_cells$row - 1
  ) %>% hot_col(col = c(1,2),
                renderer = myrenderer) %>% 
      hot_col(col = c(3),
              renderer = myrendererBoolean)})
  
})

observeEvent(input$legendInfo,{
  
  req(reactives$data_loaded)
  updateSelectInput(session,"mol2changeColor","Set colour",
                    isolate(get_legend_from_rhandTable(input$legendInfo)),isolate(input$mol2changeColor))
  
})

observeEvent(input$colorForLegend,{
  
  req(reactives$data_loaded)
  isolate({
    legends <- get_legend_from_rhandTable(input$legendInfo)
    colors  <- get_colors_from_rhandTable(input$legendInfo)
    sels    <- get_sel_from_rhandTable(input$legendInfo)
    
    idx <- which(legends == input$mol2changeColor)
    
    colors[idx] <- input$colorForLegend
    
    legendDf <- data.frame(legends = legends,color=colors,select = as.logical(sels))
    
    color_cells <- data.frame(col=2,row=1:length(colors))
    output$legendInfo <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,colHeaders=NULL,
                                                            col_highlight = color_cells$col - 1,
                                                            row_highlight = color_cells$row - 1
    ) %>% hot_col(col = c(1,2),
                  renderer = myrenderer) %>% 
        hot_col(col = c(3),
                renderer = myrendererBoolean)})
    
  })
  
})

output$counts_plot <- renderPlotly({
  
  if (is.null(modify_refeyn_data())) {return(NULL)}
  req(refeyn$massesLoaded)
  req(input$legendInfo)
  
  legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
  colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
  sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))
  
  plot <-   plotRefeynFit(refeyn,input$baseline,input$plot_width, input$plot_height, 
                          input$plot_type, input$plot_axis_size,legends,colors,sels,
                          input$show_massesLegend)
  
  if (input$runSimulation) plot <- addSimulation2plotRefeynFit(
    plot,input$positionSimulate,input$stdSimulate,input$amplitudeSimulate,input$leftLimitSimulate)

  if (input$show_massesPlot) {
    plot <- addLabels2plotRefeynFit(plot,refeyn$fit_table[,1],refeyn$fit_table[,5],
                                    sels,input$plot_axis_size)
  }
  
  return(plot)
 # defined in server_files/plot_functions.R
}
)

output$fittedParams <- renderTable({
  
  if (is.null(modify_refeyn_data())) {return(NULL)}
  req(refeyn$massesLoaded)
  table <- refeyn$fit_table
  return(table[,1:(ncol(table)-1)])
},digits = 0)

output$counts_plotNormalized <- renderPlotly({
  
  if (is.null(modify_refeyn_data())) {return(NULL)}

  req(refeyn$massesLoaded)
  req(input$legendInfo)
  
  legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
  colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
  sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))
  # see server_files/plot_functions.R
  plot <-   plotRefeynFitNormalized(refeyn,input$baseline,input$plot_width, input$plot_height, 
                          input$plot_type, input$plot_axis_size,legends,colors,sels,
                          input$show_massesLegend)
  
  if (input$show_massesPlot) {
    plot <- addLabels2plotRefeynFit(plot,refeyn$fit_table[,1],
                                    refeyn$fit_table[,5]/sum(refeyn$hist_counts),
                                    sels,input$plot_axis_size)
  }
  
  return(plot)
  
}
)

output$binding_plot <- renderPlotly({
  
  if (is.null(modify_refeyn_data())) {return(NULL)}
  req(refeyn$massesLoaded)

  # see server_files/plot_functions.R
  return(plotRefeynMassHist(refeyn,input$plot_width, input$plot_height, 
                            input$plot_type, input$plot_axis_size))
  
}
)




