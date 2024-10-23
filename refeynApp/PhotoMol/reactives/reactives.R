reactives <- reactiveValues(data_loaded=FALSE,data_loadedCalibration=FALSE,
                            calibrationMethod="calibrationFile",
                            nFiles=0)

output$nFiles <- reactive({reactives$nFiles})
outputOptions(output, "nFiles", suspendWhenHidden = FALSE)

updateInputBox <- function(){

  hist_counts <- photoMolModels$get_properties('hist_counts')
  hist_mass   <- photoMolModels$get_properties('hist_mass')

  limits <- lapply(1:length(hist_counts),function(i) get_mass_limits(hist_counts[[i]],hist_mass[[i]]))

  minValue <- min(unlist(limits))
  maxValue <- max(unlist(limits))

  updateSliderInput(session,"window_range",NULL,min = minValue,
                    max = maxValue*1.05,value = c(0,maxValue),step = 5)

  updateNumericInput(session,"leftLimitWindowRange",
                     value = minValue, min = -1e12, max = 0, step = 1)
  updateNumericInput(session,"rightLimitWindowRange",
                     value = maxValue, min = 0, max = 1e12, step = 1)

  pks_initial <- photoMolModels$get_properties('pks_initial')

  for (i in 1:length(pks_initial)) {

    if (length(pks_initial[[i]]) > 0) {
        updateTextInput(session, paste0("starting_values",i), value = paste(pks_initial[[i]],collapse=" "))
    }

  }

  if (max(unlist(pks_initial)) > 500) {

    updateSliderInput(session,"upper_limit_std",NULL,min = 5, max = 300,value = 200)
    updateSliderInput(session,"position_tolerance",NULL,min = 1, max = 300,value = 200)

  }

}

# Load example dataset
observeEvent(input$GoLoadExample,{
  
  reactives$data_loaded <- FALSE

  photoMolModels$load_models("example","www/demo.h5")
  updateInputBox()
  reactives$data_loaded <- TRUE
  Sys.sleep(1)
})

resetPlotsAndTables <- function() {

  output$counts_plot           <- NULL
  output$binding_plot          <- NULL
  output$counts_plotNormalized <- NULL
  output$fittedParams          <- NULL
  output$legendInfo            <- NULL
  output$countPlot             <- NULL

  output$counts_plot_stacked           <- NULL
  output$counts_plotNormalized_stacked <- NULL

  return(NULL)

}

observeEvent(input$massPhotometryFile,{
  
  req(input$massPhotometryFile)

  resetPlotsAndTables()

  updateCheckboxInput(session,"automaticFit",value = FALSE)

  for (i in 1:5) updateTextInput(session, paste0("starting_values",i), value = '')

  reactives$data_loaded <- FALSE

  if (length(input$massPhotometryFile$name) > 5) {
    shinyalert(title="Please select a maximum of 5 files.", type = "warning")
    return(NULL)
  }

  withBusyIndicatorServer("Go",{

    photoMolModels$load_models(input$massPhotometryFile$name,input$massPhotometryFile$datapath)

    if (photoMolModels$allMassesLoaded) {

      updateInputBox()
      reactives$nFiles      <- length(photoMolModels$models)
      reactives$data_loaded <- TRUE
    } else  {
      shinyalert(title="Masses were not found in the input file(s).
                 Please do the calibration and use the fitted intercept and slope to convert 
                 from contrasts to masses.", type = "warning")
    }
    
    Sys.sleep(0.5)
    
  })
  updateCheckboxInput(session,"automaticFit",value = TRUE)

},priority = 10)

output$dataLoaded <- reactive({reactives$data_loaded})
outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)

observeEvent(list(input$leftLimitWindowRange,input$rightLimitWindowRange),{
  
  updateSliderInput(session,"window_range",NULL,min = input$leftLimitWindowRange, 
                    max = input$rightLimitWindowRange,value = c(input$leftLimitWindowRange,input$rightLimitWindowRange),
                    step = 1)
})

createPlotsAndTables <- function() {

  lower_limit_histogram <- input$window_range[1]
  upper_limit_histogram <- input$window_range[2]
  window                <- c(lower_limit_histogram,upper_limit_histogram)

  if(lower_limit_histogram >= upper_limit_histogram) return(NULL)

  fit_tables <- list()

  i     <- 0
  cnt1  <- 0
  cnt2  <- 0

  legends_all     <- c()

  for (model in photoMolModels$models) {

    i <- i + 1

    model$create_histo(window=window,bin_width=input$bin_width)

    starting_values <- get_guess_positions(input[[paste0("starting_values",i)]])
    starting_values <- starting_values[starting_values > input$min_observed_mass]

    if (length(starting_values) == 0) {
      model$fit <- NULL
      next
    }

    model$fit_histo(guess_pos=starting_values,tol=input$position_tolerance,
                    max_std=input$upper_limit_std,
                    min_observed_mass=input$min_observed_mass,
                    baseline=input$baseline)


    Sys.sleep(0.1)

    table  <- as.data.frame(model$fit_table)

    legends <- paste0("Mol. #",cnt2+1:length(starting_values))
    if (length(legends) > 1) {legends <- c(paste0("Gaussian sum (",i,")"),legends)}

    cnt2         <- cnt2 + length(starting_values)
    legends_all  <- c(legends_all,legends)

    if (length(photoMolModels$models) >1) table$File <- model$name

    if (nrow(table) > 0) {
      cnt1 <- cnt1 + 1
      fit_tables[[cnt1]] <- table
    }

  }

  if (length(fit_tables) > 0) {

    fit_table <- do.call(rbind,fit_tables)
    # Render without the amplitude column
    output$fittedParams <- renderTable({fit_table[,-5]},digits = 0)

  } else {

    output$fittedParams <- NULL

  }

  numberOfLegends <- length(legends_all)

  if (numberOfLegends == 0) {
    legendDf <- data.frame(legends = 'No fits',color = 'NA',select = FALSE)
    color_cells <- data.frame(col=2,row=1)
  } else {
    colorPalette <- colorPalette9
    if(numberOfLegends >= 10) colorPalette <- colorPalette12
    if(numberOfLegends >= 13) colorPalette <- colorPalette40

    legendDf <- data.frame(
      legends = legends_all,color=colorPalette[1:numberOfLegends],
      select  = as.logical(rep(TRUE,numberOfLegends)))

    color_cells <- data.frame(col=2,row=1:numberOfLegends)
  }

  output$legendInfo <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,colHeaders=NULL,
                                                          col_highlight = color_cells$col - 1,
                                                          row_highlight = color_cells$row - 1
  ) %>% hot_col(col = c(1,2),
                renderer = myrenderer) %>%
      hot_col(col = c(3),
              renderer = myrendererBoolean)})

  output$counts_plot <- renderPlotly({

    req(photoMolModels$allMassesLoaded)
    req(input$legendInfo)

    legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
    colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
    sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))

    plot <-   plotRefeynFit(photoMolModels$models,input$baseline,input$plot_width, input$plot_height,
                            input$plot_type, input$plot_axis_size,legends,colors,sels,
                            input$show_massesLegend,FALSE,FALSE,input$show_massesPlot)

    if (input$runSimulation) plot <- addSimulation2plotRefeynFit(
      plot,input$positionSimulate,input$stdSimulate,input$amplitudeSimulate,input$leftLimitSimulate)

    return(plot)})

  output$counts_plot_stacked <- renderPlotly({

    req(photoMolModels$allMassesLoaded)
    req(input$legendInfo)

    legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
    colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
    sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))

    plot <-   plotRefeynFit(photoMolModels$models,input$baseline,input$plot_width, input$plot_height,
                            input$plot_type, input$plot_axis_size,legends,colors,sels,
                            input$show_massesLegend,FALSE,FALSE,input$show_massesPlot,TRUE)

    if (input$runSimulation) plot <- addSimulation2plotRefeynFit(
      plot,input$positionSimulate,input$stdSimulate,input$amplitudeSimulate,input$leftLimitSimulate)

    return(plot)})

  output$counts_plotNormalized <- renderPlotly({

    req(photoMolModels$allMassesLoaded)
    req(input$legendInfo)

    legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
    colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
    sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))

    plot <-   plotRefeynFit(photoMolModels$models,input$baseline,input$plot_width, input$plot_height,
                            input$plot_type, input$plot_axis_size,legends,colors,sels,
                            input$show_massesLegend,FALSE,TRUE,input$show_massesPlot)

    return(plot)})

  output$counts_plotNormalized_stacked <- renderPlotly({

    req(photoMolModels$allMassesLoaded)
    req(input$legendInfo)

    legends <- isolate(get_legend_from_rhandTable(input$legendInfo))
    colors  <- isolate(get_colors_from_rhandTable(input$legendInfo))
    sels    <- isolate(get_sel_from_rhandTable(input$legendInfo))

    plot <-   plotRefeynFit(photoMolModels$models,input$baseline,input$plot_width, input$plot_height,
                            input$plot_type, input$plot_axis_size,legends,colors,sels,
                            input$show_massesLegend,FALSE,TRUE,input$show_massesPlot,TRUE)

    return(plot)})

  output$binding_plot <- renderPlotly({
  # see server_files/plot_functions.R
  return(plotRefeynMassHist(
    photoMolModels$models,input$plot_width, input$plot_height,
    input$plot_type, input$plot_axis_size))
  })

  return(NULL)

}

observeEvent(list(
    input$starting_values1,
    input$starting_values2,
    input$starting_values3,
    input$starting_values4,
    input$starting_values5,
    input$window_range[1],input$window_range[2],
    input$bin_width,input$upper_limit_std,input$position_tolerance,
    input$baseline,input$min_observed_mass,
    input$automaticFit
    ),{

  req(input$automaticFit)
  req(reactives$data_loaded)
  req(photoMolModels$allMassesLoaded)

  resetPlotsAndTables()

  for (i in 1:5) {

    pks_initial <- input[[paste0("starting_values",i)]]
    pks_initial <- get_guess_positions(pks_initial)

    if (length(pks_initial) == 0) {next}

    pks_initial <- as.numeric(pks_initial)

    rightLimit <- input$window_range[2]

    if (any(pks_initial < rightLimit)) {
      pks_initial <- as.numeric(pks_initial[pks_initial < rightLimit])
      updateTextInput(session, paste0("starting_values",i), value = paste(pks_initial,collapse=" "))
    }

  }

  createPlotsAndTables()

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




