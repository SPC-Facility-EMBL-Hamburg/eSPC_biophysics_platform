observeEvent(input$massPhotometryFileCalibration,{
  
  req(input$massPhotometryFileCalibration)
  reactives$data_loadedCalibration <- FALSE
  
  withBusyIndicatorServer("Go2",{
    
    fileExtension <- getFileNameExtension(input$massPhotometryFileCalibration$datapath)
    
    if (fileExtension == "h5") {
      refeynCalib$load_data_h5(input$massPhotometryFileCalibration$datapath)
    }
    
    if (fileExtension == "csv") {
      refeynCalib$load_data_csv(input$massPhotometryFileCalibration$datapath)
    }
    
    pks_initial <- sapply(refeynCalib$pks_initial, function(x) signif(x*cstFactorForContrast,2))
    maxNpeaks   <- min(6,length(pks_initial))
    pks_initial <- pks_initial[1:maxNpeaks]
    
    updateTextInput(session, "starting_valuesContrast", 
                    value = paste(pks_initial,collapse=" "))
    
    knownMassesInitial <- c(480,148,66,200,120,600,1000,300)[1:maxNpeaks] 
    minLimit <- floor(min(refeynCalib$contrasts)*cstFactorForContrast)
    updateNumericInput(session,"leftLimitWindowRangeContrast",
                       value = minLimit, min = -1e6, max = 0, step = 1)
    updateSliderInput(session,"window_rangeContrast",NULL,min = minLimit, 
                      max = 0,value = c(minLimit,0),step = 1)
    updateTextInput(session, "knownMasses", 
                    value = paste(knownMassesInitial,collapse=" "))
    reactives$data_loadedCalibration <- TRUE
    Sys.sleep(1)
    
  })
},priority = 10)

output$data_loadedCalibration <- reactive({reactives$data_loadedCalibration})
outputOptions(output, "data_loadedCalibration", suspendWhenHidden = FALSE)


updateModels <- function (slope,intercept) {

  if (length(photoMolModels$models) > 0) {

  reactives$data_loaded <- FALSE

  for (model in photoMolModels$models) {
    if (py_has_attr(model,'contrasts')) {
      model$contrastsToMass(slope,intercept)
    }
  }

  updateInputBox()
  reactives$data_loaded <- TRUE

  }
  return(NULL)
}

observeEvent(list(input$leftLimitWindowRangeContrast,input$rightLimitWindowRangeContrast),{
  
  isolate({
    
    if (input$leftLimitWindowRangeContrast < input$rightLimitWindowRangeContrast) {
      
      updateSliderInput(session,"window_rangeContrast",NULL,min = input$leftLimitWindowRangeContrast, 
                        max = input$rightLimitWindowRangeContrast,
                        value = c(input$leftLimitWindowRangeContrast,0),
                        step = 1)
      
    }
    
  })
  
})

modify_refeynCalibration_data <- reactive({
  
  if (!(reactives$data_loadedCalibration)) {return(NULL)}
  
  Sys.sleep(0.2)
  
  lower_limit_histogram <- input$window_rangeContrast[1] / cstFactorForContrast
  upper_limit_histogram <- input$window_rangeContrast[2] / cstFactorForContrast
  window                <- c(lower_limit_histogram,upper_limit_histogram)
    
  refeynCalib$create_histo(window=window,bin_width=input$bin_widthContrast   / cstFactorForContrast)
  
  starting_values <- get_guess_positions(input$starting_valuesContrast,cstFactorForContrast)
  
  refeynCalib$fit_histo(guess_pos=starting_values, max_std=0.1)
  
  knownMasses <- get_guess_positions(input$knownMasses)
  
  if (length(knownMasses) != length(starting_values)) {
    shinyalert(title="The number of markers (known masses) does not match the number of fitted gaussians.
               ", type = "warning",closeOnClickOutside=F)
    
    Sys.sleep(3)
    return(NULL)
  } 
  
  refeynCalib$calibrate(knownMasses)

  updateModels(refeynCalib$calib_params[1],refeynCalib$calib_params[2])

  # We want to evaluate this expression everytime the user changes sth in the UI 
  return(refeynCalib) 
})

observeEvent(list(input$interceptCustom,input$slopeCustom),{
  req(input$activateCalibration)
  req(input$interceptCustom != 0)
  req(input$slopeCustom != 0)

  updateModels(input$slopeCustom/1e6,input$interceptCustom/1e6)

})

observeEvent(length(get_guess_positions(input$starting_valuesContrast)),{
  
  req(input$massPhotometryFileCalibration)
  legends <- paste0("M #",1:length(get_guess_positions(input$starting_valuesContrast)))
  
  if (length(legends) > 1) {legends <- c("Gaussian sum",legends)}
  numberOfLegends <- length(legends)
  
  colorPalette <- colorPalette9
  if(numberOfLegends >= 10) colorPalette <- colorPalette12
  
  legendDf <- data.frame(legends = legends,color=colorPalette[1:numberOfLegends],
                         select  = as.logical(rep(TRUE,numberOfLegends)))
  
  color_cells <- data.frame(col=2,row=1:numberOfLegends)
  
  output$legendInfoCalibration <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,colHeaders=NULL,
                                                          col_highlight = color_cells$col - 1,
                                                          row_highlight = color_cells$row - 1
  ) %>% hot_col(col = c(1,2),
                renderer = myrenderer) %>% 
      hot_col(col = c(3),
              renderer = myrendererBoolean)})
  
})

observeEvent(input$legendInfoCalibration,{
  
  req(input$massPhotometryFileCalibration)
  updateSelectInput(session,"mol2changeColorCalibration","Set colour",
                    isolate(get_legend_from_rhandTable(input$legendInfoCalibration)),
                    isolate(input$mol2changeColorCalibration))
  
})

observeEvent(input$colorForLegendCalibration,{
  
  req(input$massPhotometryFileCalibration)
  isolate({
    legends <- get_legend_from_rhandTable(input$legendInfoCalibration)
    colors  <- get_colors_from_rhandTable(input$legendInfoCalibration)
    sels    <- get_sel_from_rhandTable(input$legendInfoCalibration)
    
    idx <- which(legends == input$mol2changeColorCalibration)
    colors[idx] <- input$colorForLegendCalibration
    
    legendDf <- data.frame(legends = legends,color=colors,select = as.logical(sels))
    
    color_cells <- data.frame(col=2,row=1:length(colors))
    output$legendInfoCalibration <- renderRHandsontable({rhandsontable(legendDf,rowHeaders=NULL,colHeaders=NULL,
                                                            col_highlight = color_cells$col - 1,
                                                            row_highlight = color_cells$row - 1
    ) %>% hot_col(col = c(1,2),
                  renderer = myrenderer) %>% 
        hot_col(col = c(3),
                renderer = myrendererBoolean)})
    
  })
  
})

output$contrast_plot_calib <- renderPlotly({
  
  if (is.null(modify_refeynCalibration_data())) {return(NULL)}
  
  req(input$legendInfoCalibration)
  
  legends <- isolate(get_legend_from_rhandTable(input$legendInfoCalibration))
  colors  <- isolate(get_colors_from_rhandTable(input$legendInfoCalibration))
  sels    <- isolate(get_sel_from_rhandTable(input$legendInfoCalibration))
  
  plot <-   plotRefeynFit(list(refeynCalib),input$baseline,input$plot_widthCalibration,
                          input$plot_heightCalibration, input$plot_typeCalibration,
                          input$plot_axis_sizeCalibration,legends,colors,sels,
                          input$show_contrastLegend,TRUE,FALSE,input$show_contrastPlot)

  return(plot)
  # defined in server_files/plot_functions.R
})

observeEvent(list(input$calibrationMethod,input$activateCalibration),{
  reactives$calibrationMethod <- isolate(input$calibrationMethod)
})

output$calibrationMethod <- reactive({reactives$calibrationMethod})
outputOptions(output, "calibrationMethod", suspendWhenHidden = FALSE)

output$mass_vs_contrast <- renderPlotly({
  
  if (is.null(modify_refeynCalibration_data())) {return(NULL)}
  
  contrasts <- refeynCalib$calib_points
  masses    <- refeynCalib$calib_stand
  slope     <- refeynCalib$calib_params[1]
  intercept <- refeynCalib$calib_params[2]
  
  plot <- plotMass_vs_contrast(masses,contrasts,slope,intercept,
                               input$plot_widthCalibration, 
                               input$plot_heightCalibration, input$plot_typeCalibration,
                               input$plot_axis_sizeCalibration)
  return(plot)
  # defined in server_files/plot_functions.R
}
)

output$fittedParamsCalibration <- renderTable({
  
  if (is.null(modify_refeynCalibration_data())) {return(NULL)}
  
  table <- refeynCalib$fit_table
  table[,1] <- table[,1]*cstFactorForContrast
  table[,2] <- table[,2]*cstFactorForContrast
  table[,1] <- paste0(signif(table[,1],2), " / 1e3")
  table[,2] <- paste0(signif(table[,2],2), " / 1e3")
  
  return(table[,1:(ncol(table)-1)])
},digits = 0)

output$calibParams <- renderTable({

  if (is.null(modify_refeynCalibration_data())) {return(NULL)}

  slope     <- refeynCalib$calib_params[1] * 1e6
  intercept <- refeynCalib$calib_params[2] * 1e6
  r_sq      <- refeynCalib$calib_r2
  
  table <- data.frame(slope,intercept,r_sq)
  colnames(table) <- c("Slope * 1e6","Intercept * 1e6","R squared")
  return(table)
},digits = 4)

