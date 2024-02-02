## Decide if we need to render the L curve 
output$ui_contributions_plots <- renderUI({
  
  req(dlsDataUpdated())
  
  if (input$regMethod == 'fixedAlpha') {
    tabBox(
      title = "", width = 6,id = "tabset2",
      tabPanel("Intensity dist",withSpinner(plotOutput("intensityDistribution"))),
      tabPanel("Intensity dist - Collapsed",withSpinner(plotlyOutput("intensityDistributionPlotly"))),
      tabPanel("Residuals",withSpinner(plotOutput("residuals"))),
      tabPanel('Table for average',rHandsontableOutput('tableForAvg')),
      tabPanel("Intensity dist - average",withSpinner(plotOutput("intensityDistributionAvg")))
    )
  } else {
    tabBox(
      title = "", width = 6,id = "tabset2",
      tabPanel("Intensity dist",withSpinner(plotOutput("intensityDistribution"))),
      tabPanel("Intensity dist - Collapsed",withSpinner(plotlyOutput("intensityDistributionPlotly"))),
      tabPanel("Residuals",withSpinner(plotOutput("residuals"))),
      tabPanel("L curve",withSpinner(plotOutput("lCurve"))),
      tabPanel('Table for average',rHandsontableOutput('tableForAvg')),
      tabPanel("Intensity dist - average",withSpinner(plotOutput("intensityDistributionAvg")))
    )
  }
})

# Render autocorrelation versus time plot
output$autocorrelation_input_tab <- renderPlotly({
  
  req(reactives$data_loaded)
  data <- dlsData(dlsAnalyzer)
  
  if (is.null(data)) return(NULL) # Return null in the case that there is no data!
  
  fig <- plotAutocorrelation(data,input$plot_width,input$plot_height,
                             input$plot_type,input$plot_axis_size,
                             input$logScaleType,input$splitFactorPreproccess,
                             dlsAnalyzer$experimentNames)
  
  return( fig ) 
})

# Render autocorrelation versus time plot
# Filter the curves according to the residuals threshold
output$autocorrelation <- renderPlotly({
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  req(reactives$renderData)
  
  data <- dlsDataUpdated()
  
  residuals   <- data$residuals 
  dlsOridata  <- data$data
  
  dlsOridata  <- dlsOridata[dlsOridata$variable %in% selectedVariables(),]
  
  fig <- plotAutocorrelation(dlsOridata,input$plot_width,input$plot_height,
                             input$plot_type,input$plot_axis_size,
                             input$logScaleType,input$splitFactor,
                             dlsAnalyzer$experimentNames)
  
  return( fig ) 
})

output$fittedAutocorrelation1 <- renderPlotly({
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  req(reactives$renderData)
  
  data <- dlsDataUpdated()
  
  dlsOridata  <- data$data
  dlsPredData <- data$predictedData
  
  dlsOridata  <- dlsOridata[dlsOridata$variable %in% selectedVariables(),]
  dlsPredData <- dlsPredData[dlsPredData$variable %in% selectedVariables(),]
  
  fig <- plotFittedAutocorrelation1(dlsOridata,dlsPredData,
                                    input$plot_width,input$plot_height,
                                    input$plot_type,input$plot_axis_size,
                                    input$logScaleType,input$splitFactor,
                                    dlsAnalyzer$experimentNames)
  
  return( fig ) 
})

output$intensityDistribution <- renderPlot({
  
  req(estimatedContributions())
  req(reactives$renderData)

  estimatedContributions <- estimatedContributions()
  if (input$distribution_type == 'diff') {
    estimatedContributions <- estimatedContributionsDiff()
  }
  
  # Fix plotting the distribution of diff coeff!!!!!!
  fig <- plot_type_hr_distribution(input$plot_type_hr_distribution,
                                   estimatedContributions,
                                   input$plot_axis_size,
                                   input$overlapFactor,
                                   input$rowOrder,
                                   input$peak_filter
                                   )
  
  leftBounds   <- get_positions(input$peakLeftBounds)
  rightBounds  <- get_positions(input$peakRightBounds)
  
  nLeftBounds  <- length(leftBounds)
  nRightBounds <- length(rightBounds)
  
  c1 <- any(is.na(leftBounds))
  c2 <- any(is.na(rightBounds))
  c3 <- nLeftBounds != nRightBounds
  
  if (c1 | c2 | c3 ) return(fig) # Error in the input: no changes done!
  
  # Change scale factor in order to add the areas to the plot
  if (input$plot_type_hr_distribution == 'histogram') {
    
    useDiffCoeff <- 'diff' %in% colnames(estimatedContributions)
    if (useDiffCoeff) {
      logScale <- estimatedContributions$diff[2] / estimatedContributions$diff[1]
    } else {
      logScale <- estimatedContributions$hr[2] / estimatedContributions$hr[1]
    }
    
    leftBounds <- log(leftBounds,base = logScale)
    rightBounds <- log(rightBounds,base = logScale)
  }
  
  if (input$show_peaks_in_plot) {
    fig <- add_peak_areas_to_plot(fig,leftBounds,rightBounds)
  }
  
  return( fig ) 
})

output$intensityDistributionAvg <- renderPlot({
  
  req(estimatedContributions())
  req(reactives$renderData)
  req(input$tableForAvg)
  
  estimatedContributions <- estimatedContributions()
  if (input$distribution_type == 'diff') {
    estimatedContributions <- estimatedContributionsDiff()
  }
  
  # Fix plotting the distribution of diff coeff!!!!!!
  
  mapping <- hot_to_r(input$tableForAvg)

  fig <- plot_type_hr_distribution(input$plot_type_hr_distribution,
                                   estimatedContributions,
                                   input$plot_axis_size,
                                   input$overlapFactor,
                                   input$rowOrder,
                                   input$peak_filter,
                                   TRUE,mapping)
  
  leftBounds   <- get_positions(input$peakLeftBounds)
  rightBounds  <- get_positions(input$peakRightBounds)
  
  nLeftBounds  <- length(leftBounds)
  nRightBounds <- length(rightBounds)
  
  c1 <- any(is.na(leftBounds))
  c2 <- any(is.na(rightBounds))
  c3 <- nLeftBounds != nRightBounds
  
  if (c1 | c2 | c3 ) return(fig) # Error in the input: no changes done!
  
  # Change scale factor in order to add the areas to the plot
  if (input$plot_type_hr_distribution == 'histogram') {
    
    useDiffCoeff <- 'diff' %in% colnames(estimatedContributions)
    if (useDiffCoeff) {
      logScale <- estimatedContributions$diff[2] / estimatedContributions$diff[1]
    } else {
      logScale <- estimatedContributions$hr[2] / estimatedContributions$hr[1]
    }
    
    leftBounds  <- log(leftBounds,base = logScale)
    rightBounds <- log(rightBounds,base = logScale)
  }
  
  if (input$show_peaks_in_plot) {
    fig <- add_peak_areas_to_plot(fig,leftBounds,rightBounds)
  }
  
  return( fig ) 
})

output$residuals <- renderPlot({
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  req(reactives$renderData)
  
  data <- dlsDataUpdated()
  
  autocorrelationData <- data$data
  residuals           <- data$residuals  
  
  if (length(residuals) != nrow(autocorrelationData)) return(NULL)
  
  fig <- plotResiduals(autocorrelationData,residuals,
                       input$plot_axis_size,
                       input$residualsPlotSelection)
  
  return( fig ) 
})

output$intensityDistributionPlotly <- renderPlotly({

  req(estimatedContributions())
  req(req(reactives$renderData))
  
  estimatedContributions <- estimatedContributions()
  if (input$distribution_type == 'diff') {
    estimatedContributions <- estimatedContributionsDiff()
  }
  
  fig <- plothydrodynamicRadiusDistributionCollapsed(
    estimatedContributions,input$plot_width,
    input$plot_height,input$plot_type,
    input$plot_axis_size,input$peak_filter)
  
  return( fig ) 
})

output$lCurve <- renderPlot({
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  req(req(reactives$renderData))
  
  lCurveData <- dlsDataUpdated()$lCurveData
  
  if (is.null(lCurveData)) return(NULL)
  
  fig <- plotLcurve(lCurveData$lCurveDf,lCurveData$dfSel,
                    input$plot_axis_size,
                    input$residualsPlotSelection)
  
  return( fig ) 
})

## Plots from the simulation panel

output$hrDistributionNumberSim <- renderPlot({
  
  req(simulatedDistributions())
  
  estimatedContributions <- simulatedDistributions()$allDistributions
  estimatedContributions$value <- estimatedContributions$contributionsNumber
  
  yLabel <- "Number weighted contribution (%)" 
  
  # Plot type when we only have one dataset
  if (length(unique(estimatedContributions$variable))==1) {
    return(plot_simulated_distribution(estimatedContributions)+
             ylab(yLabel))
  }
  
  fig <- plot_type_hr_distribution('histogram',
    estimatedContributions,input$plot_axis_size,
    input$overlapFactor/8,input$rowOrder) + ylab(yLabel)
  
  return( fig ) 
})

output$hrDistributionVolumeSim <- renderPlot({
  
  req(simulatedDistributions())
  
  estimatedContributions       <- simulatedDistributions()$allDistributions
  estimatedContributions$value <- estimatedContributions$contributionsVolume
  
  yLabel <- "Volume weighted contribution (%)" 
  
  # Plot type when we only have one dataset
  if (length(unique(estimatedContributions$variable))==1) {
    return(plot_simulated_distribution(estimatedContributions)+
             ylab(yLabel))
  }
  
  fig <- plot_type_hr_distribution('histogram',
                                   estimatedContributions,input$plot_axis_size,
                                   input$overlapFactor/8,input$rowOrder) + 
    ylab(yLabel)
  
  return( fig ) 
})

output$hrDistributionIntensitySim <- renderPlot({
  
  req(simulatedDistributions())
  
  estimatedContributions       <- simulatedDistributions()$allDistributions
  estimatedContributions$value <- estimatedContributions$contributionsIntensity
  
  yLabel <- "Intensity weighted contribution (%)" 
  
  # Plot type when we only have one dataset
  if (length(unique(estimatedContributions$variable))==1) {
    return(plot_simulated_distribution(estimatedContributions)+
             ylab(yLabel))
  }
  
  fig <- plot_type_hr_distribution('histogram',
                                   estimatedContributions,input$plot_axis_size,
                                   input$overlapFactor/8,input$rowOrder) + 
    ylab(yLabel)
  
  return( fig ) 
})

simulatedAutocorrelation <- reactive({
  
  simData       <- simulatedDistributions()
  contributions <- simData$allDistributions
  
  # dataframe to store the autocorrelation curves
  dfs           <- list()
  
  # Iterate the simulated distributions
  i <- 0
  for (var in simData$simNames) {
    i <- i+1
    
    tempDF  <- contributions[contributions$variable == var,]
    hr      <- tempDF$hr
    contEst <- tempDF$contributionsIntensity
    
    SimAngleInRadians <- simData$simAngles[i] * 0.0174533    #Degrees to radians
    
    # Remove all spaces to avoid error in the python code when converting to complex
    # Uncomment the next 2 lines only if we change to a light absorbing medium
    # refractiveIndex <- gsub("\\s", "", simData$simRefInd[i])
    # q         <- get_q(simData$simWL[i],refractiveIndex,SimAngleInRadians)
    
    # Comment if we change to a light absorbing medium
    q         <- get_q(simData$simWL[i],simData$simRefInd[i],SimAngleInRadians)
    
    # Hydrodynamic radius in meters, temperature in kelvin & viscosity in pascal-second
    diffusion <- np_array(
      diffusion_from_hydrodynamic_radius(
        np_array(hr/1e9),simData$simTemp[i]+273,simData$simViscosity[i]*1e-4))
    
    minTime <- log(1e-8,base = 1.5)
    maxTime <- log(10,base = 1.5)
    
    time    <- 1.5^(seq(minTime,maxTime,length.out = 500))
    
    autocorrelationPredicted <- g2_finite_aproximation(
      np_array(1/s_inverse_decay_rate(diffusion,q)),np_array(time),
      input$autocorrelationIntercept,
      np_array(contEst))
    
    df <- data.frame(time,value=autocorrelationPredicted,variable=var)
    
    gaussianError <- input$autocorrelationError
    if (gaussianError > 0){
      errors <- rnorm(n = nrow(df), mean = 0, sd = gaussianError)
      df$value <- df$value + sample(errors)
    }
    
    dfs[[i]] <- df
  }
  
  df <- do.call(rbind,dfs)
  
  return(df)
})

output$autocorrelationSimulate <- renderPlotly({

  req(simulatedDistributions())
  
  return(plotAutocorrelation(
    simulatedAutocorrelation(),10,10,'png',18,'microseconds',
    'dummySplitFactor',simData$simNames)) # This two parameters are useless here
  
})



