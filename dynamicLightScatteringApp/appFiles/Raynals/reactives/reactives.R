reactives <- reactiveValues(data_loaded=FALSE,simulatedMetaData=NULL,
                            inputFileGiven=FALSE,renderData=FALSE)

# Useful to decide if we need to render the load example button or not
output$inputFileGiven             <- reactive({
  return(reactives$inputFileGiven)})

# Useful to conditionally render two ui elements:
#         ui_preprocessing_box.R and ui_signal_tab_box0.R
output$data_loaded             <- reactive({
  return(reactives$data_loaded)})

# Allow dynamic render of ui objects
outputOptions(output, "inputFileGiven", suspendWhenHidden = FALSE)
outputOptions(output, "data_loaded",    suspendWhenHidden = FALSE)

renderDlsPlateInfo <- function(name,tables) {
  
  # I don't know why a for loop didn't work to render the tables. 
  # We will have 1,2,3 or 6 tables.

  dlsPlateInfoNames <- getDlsPlateInfoNames(name,length(tables))
  
  output[[dlsPlateInfoNames[1]]] <- renderRHandsontable({rhandsontable(tables[[1]],rowHeaders=F)})
  
  if (length(tables) > 1) {
    output[[dlsPlateInfoNames[2]]] <- renderRHandsontable({rhandsontable(tables[[2]],rowHeaders=F)})
  }
  
  if (length(tables) > 2) {
    output[[dlsPlateInfoNames[3]]] <- renderRHandsontable({rhandsontable(tables[[3]],rowHeaders=F)})
  }

  if (length(tables) > 3) {
    
    output[[dlsPlateInfoNames[4]]] <- renderRHandsontable({rhandsontable(tables[[4]],rowHeaders=F)})
    output[[dlsPlateInfoNames[5]]] <- renderRHandsontable({rhandsontable(tables[[5]],rowHeaders=F)})
    output[[dlsPlateInfoNames[6]]] <- renderRHandsontable({rhandsontable(tables[[6]],rowHeaders=F)})
    
  }
}

# Load example dataset
observeEvent(input$GoLoadExample,{
  
  load <- dlsAnalyzer$loadExperiment('test.csv','test')
  dlsAnalyzer$experimentsOri[['test']]$setAutocorrelationData()

  output$dlsFilesInfo    <- NULL
  df <- generateDTtable(dlsAnalyzer)
  
  # To allow the creation of selectInput inside DT table
  session$sendCustomMessage('unbind-DT', 'dlsFilesInfo')
  output$dlsFilesInfo <- renderDTtable(df)
  
  updateSelectInput(session,"experiment2delete",choices = c("test"))
  
  nMeasurements         <- dlsAnalyzer$experimentsOri[["test"]]$nMeasurements
  numberOfDesiredTables <- get_numberOfDesiredTables(nMeasurements)
  
  df <- dlsAnalyzer$experimentsOri[["test"]]$sampleInfo
  
  tables                <- get_Table_list(df$conditions,1,1,numberOfDesiredTables)
  tables                <- split_table_into_list_of_tables(df,numberOfDesiredTables)
  
  tabP <- generateTabPanel(numberOfDesiredTables,"test")
  appendTab("tabBoxDlsWellsInfo",tabP,select=TRUE)
  
  renderDlsPlateInfo("test",tables)

  reactives$data_loaded    <- TRUE
  reactives$inputFileGiven <- TRUE
  
  Sys.sleep(1)
  
})

observeEvent(input$dlsFiles,{
  req(input$dlsFiles)
  reactives$data_loaded <- FALSE
  
  withBusyIndicatorServer("Go",{
    
    datapath  <- input$dlsFiles[[1, 'datapath']]
    name      <- input$dlsFiles[[1]]
    
    fileExt   <- getFileNameExtension(name)
    name      <- sub(paste0(".",fileExt),'',name)
    
    newFilePath   <- paste0("0.",fileExt)
    file.copy(datapath,newFilePath,overwrite=TRUE)
    
    if (fileExt == '7z')  system(paste0('7z e ',newFilePath))
    if (fileExt == 'zip') unzip(newFilePath)
    
    if (fileExt %in% c('7z','zip')) {
      
      xlsx_files <- list.files(path = ".", pattern = "\\.xlsx$", recursive = TRUE, full.names = TRUE)
      
      fileType <- detect_Excel_files_type(xlsx_files)
      
      if (fileType == 'Acquisition') {
        
        merged_df <- load_Excel_files_with_acquisitions(xlsx_files)
        
        shinyalert(text = paste("<b>Caution: The autocorrelation curves were averaged according to
                                the 'Acquisition' folders</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        
      } else {
        
        merged_df <- load_Excel_files_without_acquisitions(xlsx_files)
        
      }
      
      # Remove the files
      for (f in xlsx_files) {
        system(paste0("rm -f ",f))
      }  
      
      write.csv(merged_df,'0.csv',row.names = F)
      
    }
    ###
    
    exps      <- dlsAnalyzer$experimentNames
    
    if (name %in% exps) {
      shinyalert(text = paste("<b>File name is already being used.</b>"),
                 type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                 html=T)
      
      reactives$data_loaded <- TRUE
      return(NULL)
    } else { # Load last loaded file
      
      load <- dlsAnalyzer$loadExperiment('0.csv',name)
      
      # Change the angle of detection and wavelength (Prometheus Panta defaults)
      if (fileExt %in% c('7z','zip')) {
        
        dlsAnalyzer$experimentsOri[[name]]$lambda0         <- 405
        dlsAnalyzer$experimentsOri[[name]]$scatteringAngle <- 147 / 180 * pi
        
      }
      
      dlsAnalyzer$experimentsOri[[name]]$setAutocorrelationData()

      # Catch exception when data has wrong format
      if (load != "Data loaded successfully!!!") {
        shinyalert(text = paste("<b>The file was not loaded! Please verify the file format (User guide).</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        
        reactives$data_loaded <- TRUE
        return(NULL)
      }
      
      output$dlsFilesInfo    <- NULL
      Sys.sleep(0.1)

      df <- generateDTtable(dlsAnalyzer)
      # To allow the creation of selectInput inside DT table
      session$sendCustomMessage('unbind-DT', 'dlsFilesInfo')
      output$dlsFilesInfo <- renderDTtable(df)
      
      updateSelectInput(session,"experiment2delete",choices = c("None",dlsAnalyzer$experimentNames))

      nMeasurements <- dlsAnalyzer$experimentsOri[[name]]$nMeasurements
      
      numberOfDesiredTables <- get_numberOfDesiredTables(nMeasurements)
      
      df <- dlsAnalyzer$experimentsOri[[name]]$sampleInfo
      
      tables                <- get_Table_list(df$conditions,1,1,numberOfDesiredTables)
      tables                <- split_table_into_list_of_tables(df,numberOfDesiredTables)

      tabP <- generateTabPanel(numberOfDesiredTables,name)
      appendTab("tabBoxDlsWellsInfo",tabP,select=TRUE)
      
      renderDlsPlateInfo(name,tables)

    }
    Sys.sleep(2)
    
    shinyalert(text = paste("<b>The DLS file was successfully loaded. 
                            If you want to automatically complete the scan and read number
                            information, use the 'Experiments parameters' Box.</b>"),
               type = "success",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    
    reactives$data_loaded    <- TRUE
    reactives$inputFileGiven <- TRUE
  })
},priority = 10)

observeEvent(input$triggerDeletion,{
  req(reactives$data_loaded)

  selected <- input$experiment2delete
  
  if (selected == "None") return(NULL)
  
  reactives$data_loaded <- FALSE
  
  exps <- dlsAnalyzer$experimentNames
  
  if (selected == input$tabBoxDlsWellsInfo & length(exps) > 1) {
    shinyalert(text = paste("The experiment can't be removed, 
                            please change the active tab and try to remove it again,
                            or select another experiment."),
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    reactives$data_loaded <- TRUE
    return(NULL)
  }
  
  dlsAnalyzer$deleteExperiment(selected)
  removeTab(inputId = "tabBoxDlsWellsInfo", target = selected)
  
  if (length(exps) > 0) {
    
    updateSelectInput(session,"experiment2delete",selected = "None",
                      choices = c("None",dlsAnalyzer$experimentNames))
    df <- generateDTtable(dlsAnalyzer)
    
    # To allow the creation of selectInput inside DT table
    session$sendCustomMessage('unbind-DT', 'dlsFilesInfo')
    output$dlsFilesInfo <- renderDTtable(df)
    
  } else {
    output$dlsFilesInfo <- NULL
    updateSelectInput(session,"experiment2delete",choices = c('None'))
  }
  
  reactives$data_loaded <- length(dlsAnalyzer$experimentNames) > 0

})

# React to changes in reads, scans, initial wavelength, and scattering angle
observeEvent(input$dlsFilesInfo_cell_edit, {
  req(reactives$data_loaded)
  reactives$data_loaded <- FALSE

  info <- input$dlsFilesInfo_cell_edit
  
  columnIndex <- info$col[1] + 1 # Start index from 1
  rowIndex    <- info$row[1]
  value       <- as.numeric(info$value[1])
  
  expName <- dlsAnalyzer$experimentNames[rowIndex]
  
  if (columnIndex == 2) {
    selectedVar <- "scans"
  } else if (columnIndex == 3) {
    selectedVar <- "reads"
  } else if (columnIndex == 4) {
    selectedVar <- "lambda0"
  } else if (columnIndex == 5) {
    selectedVar <- "scatteringAngle"
    value <- value / 180 * pi # To radians!
  } else if (columnIndex == 6) {
    selectedVar <- "temperature"
    value <- value + 273 # To kelvin!
  } else if (columnIndex == 7) {
    selectedVar <- "refractiveIndex"
  } else { 
    selectedVar  <- "viscosity"
  } 
  
  dlsAnalyzer$setExperimentProperties(expName,selectedVar,value)
  
  if (columnIndex %in% c(2,3)) {
    
    idx   <- which(dlsAnalyzer$experimentNames == expName)
    nMeasurements         <- ncol(dlsAnalyzer$experimentsOri[[expName]]$autocorrelationOriginal)
    
    reads      <- dlsAnalyzer$getExperimentProperties('reads')[idx]
    scans      <- dlsAnalyzer$getExperimentProperties('scans')[idx]
    sampleInfo <- dlsAnalyzer$getExperimentProperties('sampleInfo')[[idx]]

    reads <- min(reads,nMeasurements)
    scans <- min(scans,nMeasurements)

    numberOfDesiredTables <- get_numberOfDesiredTables(nMeasurements)
    tables                <- get_Table_list(sampleInfo$conditions,reads,scans,numberOfDesiredTables)
    
    df <- getDataFrameBasedOnReadsAndScans(sampleInfo$conditions,reads,scans)
    
    dlsAnalyzer$experimentsOri[[expName]]$sampleInfo <- df
    
    renderDlsPlateInfo(expName,tables)
  }
  
  reactives$data_loaded <- TRUE

})

experimentsSamplesMetaDataFromTables <- reactive({
  
  expNames               <- dlsAnalyzer$experimentNames
  expTabData <- list()
  
  i <- 1
  for (expName in expNames) {
    
    nMeasurements         <- ncol(dlsAnalyzer$experimentsOri[[expName]]$autocorrelationOriginal)
    numberOfDesiredTables <- get_numberOfDesiredTables(nMeasurements)
    dlsPlateInfoNames     <- getDlsPlateInfoNames(expName,numberOfDesiredTables)
  
    expTables <- lapply(1:numberOfDesiredTables, function(n) input[[dlsPlateInfoNames[n]]] )
    
    expTabData[[i]] <- expTables

    i <- 1 + i
  }
  
  return(expTabData)
  
})

observeEvent(input$selectAll,{
  req(reactives$data_loaded)
  
  include <- input$selectAll %% 2 == 1
  
  reactives$data_loaded <- FALSE
  # Get metadata of samples
  experimentsSamplesMetaData  <- dlsAnalyzer$getExperimentProperties('sampleInfo')
  expNames                    <- dlsAnalyzer$experimentNames
  
  for (i in 1:length(expNames)) {
    
    expName               <- expNames[i]
    df                    <- experimentsSamplesMetaData[[i]]
    
    df$include            <- include # (Un)select all samples :)
    
    # Update dls sample metadata
    dlsAnalyzer$experimentsOri[[expName]]$sampleInfo <- df
    dlsAnalyzer$experimentsOri[[expName]]$setAutocorrelationData()
    
    numberOfDesiredTables <- get_numberOfDesiredTables(nrow(df))
    
    tables                <- split_table_into_list_of_tables(df,numberOfDesiredTables)
    
    renderDlsPlateInfo(expName,tables)

  }
  reactives$data_loaded <- TRUE
})

observeEvent(input$applyChanges,{
  
  req(reactives$data_loaded)
  expNames               <- dlsAnalyzer$experimentNames
  req(length(expNames) > 0) 
  
  reactives$data_loaded  <- FALSE                    
  autocorrelationList    <- dlsAnalyzer$getExperimentProperties('autocorrelationOriginal')
  timeList               <- dlsAnalyzer$getExperimentProperties('time')
  expTabData             <- experimentsSamplesMetaDataFromTables() 
  
  i <- 0
  for (expName in expNames) {
    i <- i+1
    df                    <- getTablesInfo(expTabData[[i]])
    
    #Change name so we can use it in the python class in in other functions
    colnames(df)          <- c("conditions","include","read","scan")
    
    # Update dls sample metadata
    dlsAnalyzer$experimentsOri[[expName]]$sampleInfo <- df
    
  }
  
  # Get metadata of samples
  experimentsSamplesMetaData  <- dlsAnalyzer$getExperimentProperties('sampleInfo')
  expNames                    <- dlsAnalyzer$experimentNames
  
  for (i in 1:length(expNames)) {
    
    autocorrelationData   <- autocorrelationList[[i]]
    expName               <- expNames[i]
    
    maxValues             <- apply(autocorrelationData,2, max)
    selValues             <- maxValues >= input$filterByInitialValue
    
    timeSel <- timeList[[i]] >  input$bumpRemovalTimeLimit/1e6 # To seconds     

    maxValuesAfterTimeLimit <- apply(autocorrelationData[timeSel,],2, max)
    selValues <- selValues & (maxValuesAfterTimeLimit < input$bumpRemovalTolerance)
    
    df                    <- experimentsSamplesMetaData[[i]]

    nMeasurements         <- ncol(dlsAnalyzer$experimentsOri[[expName]]$autocorrelationOriginal)
    numberOfDesiredTables <- get_numberOfDesiredTables(nMeasurements)
    
    df$include            <- df$include & selValues 
    
    # Update dls sample metadata
    dlsAnalyzer$experimentsOri[[expName]]$sampleInfo <- df
    dlsAnalyzer$experimentsOri[[expName]]$setAutocorrelationData()

    tables                <- split_table_into_list_of_tables(df,numberOfDesiredTables)

    renderDlsPlateInfo(expName,tables)
  }
  reactives$data_loaded <- TRUE
})

# Reactive expression to select the DLS data and format it
dlsDataUpdated <- eventReactive(input$updateInfo,{
  
  expNames               <- dlsAnalyzer$experimentNames
  timeList               <- dlsAnalyzer$getExperimentProperties('time')

  autocorrelationList  <- dlsAnalyzer$getExperimentProperties('autocorrelation')
  # Index position of experiments with at least 1 selected curve
  idx                  <- sapply(autocorrelationList, function(x) dim(x)[2] > 0 )
  
  # No data to analyse
  if (sum(idx) == 0) return(NULL)
  
  expNames               <- expNames[idx]
  timeList               <- timeList[idx]
  
  autocorrelationListPredicted <- list()

  # Store the values for optimal alpha (L-curve criteria)
  allOptimalAlpha <- c()
  
  i <- 1
  for (expName in expNames) {
    
    nMeasurements     <- dlsAnalyzer$experimentsOri[[expName]]$nMeasurements
    dlsAnalyzer$experimentsOri[[expName]]$getQ()
    dlsAnalyzer$experimentsOri[[expName]]$createFittingS_space(0.09,1e6,200) 
    
    # In min and max Rh (nm)
    dlsAnalyzer$experimentsOri[[expName]]$getBetaEstimate()
    dlsAnalyzer$experimentsOri[[expName]]$getG1correlation()
    
    if (input$regMethod == 'fixedAlpha') {
      
      alphaRegTerm <- input$alphaRegTerm
      if (alphaRegTerm < 1e-6) {
        shinyalert(text = paste("<b>Sorry! The regularization term should be higher than 1e-6.
                            We will fit using 1e-6.</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
        alphaRegTerm <- 1e-6
      }
      
      if (alphaRegTerm > 5) {
        shinyalert(text = paste("<b>Sorry! The regularization term should be lower than 5.
                            We will fit using 5.</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
        alphaRegTerm <- 5
      }
      
      dlsAnalyzer$experimentsOri[[expName]]$getInitialEstimates(alphaRegTerm,input$maxTime)
    } else {
      ## L curve criteria to select optimal reg parameter
      # Custom alpha vector which should work okayish
      alphaVec <- (5**seq(input$alphaRegTermLcurve_start,
                          input$alphaRegTermLcurve_stop,
                          input$alphaRegTermLcurve_step))**2
      
      # Check that the user doesn't want to perform more than 2000 fittings
      # Otherwise it takes too long...
      nFittings <- length(alphaVec)*nMeasurements
      if (nFittings > 2000) {
        
        shinyalert(text = paste("<b>We can't perform more than 2000 fittings.
                                Based on the number of samples and regularisation 
                                parameters to be tested, you require ",nFittings,
                                " fittings. To reduce this number, select less samples
                                or edit the L-curve criteria parameters from the '
                                Analysis' Box.)</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        
        return(NULL)
      }
      
      dlsAnalyzer$experimentsOri[[expName]]$getInitialEstimatesManyAlpha(alphaVec,input$maxTime)
      dlsAnalyzer$experimentsOri[[expName]]$getOptimalAlphaLcurve()
      dlsAnalyzer$experimentsOri[[expName]]$getInitialEstimatesOptimalAlphaLcurve(input$maxTime)
      
      allOptimalAlpha <- c(allOptimalAlpha,dlsAnalyzer$experimentsOri[[expName]]$optimalAlpha)
      
    }

    #Give time to the fitting to finish
    Sys.sleep(0.02*nMeasurements)
    
    dlsAnalyzer$experimentsOri[[expName]]$predictAutocorrelationCurves()
    dlsAnalyzer$experimentsOri[[expName]]$getMassWeightedContributions()

    autocorrelationListPredicted[[expName]] <- dlsAnalyzer$experimentsOri[[expName]]$autocorrelationPredicted
    
    }
  
  # Retrieve the 2nd order autocorrelation data in an useful format
  data <- dlsData(dlsAnalyzer) 
  
  experimentsSamplesMetaData  <- dlsAnalyzer$getExperimentProperties('sampleInfoRelevant')[idx]
  
  predictedData <- formatDlsPredictedInfoForPlotting(autocorrelationListPredicted,timeList,
                                           expNames,experimentsSamplesMetaData)
  
  estimatedContributions <- dlsAnalyzer$getExperimentProperties('contributionsGuess')[idx]
  
  estimatedContributionsMassWeighted <- dlsAnalyzer$getExperimentProperties(
    'contributionsGuessMassWeighted')[idx]
  
  contributions          <- formatContributions(estimatedContributions,
                                                dlsAnalyzer$experimentsOri[[expNames[1]]]$hrs,
                                                expNames,experimentsSamplesMetaData)
  
  # Contributions as a function of the diffusion coefficients
  contributionsDiff          <- formatContributions(estimatedContributions,
                                                dlsAnalyzer$experimentsOri[[expNames[1]]]$ds,
                                                expNames,experimentsSamplesMetaData,
                                                'diff')
  
  # Mass weighted contributions
  contributionsMass          <- formatContributions(estimatedContributionsMassWeighted,
                                                dlsAnalyzer$experimentsOri[[expNames[1]]]$hrs,
                                                expNames,experimentsSamplesMetaData)

  nChoices <- length(unique(data$variable))
  
  choices <- get_choices_residuals_plot(nChoices)
  updateSelectInput(session,"residualsPlotSelection",choices = choices)
  
  residuals           <- predictedData$value - data$value
  maxR                <- max(abs(residuals))
  minR                <- min(abs(residuals))
  step                <- signif((maxR-minR)/20,2)
  
  updateNumericInput(session, "filterByResidualsValue",
                     label = NULL,
                     value = signif(maxR,2)+step, min = 0, max = 1e6, step = step)
  
  shinyalert(text = 
  paste("<b>Fitting done! Use the estimated standard deviation
  only for comparison purposes.</b>"),
  type = "success",closeOnEsc = T,closeOnClickOutside = T,
  html=T)
  
  # If we used the l criteria
  if (input$regMethod != 'fixedAlpha') {
    allSampleNames  <- do.call(rbind,experimentsSamplesMetaData)$conditions
    allOptimalAlpha <- data.frame(allOptimalAlpha,allSampleNames)
    
    # Get norms to plot the l curve
    expN   <- dlsAnalyzer$experimentNames[idx]
    sNames <- dlsAnalyzer$getExperimentProperties("sampleInfoRelevant")[idx]
    rn     <- dlsAnalyzer$getExperimentProperties("curvesResidualNorm")[idx]
    pn     <- dlsAnalyzer$getExperimentProperties("curvesPenaltyNorm")[idx]
    idSel  <- dlsAnalyzer$getExperimentProperties("alphaOptIdx")[idx]
    
    lCurveData  <- formatNorms(expN,sNames,rn,pn,idSel,alphaVec)
  } else{
    allOptimalAlpha <- NULL
    lCurveData      <- NULL
  }
  
  reactives$renderData <- TRUE
    
  return(list("data"=data,"predictedData"=predictedData,
              "estimatedContributions"=contributions,
              "estimatedContributionsDiff"=contributionsDiff,
              'estimatedContributionsMass'=contributionsMass,
              "residuals"=residuals,
              "allOptimalAlpha"=allOptimalAlpha,
              'lCurveData'=lCurveData))
})

#Get selected variables based on residuals threshold 
selectedVariables <- reactive({
  
  req(dlsDataUpdated())
  
  data <- dlsDataUpdated()
  
  residuals   <- data$residuals 
  dlsOridata  <- data$data
  
  selVariables <- get_selected_variables_based_on_residuals(dlsOridata,residuals,input$filterByResidualsValue)

  return(selVariables)
  
})

estimatedContributions <- reactive({
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  estimatedContributions  <- dlsDataUpdated()$estimatedContributions
  estimatedContributions  <- estimatedContributions[estimatedContributions$variable %in% selectedVariables(),]
  
  return(estimatedContributions)
})

# using diffusion coefficients
estimatedContributionsDiff <- reactive({
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  estimatedContributions  <- dlsDataUpdated()$estimatedContributionsDiff
  estimatedContributions  <- estimatedContributions[estimatedContributions$variable %in% selectedVariables(),]
  
  return(estimatedContributions)
})

# Using mass weighted contributions
estimatedContributionsMass <- reactive({
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  estimatedContributions  <- dlsDataUpdated()$estimatedContributionsMass
  estimatedContributions  <- estimatedContributions[estimatedContributions$variable %in% selectedVariables(),]
  
  return(estimatedContributions)
})

## Render only the peak table info
## or also (in case we applied the l-curve criteria) the values of the optimal alpha

output$ui_regularisationTerm_table <- renderUI({
  
  req(dlsDataUpdated())

  if (input$regMethod == 'fixedAlpha') {
    tabBox(title = "", width = 10,id = "tabsetParamsTable",
           tabPanel("Peak parameters",tableOutput("params_table"))
    )
  } else {
    
    tabBox(title = "", width = 10,id = "tabsetParamsTable",
           tabPanel("Peak parameters",tableOutput("params_table")),
           tabPanel("Regularisation terms",tableOutput("regularisationTerm_table"))
    )
  }
  
})

output$params_table <- renderTable({
  req(estimatedContributions())
  req(reactives$renderData)
  
  estimatedContributions <- estimatedContributions()
  if (input$distribution_type == 'diff') {
    estimatedContributions <- estimatedContributionsDiff()
  }

  return(getPeakInfo(estimatedContributions,
                     estimatedContributionsMass(),
                      input$peakLeftBounds,input$peakRightBounds,
                      input$hrEstimationMethod))
})

# Automatically change peak borders when selecting diffusion coefficients
# or hydrodynamic radii
observeEvent(input$distribution_type,{
  
  req(reactives$data_loaded)
  
  reactives$renderData <- FALSE
  
  leftBounds  <- get_positions(input$peakLeftBounds)
  rightBounds <- get_positions(input$peakRightBounds)
  
  c1 <- any(is.na(leftBounds))
  c2 <- any(is.na(rightBounds))
  c3 <- length(leftBounds) != length(rightBounds)
  
  # Return null if the bounds are not well defined...
  if (c1 | c2 | c3) return(NULL)
  
  # Otherwise, detect if the user wants to analyse 
  # diff coefficients or hydrodynamic radii
  
  # Case 1 - we suppose the user has input hr in the peak box
  # and wants diff coefficients
  if (rightBounds[1] > 1e-4 & input$distribution_type == "diff") {
    
    # Convert to meters
    leftBounds  <- leftBounds  / 1e9
    rightBounds <- rightBounds / 1e9
    
    # Use the temperature of the first experiment
    rightBounds <- diffusion_from_hydrodynamic_radius(
      rightBounds,dlsAnalyzer$experimentsOri[[1]]$temperature,
      dlsAnalyzer$experimentsOri[[1]]$viscosity)
    
    leftBounds <- diffusion_from_hydrodynamic_radius(
      leftBounds,dlsAnalyzer$experimentsOri[[1]]$temperature,
      dlsAnalyzer$experimentsOri[[1]]$viscosity)    
    
    rightBounds <- paste0(signif(rightBounds,2),collapse = " ")
    leftBounds  <- paste0(signif(leftBounds,2), collapse = " ")
    
    # Exchange bounds - remember that diff ~ 1/hr
    updateTextInput(session, "peakRightBounds", label = NULL, value = leftBounds)
    updateTextInput(session, "peakLeftBounds", label = NULL,  value = rightBounds)
    
    shinyalert(
      text = 
      paste("<b>Peak boundaries were automatically 
      adjusted using the temperature and 
      viscosity of the first loaded experiment.</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
    
  }
  
  # Case 2 - we suppose the user has input diff coef values in the peak box
  # and wants hr 
  if (rightBounds[1] < 1e-4 & input$distribution_type == "hr") {
    
    # Use the temperature and viscosity of the first experiment
    rightBounds <- hydrodynamic_radius(
      rightBounds,dlsAnalyzer$experimentsOri[[1]]$temperature,
      dlsAnalyzer$experimentsOri[[1]]$viscosity)
    
    leftBounds <- hydrodynamic_radius(
      leftBounds,dlsAnalyzer$experimentsOri[[1]]$temperature,
      dlsAnalyzer$experimentsOri[[1]]$viscosity)    
    
    # Convert to nanometers
    leftBounds  <- leftBounds  * 1e9
    rightBounds <- rightBounds * 1e9
    
    rightBounds <- paste0(signif(rightBounds,2),collapse = " ")
    leftBounds  <- paste0(signif(leftBounds,2), collapse = " ")
    
    # Exchange bounds - remember that diff ~ 1/hr
    updateTextInput(session, "peakRightBounds", label = NULL, value = leftBounds)
    updateTextInput(session, "peakLeftBounds", label = NULL,  value = rightBounds)
    
    shinyalert(
      text = 
        paste("<b>Peak boundaries were automatically 
      adjusted using the temperature and 
      viscosity of the first loaded experiment.</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
    
  }
  
  Sys.sleep(2)
  reactives$renderData <- TRUE
  
})

output$regularisationTerm_table <- renderTable({
  
  req(dlsDataUpdated()) 
  req(reactives$renderData)
  
  data <- dlsDataUpdated()$allOptimalAlpha
  
  if (is.null(data)) return(NULL) # If it is null, don't proceed
  
  colnames(data) <- c("Selected alpha","Measurement name")
  
  return(data)
},digits = 8)

output$tableForAvg <- renderRHandsontable({
  
  req(dlsDataUpdated())
  req(reactives$renderData)
  
  data <- dlsData(dlsAnalyzer)
  df   <- data.frame('Condition' = unique(data$variable), 'Group'='A')
  
  rTable <- rhandsontable(df,rowHeaders = NULL) %>% 
    hot_col(c(1),readOnly = TRUE) %>% 
    hot_table(stretchH='all') 
  
  return(rTable)
  
})  
  



