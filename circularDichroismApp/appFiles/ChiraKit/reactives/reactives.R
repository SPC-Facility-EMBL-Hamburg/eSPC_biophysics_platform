append_record_to_logbook <- function(record_str,include_time=FALSE) {
  
  if (include_time) {
    record_str <- paste0(as.character(format(Sys.time(),usetz = TRUE)),' ',record_str)
  }
  # Append '' to print in the output a new empty line after the record
  record_str <- c(record_str,'')
    
  reactives$logbook <- append(reactives$logbook, record_str)
  
  return(NULL) 
}

updateProcessingTable <- function(
    operation = 'Sum',operationUnits='millidegrees') {
  
  # Generate table to preprocess spectra (subtract, sum, smooth, average, etc.)
  df <- generateDTtableProcessing(cdAnalyzer,operation,operationUnits)
  
  # To allow the creation of selectInput inside DT table
  session$sendCustomMessage('unbind-DT', 'proccesingInfo')
  output$proccesingInfo <- renderDTtable(df)
  
}

updateCDFilesInfoTable <- function() {
  
  df <- generateDTtable(cdAnalyzer)
  
  # To allow the creation of selectInput inside DT table
  session$sendCustomMessage('unbind-DT', 'cdFilesInfo')
  output$cdFilesInfo <- renderDTtable(df)
  
  Sys.sleep(0.5)
  observeEvent(lapply(1:length(cdAnalyzer$experimentNames), function(i) input[[paste0("inputUnits", i)]] ), {
    
    reactives$data_loaded <- NULL
    convertExperimentToWorkingUnits()
    reactives$data_loaded <- TRUE
    
  },ignoreInit = T,ignoreNULL = T)
  
}

updateMaxVoltageValue <- function() {
  allVoltageData <- unlist(cdAnalyzer$get_experiment_properties('signalHT'))
  maxVoltage     <- max(c(allVoltageData,0),na.rm = T)
  
  reactives$showVoltageThreshold <- maxVoltage > 0
  
  updateNumericInput(session,'maxHTvalue',NULL,value = maxVoltage)
}

output$showVoltageThreshold   <- reactive( { return( reactives$showVoltageThreshold  ) } )
outputOptions(output, "showVoltageThreshold" , suspendWhenHidden = FALSE)

load_one_experiment <- function(cd_data_file,name,inputUnits = 'millidegrees') {
  
  nameOri <- name
  name    <- remove_file_extension(name)

  # Remove ':' which will break our app
  name <- gsub(':','',name)
  
  withBusyIndicatorServer("Go",{
  
    exps      <- cdAnalyzer$experimentNames
    
    if (name %in% exps) {
      shinyalert(text = paste("<b>File name is already being used.</b>"),
                 type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                 html=T)
      return(NULL)
    } else { # Load last loaded file
      
      load <- cdAnalyzer$load_experiment(cd_data_file,name)
      
      # Catch exception when data has wrong format
      if (load != "Data loaded successfully!!!") {
        shinyalert(text = paste("<b>The file was not loaded! Please verify the file format (User guide).</b>"),
                   type = "warning",closeOnEsc = T,closeOnClickOutside = T,
                   html=T)
        return(NULL)
      }
      
      l1        <- list()
      
      # Case 1 - Override default units
      if (inputUnits != 'millidegrees') {
        l1[[1]]   <- inputUnits
        # Case 2 - Use experiments units
      } else {
        l1[[1]]   <- cdAnalyzer$experimentsOri[[name]]$units
      }
      
      names(l1) <- c(name)
      
      # Convert to absorbance
      cdAnalyzer$experiments_to_absorbance_units(l1)
      
      # Convert to desired units
      l1[[1]]   <- input$workingUnits
      cdAnalyzer$experiments_absorbance_units_to_other_units(l1)
      
      metadata_info     <- cdAnalyzer$experimentsOri[[name]]$metadata
      
      tabP <- tabPanel(name,value=name,fluidRow(
        column(12,tableOutput(paste0('metadata_',name)))))
      appendTab("metadata",tabP,select=TRUE)
      
      if (length(metadata_info) > 0) {

        metadataFeature <- (names(metadata_info))
        metadataValue   <- unlist(metadata_info)
        
        metadata_df <- data.frame(metadataFeature,metadataValue)
        colnames(metadata_df) <- c('Metadata feature (read from file)','Value') 
        
      } else {
        metadata_df <- data.frame('Metadata'='No information available')
      }
      
      output[[paste0('metadata_',name)]] <- renderTable({metadata_df})
      Sys.sleep(0.02)
      updateProcessingTable()
      
      Sys.sleep(0.02)
      
    }
    
  })
  
  append_record_to_logbook(paste0('File loaded: ',nameOri,
                                  ' | Input units: ',inputUnits),
                           include_time = TRUE)
  
}

observeEvent(input[["inputUnitsAll"]],{
  req(reactives$data_loaded)
  
  reactives$data_loaded <- NULL
  convertExperimentToWorkingUnits()
  reactives$data_loaded <- TRUE
},ignoreInit = T,ignoreNULL = T)


# Generate and render
# 1) The Input parameters Table,
# 2) The legends Table
# 3) The wavelength range slider
renderInputData <- function() {
  
  output$cdFilesInfo    <- NULL
  
  wls  <- cdAnalyzer$get_experiment_properties('wavelength')
  
  if (length(wls) == 0) return(NULL)
  
  Sys.sleep(length(wls)*0.015)
  
  minV <- min(unlist(wls))
  maxV <- max(unlist(wls))
  
  updateSliderInput(session,inputId="wavelengthRange",label = NULL,min=minV,max=maxV,value = c(minV,maxV))
  updateSelectInput(session,"experiment2delete",choices     = c("ALL",cdAnalyzer$experimentNames))
  updateSelectInput(session,"selected_cd_exp",choices       = cdAnalyzer$experimentNames)
  
  updateCDFilesInfoTable()

  cdAnalyzer$sharedParameters <- FALSE
  updateCheckboxInput(session,'sharedExperimentParameters',label = NULL,FALSE)
  
  legendDf              <- getPlottingDF(cdAnalyzer)
  
  output$legendInfo     <- helperRenderRHandsontable(legendDf)
  
  Sys.sleep(length(wls)*0.015)
  
}

observeEvent(input$cdFiles,{
  
  req(input$cdFiles)
  
  reactives$data_loaded <- NULL
  output$legendInfo     <- NULL
  
  cd_data_files   <- input$cdFiles$datapath
  names           <- input$cdFiles[[1]]
  
  sorted_indices  <- order(names)
  
  cd_data_files   <- cd_data_files[sorted_indices]
  names           <- names[sorted_indices]
  
  i <- 0
  
  # iterate over the files
  for (name in names) {
    i              <- i + 1
    cd_data_file   <- cd_data_files[i]
    
    fileExtension <- getFileNameExtension(name)
    
    if (fileExtension == 'zip') {
      
      # Create a temporary directory to unzip files
      unzipDir <- tempfile()
      dir.create(unzipDir)
      
      unzip(cd_data_file,exdir = unzipDir)
      
      # List all possible types of files
      file_list <- list.files(unzipDir, 
                              pattern = "\\.xlsx$|\\.gen$|\\.dat$|\\.txt$|\\.csv$|d\\d+$",
                              ignore.case = TRUE,
                              recursive = TRUE,
                              full.names = TRUE)
      
      for (file in file_list) {
      
        load_one_experiment(file,basename(file))
        
      }

    } else {
      
      load_one_experiment(cd_data_file,name)
      
    }
  }
  
  # Copy experiments to allow modifying the wavelength range
  cdAnalyzer$initialize_experiment_modif()
  
  updateMaxVoltageValue()
  
  renderInputData()
  
  Sys.sleep(0.5)
  reactives$data_loaded <- TRUE
  
},priority = 10)

# Modal dialog to ask for the input and desired units
observeEvent(input$automaticProcess,{
  
  req(input$cdFilesSample)
  req(input$cdFilesBaseline)
  
  showModal(modalDialog(
    
    tags$h3('Please enter the units of the CD data files,
     the desired units of the final processed spectrum, 
    and if the final spectrum should be zeroed:'),
    selectInput("inputUnitsAutomatic"  , 'CD input files unit',
                global_cd_units_choice_short),
    selectInput("workingUnitsAutomatic", 'Desired units of the final spectrum',
                global_cd_units_choice),
    checkboxInput('zeroFinalSpectrum',
    'Subtract the average signal of the highest 
    10 nm wavelength interval from the final CD spectrum?',
    FALSE),
    
    footer=tagList(
      actionButton('submitUnits', 'Submit'),
      modalButton('Cancel')
    )
  ))
  
})

# Modal dialog to ask for the molecular weight, concentration, path length, and
# number of chromophore units, if required.
observeEvent(input$submitUnits, {
  
  removeModal()
  
  c1 <- grepl('molar'  , input$workingUnitsAutomatic, ignore.case = TRUE)
  c2 <- grepl('unit'   , input$workingUnitsAutomatic, ignore.case = TRUE)
  
  if (c1 & c2) {
    
    showModal(modalDialog(
      
      tags$h3('Please enter the following experimental details:'),
      numericInput('pathLength'   , 'Path length (mm)',     1,min = 0,max = 100,step = 1),
      numericInput('molWeight'    , 'Molecular weight (Da)',1,min = 0,max = 1e6,step = 1),
      numericInput('concentration', 'Concentration (mg/ml)',1,min = 0,max = 100,step = 0.2),
      numericInput('numberOfC'   ,  'Number of chromophore units'   ,1,min = 0,max = 1e6,step = 1),
      footer=tagList(
        actionButton('submitExperimentalDetails', 'Submit'),
        modalButton('Cancel')
      )
    ))
    
  } else if (c1 & !c2) {
    
    showModal(modalDialog(
      
      tags$h3('Please enter the following experimental details:'),
      numericInput('pathLength'   , 'Path length (mm)',     1,min = 0,max = 100,step = 1),
      numericInput('molWeight'    , 'Molecular weight (Da)',1,min = 0,max = 1e6,step = 1),
      numericInput('concentration', 'Concentration (mg/ml)',1,min = 0,max = 100,step = 0.2),
      footer=tagList(
        actionButton('submitExperimentalDetails', 'Submit'),
        modalButton('Cancel')
      )
    ))
    
  } else {
    
    selectedIDx <- which(global_cd_units_choice == input$workingUnitsAutomatic)
    
    showModal(modalDialog(
      
      tags$h3(paste('Selected units for the final spectrum are: ' ,names(global_cd_units_choice)[selectedIDx])),
      footer=tagList(
        actionButton('submitExperimentalDetails', 'Okay'),
        modalButton('Cancel')
      )
    ))    
    
  } 
  
})

# Generate programatically one CD experiment,
# Add this experiment to the python class cdAnalyzer
append_cd_experiment <- function(expName,spectraNames,wlNew,signalNew,inputUnits,
                                 metadata_list,temperature=NULL,
                                 signalHT = NULL,isFakeExperiment=FALSE,
                                 fakeExperimentSignal=NULL) {
  
  signalNew <- as.matrix(signalNew)

  cdAnalyzer$experimentsOri[[expName]]                   <- CdExperimentGeneral()
  cdAnalyzer$experimentsOri[[expName]]$wavelength        <- np_array(wlNew)
  cdAnalyzer$experimentsOri[[expName]]$spectraNames      <- np_array(c(spectraNames))
  cdAnalyzer$experimentsOri[[expName]]$internalID        <- np_array(c(spectraNames))
  cdAnalyzer$experimentsOri[[expName]]$signalInput       <- np_array(signalNew)

  if (!is.null(temperature)) {
    cdAnalyzer$experimentsOri[[expName]]$temperature     <- np_array(temperature)
  }
  
  na_matrix <- matrix(NA, nrow = nrow(signalNew), ncol = ncol(signalNew))
  if (is.null(signalHT)) {
    cdAnalyzer$experimentsOri[[expName]]$signalHT          <- na_matrix
  } else {
    cdAnalyzer$experimentsOri[[expName]]$signalHT          <- np_array(signalHT)
  }
  
  # Assign the desired signal
  # For fake experiments, it can't be changed once it's created!!
  if (isFakeExperiment) {
    cdAnalyzer$experimentsOri[[expName]]$signalDesiredUnit    <- np_array(signalNew)
    cdAnalyzer$experimentsOri[[expName]]$fakeExperimentSignal <- fakeExperimentSignal
    cdAnalyzer$experimentsOri[[expName]]$signalAbs            <- na_matrix
    # Trick to avoid further processing with other spectra using
    # other operation units different from 'fakeExperimentSignal'
    cdAnalyzer$experimentsOri[[expName]]$concentration        <- -1
  } else {
    cdAnalyzer$experimentsOri[[expName]]$experiment_to_absorbance_units(inputUnits)
  }
  
  cdAnalyzer$experimentsOri[[expName]]$isGenerated       <- TRUE
  
  newExperimentNames <- c(expName)
  
  cdAnalyzer$experimentsOri[[expName]]$isFakeExperiment  <- isFakeExperiment
  
  # Start the lines vector, to fill the output file
  
  lines <- c(paste0('New generated experiment: ',expName))
  lines <- c(lines,
             paste0('Number of generated CD curves:  ',length(spectraNames)))
  
  if (!is.null(metadata_list)) {
    
    metadataFeature   <- names(metadata_list)
    metadataFeature   <- paste0('#',metadataFeature,' :')
    metadataFeature   <- identate(metadataFeature,48)
    metadataValue     <- unlist(metadata_list)
    metadata_df       <- data.frame(metadataFeature,metadataValue)
    
    for (i in 1:nrow(metadata_df)) {
      lines <- c(lines,paste0(metadata_df[i,],collapse = ''))
    }
  }
  
  append_record_to_logbook(lines)

  # end of add new experiment to the logbook
  
  metadata_list[['Timestamp']]       <- format(Sys.time(),usetz = TRUE)
  metadata_list[['CD curve origin']] <- 'ChiraKit (https://spc.embl-hamburg.de/)'
  
  cdAnalyzer$experimentsOri[[expName]]$metadata          <- metadata_list
  
  cdAnalyzer$experimentNames <- c(cdAnalyzer$experimentNames,expName)
  
  return(NULL)
}

observeEvent(input$submitExperimentalDetails,{
  
  req(input$submitExperimentalDetails)
  
  removeModal()
  
  sample_data_files     <- input$cdFilesSample$datapath
  scans_sample_names    <- sapply(input$cdFilesSample[[1]],remove_file_extension)
  
  baseline_data_files   <- input$cdFilesBaseline$datapath
  scans_bl_names        <- sapply(input$cdFilesBaseline[[1]],remove_file_extension)
  
  # Check that all the files are different
  intersectFiles <- intersect(scans_sample_names,scans_bl_names)
  
  if (length(intersectFiles) > 0) {
    
    shinyalert(text = paste("<b>Please use different files for the sample and baseline</b>"),
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    
    return(NULL)
  }
  
  reactives$data_loaded      <- NULL

  names <- c(scans_sample_names,scans_bl_names)
  files <- c(sample_data_files,baseline_data_files)
  
  i <- 0
  # Iterate over the files
  for (name in names) {
    
    i <- i + 1
    cd_data_file   <- files[i]

    load_one_experiment(cd_data_file,name,input$inputUnitsAutomatic)

  }

  updateMaxVoltageValue()
  ## Average sample files

  # List of signal dataframes 
  signalSel <- lapply(scans_sample_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalInput})
  
  htSel <- lapply(scans_sample_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalHT})
  
  # List of temperatures
  temperatureSel <- lapply(scans_sample_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$temperature
  })
  
  # Assign the mean of the sample scans temperature to the rest of the curves
  temps          <- unlist(temperatureSel)
  
  if (all(is.null(temps))) {
    newTemp <- NA
  } else {
    newTemp        <- mean(temps,na.rm = T)
  }
  
  signalSelVectorized <- list()
  htSelVectorized     <- list()
  wlSelVectorized     <- list()

  spectra_counter    <- 0
  experiment_counter <- 0
  
  # Transform list of signal dataframes into a list of signal vectors
  for (signalDF in signalSel) {
    
    experiment_counter <- experiment_counter + 1
    expName            <- names[experiment_counter]
    htDF               <- htSel[[experiment_counter]]
      
    for (columnID in 1:ncol(signalDF)) {
      spectra_counter                        <- spectra_counter + 1
      signalSelVectorized[[spectra_counter]] <- signalDF[,columnID]
      htSelVectorized[[spectra_counter]]     <- htDF[,columnID]
      wlSelVectorized[[spectra_counter]]     <- cdAnalyzer$experimentsOri[[expName]]$wavelength
    }
    
  }
  
  interpolatedAverage <- averageListOfVectors(wlSelVectorized,signalSelVectorized)
  signalNewSample     <- matrix(interpolatedAverage$average,ncol = 1)
  
  wlNewSample         <- interpolatedAverage$x
  htNewSample         <- get_ht_average(wlSelVectorized,htSelVectorized)
  
  # Find if we have already used the spectra names we want to use 
  currentSpectraNames <- getSpectraNames(cdAnalyzer)
    
  time_str <- get_hour_minute_sec()
  
  sampleAvgName <- 'Sample average'
  if ('Sample average' %in% currentSpectraNames) {
    sampleAvgName <- paste0('Sample average ',time_str)
  } 
  
  processingUnitsString <- clean_html_sup_tag(
    workingUnits2ProperLabel(input$inputUnitsAutomatic)
  )
  
  metadata_list <- list(
    'Processing'           = 'Generated by averaging many CD curves',
    'Original CD files'    = paste0(scans_sample_names,collapse = ' '),
    'Units for processing' = processingUnitsString
  )
  
  append_cd_experiment(sampleAvgName,sampleAvgName,wlNewSample,
                       signalNewSample,input$inputUnitsAutomatic,
                       metadata_list,newTemp,htNewSample)
  
  ## End of average sample files
  
  ## Average baseline files
  
  # List of signal dataframes 
  signalSel <- lapply(scans_bl_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalInput})
  
  htSel    <- lapply(scans_bl_names, function(exp) {
    cdAnalyzer$experimentsOri[[exp]]$signalHT})
  
  signalSelVectorized <- list()
  htSelVectorized     <- list()
  wlSelVectorized     <- list()
  
  spectra_counter     <- 0

  n_loaded_experiments <- experiment_counter
  
  # Transform list of signal dataframes into a list of signal vectors
  for (signalDF in signalSel) {
    
    experiment_counter <- experiment_counter + 1
    expName            <- names[experiment_counter]
    htDF               <- htSel[[experiment_counter - n_loaded_experiments]]
    
    for (columnID in 1:ncol(signalDF)) {
      spectra_counter                        <- spectra_counter + 1
      signalSelVectorized[[spectra_counter]] <- signalDF[,columnID]
      htSelVectorized[[spectra_counter]]     <- htDF[,columnID]
      wlSelVectorized[[spectra_counter]]     <- cdAnalyzer$experimentsOri[[expName]]$wavelength
    }
    
  }
  
  interpolatedAverage <- averageListOfVectors(wlSelVectorized,signalSelVectorized)
  signalNewBaseline   <- matrix(interpolatedAverage$average,ncol = 1)
  wlNewBaseline       <- interpolatedAverage$x
  
  htNewBaseline       <- get_ht_average(wlSelVectorized,htSelVectorized)
  
  baselineAvgName <- 'Baseline average'
  if ('Baseline average' %in% currentSpectraNames) {
    baselineAvgName <- paste0('Baseline average ',time_str)
  }
  
  metadata_list <- list(
    'Processing'           = 'Generated by averaging many CD curves',
    'Original CD files'    = paste0(scans_sample_names,collapse = ' '),
    'Units for processing' = processingUnitsString
  )
  
  append_cd_experiment(baselineAvgName,baselineAvgName,wlNewBaseline,
                       signalNewBaseline,input$inputUnitsAutomatic,
                       metadata_list,newTemp,htNewBaseline)
  
  ## End of average baseline files
  
  ## Subtract the baseline
  
  interpolatedSignals <- interpolateVectors(wlNewSample,signalNewSample,wlNewBaseline,signalNewBaseline)
  
  signalNew <- interpolatedSignals$y1 - interpolatedSignals$y2
  wlNew     <- interpolatedSignals$x1
  
  metadata_list <- list(
    'Processing'           = 'Generated by subtracting the sample average from the baseline average',
    'Original CD files_sample'    = paste0(scans_sample_names,collapse = ' '),
    'Original CD files_baseline'  = paste0(scans_bl_names,collapse = ' '),
    'Units for processing'        = processingUnitsString
  )
  
  baselineSubName <- 'Baseline subtracted'
  if ('Baseline subtracted' %in% currentSpectraNames) {
    baselineSubName <- paste0('Baseline subtracted ',time_str)
  }
  
  append_cd_experiment(baselineSubName,baselineSubName,wlNew,
                       signalNew,input$inputUnitsAutomatic,
                       metadata_list,newTemp,htNewSample)
  
  ## End of subtract the baseline
  
  ## Zeroing to obtain the final spectra
  ## TO DO - add modal dialog!!
  zeroFinalSpectrum <- input$zeroFinalSpectrum
  if (zeroFinalSpectrum) {
    
    maxV <- max(unlist(wlSelVectorized))
    
    signalNew <- signalNew - mean(signalNew[wlNew > (maxV - 10)])
    
    wlRange <- paste0(maxV - 10,' - ',maxV, ' nm')
    
    metadata_list <- list(
      'Processing' = 'Generated by 1) subtracting the sample average from the baseline average and 2) subtracting the mean signal from the last 10 nm',
      'Original CD files of the sample'    = paste0(scans_sample_names,collapse = ' '),
      'Original CD files of the baseline'  = paste0(scans_bl_names,collapse = ' '),
      'Units for processing'               = processingUnitsString,
      'Wavelength range to calculate offset' = wlRange
    )
    
    baselineZeroedName <- 'Baseline subtracted with offset'
    if ('Baseline subtracted with offset' %in% currentSpectraNames) {
      baselineZeroedName <- paste0('Baseline subtracted with offset ',time_str)
    }
    
    append_cd_experiment(baselineZeroedName,baselineZeroedName,wlNew,
                         signalNew,input$inputUnitsAutomatic,
                         metadata_list,newTemp,htNewSample)
    
  }
  
  ## End of zeroing to obtain the final spectra
  
  newExperimentNames <- c(baselineAvgName,sampleAvgName,baselineSubName)
  
  if (zeroFinalSpectrum) {
    newExperimentNames <- c(newExperimentNames,baselineZeroedName)
  }

  # We need to use lapply instead of a for-loop to avoid R lazy evaluation!
  lapply(newExperimentNames, function (exp) {
    
    saneExpName <- gsub(':','',exp)
    
    tabP <- tabPanel(saneExpName,value=saneExpName,fluidRow(column(12,tableOutput(paste0('metadata_',saneExpName)))))
    
    appendTab("metadata",tabP,select=TRUE)
    
    metadata_info     <- cdAnalyzer$experimentsOri[[exp]]$metadata
    metadataFeature   <- (names(metadata_info))
    metadataValue     <- unlist(metadata_info)
    
    metadata_df <- data.frame(metadataFeature,metadataValue)
    colnames(metadata_df) <- c('Metadata feature','Value') 
    
    output[[paste0('metadata_',saneExpName)]] <- renderTable({metadata_df})
    return(NULL)
    
  })
  
  # Add the information about the molecular weight, concentration, path length
  # and number of chromophore units, if required.
  
  c1 <- grepl('molar'  , input$workingUnitsAutomatic, ignore.case = TRUE)
  c2 <- grepl('unit'   , input$workingUnitsAutomatic, ignore.case = TRUE)
  
  cdAnalyzer$initialize_experiment_modif()
  
  # Spectra to convert into molar ellipticity / extinction
  experimentsToModifyParams    <- c(baselineSubName)
  
  if (zeroFinalSpectrum) {
    experimentsToModifyParams    <- c(experimentsToModifyParams,baselineZeroedName)
  }
 
  if (c1) {
    
    for (name in experimentsToModifyParams) {
      
      cdAnalyzer$set_experiment_properties(name,'concentration',   input$concentration)
      cdAnalyzer$set_experiment_properties(name,'molecularWeight', input$molWeight)
      cdAnalyzer$set_experiment_properties(name,'pathLength',      input$pathLength / 10) # Divide by 10 (from mm to cm)
      
      append_record_to_logbook(paste0('Setting concentration (mg/ml) of '    ,name,' to ',input$concentration))
      append_record_to_logbook(paste0('Setting molecular weight (Dalton) of ',name,' to ',input$molWeight))
      append_record_to_logbook(paste0('Setting path length (cm) of '         ,name,' to ',input$pathLength / 10))
      
    }
    
  }
  
  if (c2) {
    
    for (name in experimentsToModifyParams) {
      
      cdAnalyzer$set_experiment_properties(name,'numberOfCroms', input$numberOfC)
      append_record_to_logbook(paste0('Setting chromophore units of '    ,name,' to ',input$numberOfC))
    }
  }
  
  # Convert to selected units
  
  updateSelectInput(session,'workingUnits',NULL,
                    choices = getChoices(input$workingUnitsAutomatic))
  
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input$workingUnitsAutomatic)
  
  Sys.sleep(0.5)
  renderInputData()
  
  Sys.sleep(0.5)
  reactives$data_loaded     <- TRUE
  
})


observeEvent(input$triggerDeletion,{
  req(reactives$data_loaded)

  selected <- input$experiment2delete
  
  # Do not delete anything
  if (selected == "None") return(NULL)
  # Delete all experiments
  if (selected == "ALL") selected <- cdAnalyzer$experimentNames
  
  # Delete the legend table
  output$legendInfo <- NULL
  cdAnalyzer$delete_experiment(selected)
  
  # Remove the metadata Tabs
  for (sel in c(selected)) {
    
    target <-  gsub(':','',sel) # Remove the ':' character which can't be used in the Shiny Tabs names
    removeTab(inputId = "metadata", target = target)
    
    append_record_to_logbook(paste0('Deleting experiment ',sel))
    
  }
  
  exps <- cdAnalyzer$experimentNames

  if (length(exps) > 0) {
    
    updateSelectInput(session,"experiment2delete",selected = "ALL",
                      choices = c("ALL",cdAnalyzer$experimentNames))
    
    updateSelectInput(session,"selected_cd_exp",choices     = cdAnalyzer$experimentNames)
    
    # Create the legend table
    legendDf              <- getPlottingDF(cdAnalyzer)
    output$legendInfo     <- helperRenderRHandsontable(legendDf)
    
    # If the parameters are not shared, then
    # update the Table with the parameters MW, #AA, etc.
    if (!cdAnalyzer$sharedParameters) {
      
      updateCDFilesInfoTable()
      
    }
      
    # Update the Table with the possible spectra to be processed
    updateProcessingTable()
    
  } else {
    output$cdFilesInfo    <- NULL
    output$proccesingInfo <- NULL
    
    updateSelectInput(session,"experiment2delete",choices = c('None'))
    reactives$data_loaded <- FALSE
  }
  
})

convertExperimentToWorkingUnits <- function() {
  
  # Convert everything (except 'fake' experiments) to absorbance units (from input units)
  namesAll  <- cdAnalyzer$experimentNames
  
  trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
  names             <- namesAll[trueExperiments]
  namesN            <- length(names)
  
  # Check if the input units are shared and convert to absorbance
  if (cdAnalyzer$sharedParameters) {
    
    cdAnalyzer$all_experiments_to_absorbance_units(input[["inputUnitsAll"]])
    append_record_to_logbook(paste0('Setting input units to ',input[["inputUnitsAll"]]))
    
  } else {
    
    rList <- list()
    for(i in 1:namesN) {
      selectedInputUnits <- input[[(paste0("inputUnits", i))]]
      req(selectedInputUnits)
      rList[[i]]         <- selectedInputUnits
    }
    if (length(rList) != namesN) return(NULL)
    
    names(rList) <- names
    
    modifiedExperiments <- cdAnalyzer$experiments_to_absorbance_units(rList)
    append_record_to_logbook(paste0('Setting input units of ',modifiedExperiments))
  }
  
  # Convert everything from absorbance units to selected units
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input[["workingUnits"]])
  
  append_record_to_logbook(paste0('Converting CD units to: ',input[["workingUnits"]]))
  
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])

}

# React to changes in molecular weight, path length or concentration 
observeEvent(input$cdFilesInfo_cell_edit, {

  req(reactives$data_loaded)
  
  reactives$data_loaded <- NULL
  info <- input$cdFilesInfo_cell_edit

  columnIndex <- info$col[1] + 1 # Start index from 1
  
  rowIndex    <- info$row[1]
  value       <- as.numeric(info$value[1])

  if (cdAnalyzer$sharedParameters) {
    
    expNames          <- cdAnalyzer$experimentNames
    trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
    expNames          <- expNames[trueExperiments]
  } else {
    expNames <- c(cdAnalyzer$experimentNames[rowIndex])
  }
  
  for (expName in expNames) {
    
    valueNew <- value
    if (columnIndex == 2) {
      selectedVar <- "molecularWeight"
      selectedVar2print <- 'molecular weight (Da)'
    } else if (columnIndex == 3) {
      selectedVar <- "numberOfCroms"
      selectedVar2print <- 'chromophore units number'
    } else if (columnIndex == 4) {
      selectedVar <- "concentration"
      selectedVar2print <- 'concentration (mg/ml)'
    } else { 
      selectedVar  <- "pathLength"
      selectedVar2print <- 'path length (cm)'
      valueNew <- valueNew / 10 # To cm!
    } 
    
    cdAnalyzer$set_experiment_properties(expName,selectedVar,valueNew)
    
    append_record_to_logbook(paste0('Setting ',selectedVar2print,' of ',
                                    expName,' to ',valueNew))
    
  }
  
  reactives$data_loaded <- NULL
  convertExperimentToWorkingUnits()
  reactives$data_loaded <- TRUE
})

observeEvent(input$workingUnits,{
  
  req(reactives$data_loaded)
  reactives$data_loaded <- NULL
  # Convert everything from absorbance units to selected units
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input[["workingUnits"]])
  
  append_record_to_logbook(paste0('Converting CD units to: ',input[["workingUnits"]]))
  
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])

  reactives$data_loaded <- TRUE
  
})

observeEvent(input$sharedExperimentParameters,{
  
  req(reactives$data_loaded)
  reactives$data_loaded <- FALSE
  
  cdAnalyzer$sharedParameters <- input$sharedExperimentParameters
  
  # Force CD samples to share the same molecular weight, concentration, path length, input units and concentration
  if (input$sharedExperimentParameters) {
    
    df   <- data.frame('All experiments',1,1,1,
                       as.character(selectInput('inputUnitsAll',label=NULL, choices = getChoices('Millidegrees'))),
                       1)
    
    colnames(df) <- c('File name',"Mol. weight (Dalton)","#Chromophore units","Conc. (mg/ml)",
                      "Input units","Path length (mm)")
    
    # Set default values
    exps              <- cdAnalyzer$experimentNames
    trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
    exps              <- exps[trueExperiments]
    
    for (exp in exps) {
      
      cdAnalyzer$set_experiment_properties(exp,'units'            , 'millidegrees')
      cdAnalyzer$set_experiment_properties(exp,'numberOfCroms'    , 1)
      cdAnalyzer$set_experiment_properties(exp,'concentration'    , 1)
      cdAnalyzer$set_experiment_properties(exp,'pathLength'       , 0.1)
      cdAnalyzer$set_experiment_properties(exp,'molecularWeight'  , 1)
      
      append_record_to_logbook(paste0('Setting number of chromophore units / concentration (mg/ml)/ molecular weight (Da) of ',
                                      exp,' to 1'))
      
      append_record_to_logbook(paste0('Setting path length (cm) of: ',
                                      exp,' to 0.1'))
      
    }
    
    cdAnalyzer$all_experiments_to_absorbance_units('millidegrees')
    
    append_record_to_logbook('Setting all experiment input units to millidegrees')
    
                       
  } else {
    
    df <- generateDTtable(cdAnalyzer)
  
  }
  
  cdAnalyzer$all_experiments_absorbance_units_to_other_units(input$workingUnits)

  append_record_to_logbook(paste0('Converting CD units to: ',input[["workingUnits"]]))
  
  # To allow the creation of selectInput inside DT table
  session$sendCustomMessage('unbind-DT', 'cdFilesInfo')
  output$cdFilesInfo <- renderDTtable(df)
  reactives$data_loaded <- TRUE
  
})

observeEvent(input$legendInfo,{
  
  req(reactives$data_loaded)

  updateSelectInput(session,"mol2changeColor",NULL,
                    get_legend_from_rhandTable(input$legendInfo),input$mol2changeColor)
  

})

observeEvent(input$legendInfo$changes$changes, {
  
  # Get the changes made to the table
  changes <- input$legendInfo$changes$changes[[1]]
  
  # Check if any changes were made in the column with the legends
  if (any(changes[[2]] == 2)) {

    names          <- cdAnalyzer$experimentNames
    internalIDs    <- cdAnalyzer$get_experiment_properties('internalID')
    spectraID      <- get_ID_from_rhandTable(input$legendInfo)[changes[[1]]+1]
    newSpectraName <- get_legend_from_rhandTable(input$legendInfo)[changes[[1]]+1]
    
    id                 <- which(unlist(internalIDs) == spectraID)
    names(internalIDs) <- names
    id                 <- found_ids(internalIDs,c(id))
    
    append_record_to_logbook(paste0("Changing spectrum name from '",
                                    cdAnalyzer$experimentsOri[[names(id)[1]]]$spectraNames[id[1]],
                                    "' to '",newSpectraName,"'"))
    
    cdAnalyzer$experimentsOri[[names(id)[1]]]$spectraNames[id[1]]   <- newSpectraName
    cdAnalyzer$experimentsModif[[names(id)[1]]]$spectraNames[id[1]] <- newSpectraName
    
  }
  
})

observeEvent(input$colorForLegend,{
  
  req(reactives$data_loaded)
  legendDf <- getLegendDF(input$legendInfo)
  
  idx                  <- which(legendDf$Legend == input$mol2changeColor)
  legendDf$Color[idx]  <- input$colorForLegend
  output$legendInfo    <- helperRenderRHandsontable(legendDf)

})

observeEvent(input$triggerProcessing,{
  
  req(reactives$data_loaded)
  # Retrieve the experiment names 
  names          <- cdAnalyzer$experimentNames
  
  # Retrieve the internal IDs
  internalIDs        <- cdAnalyzer$get_experiment_properties('internalID')
  names(internalIDs) <- names
  
  # Retrieve the user selected label for the new spectrum / spectra
  newSpectraName <- input$newSpectraName
  
  if (newSpectraName %in% names) {
    newSpectraName <- paste0(newSpectraName,' ',get_hour_minute_sec())
  }
  
  # Find the second spectrum label, the desired operation and the units of operation
  spectra2       <- input$inputSpectra2
  operation      <- input$operation
  operationUnits <- input$operationUnits
  
  legendDf <- getLegendDF(input$legendInfo)
  # Find the spectrum / spectra to be processed
  if (input$inputSpectra1 == "Selected in 'Show' Column") {
    
    # Sanity check - Stop if we haves less than 2 spectra and we selected the option "Selected in 'Show' Column"
    if (sum(legendDf$Show) < 2) {
      
      shinyalert(text = paste(
        "<b>Please select at least two spectra to use this option.</b>"),
        type = "info",closeOnEsc = T,closeOnClickOutside = T,
        html=T)
      
      updateCDFilesInfoTable()
      updateProcessingTable()
      
      return(NULL)
      
    }

    # Spectra from 'Show' legend
    spectra1 <- legendDf$Internal.ID[legendDf$Show]
    
  } else {
    
    spectra1 <- c(input$inputSpectra1)
    
  } 
  
  # Start of -- 1st step Find the first and second spectrum IDs
  
  id1s <- c()
  
  for (spectrum in spectra1) {
    id1s <- c(id1s,which(unlist(internalIDs) == spectrum))
  }
  
  # Find the second spectrum, if required
  
  if (operation %in% c('Sum','Subtract')) {
    id2                <- which(unlist(internalIDs) == spectra2)
    ids                <- found_ids(internalIDs,c(id1s,id2))
  } else {
    ids                <- found_ids(internalIDs,c(id1s))
  }
  
  # Find the experiment names that are required for the processing
  selected_id <- names %in% (unique(names(ids)))  

  # End of -- 1st step Find the first and second spectrum IDs

  # Decide if the generated experiments are of the type 'fake'
  # In other words, the calculated spectra will not be available for further processing
  operationIsInMolarUnits       <- grepl('molar', operationUnits, ignore.case = TRUE)
  operationIsInMeanUnits        <- grepl('mean' , operationUnits, ignore.case = TRUE)
  
  molWeights   <- (unlist(cdAnalyzer$get_experiment_properties('molecularWeight'  )))[selected_id]
  pathLengths  <- (unlist(cdAnalyzer$get_experiment_properties('pathLength'       )))[selected_id]
  concs        <- (unlist(cdAnalyzer$get_experiment_properties('concentration'    )))[selected_id]
  numberOfCU   <- (unlist(cdAnalyzer$get_experiment_properties('numberOfCroms' )))[selected_id]
  
  ##  Check that we can process normalized CD units
  
  mw_is_not_zero <- molWeights   != 0 & !is.na(molWeights)  
  pl_is_not_zero <- pathLengths  != 0 & !is.na(pathLengths) 
  cn_is_not_zero <- concs        != 0 & !is.na(concs)       
  nr_is_not_zero <- numberOfCU   != 0 & !is.na(numberOfCU) 
  
  check_list <- c()
  
  if (operationIsInMolarUnits) {
    
    check_list <- c(check_list,mw_is_not_zero,pl_is_not_zero,cn_is_not_zero)
    
  }
  
  if (operationIsInMeanUnits) {
    
    should_stop_process <- c(check_list,nr_is_not_zero)  
    
  }
  
  should_stop_process <- !all(check_list)
  
  if (should_stop_process) {
    
    shinyalert(text = 
    paste("<b>Please add information about the experimental parameters.
          For example, the path length can not be zero.</b>"),
    type = "warning",closeOnEsc = T,closeOnClickOutside = T,
    html=T)
    
    updateCDFilesInfoTable()
    updateProcessingTable()
    
    return(NULL)
    
  }

  # Special check for fake experiments which can only be processed using the same units
  # that were used to generate them
  areSelSpecFake <- (unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment')))[selected_id]
  fakeSignals     <- (unlist(cdAnalyzer$get_experiment_properties('fakeExperimentSignal')))[selected_id]
  
  if (any(areSelSpecFake)) {
    
    fakeSignals <- fakeSignals[fakeSignals != 'unknown']
    
    if (any(fakeSignals != operationUnits)) {
      
      shinyalert(text = 
      paste("<b>The desired operation units can't be used
      with the selected spectra. Please check the User Guide 
            ('Processing' Section).</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
      
      updateCDFilesInfoTable()
      updateProcessingTable()
      
      return(NULL)
    }
  }
  
  append_record_to_logbook(paste0('A processing step was run with the following options: [',
                                  paste0(spectra1,collapse = ' ; '),'] [',
                                  operation,'] [',
                                  operationUnits,'] [',
                                  spectra2,']',
                                  '[Linear interpolation mode: ',input$allowLinearInterpolation,']'))
  
  # Check if they are equal and decide if we can further process 
  # the generated spectra or not ... 
  
  isFake <- length(c(unique(molWeights),unique(pathLengths),
                      unique(concs),unique(numberOfCU))) != 4    
  
  ## End of check that we can process normalized CD units
  
  # 2nd step - convert to the appropriate operation units
  rList2        <- as.list(rep(operationUnits,length(unique(names(ids)))))
  names(rList2) <- unique(names(ids))
  
  cdAnalyzer$experiments_absorbance_units_to_other_units(rList2)

  # 3rd step - Retrieve the desired signals and process them
  
  if (operation %in% c('Sum','Subtract')) {
    
    id1s  <- head(ids, -1)        # Remove the last element (spectra #2), if required
    exp1s <- head(names(ids), -1) # Remove the last element (spectra #2), if required
    
  } else {
    id1s  <- ids
    exp1s <- names(ids)
  }
  
  # List of wavelength vectors - spectra 1
  wlSel     <- lapply(exp1s, function(exp) cdAnalyzer$experimentsOri[[exp]]$wavelength)
  
  # List of signal vectors - spectra 1
  signalSel <- lapply(1:length(exp1s), function(i) {
    cdAnalyzer$experimentsOri[[exp1s[i]]]$signalDesiredUnit[,id1s[i]]
  })
  
  htSel  <- lapply(1:length(exp1s), function(i) {
    cdAnalyzer$experimentsOri[[exp1s[i]]]$signalHT[,id1s[i]]
  })
  
  # List of temperature - spectra 1
  temperatureSel <- lapply(1:length(exp1s), function(i) {
    cdAnalyzer$experimentsOri[[exp1s[i]]]$temperature[id1s[i]]
  })
  
  # Check if all the wavelength vectors are equal
  if (!input$allowLinearInterpolation) {
    
    wlToCheck <- wlSel
    
    if (operation %in% c('Sum','Subtract')) {
      
      # Add  the second spectrum to the wavelength vectors list which is going to be checked
      wlToCheck[[length(wlToCheck)+1]] <- cdAnalyzer$experimentsOri[[last(names(ids))]]$wavelength
      
    }
    
    all_equal <- all(sapply(wlToCheck[-1], function(vec) identical(wlToCheck[[1]], vec)))
    
    # Abort processing if the wavelength vectors are not equal
    if (!all_equal) {
      
      shinyalert(text = 
      paste("<b>Processing aborted: The selected spectra were measured 
      at different wavelengths. 
      To continue you need to allow linear interpolation 
            (not recommended for baseline/reference subtraction).</b>"),
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,
      html=T)
      
      updateCDFilesInfoTable()
      updateProcessingTable()
      return(NULL)
    }
      
  }
    
  # Change data_loaded to NULL to stop other reactives
  reactives$data_loaded <- NULL   
  
  processingUnitsString <- clean_html_sup_tag(
    workingUnits2ProperLabel(operationUnits)
  )
  
  # Case 1 - we generated an average spectrum
  
  if (operation == 'Average') {
    
    interpolatedAverage <- averageListOfVectors(wlSel,signalSel)
    signalNew           <- matrix(interpolatedAverage$average,ncol = 1)
    wlNew               <- interpolatedAverage$x
    
    htNew <- get_ht_average(wlSel,htSel)
    
    # 5th step (Case 1) - Create one new experiment per selected experiment
    newExperimentNames <- c(newSpectraName)
    legendDf <- appendRowsToLegendDataFrame(legendDf,c(newSpectraName))
    
    metadata_list <- list(
      'Processing'           = 'Generated by averaging many CD curves',
      'Original CD curves'   = paste0(spectra1,collapse = ' '),
      'Units for processing' = processingUnitsString
      )
    
    newTemp              <- mean(unlist(temperatureSel))
    
    append_cd_experiment(newSpectraName,newSpectraName,wlNew,signalNew,
                         operationUnits,metadata_list,newTemp,
                         htNew,isFake,operationUnits)
    
    # Case 2 - Batch average
  } else if (operation == 'Batch average') {
    
    # In this case, spectra 2 should be the string 'N = ' followed by a number
    n <- as.numeric(strsplit(spectra2,'=')[[1]][2])
    
    total_batches <- floor(length(spectra1) / n)
    newExperimentNames <- c()
    
    for (i in 1:total_batches) {
      
      newSpectraNameBatch <- paste0('Batch ',i,' ',newSpectraName)

      if (newSpectraNameBatch %in% names) {
        newSpectraName <- paste0(newSpectraNameBatch,' ',get_hour_minute_sec())
      }
      
      wlSelSubset          <- wlSel[(1+n*(i-1)):(n*i)]
      signalSelSubset      <- signalSel[(1+n*(i-1)):(n*i)]
      htSelSubset          <- htSel[(1+n*(i-1)):(n*i)]
      temperatureSelSubset <- temperatureSel[(1+n*(i-1)):(n*i)]
      newTemp              <- mean(unlist(temperatureSelSubset))
        
      interpolatedAverage  <- averageListOfVectors(wlSelSubset,signalSelSubset)
      signalNew            <- matrix(interpolatedAverage$average,ncol = 1)
      wlNew                <- interpolatedAverage$x
      
      htNew <- get_ht_average(wlSelSubset,htSelSubset)
      
      # 5th step (Case 1) - Create one new experiment per selected experiment
      metadata_list <- list(
        'Processing'           = 'Generated by averaging many CD curves',
        'Original CD curves'   = paste0(spectra1[(1+n*(i-1)):(n*i)],collapse = ' '),
        'Units for processing' = processingUnitsString
      )
      
      append_cd_experiment(newSpectraNameBatch,newSpectraNameBatch,wlNew,
                           signalNew,operationUnits,metadata_list,newTemp,
                           htNew,isFake,operationUnits)
      
      newExperimentNames <- c(newExperimentNames,newSpectraNameBatch)
      legendDf           <- appendRowsToLegendDataFrame(legendDf,c(newSpectraNameBatch))
      
    }
    
  } else {
    
    signalNew  <- list()
    wlNew      <- list()
    htNew      <- list()

    if (operation == 'Smooth') {
      
      windowLength <- as.numeric(spectra2) 
      
      # smooth each of the spectra
      for (i in 1:length(exp1s)) {
        
        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]
        
        deltaWL      <- ( max(wlSel1) - min(wlSel1) ) / (length(wlSel1) - 1)
        odd_n_data_points_window_len <- ceiling(windowLength / deltaWL) %/% 2 * 2 + 1
        
        signalNew[[i]]  <- sgolayfilt(signal1, p = 2, n = odd_n_data_points_window_len)
        wlNew[[i]]      <- wlSel1
        htNew[[i]]      <- ht1
        
      }
      
    } else if (operation == 'Zero') {
      
      wlRange <- as.numeric(spectra2) - 1

      # Apply an offset to each of the spectra. 
      # Mathematically, we subtract from each spectra the mean signal 
      # within a selected wavelength range
      for (i in 1:length(exp1s)) {
        
        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]
        
        wlSel_temp      <- wlSel1 >= (input$wavelengthRange[2] - wlRange) & wlSel1 <= input$wavelengthRange[2]
        
        offsetValue <- mean(signal1[wlSel_temp])

        signalNew[[i]]  <- signal1 - offsetValue
        
        wlNew[[i]]      <- wlSel1
        htNew[[i]]      <- ht1
        
      }
      
    } else if (operation %in% c('Subtract','Sum')) {
      
      # Retrieve the second spectra
      wlSel2      <- cdAnalyzer$experimentsOri[[last(names(ids))]]$wavelength
      signalSel2  <- cdAnalyzer$experimentsOri[[last(names(ids))]]$signalDesiredUnit[,last(ids)]
    
      for (i in 1:length(exp1s)) {
        
        wlSel1     <- wlSel[[i]]
        signal1    <- signalSel[[i]]
        ht1        <- htSel[[i]]
        
        if (input$allowLinearInterpolation) {
         
          interpolatedSignals <- interpolateVectors(wlSel1,signal1,wlSel2,signalSel2)
          y1                  <- interpolatedSignals$y1
          y2                  <- interpolatedSignals$y2
          wlNew[[i]]          <- interpolatedSignals$x1
          
          if (all(is.na(ht1))) {
            
            htNew[[i]] <- rep(NA,length(interpolatedSignals$x1))
            
          } else {
            
            # We do not sum the HT signals!
            # We return the HT signal of the 'left' selected experiment
            # The HT signal data points are interpolated (if required)
            # to match the wavelength of the new generated experiment
            htNew[[i]]      <- approx(x = wlSel1, 
                                      y = ht1, 
                                      xout = interpolatedSignals$x1 , method = "linear")$y
            
          }
          
        } else {
          
          y1                  <- signal1
          y2                  <- signalSel2
          wlNew[[i]]          <- wlSel1
          htNew[[i]]          <- ht1
        }
        
        if (operation == 'Subtract') {
          signalNew[[i]] <- y1 - y2
        } else {
          signalNew[[i]] <- y1 + y2
        }
      
      }
    } 
    
    newExperimentNames <- c()
    
    for (exp in unique(exp1s)) {
      
      correspondingIDs <- which(exp1s == exp)
      
      # All spectra from the same experiment should share the same wavelength data points
      wlExperiment     <- wlNew[[correspondingIDs[1]]]
      
      c1 <- length(unique(exp1s))    == 1
      c2 <- length(correspondingIDs) == 1
      
      # Assign the new experiment name
      if (c1) {
        newExperimentName <- newSpectraName
      } else {
        newExperimentName <- paste0(exp,' ',newSpectraName)
        
        if (newExperimentName %in% names) {
          newExperimentName <- paste0(newExperimentName,' ',get_hour_minute_sec())
        }
      }
      
      # Select the signal and assign spectra names 
      if (c2) {
        signalExperiment  <- matrix(signalNew[[correspondingIDs]],ncol = 1)
        htExperiment      <- matrix(htNew[[correspondingIDs]],ncol = 1)
        newSpectraNames   <- c(newExperimentName)
      } else {
        signalExperiment  <- as.matrix(do.call(cbind,signalNew[correspondingIDs]))
        htExperiment      <- as.matrix(do.call(cbind,htNew[correspondingIDs]))
        selSpectraNames   <- cdAnalyzer$experimentsOri[[exp]]$spectraNames[id1s[correspondingIDs]]
        newSpectraNames   <- paste0(selSpectraNames,' ',newSpectraName)
      }
     
      # Assign the temperature from the spectrum/a #1
      
      expTemperature    <- unlist(temperatureSel[correspondingIDs])
      
      # 5th step (Case 2) - Create one new experiment per selected experiment
      metadata_list <- list(
        'Processing'           = paste0('Generated by applying the ',operation,' operation'),
        'Original CD curves'   = paste0(cdAnalyzer$experimentsOri[[exp]]$spectraNames[id1s[correspondingIDs]],
                                        collapse = ' '),
        'Units for processing' = processingUnitsString
      )
      
      # Assign the temperature from the spectrum/a #1

      newExperimentNames <- c(newExperimentNames,newExperimentName)
      
      legendDf <- appendRowsToLegendDataFrame(legendDf,newSpectraNames)
      
      # Add the parameters for the smoothing
      if (operation == 'Smooth') {
        
        metadata_list[['Savitzky-Golay smoothing filter window']] = paste0(windowLength,' nm') 
        metadata_list[['Savitzky-Golay smoothing filter order']]  = '2'
        
        deltaWL      <- ( max(wlExperiment) - min(wlExperiment) ) / (length(wlExperiment) - 1)

        metadata_list[['Savitzky-Golay smoothing filter length']] = paste0(ceiling(windowLength / deltaWL) %/% 2 * 2 + 1)
      }
      
      # Add the parameters for the sum/subtract
      if (operation %in% c('Sum','Subtract')) {
        
        metadata_list[['Second CD curve (for addition / subtraction)']] = spectra2

      }
      
      # Add the parameters for the smoothing
      if (operation == 'Zero') {
        
        wlRangeStr <- paste0(input$wavelengthRange[2] - wlRange,' - ',input$wavelengthRange[2], ' nm')
        metadata_list[['Wavelength range to calculate offset']] = wlRangeStr
        
      }
      
      append_cd_experiment(newExperimentName,newSpectraNames,np_array(wlExperiment),
                           signalExperiment,operationUnits,metadata_list,
                           expTemperature,htExperiment,isFake,operationUnits)
      
    }
    
  } 
  
  ## If shared parameters, assign known values of MW, #AA, ...
  if (cdAnalyzer$sharedParameters & !isFake) {
    
    # Use first experiment as proxy for the parameters
    mol_weight  <- cdAnalyzer$experimentsOri[[1]]$molecularWeight
    path_length <- cdAnalyzer$experimentsOri[[1]]$pathLength
    cro_num     <- cdAnalyzer$experimentsOri[[1]]$numberOfCroms
    conc        <- cdAnalyzer$experimentsOri[[1]]$concentration
    
    # Assign the same values to the generated experiments
    for (exp in newExperimentNames) {
      
      cdAnalyzer$experimentsOri[[exp]]$molecularWeight  <- mol_weight
      cdAnalyzer$experimentsOri[[exp]]$pathLength       <- path_length
      cdAnalyzer$experimentsOri[[exp]]$numberOfCroms    <- cro_num
      cdAnalyzer$experimentsOri[[exp]]$concentration    <- conc
      
    }
    
  } else {
    
    # Assign the same value as the experiments used to generate the new CD curves
    if (operation %in% c('Smooth','Subtract','Sum','Zero')) {
      
      i <- 0
      for (expNew in newExperimentNames) {
        i <- i + 1
        
        mol_weight  <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$molecularWeight
        path_length <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$pathLength
        cro_num     <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$numberOfCroms
        conc        <- cdAnalyzer$experimentsOri[[unique(exp1s)[i]]]$concentration
        
        cdAnalyzer$experimentsOri[[expNew]]$molecularWeight  <- mol_weight
        cdAnalyzer$experimentsOri[[expNew]]$pathLength       <- path_length
        cdAnalyzer$experimentsOri[[expNew]]$numberOfCroms    <- cro_num
        cdAnalyzer$experimentsOri[[expNew]]$concentration    <- conc
        
      }
      
    }
    
    # Modify the input parameters Table
    updateCDFilesInfoTable()
  }
  
  # 6th step - convert back to the appropriate units
  rList2        <- as.list(rep(input$workingUnits,length(unique(names(ids)))+length(newExperimentNames)))
  names(rList2) <- c(unique(names(ids)),newExperimentNames)
  
  cdAnalyzer$experiments_absorbance_units_to_other_units(rList2)
  
  # Filter the new experiments
  wlRange <- input$wavelengthRange
  cdAnalyzer$filter_data_by_wavelength(wlRange[1], wlRange[2])
  
  updateSelectInput(session,"experiment2delete",selected = "ALL",
                    choices = c('ALL',cdAnalyzer$experimentNames))

  updateSelectInput(session,"selected_cd_exp",choices     = cdAnalyzer$experimentNames)
  
  output$legendInfo           <- helperRenderRHandsontable(legendDf)
  
  updateProcessingTable(operation,operationUnits)
  
  # We need to use lapply instead of a for-loop to avoid R lazy evaluation!
  lapply(newExperimentNames, function (exp) {
    
    saneExpName <- gsub(':','',exp)
    
    tabP <- tabPanel(saneExpName,value=saneExpName,fluidRow(column(12,tableOutput(paste0('metadata_',saneExpName)))))
    
    appendTab("metadata",tabP,select=TRUE)
    
    metadata_info     <- cdAnalyzer$experimentsOri[[exp]]$metadata
    metadataFeature   <- names(metadata_info)
    metadataValue     <- unlist(metadata_info)
    
    metadata_df <- data.frame(metadataFeature,metadataValue)
    colnames(metadata_df) <- c('Metadata feature','Value') 
    
    output[[paste0('metadata_',saneExpName)]] <- renderTable({metadata_df})
    return(NULL)
    
  })
  
  Sys.sleep(2)
  reactives$data_loaded <- TRUE
  
})

# To update which the selectInput options
# For example, if we select 'Average' as the operation, there is no need to allow selecting
# a second spectrum
observeEvent(input$operation,{
  
  req(input$operation)
  updateProcessingTable(input$operation,input$operationUnits)
  
})

# Update the slider range based on the maximum accepted voltage value
# For this filter to work, the HT values of the experiments should have similar ranges
# E.g., it won't work if some CD experiments have HT values ranging from 400 to 600
# and other CD experiments have HT values from 2 to 6
observeEvent(input$maxHTvalue,{
  
  req(reactives$data_loaded)
  req(is.numeric(input$maxHTvalue))
  
  allVoltageData <- cdAnalyzer$get_experiment_properties('signalHT')
  allWLData      <- cdAnalyzer$get_experiment_properties('wavelength')
  
  # current values of the wavelength range slider
  minWL <- input$wavelengthRange[1]
  maxWL <- input$wavelengthRange[2]
  
  setNewWL <- FALSE
  
  minWL_news <- c()
  maxWL_news <- c()
  
  for (i in length(allVoltageData)) {
    voltage_temp <- allVoltageData[[i]]
    
    if (sum(!is.na(voltage_temp)) != 0 ) {
      wlTemp    <- allWLData[[i]]
      
      # Find the maximum wavelength that has associated HT values
      # under the HT threshold 
      filteredA <- apply(voltage_temp, 2, function(x) {
        sel <- which(x <= input$maxHTvalue)
        return(ifelse(length(sel)>0,max(sel),NA))
      })
      
      # Find the minimum wavelength that has associated HT values
      # under the HT threshold 
      filteredB <- apply(voltage_temp, 2, function(x) {
        sel <- which(x <= input$maxHTvalue)
        return(ifelse(length(sel)>0,min(sel),NA))
      })
      
      # Find the overall minimum and maximum wavelength for this CD experiment
      # Append these values to minWL_news and maxWL_news
      if (!all(is.na(filteredA)) & !all(is.na(filteredB))) {
        setNewWL    <- TRUE
        minWL_new_i <- min(c(wlTemp[filteredA],wlTemp[filteredB]),na.rm = T)
        maxWL_new_i <- max(c(wlTemp[filteredA],wlTemp[filteredB]),na.rm = T)
        
        minWL_news <- c(minWL_news,minWL_new_i)
        maxWL_news <- c(maxWL_news,maxWL_new_i)
        
      }
      
    }
    
  }
  
  
  if (setNewWL) {
    
    # Set the new wavelength limits

    minWL <- max(minWL_news)
    maxWL <- min(maxWL_news)
    
    updateSliderInput(session,'wavelengthRange',NULL,value = c(minWL,maxWL))
      
  }
  
})

