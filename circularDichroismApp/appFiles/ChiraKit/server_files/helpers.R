## Welcome message
welcomeMessage <- function() {
  shinyalert(paste("Welcome to", appName," <br><small>By clicking the 'I accept' button and using the eSPC Software, 
  you agree to be bound by the
  terms of this <a href='eSPC_academicSoftwareLicenseAgreement _EMBLEM.pdf' target='_blank' rel='noopener noreferrer'>Academic Software License Agreement</a>. You acknowledge that this Agreement
  is enforceable like any written agreement negotiated and signed by you. If you do not agree,
  please click the 'I decline' button and you will not have access to the Software. If you
  are not a member of a public funded academic and/or education and/or research institution,
  please contact <a href='https://embl-em.de/company/contact/' target='_blank' rel='noopener noreferrer'>EMBLEM</a>. </small>"),#imageUrl="embl_logo.svg",
             imageWidth = 180,imageHeight = 180,closeOnClickOutside=FALSE,closeOnEsc=FALSE,
             confirmButtonText="I accept",size = "m", 
             showCancelButton=TRUE,cancelButtonText="I decline",html=TRUE,
             confirmButtonCol="#8bb8e8",
             callbackR = function(x) {
               if (!x) welcomeMessage()
             })
}

# Return a string with the pattern %H:%M
get_hour_minute_sec <- function() {
  time_str <- as.character(Sys.time())
  time_str <- strsplit(time_str,' ')[[1]][2]
  return(time_str)
}

popUpWarning <- function(string) shinyalert(text = string,type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
popUpInfo    <- function(string) shinyalert(text = string,type = "info",   closeOnEsc = T,closeOnClickOutside = T,html=T)
popUpSuccess <- function(string) shinyalert(text = string,type = "success",closeOnEsc = T,closeOnClickOutside = T,html=T)

## Count folders in the current directory
count_folders <- function(dir) {
  
  total_folders <- length(list.files(dir))
  
  return(total_folders)
}

# Specify the numbers of decimal (k) for a certain number (x)
specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# Clean html sup tags and symbols that may not be read by computers...
# Replace the html like <sup> tag with the character '^' 
clean_html_sup_tag <- function(string) {

  string <- gsub('<sup>','^',        string)
  string <- gsub('</sup>','',        string)
  string <- gsub('\u2219','*',       string)
  string <- gsub('\u0394','Delta',   string)
  string <- gsub('\u03B5',' epsilon',string)
  
  string <- gsub('(\u0394A)','', string)
  string <- gsub('(\u0394mA)','',string)
  string <- gsub('\\s+',' ',     string)
  
  return(string)
}

# Fix indentation to n characters
# The same could be achieved with the function 'sprintf'
# but it does not work with same characters like '°'
identate <- function(string_vector,nChar) {
  
  for (i in 1:length(string_vector)) {
    nCharDiff <- nChar - nchar(string_vector[i])
    
    if (nCharDiff > 0 ) {
      toAppend          <- paste(rep(" ", nCharDiff), collapse = "")
      string_vector[i] <- paste0(string_vector[i],toAppend)
    }
  }
  
  return(string_vector)
}

## Splits a vector into a list of n elements 
split_vec <- function(vector,chunck_n) {
  sels    <- list()
  chuncks <- ceiling( length(vector) /chunck_n )
  
  for (i in 1:chuncks) {
    idx <- seq(1+chunck_n*(i-1),chunck_n*i)
    sels[[i]] <- na.omit(vector[idx])
  }
  
  return(sels)
}

getFileNameExtension <- function (fn) {
  # remove a path
  splitted    <- strsplit(x=fn, split='/')[[1]]   
  # or use .Platform$file.sep in stead of '/'
  fn          <- splitted [length(splitted)]
  ext         <- ''
  splitted    <- strsplit(x=fn, split='\\.')[[1]]
  l           <-length (splitted)
  if (l > 1 && sum(splitted[1:(l-1)] != ''))  ext <-splitted [l] 
  # the extention must be the suffix of a non-empty name    
  return(ext)
}


remove_file_extension <- function(name) {
  
  # Heuristic rule - if no '.' in the last 6 characters, return the same input
  
  last_six_characters <- substr(name, nchar(name) - 6, nchar(name))
  
  if ( !grepl('\\.',last_six_characters) ) return( name )
  
  fileExtension <- getFileNameExtension(name)
  
  # Remove extension, if it does not end with .dx where x is a number
  if (!grepl(".*\\.d[0-9]+$", name)) {
    name      <- sub(paste0(".",fileExtension,"$"), "", name)
  }
  return(name)
}

# Move one element to the beginning of the vector
# Based on the elements names
orderChoices <- function(selected,choices) {
  
  levels  <- c(names(choices[choices==selected]),names(choices[choices!=selected]))
  choices <- choices[order(factor(names(choices), levels = levels))]
  
  return(choices)
  
}

getChoices <- function(selected) {
  
  return(orderChoices(selected,global_cd_units_choice))
  
}

replace0valuesWithNAs <- function(vector) {
  idx <- vector == 0
  vector[idx] <- NA
  return(vector)
}

# Generate the DT table (user editable) with the experiments information about
# molecular, path length, number of chromophore units, etc. 
# This table is shown in the Tab 1. Load Input
# Requires:
#   - cdAnalyzer, a python object used to handle CD data 
generateDTtable <- function(cdAnalyzer) {
  
  exps              <- cdAnalyzer$experimentNames
  trueExperiments   <- !unlist(cdAnalyzer$get_experiment_properties('isFakeExperiment'))
  
  if (sum(trueExperiments) == 0) {
    df   <- data.frame(col1= character(),col2= numeric()  ,col3= numeric(),
                       col4= numeric()  ,col5= character(),col6= numeric()) 
  } else {
    
    # Generate table data
    molWeights       <- (unlist(cdAnalyzer$get_experiment_properties('molecularWeight')))
    pathLengths      <- (unlist(cdAnalyzer$get_experiment_properties('pathLength')))*10
    concs            <- (unlist(cdAnalyzer$get_experiment_properties('concentration')))
    numberOfCroms <- (unlist(cdAnalyzer$get_experiment_properties('numberOfCroms')))
    
    iUnits           <- unlist(cdAnalyzer$get_experiment_properties('units'))
    
    # Subset data
    exps                <- exps[trueExperiments]
    molWeights          <- molWeights[trueExperiments]
    pathLengths         <- pathLengths[trueExperiments]
    concs               <- concs[trueExperiments]
    numberOfCroms       <- numberOfCroms[trueExperiments]
    iUnits              <- iUnits[trueExperiments]
    currentN            <- length(exps)
    
    inputUnits  <- sapply(1:currentN, function(i) {
      as.character(selectInput(paste0('inputUnits',i), label=NULL, 
                               choices = getChoices(iUnits[i]), selectize=FALSE))
    })
    
    df   <- data.frame(exps,molWeights,numberOfCroms,concs,inputUnits,pathLengths) 
    
  }

  colnames(df) <- c('File name',"Mol. weight (Dalton)","#Chromophore units","Conc. (mg/ml)",
                    "Input units","Path length (mm)")
  
  return(df)
}

# Auxiliary function to render the DT Table on the server
renderDTtable <- function(df) {
  
  scrollY <- FALSE

  columns2disable <- which(colnames(df) %in% c("File name","Input units")) - 1
  
  DT::renderDataTable({
    df},editable = list(target = "cell", disable = list(columns = columns2disable)),
    escape=FALSE,rownames=FALSE,
    options = list(info = FALSE, dom="t",autoWidth = TRUE,ordering=FALSE,
                   scrollX = FALSE,scrollY = scrollY,pageLength = 1000,fillContainer = TRUE,
                   drawCallback    = JS('function() { Shiny.bindAll(this.api().table().node()); } ')
    ))
  
}

## Add replicate vector to dataframe according to variable and legend columns
## variable differs for each replicate
## The column legend will be replacted with df$legend rep X where X is the replicate number
add_rep_vec <- function(df) {
  
  df <- df %>% 
    group_by(variable) %>% 
    nest() 
  
  nRows <- nrow(df)
  if (nRows > 1) {
    df$rep <- paste0(" Rep ", 1:nRows)
  } else {
    df$rep <- ""
  }
  
  df        <- df %>% unnest(cols = c(data)) %>% ungroup()
  df$legend <- paste0(df$legend,df$rep)
  
  return(df)
}

counter2letter <- function(n) {
  vec <- c("A","B","C","D","E","F","G","H","I","J","K","L",
    "M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  
  temp <- expand.grid(vec,vec)
  
  return(c(vec,paste0(temp$Var1,temp$Var2))[n])
  
}

# Obtain a set of n colours (nColors) based on the Set 3 from the Brewer palette
getPalette <- function(nColors) {
  
  if (nColors <= 9) return(RColorBrewer::brewer.pal(9, "Set1")[1:nColors])
  if (nColors <= 12) return(RColorBrewer::brewer.pal(12, "Set3")[1:nColors])
  
  return( colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nColors) )
  
}

## Generate Internal IDs
getInternalIDs <- function(cdAnalyzer) {

  ID <- cdAnalyzer$get_experiment_properties('internalID')
  ID <- unlist(ID)
  
  return(ID)
  
}

## Generate a boolean vector to known which are 'fake' experiments
getFakeIDs <- function(cdAnalyzer) {
  
  fakeIDs <- cdAnalyzer$get_experiment_properties('isFakeExperiment')
  fakeIDs <- unlist(fakeIDs)
  
  return(fakeIDs)
  
}

## Generate spectra names
getSpectraNames <- function(cdAnalyzer) {
  
  spectraNames <- cdAnalyzer$get_experiment_properties('spectraNames')
  spectraNames <- unlist(spectraNames)
  
  return(spectraNames)
  
}

## Get plotting dataframe
## Takes as input the python object 
## Returns a dataframe  
getPlottingDF <- function(cdAnalyzer) {
  
  ID           <- getInternalIDs(cdAnalyzer)
  names        <- getSpectraNames(cdAnalyzer)
  total        <- length(ID)
  legend       <- gsub(pattern = "_",replacement = " ",names)
  colorPalette <- getPalette(total)
  
  df <- data.frame("Internal ID"=ID,
                   "Color"=colorPalette,
                   "Legend" = legend,
                   "Show"  = as.logical(rep(TRUE,total))
                   )
  return(df)
  
}

get_legend_from_rhandTable   <- function(table) return(as.character(hot_to_r(table)$Legend))
get_colors_from_rhandTable   <- function(table) return(hot_to_r(table)$Color)
get_sel_from_rhandTable      <- function(table) return(hot_to_r(table)$Show)
get_ID_from_rhandTable       <- function(table) return(hot_to_r(table)$Internal.ID)

# Convert the rhandsontable into a dataframe
# Requires
# - the rhandsontable table
# Output
# - a dataframe that can be used for further processing
getLegendDF <- function(legendInfo) {
  
  Legend  <- get_legend_from_rhandTable(legendInfo)
  Color   <- get_colors_from_rhandTable(legendInfo)
  Show    <- get_sel_from_rhandTable(legendInfo)

  Internal.ID   <- get_ID_from_rhandTable(legendInfo)

  legendDf          <- data.frame(Internal.ID,Color,Legend,Show)
  
  return(legendDf)
}

# Auxiliary function to render the legend dataframe
# We colour cells from the second column using the hex code inside its cells,
# we fix the first column
# we render the fourth column as a boolean 
# Requires
# - The dataframe
# Output
# - The rendered table to show in the server side
helperRenderRHandsontable <- function(legendDf) {
  
  color_cells <- data.frame(col=2,row=1:nrow(legendDf))
  
  rtable <- rhandsontable(legendDf,rowHeaders=NULL,
                          col_highlight = color_cells$col - 1,
                          row_highlight = color_cells$row - 1
  )
  
  rtable <- rtable %>% hot_col(col = c(1), readOnly=TRUE,
                               renderer = myrenderer) %>% 
    hot_col(col = c(2,3),renderer = myrenderer) %>%   
    hot_col(col = c(4),renderer = myrendererBoolean) %>% 
    hot_table(stretchH='all')

  return(renderRHandsontable({rtable}))
  
}

# Create a nice palette (viridis) based on numeric values
get_colors_from_numeric_values <- function(values,useLogScale=F) {
  
  viridis <- c("#440154","#471063","#481d6f","#472a7a",
               "#414487","#3c4f8a","#375a8c",
               "#32648e","#2a788e","#26828e",
               "#228b8d","#1f958b","#22a884",
               "#2cb17e","#3bbb75","#4ec36b",
               "#7ad151","#95d840","#b0dd2f","#cae11f",
               "#fde725")
  
  if (useLogScale) {
    maxL <- max(log10(values))
    minL <- min(log10(values))
  } else {
    maxL <- max(values)
    minL <- min(values)
  }
 
  seq <- seq(minL,maxL,length.out = 21)
  
  if (useLogScale) {
    idx <- sapply(values,function(v) which.min(abs(log10(v) - seq)))
  } else {
    idx <- sapply(values,function(v) which.min(abs(v - seq)))
  }
  
  return(viridis[idx])
}

# Generate the DT table for the processing step
# Requires:
# - cdAnalyzer, the python object to handle CD data
# - processingOption, a string to decide how and which checkbox will be created
# Output
# - the dataframe with the preprocessing options
generateDTtableProcessing <- function(cdAnalyzer,
                                      operation='Sum',operationUnits='millidegrees') {
  
  internal.IDs <- getInternalIDs(cdAnalyzer) # Internal ID

  if (length(internal.IDs) > 1) {
    possibleSpectra1choices <- c("Selected in 'Show' Column")
  } else  {
    possibleSpectra1choices <- c()
  }
  
  # Set Selected in 'Show' Column as the only option for the average or batch average operation
  if (!grepl('average',operation,ignore.case = TRUE)) {
    possibleSpectra1choices <- c(possibleSpectra1choices,internal.IDs)
  }
  
  possibleSpectra1  <- as.character(selectInput('inputSpectra1', label=NULL, 
                                                choices = possibleSpectra1choices,
                                                selectize=FALSE))
  
  lastNforBatchAverage <- min(10,length(internal.IDs)-1)
  
  if (operation == 'Batch average') {
    possibleSpectra2choices <- paste0('N = ',2:lastNforBatchAverage)
  } else if (operation == 'Zero') {
    possibleSpectra2choices <- c(3,5,10,20,40)
  } else if (operation == 'Smooth') {
    possibleSpectra2choices <- c(3,6,8,10)
  } else {
    possibleSpectra2choices <- c(internal.IDs)
  }
  
  if (!(operation %in% c('Sum','Subtract','Batch average','Zero','Smooth'))) {
    second_spectrum <- c('Not required for this operation')
  } else {
    second_spectrum <- possibleSpectra2choices
  }
  
  possibleSpectra2  <- as.character(selectInput('inputSpectra2', label=NULL, 
                                                choices = second_spectrum,
                                                selectize=FALSE))
  
  operationsChoices <- c('Subtract','Sum','Smooth','Average','Zero')
  
  # Add batch average if we have more than 1 curve
  if (length(internal.IDs) > 2) {
    operationsChoices <- c(operationsChoices,'Batch average')
  } 
  
  operationsChoices <- c(operation,operationsChoices[operationsChoices!=operation])
  
  operations        <- as.character(
    selectInput('operation', label=NULL, choices = operationsChoices,
                selectize=FALSE))
  
  operationUnitsChoices <- getChoices(operationUnits)
  
  operationUnits    <- as.character(
    selectInput('operationUnits',label=NULL,
                choices = operationUnitsChoices,
                selectize=FALSE))
  
  newSpectraName <- as.character(textInput('newSpectraName',label = NULL,'test spectra'))
  
  df   <- data.frame(operations,possibleSpectra1,operationUnits,possibleSpectra2,newSpectraName) 
  
  if (operation == 'Batch average') {
    third_col_name <- 'N (for batch average)'
  } else if (operation == 'Zero') {
    third_col_name <- 'WL range (for zeroing)'
  } else if (operation %in% c('Sum','Subtract')) {
    third_col_name <- 'Spectrum 2 (for + or -)'
  } else if (operation  == 'Smooth') {
    third_col_name <- 'Smoothing window size (nm)'
  } else {
    third_col_name <- 'Not required for this operation'
  }
  
  colnames(df) <- c('Operation','Spectrum/a 1','Operation Units',
                    third_col_name,'Output name')
  
  return(df)
}

# Given a list of vectors and a list of integers (n1, n2, ...)
# Let's call each vector an item list, and the elements of the vectors items
# We  
#   1) Flatten the list of vectors
#   2) Find the n-th item (of the flattened list)
#   3) Return     a) the name of the item list containing the item
#             and b) the subindex of the item inside the corresponding item list

# For example, the input 
#   'listOfVectors'       : list(c('a','b'),c('c','d')) ; names('listOfVectors') gives c('hola','chau')
#   'vectorOfIntegers'    : c(3,2)

# Returns: the named vector c(1,2) with names c('chau','hola')

found_ids <- function(listOfVectors,vectorOfIntegers) {
  
  c <- 0 # counter of ALL items
  t <- 0 # counter of list of items
  
  selectedItemIndex      <- c()
  itemListName           <- c()  
  vectorOfIntegersFound  <- c()
  vectorOfIntegersOri    <- vectorOfIntegers
  
  for (name in names(listOfVectors)) {
    
    t <- t + 1
    totalItemsOfItemlist <- length(listOfVectors[[t]])
    
    for (integerID in vectorOfIntegers) {

      if (integerID > c &  integerID <= (c + totalItemsOfItemlist)) {
        
        selectedItemIndex     <- c(selectedItemIndex,integerID - c) 
        itemListName          <- c(itemListName,name)
        vectorOfIntegers      <- vectorOfIntegers[vectorOfIntegers != integerID]
        vectorOfIntegersFound <- c(vectorOfIntegersFound,integerID)
      }
      
    }

    c <- c + totalItemsOfItemlist
    
  }

  matching_indexes  <- sapply(vectorOfIntegersOri, function (x) which(vectorOfIntegersFound == x)[1])

  selectedItemIndex   <- selectedItemIndex[matching_indexes]
  itemListName        <- itemListName[matching_indexes]
  
  names(selectedItemIndex) <- itemListName
  return(selectedItemIndex) 
}

# Combine the two x-axis of two vectors and interpolate the y-values
interpolateVectors <- function(x1,y1,x2,y2) {
  
  # Restrict the vectors to the range of analysis
  combinedX <- unique(c(x1,x2))
  
  # Y1 and Y2 linearly approximated on the combined X-axis of Y1 and Y2
  y1_approx <- approx(x = x1, y = y1, xout = combinedX, method = "linear")$y
  y2_approx <- approx(x = x2, y = y2, xout = combinedX, method = "linear")$y
  
  toDeleteId     <- is.na(y1_approx) | is.na(y2_approx)
  
  y1_approx <- y1_approx[!toDeleteId]
  y2_approx <- y2_approx[!toDeleteId]
  combinedX <- combinedX[!toDeleteId]
  
  return(list('x1'=combinedX,'y1'=y1_approx,'y2'=y2_approx))
}

averageListOfVectors <- function(listXaxis,listYaxis) {
  
  # Restrict the vectors to the range of analysis
  combinedX <- unique(unlist(listXaxis))
  
  # Y linearly approximated on the combined X-axis 
  y_approx <- lapply(1:length(listYaxis), function (i) {
    
    approx(x = listXaxis[[i]], y = listYaxis[[i]], xout = combinedX)$y
    
  })
  
  total      <- Reduce('+',y_approx)
  toDeleteId <- is.na(total) # Find values outside the overlapping x-axis range
  
  y_approx  <- total[!toDeleteId] / length(listXaxis)
  combinedX <- combinedX[!toDeleteId]
  
  return(list('x'=combinedX,'average'=y_approx))
}

get_ht_average <- function(wlSel,htSel) {
  
  if (any(is.na(unlist(htSel)))) {
    htNew <- NULL
  } else {
    interpolatedAverageHT <- averageListOfVectors(wlSel,htSel)
    htNew                 <- matrix(interpolatedAverageHT$average,ncol = 1)
  }
  
  return(htNew)
}

appendRowsToLegendDataFrame <- function(legendDf,newSpectraNames) {
  
  # Create fake palette
  fakePalette <- rep('#000000',length(newSpectraNames))

  # Append information to the legends data frame
  newRow   <-   data.frame("Internal ID" = newSpectraNames,
                           "Color"       = fakePalette,
                           "Legend"      = newSpectraNames,
                           "Show"        = as.logical(TRUE))
  
  legendDf                    <- rbind(legendDf,newRow)
  
  # Replace the palette (of the whole dataframe!)
  legendDf$Color <- getPalette(nrow(legendDf))
  
  return(legendDf)
}

# Full join list of dataframes
full_join_df_lst <- function(dfs,by_column = 'wavelength') {
  
  number_of_dfs <- length(dfs)
  
  # Case 1 - only one dataframe
  if (number_of_dfs == 1) {
    return(dfs[[1]])
  } else {
    df <- full_join(dfs[[1]],dfs[[2]],by = by_column)
    
    if (number_of_dfs > 2) {
      for (i in 3:number_of_dfs) {
        df <- full_join(df,dfs[[i]],by = by_column)
      }
    }
    
  }
  return(df)
}

## Find experiments where the 'signalDesiredUnit' matrix units do not match the selected working units
## In other words, experiments of the type 'fakeExperiment' which were created by preprocessing spectra
## using (mean unit) molar extinction / ellipticity as units

## Input
##  - cdAnalyzer, python object to handle CD data
## Output
##  - the ids of the non matching experiments
find_non_matching_units_experiments <- function(cdAnalyzer,workingUnits) {

  signalsAll             <- cdAnalyzer$get_experiment_properties_modif('signalDesiredUnit')
  isFake                 <- unlist(cdAnalyzer$get_experiment_properties_modif('isFakeExperiment'))
  fakeExperimentSignal   <- unlist(cdAnalyzer$get_experiment_properties_modif('fakeExperimentSignal'))

  ids <- c()
  i   <- 0
  for (expName in cdAnalyzer$experimentNames) {

    i   <- i + 1
    signals     <- signalsAll[[i]]

    c1 <- all(is.na(signals))
    c2 <- isFake[i] & (fakeExperimentSignal[i] != workingUnits)

    ids <- c(ids,c1 | c2)

  }

  return(ids)
}

## Convert dataframe to lines
df_to_lines <- function(df) {
  
  numeric_columns <- sapply(df, is.numeric)
  
  original_colnames <- colnames(df)
  
  if (sum(numeric_columns) > 0 ){
    
    df[, numeric_columns] <- signif(df[, numeric_columns], 5)
    
    # Sort from high to low using the first column
    df   <- as.data.frame(df)
    df <- df[rev(order(df[,which(numeric_columns)[1]])),]
    
  }
  
  df <- data.frame(lapply(df, as.character), stringsAsFactors=FALSE)
  
  # Assign original column names
  colnames(df) <- original_colnames
  
  # Create a character vector with custom text
  lines <- c(paste0(colnames(df),collapse = ','))
  
  for (i in 1:nrow(df)) {
    lines <- c(lines,paste0(df[i,],collapse = ','))
  }
  
  # Specify the file and use the cat() function to create the text file
  return(lines)
  
}

## Add custom comments to the exported file
df_usual_comments <- function(workingUnits,
                              concentration_vector,pathLength_vector,
                              molWeight_vector,nPepBond_vector) {
  
  # Convert to a nice format
  workingUnits <- workingUnits2ProperLabel(workingUnits)
  
  # Replace the html like <sup> tag with the character '^' 
  workingUnits <- clean_html_sup_tag(workingUnits)
  
  commentsLeft  <- identate(c('#File origin :',
                              '#Time stamp :',
                              '#CD units :',
                              '#Sample concentration (mg/ml) :',
                              '#Path length (cm) :',
                              '#Molecular weight (Da) :',
                              '#Number of chromophore units :',
                              '#Wavelength units :'),48)
  
  if (length(unique(concentration_vector)) == 1) concentration_vector  <- concentration_vector[1]
  if (length(unique(pathLength_vector))    == 1) pathLength_vector     <- pathLength_vector[1]
  if (length(unique(molWeight_vector))     == 1) molWeight_vector      <- molWeight_vector[1]
  if (length(unique(nPepBond_vector))     == 1)  nPepBond_vector       <- nPepBond_vector[1]
  
  commentsRight <- c('Exported with ChiraKit',
                     format(Sys.time(),usetz = TRUE),
                     workingUnits,
                     paste0(concentration_vector,collapse = ' '),
                     paste0(pathLength_vector,collapse = ' '),
                     paste0(molWeight_vector,collapse = ' '),
                     paste0(nPepBond_vector,collapse = ' '),
                     'nanometers (nm)')

  
  
  comments <- paste0(commentsLeft,commentsRight)
  
  return(comments)
}

# Function to delete elements equal to their previous element
delete_duplicates <- function(vector) {
  keep <- c(TRUE, vector[-length(vector)] != vector[-1])
  vector[keep]
}

# For each ith line, delete them if the ith+2 line contains the same pattern
delete_duplicate_pattern <- function(lines,pattern) {
  
  id2keep <- c()
  
  for (i in 1:(length(lines)-2)) {
    
    c <- grepl(pattern,lines[c(i,i+2)])
    
    id2keep <- c(id2keep,!all(c))
    
  }
  
  lines <- lines[id2keep]
  
  return(lines)
}

# Delete unwanted lines
purge_logbook_lines <- function(lines) {
  
  # Remove empty line
  toRemove <- lines == 'Setting input units of '

  lines <- lines[!toRemove]

  # Remove non-informative entries of units conversions
  # For example, lines with the text 'Converting CD units to: millidegrees'
  # if the previous units conversion was also the same
  unitsConversions <- grepl("Converting CD units to:",lines)
  
  if (any(unitsConversions)) {
    
    ids           <- (1:length(lines))[unitsConversions]
    selectedUnits <- sapply(lines[unitsConversions], function (x) {
      strsplit(x,'Converting CD units to: ')[[1]][2]
    })

    # Find indices where the element is equal to the previous one
    indices <- which(selectedUnits[-1] == selectedUnits[-length(selectedUnits)])
    
    ids2delete <- c()
    
    if (length(indices) > 0) {
      
      # Adjust indices by adding 1 to account for the offset introduced by [-1]
      indices    <- indices + 1
      
      # Select those indices for deletion
      ids2delete <- c(ids2delete,ids[indices])
      
      lines <- lines[-ids2delete]
      
    }
    
  }
  
  # Remove non-informative entries 
  # For example, if we have consecutive lines with the text
  #  'Setting wavelength range (nm) to ...
  # we only leave the last line
  lines <- delete_duplicate_pattern(lines,'Setting wavelength range')
  lines <- delete_duplicate_pattern(lines,'Setting the signal to analyse according to the wavelengths:')
  
  # remove duplicates entries, useful to remove extra white lines
  lines <- delete_duplicates(lines)
  
  return(lines)
}

# Create a dataframe useful for exporting the data
# Requires the experimental signal in delta epsilon units and the fitted signal

# Input: 
#   - list_of_signals,        list with sublists, each sublist contains vectors
#   - list_of_fitted_signals, list with sublists, each sublist contains vectors
# Output:
#   - nicely formatted dataframe
lst_of_cd_spectra_to_df <- function(list_of_signals,list_of_fitted_signals) {
  
  dfs_list <- list()
  
  for (i in 1:length(list_of_signals)) {
    
    exp_signals  <- list_of_signals[[c(i)]]
    exp_fittings <- list_of_fitted_signals[[c(i)]]
    
    for (ii in 1:length(exp_signals)) {
      
      exp_signal  <- exp_signals[[c(ii)]]
      exp_fitting <- exp_fittings[[c(ii)]]
      name        <- names(exp_signals)[c(ii)]
      
      df  <- data.frame('wavelength_nm'      = seq(240,240-length(exp_fitting)+1),
                        'query_deltaEpsilon' = exp_signal,
                        'fitting'            = exp_fitting,
                        'name'               = name) 
      
      dfs_list <- append(dfs_list,list(df))
      
    }
  }
  
  all_dfs <- do.call(rbind,dfs_list)
  
  return(all_dfs)
}


