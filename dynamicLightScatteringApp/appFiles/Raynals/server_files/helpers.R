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

## Count folders in the current directory
count_folders <- function(dir) length(list.files(dir))

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

# Remove some specific strings 
cleanSampleName <- function(sampleName) {
  
  sampleName <- gsub('Replicate','Rep',sampleName,ignore.case = TRUE)
  sampleName <- gsub('Acquisition','Acq',sampleName,ignore.case = TRUE)
  
  patterns <- c('.xlsx','autocorrelation','acf','capillary','-','_','\\./',
                '/','export','panta','\\s+')
  
  for (p in patterns) {
    sampleName <- gsub(p,' ',sampleName,ignore.case = TRUE)
  }
  
  sampleName <- trimws(sampleName)
  
  return(sampleName)
}

# Detect if the excel files belong to a folder with the name Acquisition X, where X is any number
# In that case, return 'Acquisition', otherwise return 'Independent'
detect_Excel_files_type <- function(xlsx_files) {
  
  c1 <- sum(grepl('cquisit',xlsx_files,ignore.case = TRUE)) > (length(xlsx_files) / 2)
  
  if (c1) return('Acquisition')
  
  return('Independent')
}

# Try to find the relevant sheet and read the first two columns (time / autocorrelation function)
read_excel_file <- function(xlsx_file,sampleName=NULL) {
  
  # replace sample name with the file name and path
  if (is.null(sampleName)) sampleName <- xlsx_file
  
  result <- tryCatch({
    
    sheet_names <- openxlsx::getSheetNames(xlsx_file)
    
    for (sheet_name in sheet_names) {
      
      data        <- openxlsx::read.xlsx(xlsx_file,sheet = sheet_name)
      
      if (is.null(data)) next
      
      if (grepl('tau',colnames(data)[1],ignore.case = TRUE) | grepl('time',colnames(data)[1],ignore.case = TRUE)) {
        
        data           <- data.frame(data[,c(1:2)])
        
        colnames(data) <- c('time','ac')
        # Remove duplicates ... 
        
        data           <- data %>% group_by(time) %>% 
          summarise(ac = mean(ac))
        
        # From seconds to microseconds
        if (max(data$time) < 1e4) data$time <- data$time*1e6
        
        data$time <- signif(data$time,3)

        colnames(data) <- c('time',cleanSampleName(sampleName))
        
        return(data)
      }
      return(NULL)
    }
  }, error = function(err) {
    return(NULL)
  })
  return(result)
}

# Given a list of dataframes, find
# the maximum number of rows, and remove
# the dataframes with less rows than half of that number
filter_corrupt_datasets <- function(allData) {
  
  # Remove corrupt data
  nrows     <- sapply(allData, nrow)

  idx2keep  <- sapply(nrows, function(x) x > max(nrows) / 2)
  toDelete  <- length(allData) - sum(idx2keep)
  
  if (toDelete > 0) {
    
    shinyalert(text = paste("<b>Caution: ",toDelete," datasets were removed from the analysis
                            because of insufficient data points.</b>"),
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    
  }
  
  allData   <- allData[idx2keep]
  
  return(allData)
}

load_Excel_files_without_acquisitions <- function(xlsx_files) {
  
  allData      <- list()
  
  i <- 0
  for (file in xlsx_files) {
    
    # Load the first sheet into a dataframe
    data <- read_excel_file(file)
    
    # Delete the file
    if (!is.null(data)) {
      i <- i + 1
      allData[[i]] <- data
      
    }
    
  }
  
  allData <- filter_corrupt_datasets(allData)
  
  return( Reduce(merge, allData) )
}

load_Excel_files_with_acquisitions <- function(xlsx_files) {
  
  allData <- list()
  
  i <- 0
  for (f in xlsx_files) {
    
    splitted     <- strsplit(f, '/', fixed = FALSE)[[1]]
    len_splitted <- length(splitted)
    
    if (grepl('cquisit',splitted[len_splitted - 1],ignore.case = TRUE)) {
      sampleName  <- paste0(splitted[(len_splitted - 3):(len_splitted - 2)],collapse = ' ')
      
      data <- read_excel_file(f,sampleName)
      
      if (!is.null(data)) {
        
        i <- i + 1
        allData[[i]]   <- data

      }
      
    }
    
  }
  
  allData <- filter_corrupt_datasets(allData)
  
  df_names   <- sapply(allData, function(x) colnames(x)[2])
  index_list <- split(seq_along(df_names), df_names)
  
  allDataAvg <- list()
  
  i <- 0
  
  for (idxs in index_list){
    
    i  <- i + 1
    df <- do.call(rbind,allData[idxs])
    
    colnamesOri    <- colnames(df)
    colnames(df)   <- c('time','ac')
    
    df <- df %>% group_by(time) %>% 
      summarise(ac = mean(ac))
    
    colnames(df) <- colnamesOri
    
    allDataAvg[[i]] <- df
    
  }
  
  return( Reduce(merge, allDataAvg) )
}

specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))

# Return the number of tables we want to show according to the number of
# measurements that were done
get_numberOfDesiredTables <- function(numberOfConditionds) {
  if (numberOfConditionds < 16 ) return(1) 
  if (numberOfConditionds < 31 ) return(2) 
  if (numberOfConditionds < 46 ) return(3) 
  
  return(6)
}

# dlsAnalyzer is a python object!
generateDTtable <- function(dlsAnalyzer) {
  
  exps      <- c(dlsAnalyzer$experimentNames)
  currentN  <- length(exps)
  
  if (currentN == 0) return(NULL)
  # Generate table data
  lambda0          <- unlist(dlsAnalyzer$getExperimentProperties("lambda0"))
  scatteringAngle  <- unlist(dlsAnalyzer$getExperimentProperties("scatteringAngle")) * 180 / pi # to degree
  scans            <- unlist(dlsAnalyzer$getExperimentProperties("scans"))
  reads            <- unlist(dlsAnalyzer$getExperimentProperties("reads"))
  
  viscosity          <- unlist(dlsAnalyzer$getExperimentProperties("viscosity"))
  refractiveIndex    <- unlist(dlsAnalyzer$getExperimentProperties("refractiveIndex"))
  temperature        <- unlist(dlsAnalyzer$getExperimentProperties("temperature")) - 273
  
  scatteringAngle <- round(scatteringAngle,0)
  df   <- data.frame(exps,scans,reads,lambda0,scatteringAngle,
                     temperature,refractiveIndex,viscosity) 

  colnames(df) <- c('File name',"#Scans","#Reads","Wavelength (nm)","Scattering angle (°)",
                    "Temperature (°C)","Refractive Index","Viscosity (pascal-second)")
  
  return(df)
}

counter2letter <- function(n) {
  vec1 <- c("A","B","C","D","E","F","G","H","I","J","K","L",
            "M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z")
  
  vec2 <- 1:24
  temp <- expand.grid(vec2,vec1)
  
  return(paste0(temp$Var2,temp$Var1)[n])
  
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

## Generate the tabPanel including the rHandson tables
## This tables contain the conditions data, including read and scan number
generateTabPanel<- function(numberOfDesiredTables,name) {
  
  # I don't know why a for loop didn't work to create the tables :)
  
  dlsPlateInfoNames <- getDlsPlateInfoNames(name,numberOfDesiredTables)
  
  if (numberOfDesiredTables == 1) {
    tabP <- tabPanel(name,value=name,fluidRow(
      column(3,rHandsontableOutput(dlsPlateInfoNames[1],height="400px"))))
  }
  
  if (numberOfDesiredTables == 2) {
    tabP <- tabPanel(name,value=name,fluidRow(
      column(3,rHandsontableOutput(dlsPlateInfoNames[1],height="400px")),
      column(3,rHandsontableOutput(dlsPlateInfoNames[2],height="400px"))))
  }
  
  if (numberOfDesiredTables == 3) {
    tabP <- tabPanel(name,value=name,fluidRow(
      column(3,rHandsontableOutput(dlsPlateInfoNames[1],height="400px")),
      column(3,rHandsontableOutput(dlsPlateInfoNames[2],height="400px")),
      column(3,rHandsontableOutput(dlsPlateInfoNames[3],height="400px"))))
  }
  
  if (numberOfDesiredTables == 6) {

    tabP <- tabPanel(name,value=name,fluidRow(
      column(2,rHandsontableOutput(dlsPlateInfoNames[1],height="400px")), 
      column(2,rHandsontableOutput(dlsPlateInfoNames[2],height="400px")), 
      column(2,rHandsontableOutput(dlsPlateInfoNames[3],height="400px")), 
      column(2,rHandsontableOutput(dlsPlateInfoNames[4],height="400px")), 
      column(2,rHandsontableOutput(dlsPlateInfoNames[5],height="400px")), 
      column(2,rHandsontableOutput(dlsPlateInfoNames[6],height="400px")))) 
    
  }
  
  return(tabP)
  
}

# Generate measurement Data frame based on number of scans and reads
getDataFrameBasedOnReadsAndScans <- function(conditions,reads,scans) {
  
  if (reads < 1) reads <- 1
  if (scans < 1) scans <- 1
  
  nConditions <- length(conditions)
  nWells   <- nConditions / scans / reads
  readCol  <- rep(1:reads,nConditions/reads)
  scansCol <- c(sapply(1:scans, function(x) rep(x,nWells*reads)))
  condNew  <- counter2letter(1:nWells)

  condNew2 <- c()
  
  for (i in 1:scans) {
    for (cc in condNew) {
      condNew2 <- c(condNew2,rep(cc,reads))
    }
  }
  
  maxCol <- min(unlist(lapply(list(condNew2,readCol,scansCol), length)))
  missingData <- rep('Unkown',nConditions-maxCol)
  
  df <- data.frame("conditions"=conditions,
                   "include"=c(rep(T,maxCol),rep(F,length(missingData))),
                   "read"=c(readCol[1:maxCol],rep(99,length(missingData))),
                   "scan"=c(scansCol[1:maxCol],rep(99,length(missingData)))
                   )
  
  return(df)
}

# Split the conditions table
split_table_into_list_of_tables  <- function(df,numberOfDesiredTables) {
  
  numberOfconditions    <- nrow(df)
  elementsPerTable      <- ceiling(numberOfconditions / numberOfDesiredTables)
  conditionsSplitted    <- split_vec(df$conditions,elementsPerTable)
  readsSplitted         <- split_vec(df$read,elementsPerTable)
  scansSplitted         <- split_vec(df$scan,elementsPerTable)
  includeSplitted       <- split_vec(df$include,elementsPerTable)
  
  tables <- list()
  for (i in 1:numberOfDesiredTables) {
    
    data <- data.frame(Condition = as.character(conditionsSplitted[[i]]),
                       Include   = includeSplitted[[i]],
                       Read      = as.integer(readsSplitted[[i]]),
                       Scan      = as.integer(scansSplitted[[i]])
    )
                       
    
    tables[[i]] <- data
    
  }
  
  return(tables)
}

## Get the Tables that contain the measurements Conditions, Reads, Scans
get_Table_list <- function(conditions,reads,scans,numberOfDesiredTables) {
  
  df     <- getDataFrameBasedOnReadsAndScans(conditions,reads,scans)
  tables <- split_table_into_list_of_tables(df,numberOfDesiredTables)
  
  return(tables)
}

renderDTtable <- function(df,columns2disable=c("File name"),
                          autoWidth=T,scrollYvalue='400px') {
  
  if (is.null(df)) return(NULL)
  
  scrollY = T
  if (nrow(df) > 10) {
    scrollY = scrollYvalue
  }
  
  columns2disable <- which(colnames(df) %in% columns2disable) - 1

  DT::renderDataTable({
    df},editable = list(target = "cell", disable = list(columns = columns2disable)),
    escape=FALSE,rownames=FALSE,
    options = list(info = FALSE, dom="t",autoWidth = autoWidth,ordering=FALSE,
                   fillContainer = TRUE,
                   scrollX = FALSE,scrollY = scrollY,
                   drawCallback    = JS('function() { Shiny.bindAll(this.api().table().node()); } ')
    ))
}

getDlsPlateInfoNames <- function(name,numberOfTables=6) {
  
  return(paste0("dlsPlateDatasetInfo_",name,1:numberOfTables))

}

getTablesInfo <- function(expData) {
  
  dfs         <- lapply(expData, hot_to_r)
  df          <- do.call(rbind,dfs)
  
  return(df)
}

getConditionNamesFromSamplesMetaData <- function(expName,df) {
  
  if (!('conditions' %in% colnames(df))) return(NULL)
  
  diffReads <- length(unique(df$read))
  difScans <- length(unique(df$scan))
  
  NAsIdx <- df$read == 'NA' | df$scan == 'NA' | df$conditions == 'NA' | (!df$include)
  
  if (diffReads > 1){
    df$read <- paste0('Read ',df$read,' ')
  } else {
    df$read <- ''
  }
  
  if (difScans > 1){
    df$scan <- paste0('Scan ',df$scan,' ')
  } else {
    df$scan <- ''
  }
  
  sampleNames <- paste0(df$conditions," ",df$read,df$scan,expName)
  
  for (s in unique(sampleNames)) {
    idxs <- which(sampleNames == s)
    if (length(idxs) > 1) {
      c <- 1
      for (ii in idxs) {
        sampleNames[ii] <- paste0(sampleNames[ii],' Rep ',c)
        c               <- c + 1
      }
    }
  }
  
  sampleNames[NAsIdx] <- 'NA'
  
  return(sampleNames)
}

filterDataBasedOnLogTime <- function(data,length.out = 80) {
  
  # Data should have one column called 'time'
  # Returns the rows of data where the time is evenly spaced in a log scale
  
  base2filter <- 1.5
  minTime <- log(min(data$time),base = base2filter)
  maxTime <- log(max(data$time),base = base2filter)
  
  time    <- base2filter^(seq(minTime,maxTime,length.out = length.out))
  selIdx  <- sapply(time, function(t) which.min(abs(data$time - t)))

  data    <- data[unique(selIdx), ]
  
  return(data)
}

reshapeAutocorrelationDFAndFilterByTime <- function(data,sampleNames) {
  
  # Check that more at least one condition was selected
  if (sum(sampleNames=="NA") < (ncol(data)-1)) {
    data <- data[,colnames(data) != 'NA']
    lengthOut <- 80
    if (ncol(data) > 200) lengthOut <- 40  
    if (ncol(data) > 400) lengthOut <- 30
    
    data <- filterDataBasedOnLogTime(data,lengthOut)
    data <- reshape2::melt(data,id.vars="time")
  } else {
    data <- data.frame(time=double(),variable=character(),value=double())
  }
  return(data)
}

formatDlsInfoForPlotting <- function(autocorrelationList,timeList,
                                     expNames,expTabData,threshold) {
  
  allDFs <- list()
  
  for (i in 1:length(expNames)) {
    
    expName     <- ifelse(length(expNames) == 1,'',expNames[i])
    
    sampleNames <- getConditionNamesFromSamplesMetaData(expName,expTabData[[i]])
    
    autocorrelation     <- autocorrelationList[[i]]
    time                <- timeList[[i]]
    
    data <- data.frame(time,autocorrelation)
    colnames(data) <- c('time',sampleNames)
    data           <- reshapeAutocorrelationDFAndFilterByTime(data,sampleNames)
    allDFs[[i]] <- data
  }
  
  df <- do.call(rbind,allDFs)

  return(df)
  
}

formatDlsPredictedInfoForPlotting <- function(autocorrelationListPredicted,timeList,
                                     expNames,expTabData) {
  
  allDFs <- list()
  
  for (i in 1:length(expNames)) {
    
    expName     <- ifelse(length(expNames) == 1,'',expNames[i])
    
    sampleNames <- getConditionNamesFromSamplesMetaData(expName,expTabData[[i]])
    
    autocorrelation     <- autocorrelationListPredicted[[i]]
    time                <- timeList[[i]]
    
    data <- data.frame(time,autocorrelation)
    colnames(data) <- c('time',sampleNames)
    data           <- reshapeAutocorrelationDFAndFilterByTime(data,sampleNames)
    
    allDFs[[i]] <- data
  }
  
  df <- do.call(rbind,allDFs)
  
  return(df)
  
}

formatContributions <- function(estimatedContributions,contributionsAxis,
                                expNames,expTabData,contributionsAxisType='hr') {
  
  # Obtain a nice formatted dataframe from the estimated relative contributions
  # Input:
  #  'estimatedContributions' - list (one element per experiment) 
  #     of lists (one element per curve: vectors of length n)
  
  #  'contributionsAxis' vector of length n (same n as the one from the 
  #     sublist elements of 'estimatedContributions')
  
  # 'expNames'   - names of the experiments, str vector
  
  # 'expTabData' - list of dataframes, each df  has the following columns:
  #     'conditions' 'read' 'scan' 'include'
  
  # 'contributionsAxis' - should be either 'hr' or 'diff' 
  #     Use 'hr' to return the results as a function of the hydrodynamic radius
  #     and 'diff' to return them as a function of diffusion coefficients
  
  allDFs <- list()
  
  for (i in 1:length(expNames)) {
    
    expName     <- ifelse(length(expNames) == 1,'',expNames[i])
    
    sampleNames <- getConditionNamesFromSamplesMetaData(expName,expTabData[[i]])
    
    estimatedContributionsTemp     <- estimatedContributions[[i]]
    estimatedContributionsTemp     <- do.call(cbind,estimatedContributionsTemp)

    data           <- data.frame(contributionsAxis,estimatedContributionsTemp)
    
    colnames(data) <- c(contributionsAxisType,sampleNames)
    
    data <- data[,colnames(data) != 'NA']
    data <- reshape2::melt(data,id.vars=contributionsAxisType)
    
    allDFs[[i]] <- data
  }
  
  df <- do.call(rbind,allDFs)
  
  return(df)
  
}

## Get choices to plot residuals
## We can't display more than 25 conditions at a time...
get_choices_residuals_plot <- function(nExperiments) {
    nChoices <- ceiling(nExperiments/25)
    choices  <- c('Conditions 1-25')
    if (nChoices==1) return(choices)
    for (i in 1:(nChoices-1)) {
      choices <- c(choices,paste0('Conditions ',i*25-1,'-',(i+1)*25))
    }
    return(choices)
}

## Remove small peaks
## Use intensity threshold to remove 'fake' peaks
removeSmallPeaks <- function(x,threshold=0.025) {
  # x is a data frame of two columns
  # 'hr' (hydrodynamic radius) 
  # 'value' (contribution of the corresponding Hrs)
  
  peaks <- findpeaks(x$value,nups = 2,
                     minpeakheight=threshold)
  hr <- x$hr
  hrs2keep <- c()
  
  if (is.null(peaks)) return(x)
  
  for (i in 1:nrow(peaks)) {
    hrs <- hr[hr > (hr[peaks[i,3]]) & hr < (hr[peaks[i,4]])]
    hrs2keep <- c(hrs2keep,hrs)
  }
  
  x$value[!(hr %in% hrs2keep)] <- 0
  
  return(x)
}

### Return the substring (from a vector of substrings) that is present in a longer string 
### Used to get only the experiment name 
### Used to plot the autocorrelation curves using as factor the experiment names
### Input example - substrings = c("CA 5 replicates","BSA 3 replicates")
###                 longString = 'A1 Scan 5 CA 5 replicates'
### Returns - CA 5 replicates
get_exp_name <- function(substrings,longString) {
  
  for (s in substrings) {
    if (grepl(s,longString)) {
      return(s)
    }
  }
  return('x') # If none of the substrings is present in the long string, return 'x'
  
}

## Get the condition name from a string that may contain information about
##    the scan number, read number, and/or experiment name
get_condition_from_string <- function(expNames,string) {
  
  string <- strsplit(string,'Scan')[[1]][1]
  string <- strsplit(string,'Read')[[1]][1]
  string <- strsplit(string,' Rep ')[[1]][1]
  
  for (expN in expNames) {
    string <- strsplit(string,expN)[[1]][1]
  }
  
  return(string)
}

## Get the relevant factor to split the autocorrelation curves into groups
## splitFactor is defined in ui_analysis_box.R, it will be selected by the user
## Input example - splitFactor = 'Experiment'
##               - autocorrelationData should be a
##                 dataframe with at least 3 columns - 'variable', 'value', 'time')
##               - expNames = c("CA 5 replicates","BSA 3 replicates")
## Returns the autocorrelationData dataset with a new column named 'factor'
add_factor_column <- function(autocorrelationData,splitFactor,expNames) {
  
  current_factor <- as.character(autocorrelationData$variable) 
  
  if (splitFactor %in% c('Scan','Read')) {
    new_factor <- sapply(current_factor,function(x) {
      
        xSplit <- strsplit(x,splitFactor)[[1]][2]
        
        if (is.na(xSplit)) {
          return('x') # colour will be the same for all, i.e., only one experiment was loaded
        } else {
          return(strsplit(xSplit," ")[[1]][2])
        }
        
      })
  } else if (splitFactor == 'Experiment') {
    
    new_factor <- sapply(current_factor,function(x) get_exp_name(expNames,x))
    
  } else if (splitFactor == 'Condition') {
    
    new_factor <- sapply(current_factor,function(x) get_condition_from_string(expNames,x))
    
  } else {
    
    new_factor <- current_factor
    
  }
  
  autocorrelationData$factor <- new_factor
  
  return(autocorrelationData)
}

# Select only the samples where the maximum value of the residuals is below a certain threshold
# acData is a dataframe with at least 3 columns - 'time', 'variable', 'value'
# residuals is a vector of the same length as the number of rows in acData
# threshold is the max acepted value of residuals, if one autocorrelation curve contains 
# one residual (in absolut values) above this value, it will be removed
get_selected_variables_based_on_residuals <- function(acData,residuals,threshold) {
  
  acData$residuals <- abs(residuals)
  
  selVariables <- (acData %>% group_by(variable) %>% 
                     summarise(maxResiduals = max(residuals)) %>% 
                     filter(maxResiduals <= threshold))$variable
  
  return(selVariables)
  
}

# Input:  a string of numbers separated by spaces, i.e.  "1 2 3"
# Output: a vector with the corresponding numbers, i.e. c(1,2,3)
get_positions <- function(guess_positions_list) {
  
  guess_positions_list <- trimws(guess_positions_list)
  return(as.numeric(strsplit(guess_positions_list,"\\s+")[[1]]))
  
}

# Return the standard deviation of a probability distribution
# w is the vector of probabilities, m the values of the variable we want to analyse
getSD <- function(w,m) {
  
  mu <- sum(w*m) / sum(w)
  sqrt(sum(w*((m-mu)**2)))
  
}

# Used to estimate the peaks and return the value of the hydrodynamic radius and the
# peak contribution
# We will detect as many peaks as selected by the user, the minimum intensity should be > 2 %
getPeakInfo <- function(estimatedContributions,
                        estimatedContributionsMassWeighted,
                        leftBounds="1 1e2 1e4",
                        rightBounds="1e2 1e4 1e6",
                        estimationMethod="peakMax") {
  
  # We will detect peak present in different regions 
  # defined in the vectors leftBounds and rightBounds
  # estimationMethod should be 'harmonicMean'  or 'peakMax'
  
  leftBounds  <- get_positions(leftBounds)
  rightBounds <- get_positions(rightBounds)
  
  # Quality check of the input - return empty dataframe
  emptyDF1 <- data.frame("Error"="Please give only numeric values in the 'Peak selection box'.") 
  emptyDF2 <- data.frame("Error"="The number of left and right bounds is different.") 
  
  c1 <- any(is.na(leftBounds))
  c2 <- any(is.na(rightBounds))
  c3 <- length(leftBounds) != length(rightBounds)
  
  if (c1 | c2) return(emptyDF1)
  if (c3)      return(emptyDF2)
  
  samples  <- unique(estimatedContributions$variable)
  nRegions <- length(leftBounds)
  
  colNames1 <- paste0('#',1:nRegions,' ([',leftBounds,' - ',rightBounds,'] nm) - Hr')
  
  # Small hack to use instead the diffusion coefficients!!!!
  if (colnames(estimatedContributions)[1] == "diff") {
    estimatedContributions$hr <- estimatedContributions$diff
    colNames1 <- paste0('#',1:nRegions,' ([',leftBounds,' - ',rightBounds,'] m^2/s) - Diff. coef.')
  }
  
  colNames2 <- paste0('#',1:nRegions,' - Intensity contribution (%)')
  colNames3 <- paste0('#',1:nRegions,' - Standard deviation (nm)')
  colNames4 <- paste0('#',1:nRegions,' - Mass contribution (%)')
  
  colNames  <- c(sapply(1:nRegions, function(i) c(colNames1[i],colNames2[i],
                                                  colNames3[i],colNames4[i])))
  colNames  <- c('Sample',colNames)
  
  df <- data.frame(matrix(nrow = 0, ncol = 1+nRegions*4)) 
  colnames(df) <- colNames
  
  for (s in samples) {
    
    # Work with the desired sample
    x     <- estimatedContributions[estimatedContributions$variable == s,]
    
    # Same order and length as x, mass weighted (instead of intensity weighted)
    x2    <- estimatedContributionsMassWeighted[estimatedContributionsMassWeighted$variable == s,]
     
    # Get a list of dataframes according to the bounds
    xs    <- lapply(1:nRegions, function(i) x %>% filter(hr >= leftBounds[i] & hr <= rightBounds[i]) )
    peaks <- lapply(xs,         function(x) findpeaks(x$value,nups = 1,
                                                      minpeakheight=0.02) )
    
    xs2 <- lapply(1:nRegions, function(i) x2 %>% filter(x$hr >= leftBounds[i] & x$hr <= rightBounds[i]) )
    
    i      <- 0
    hrs    <- c()
    contr  <- c()
    contM  <- c()
    sds    <- c()
    
    for (p in peaks) {
      i <-  i + 1
      if (!is.null(p)) {
        
        p           <- data.frame(p)
        colnames(p) <- c("height","index","begin","end")
        
        # Select only the highest peak
        if (nrow(p) > 1) p <- p %>% arrange(-height)
        
        intensity <- p[1,1]
        idx       <- p[1,2]
        begin     <- p[1,3]
        end       <- p[1,4]
        
        # Get peak max Hr
        hrEst     <- xs[[i]]$hr[idx]
        sel       <- (begin-1):(end+1)
        xsSubset  <- xs[[i]][sel,]
        xs2Subset <- xs2[[i]][sel,]
        sd   <- getSD(xsSubset$value,xsSubset$hr)
        
        cont     <- signif(sum(xsSubset$value)  * 100,3)
        contMass <- signif(sum(xs2Subset$value) * 100,3)
        
        # Calculate and replace Hr based on the harmonic mean
        if (estimationMethod == "harmonicMean") {
          weightSum   <- sum(xsSubset$value)
          reciprocals <- sum(xsSubset$value / xsSubset$hr)
          hrEst       <- weightSum / reciprocals
          
        }
        
        sd <- signif(sd, 3)
        hr <- format(hrEst,scientific = F,digits = 3)
        
      } else {
        hr   <- "Max intensity < 2 %"
        cont <- sd <- contMass <- "" 
      }
      
      hrs    <- c(hrs,hr)
      contr  <- c(contr,cont)
      contM  <- c(contM,contMass)
      sds    <- c(sds,sd)
      
    }
    
    colValues  <- c(sapply(1:nRegions, function(i) c(hrs[i],contr[i],sds[i],contM[i])))
    colValues  <- c(s,colValues)
    df[nrow(df)+1,] <- colValues
    
  }
  
  return(df)
}

#compute the dispersity index Mw / Mn
# where Mw is the weight average molecular weight and
# Mn is the number average molecular weight
# in our case, we replace masses with Hrs!
## Input example - freq      = c(5,6)
##               - property  = c(1,10)
## Returns the dispersity index
dispersity <- function(freq,property) {
# check if this function is useful!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  n   <- length(freq)
  nm  <- freq*property
  
  #Number average molecular weight
  nAvg <- sum(nm) / sum(freq)
  
  wf  <- nm / sum(nm) 
  wm  <- property * wf 
  
  #Weight Average Molecular Weight 
  wAvg <- sum(wm)

  return(wAvg / nAvg)
  
}

dlsData <- function(dlsAnalyzer) {
  
  expNames                    <- dlsAnalyzer$experimentNames
  autocorrelationList         <- dlsAnalyzer$getExperimentProperties('autocorrelation')
  timeList                    <- dlsAnalyzer$getExperimentProperties('time')
  experimentsSamplesMetaData  <- dlsAnalyzer$getExperimentProperties('sampleInfoRelevant')
  
  idx <- sapply(autocorrelationList, function(x) dim(x)[2] > 0 )
  
  if (sum(idx)==0) return(NULL) # return null in the case that there is no data! 
  
  data                   <- formatDlsInfoForPlotting(autocorrelationList[idx],timeList[idx],
                                                     expNames[idx],experimentsSamplesMetaData[idx])
  
  return(data)
}

getPalette <- function(nColors) {
  
  if (nColors <= 9) return(RColorBrewer::brewer.pal(9, "Set1"))
  if (nColors <= 12) return(RColorBrewer::brewer.pal(12, "Set3"))
  
  return( colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nColors) )
  
}

formatNorms <- function(expN,sNames,rn,pn,idSel,alphaVec) {
  
  # Obtain a dataframe containing the penalty & fidelity terms
  # Input: 
  #        'sNames' - sample names 
  #        'expN'   - experiment names
  #        'rn'     - list of matrixes, fidelity term n*m, m=number of curves, n=tested reg terms 
  #        'pn'     - list of matrixes, penalty term n*m, m=number of curves, n=tested reg terms
  #        'idSel'  - list of lists containing the selected corners
  #        'alphaVec' - vector of the tested regularisation terms
  # Output:
  #         a list with two dataframes
  #         'lCurveDf'  - dataframe containing the penalty & fidelity terms (and alpha values)
  #         'dfSel'     - dataframe containing the penalty & fidelity terms (only the corner!)
  
  dfs    <- list()
  dfsSel <- list()
  
  for (i in 1:length(expN)) {
    
    rnI <- rn[[i]]
    pnI <- pn[[i]]
    
    colnames(rnI) <- sNames[[i]]$conditions
    colnames(pnI) <- sNames[[i]]$conditions
    
    if (length(expN) > 1){
      colnames(rnI) <- paste0(colnames(rnI),' ',expN[i])
      colnames(pnI) <- paste0(colnames(pnI),' ',expN[i])
    }
    
    selPointsRn <- sapply(1:ncol(rnI), function(x) {rnI[idSel[[i]][x],x]})
    selPointsPn <- sapply(1:ncol(pnI), function(x) {pnI[idSel[[i]][x],x]})
    rnI         <- melt(rnI,value.name = "fidelity") 
    pnI         <- melt(pnI,value.name = "penalty") 
    
    dfs[[i]]         <- inner_join(rnI,pnI) %>%  rename(variable = Var2)
    dfs[[i]]$alpha   <- rep(alphaVec,length(unique( dfs[[i]]$variable)))
    dfsSel[[i]] <- data.frame(variable = names(selPointsRn),
                              fidelity = selPointsRn,
                              penalty  = selPointsPn
                              )
    
  }
  
  df    <- do.call(rbind,dfs)
  dfSel <- do.call(rbind,dfsSel)
  rownames(dfSel) <- NULL
  
  return(list('lCurveDf'=df,'dfSel'=dfSel))
}





