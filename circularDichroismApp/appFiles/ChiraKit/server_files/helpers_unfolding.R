# Useful scripts to process the unfolding datasets (thermal and chemical)
# and also to process the 'custom' datasets

# Parse selected wavelengths
# Expected input formats numbers separated by spaces and/or dashes
# E.g., '220 224', '220-224', and '220-224 228'

parse_selected_wavelengths <- function(input_string) {
  parts <- unlist(strsplit(input_string, " "))  # Split by spaces
  
  parsed_values <- vector("numeric", length = 0)
  
  for (part in parts) {
    if (grepl("-", part)) {
      range_parts <- unlist(strsplit(part, "-"))  # Split by dash
      start <- as.numeric(range_parts[1])
      end   <- as.numeric(range_parts[2])
      
      if (!is.na(start) && !is.na(end)) {
        parsed_values <- c(parsed_values, seq(start, end))
      }
    } else {
      value <- as.numeric(part)
      if (!is.na(value)) {
        parsed_values <- c(parsed_values, value)
      }
    }
  }
  
  parsed_values <- sort(unique(parsed_values))
  return(parsed_values)
}

## Use spectra names to merge CD datasets
## Requires:
##  - relevantSpectra:  vector of strings, each string represents one CD curve
##  - cdAnalyzer:       python object to handle CD data
## Returns:
##  A merged dataset with (n + 1) columns where n
##  is the number of selected spectra and '+ 1' comes from the shared wavelengths
get_signal_dfs_from_selected_spectra <- function(relevantSpectra,cdAnalyzer) {
  
  # Retrieve a list of internalIDs
  names              <- cdAnalyzer$experimentNames
  internalIDs        <- cdAnalyzer$getExperimentProperties('internalID')
  # Assign experiment names as names for the list objects
  names(internalIDs) <- names
  
  # Find the index of the relevant spectra using as proxy the flattened list of internalIDs 
  ids2find <- c()
  
  for (spectra in relevantSpectra) {
    ids2find <- c(ids2find,which(unlist(internalIDs) == spectra))
  }
  
  # Match the indexes with the indexes of the non-flattened list
  ids                <- found_ids(internalIDs,ids2find)
  
  exps <- names(ids)
  
  # Retrieve the list of wavelength vectors 
  wlSel     <- lapply(exps, function(exp) cdAnalyzer$experimentsModif[[exp]]$wavelength)
  
  # Retrieve the list of signal vectors 
  signalSel <- lapply(1:length(exps), function(i) cdAnalyzer$experimentsModif[[exps[i]]]$signalDesiredUnit[,ids[i]])
  
  # Create the dataframes of wavelength versus signal
  signal_dfs <- lapply(1:length(exps), function (i) {
    
    df           <- data.frame(wlSel[[i]],signalSel[[i]])
    colnames(df) <- c('Wavelength',paste0('curve',i))
    
    return(df)
  })
  
  totalDFs <- length(signal_dfs)
  
  # Merge the dataframes
  if (totalDFs == 1) {
    merged <- signal_dfs[[1]]
  } else {
    
    merged <- merge(signal_dfs[[1]],signal_dfs[[2]],by='Wavelength')
    
    if (totalDFs > 2 ) {
      
      for (i in 3:totalDFs) {
        
        merged <- merge(merged,signal_dfs[[i]],by='Wavelength')
        
      }
    }
  }
  
  return( merged )
  
}

# Create a dataframe to plot the thermal curves
# Requires:
# - the python object to handle CD data, cdAnalyzer
# - a string to decide if we want to retrieve experimental data (and which wavelength range) or fitted data
#   should be 'signal_useful', 'signalDesiredUnit', 'signal_predicted', 'fitted_spectra'
# Returns:
# - a merged dataframe with 4 columns 
# The columns are 'wavelength', 'temperature' , 'value' and 'legend'
generate_thermal_ramp_df <- function(cdAnalyzer,signal_type='signal_useful') {
  
  # Initialize an empty list to hold individual dataframes
  dfs <- list()
  
  exps <- cdAnalyzer$experimentNamesThermal
  
  counter   <- 0
  
  # Loop through each thermal experiment in cdAnalyzer
  for (exp in exps) {
    
    # Increment the counter for thermal ramps
    counter <- counter + 1
    
    # Extract relevant data from the Python object
    signals     <- cdAnalyzer$experimentsThermal[[exp]][[signal_type]]
    
    if (signal_type %in% c('signal_useful','signal_predicted')) {
      wl          <- cdAnalyzer$experimentsThermal[[exp]]$wavelength_useful
    } else {
      wl          <- cdAnalyzer$experimentsThermal[[exp]]$wavelength
    }
    
    temperature <- cdAnalyzer$experimentsThermal[[exp]]$temperature
    
    # Create a dataframe with wavelength and signal data
    df <- data.frame(wl, signals)
    colnames(df) <- c('wavelength', paste0(1:length(temperature),'_',temperature))
    
    # Reshape the dataframe using 'melt' to make it suitable for plotting
    df           <- melt(df, id = c("wavelength"))
    colnames(df) <- c('wavelength', 'temperature', 'value')
    
    df[,'temperature'] <- sapply(as.character(df[,'temperature']), 
                                 function (x) strsplit(x,'_')[[1]][2])
    
    df[,'temperature'] <- as.numeric(df[,'temperature'])
    
    # Add a 'legend' column to identify the experiment
    df[,'legend'] <- exp 
    
    # Store the dataframe in the 'dfs' list
    dfs[[counter]] <- df
    
  }
  
  # Combine individual dataframes into a single dataframe
  df <- NULL
  if (length(dfs) > 0 ) df <- do.call(rbind, dfs)
  
  # Return the final dataframe for plotting
  return(df)
}

# Create a dataframe to plot the chemical unfolding ramps
# Requires:
# - the python object to handle CD data, cdAnalyzer
# - a string to decide if we want to retrieve experimental data (and which wavelength range) or fitted data
#   should be 'signal_useful', 'signalDesiredUnit', 'signal_predicted', 'fitted_spectra'
# Returns:
# - a merged dataframe with 4 columns 
# The columns are 'wavelength', 'chem_conc' , 'value' and 'legend'
generate_chemical_unfolding_df <- function(cdAnalyzer,signal_type='signal_useful') {
  
  # Initialize an empty list to hold individual dataframes
  dfs <- list()
  
  exps <- cdAnalyzer$experimentNamesChemical
  
  counter   <- 0
  # Loop through each thermal experiment in cdAnalyzer
  for (exp in exps) {
    
    # Increment the counter for thermal ramps
    counter <- counter + 1
    
    # Extract relevant data from the Python object
    signals     <- cdAnalyzer$experimentsChemical[[exp]][[signal_type]]
    
    if (signal_type %in% c('signal_useful','signal_predicted')) {
      wl          <- cdAnalyzer$experimentsChemical[[exp]]$wavelength_useful
    } else {
      wl          <- cdAnalyzer$experimentsChemical[[exp]]$wavelength
    }
    
    chem_conc   <- cdAnalyzer$experimentsChemical[[exp]]$chem_concentration
    
    # Create a dataframe with wavelength and signal data
    df <- data.frame(wl, signals)
    
    colnames(df) <- c('wavelength', paste0(1:length(chem_conc),'_',chem_conc))
    
    # Reshape the dataframe using 'melt' to make it suitable for plotting
    df           <- melt(df, id = c("wavelength"))
    colnames(df) <- c('wavelength', 'chem_conc', 'value')
    
    df[,'chem_conc'] <- sapply(as.character(df[,'chem_conc']), 
                               function (x) strsplit(x,'_')[[1]][2])
    
    df[,'chem_conc'] <- as.numeric(df[,'chem_conc'])
    
    # Add a 'legend' column to identify the experiment
    df[,'legend'] <- exp 
    
    # Store the dataframe in the 'dfs' list
    dfs[[counter]] <- df
    
  }
  
  # Combine individual dataframes into a single dataframe
  df <- NULL
  if (length(dfs) > 0 ) df <- do.call(rbind, dfs)
  
  # Return the final dataframe for plotting
  return(df)
}

# To do - Adapt this function to include the second experimental parameter

# Create a dataframe to plot the custom analysis curves
# Requires:
# - the python object to handle CD data, cdAnalyzer
# - a string to decide if we want to retrieve experimental data (and which wavelength range) or fitted data
#   should be 'signal_useful', 'signalDesiredUnit', 'signal_predicted', 'fitted_spectra'
# Returns:
# - a merged dataframe with 4 columns 
# The columns are 'wavelength', firstExperimentalParam ,
# 'value' and 'legend'

generate_custom_df <- function(cdAnalyzer,signal_type='signal_useful') {
  
  # Initialize an empty list to hold individual dataframes
  dfs <- list()
  
  exps <- cdAnalyzer$experimentNamesCustom
  
  counter   <- 0
  
  # Loop through each thermal experiment in cdAnalyzer
  for (exp in exps) {
    
    # Increment the counter for thermal ramps
    counter <- counter + 1
    
    # Extract relevant data from the Python object
    signals     <- cdAnalyzer$experimentsCustom[[exp]][[signal_type]]
    
    if (signal_type %in% c('signal_useful','signal_predicted')) {
      wl          <- cdAnalyzer$experimentsCustom[[exp]]$wavelength_useful
    } else {
      wl          <- cdAnalyzer$experimentsCustom[[exp]]$wavelength
    }
    
    firstExpParam  <- cdAnalyzer$experimentsCustom[[exp]]$first_measurement_dimension
    
    firstExpParamName  <- cdAnalyzer$experimentsCustom[[exp]]$first_exp_param_name
    
    #secondExpParam <- cdAnalyzer$experimentsCustom[[exp]]$second_measurement_dimension
    
    # Create a dataframe with wavelength and signal data
    df <- data.frame(wl, signals)
    colnames(df) <- c('wavelength', paste0(1:length(firstExpParam),'_',firstExpParam))
    
    # Reshape the dataframe using 'melt' to make it suitable for plotting
    df           <- melt(df, id = c("wavelength"))
    colnames(df) <- c('wavelength', firstExpParamName, 'value')
    
    df[,firstExpParamName] <- sapply(as.character(df[,firstExpParamName]), 
                                 function (x) strsplit(x,'_')[[1]][2])
    
    df[,firstExpParamName] <- as.numeric(df[,firstExpParamName])
    
    # Add a 'legend' column to identify the experiment
    df[,'legend'] <- exp 
    
    # Store the dataframe in the 'dfs' list
    dfs[[counter]] <- df
    
  }
  
  # Combine individual dataframes into a single dataframe
  df <- NULL
  if (length(dfs) > 0 ) df <- do.call(rbind, dfs)
  
  # Return the final dataframe for plotting
  return(df)
}

## Retrieve the adequate python class and experiment names based on the analysis type
get_py_class_and_exp_names <- function(cdAnalyzer,type) {
  
  if (type == 'Thermal') {
    
    exps      <- cdAnalyzer$experimentNamesThermal
    py_object <- cdAnalyzer$experimentsThermal
    
  } else if (type == 'Chemical') {
    
    exps      <- cdAnalyzer$experimentNamesChemical
    py_object <- cdAnalyzer$experimentsChemical
    
  } else{
    
    exps      <- cdAnalyzer$experimentNamesCustom
    py_object <- cdAnalyzer$experimentsCustom
    
  } 
  
  return(list('exps'=exps,'py_object'=py_object))
}

# Generate the dataframe with the fitted parameters or relative errors
# Requires:
# - the python object to handle CD data, cdAnalyzer
# - the type of unfolding model, either 'Thermal' or 'Chemical'
# - a boolean to decide if we want to retrieve the parameters or the relative errors
# Returns:
# - a dataframe with the errors 
get_fitted_params_unfolding <- function(cdAnalyzer,type='Thermal',errors=FALSE) {
  
  res       <- get_py_class_and_exp_names(cdAnalyzer,type)
  exps      <- res$exps
  py_object <- res$py_object
  
  dfs <- list()
  
  i <- 0
  for (exp in exps) {
    i  <- i + 1
    if (errors) {
      
      df <- py_object[[exp]]$fit_rel_errors
    } else {
      df <- py_object[[exp]]$fit_params
    }
    
    dfs[[i]] <- df
    
  }
  
  df <- do.call(rbind,dfs)
  
  numeric_cols      <- colnames(df)[-length(colnames(df))]
  for (numeric_col in numeric_cols) {
    df[,numeric_col] <- signif(as.numeric(df[,numeric_col]),5)
  }
  return(df)
}

## Generate dataframe to plot the basis spectra
get_basis_spectra_df <- function(cdAnalyzer,type='Thermal') {
  
  res       <- get_py_class_and_exp_names(cdAnalyzer,type)
  exps      <- res$exps
  py_object <- res$py_object
  
  dfs <- list()
  
  i <- 0
  for (exp in exps) {
    
    i        <- i + 1
    signal   <- py_object[[exp]]$basis_spectra
    wl       <- py_object[[exp]]$wavelength
    
    # Create a dataframe with wavelength and signal data
    df <- data.frame(wl, signal)
    
    colnames(df) <- c('wavelength', paste0('k = ',1:ncol(signal)))
    # Reshape the dataframe using 'melt' to make it suitable for plotting
    df           <- melt(df, id = c("wavelength"))
    colnames(df) <- c('wavelength', 'k', 'value')
    
    # Add a 'legend' column to identify the experiment
    df[,'legend'] <- exp 
    
    dfs[[i]] <- df
    
  }
  
  df <- do.call(rbind,dfs)
  
  return(df)
}

## Generate dataframe to plot the coefficients change
get_coefficients_df <- function(cdAnalyzer,type='Thermal') {
  
  res       <- get_py_class_and_exp_names(cdAnalyzer,type)
  exps      <- res$exps
  py_object <- res$py_object
  
  dfs <- list()
  
  i <- 0
  for (exp in exps) {
    
    i        <- i + 1
    coeff    <- py_object[[exp]]$coefficients
    
    if (type == 'Thermal')  {
      mment_factor <- py_object[[exp]]$temperature
      colname2     <- 'Temperature'
    } else if (type == 'Chemical') {
      mment_factor <- py_object[[exp]]$chem_concentration
      colname2     <- 'chem_conc'
    } else {
      mment_factor <- py_object[[exp]]$first_measurement_dimension
      colname2     <- py_object[[exp]]$first_exp_param_name
    }
    
    k        <- paste0('k = ',1:nrow(coeff))
    
    # Create a dataframe with wavelength and signal data
    df <- data.frame(k,coeff)
    
    colnames(df) <- c('k', paste0(1:length(mment_factor),'_',mment_factor))
  
    # Reshape the dataframe using 'melt' to make it suitable for plotting
    df           <- melt(df, id = c("k"))
    
    colnames(df) <- c('k', colname2, 'value')
    
    df[,2] <- sapply(as.character(df[,2]), function (x) strsplit(x,'_')[[1]][2])
    
    df[,2] <- as.numeric(df[,2])
    
    # Add a 'legend' column to identify the experiment
    df[,'legend'] <- exp 
    
    # Store the dataframe in the 'dfs' list
    dfs[[i]] <- df
    
  }
  
  df <- do.call(rbind,dfs)
  
  return(df)
}

get_explained_variance_df <- function(cdAnalyzer,type) {
  
  res       <- get_py_class_and_exp_names(cdAnalyzer,type)
  exps      <- res$exps
  py_object <- res$py_object
  
  dfs <- list()
  
  i <- 0
  for (exp in exps) {
    
    i                      <- i + 1
    explained_variance     <- py_object[[exp]]$explained_variance
    k                      <- 1:length(explained_variance)
    df                     <- data.frame(k,explained_variance)
    df[,'legend']          <- exp 
    dfs[[i]]               <- df
  }
  
  df <- do.call(rbind, dfs)
  
  return(df)
}

# Fix names before exporting melting data
# Requires -  a dataframe with four columns:
# 'wavelength', 'temperature / chem_conc', 'value', 'legend'
# Output - same dataframe with new column names
reassign_unfolding_df_colnames <- function(df,method = 'fixedWL') {
  
  svd_based <- method == 'SVD'
  pca_based <- method == 'PCA'
  
  if (svd_based) {
    colnames(df)[1] <- 'SVD_basis_spectra'
    colnames(df)[3] <- 'SVD_coefficient'
  } else if (pca_based) {
    colnames(df)[1] <- 'PCA_basis_spectra'
    colnames(df)[3] <- 'PCA_coefficient'
  } else {
    colnames(df)[1] <- 'wavelength_(nm)'
    colnames(df)[3] <- 'CD_signal_value'
  }
  
  if (grepl('tempera',colnames(df)[2],ignore.case = T)) {
    colnames(df)[2] <- 'temperature_(°C)'
  } 
  
  if (grepl('chem_conc',colnames(df)[2],ignore.case = T)) {
    colnames(df)[2] <- 'denaturant_concentration_(M)'
  } 
  
  return(df)
}

generate_fractions_df <- function(cdAnalyzer,type='Thermal') {
  
  res       <- get_py_class_and_exp_names(cdAnalyzer,type)
  exps      <- res$exps
  py_object <- res$py_object
  
  dfs <- list()
  
  i <- 0
  for (exp in exps) {
    
    i                      <- i + 1
    fractions              <- py_object[[exp]]$fractions
    df                     <- data.frame(fractions)
    
    if (type == 'Thermal')  {
      mment_factor <- py_object[[exp]]$temperature
      colname     <- 'temperature'
    } else {
      mment_factor <- py_object[[exp]]$chem_concentration
      colname     <- 'chem_conc'
    } 
    
    df[,colname]      <- mment_factor
    df                <- reshape2::melt(df,id = c(colname))
    df[,'legend']     <- exp 
    dfs[[i]]          <- df
  }
  
  df <- do.call(rbind, dfs)
  
  return(df)
  
}


