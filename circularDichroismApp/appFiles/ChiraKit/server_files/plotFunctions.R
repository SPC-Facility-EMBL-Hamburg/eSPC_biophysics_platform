plotCDexperiments <- function(cdAnalyzer,workingUnits,
                              plot_width,plot_height,plot_type,axis_size,
                              legends,colorPalette,sels,useMilliDeg = FALSE,
                              fig=NULL) {
  
  # Initialize the plot
  if (is.null(fig)) fig <- plot_ly()
  
  if (useMilliDeg) {
    # Fetch the differential absorbance signal
    signalsAll      <- cdAnalyzer$get_experiment_properties_modif('signalAbs')
    # Replace working units
    workingUnits    <- 'millidegrees'
  } else {
    signalsAll      <- cdAnalyzer$get_experiment_properties_modif('signalDesiredUnit')
  }
  
  wlsAll          <- cdAnalyzer$get_experiment_properties_modif('wavelength')
  spectraNamesAll <- cdAnalyzer$get_experiment_properties_modif('spectraNames')

  isFake                 <- unlist(cdAnalyzer$get_experiment_properties_modif('isFakeExperiment'))
  fakeExperimentSignal   <- unlist(cdAnalyzer$get_experiment_properties_modif('fakeExperimentSignal'))
  
  # Detect if we need to reduce data points
  nPointsPerExp   <- unlist(lapply(signalsAll, function(x) ncol(x) * nrow(x)))
  nPoints         <- sum(nPointsPerExp,na.rm = T) 
  reducePoints    <- nPoints > 4000 
  
  counter <- 0
  i       <- 0
  
  for (expName in cdAnalyzer$experimentNames) {
   
    i <- i+1
    
    signals     <- signalsAll[[i]]
    
    # Check that 1 ) the signal in the desired unit was computed. 
    # E.g, setting concentration to zero will give NA for the molar ellipticity
    
    # and 2) if the experiment is of type 'fake', the 'fakeExperimentSignal'
    # matches the desired unit to plot
    # Otherwise, skip to next iteration
    
    c1 <- all(is.na(signals))
    c2 <- isFake[i] & (fakeExperimentSignal[i] != workingUnits)
      
    if (c1 | c2) {
      counter <- counter + length(spectraNamesAll[[i]])
      next
    } 
    
    wavelength  <- wlsAll[[i]]
    
    for (ii in 1:ncol(signals)) {
      counter <- counter + 1
      signal  <- signals[,ii]
      
      # in place conversion - absorbance to millidegrees
      if (useMilliDeg) signal <- signal * 32980
      
      if (sels[counter]) {
        
        df <- data.frame('wavelength'=wavelength,signal)
        df <- df[order(df$wavelength), ]
        
        if (reducePoints) {
          df <- df[seq(1, nrow(df), by = 2), ]
        }
        
        # plot lines if we have at least two wavelengths
        if (nrow(df) > 1) {
          fig <- fig %>% add_trace(data=df,color=I(colorPalette[counter]),x=~wavelength,y=~signal,
                                   type = 'scatter', mode = 'lines',
                                   name = legends[counter],list(width = 2))
        } else {
          fig <- fig %>% add_trace(data=df,color=I(colorPalette[counter]),x=~wavelength,y=~signal,
                                   type = 'scatter', mode = 'markers',
                                   name = legends[counter])
        }
      }
    }
  }
  
  minWL <- min(sapply(wlsAll, min)) - 5
  maxWL <- max(sapply(wlsAll, max)) + 5
  
  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  y <- list(title = workingUnits2ProperLabel(workingUnits),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),showgrid = F)
   
  fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))
  
  fig <- configFig(fig,paste0("CDspectra_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  
  return(fig)
}

sesca_plot <- function(sescaPyClass,plot_width,plot_height,plot_type,axis_size,
    ref_wl=NULL,ref_signal=NULL,ref_name=NULL,average_ensemble=FALSE) {

  fig <- plot_ly()

  spectraNames <- as.character(sescaPyClass$spectraNames)
  wl           <- sescaPyClass$wavelength
  signals      <- sescaPyClass$CD_Spectra

  idxContr <- grepl("BB contr.", spectraNames) | grepl("SC contr.", spectraNames)

  # Remove names containing the pattern 'BB' and 'SC'
  spectraNames_filt <- spectraNames[!idxContr]

  if (average_ensemble) {

    spectraNames <- c('Ensemble average')
    signals_filtered <- signals[, !idxContr]

    if (is.vector(signals_filtered)) {
      signals_filtered <- matrix(signals_filtered, ncol = 1)
    }

    signals      <- rowMeans(signals_filtered,na.rm = T)
    signals      <- data.frame(signals)
  }

  colorPalette <- getPalette(length(spectraNames_filt))

  if (length(colorPalette) < length(spectraNames)) {
    colorPalette <- rep(colorPalette, each = 3)
  } else {
    colorPalette <- colorPalette
  }

  for (i in 1:length(spectraNames)) {

    signal <- signals[,i]

    df <- data.frame('wl'=wl,'signal'=signal)

    if (!is.null(ref_wl)) {
      df <- df[df$wl >= min(ref_wl) & df$wl <= max(ref_wl),]
    }

    linetype <- ifelse(grepl("BB contr.|SC contr.", spectraNames[i]), "dash", "solid")

    fig <- fig %>% add_trace(data=df,x=~wl,y=~signal,
                             type = 'scatter', mode = 'lines',
                             name = spectraNames[i],
                             line = list(color = I(colorPalette[i]),opacity = 0.6,dash = linetype),
                             showlegend = TRUE)

  }

  if (!is.null(ref_wl) & !is.null(ref_signal)) {
     fig <- fig %>% add_trace(x=~ref_wl,y=~ref_signal,
                             type = 'scatter', mode = 'markers',
                             name = ref_name,
                             showlegend = TRUE)

    minWL <- min(ref_wl) - 5
    maxWL <- max(ref_wl) + 5

  } else {
    minWL <- min(sescaPyClass$wavelength) - 5
    maxWL <- max(sescaPyClass$wavelength) + 5
  }

  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size),
            tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)

  y <- list(title = workingUnits2ProperLabel('meanUnitMolarExtinction'),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),showgrid = F)

  fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))

  fig <- configFig(fig,paste0("CDspectra_sesca_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)

  return(fig)
}



# Plot the voltage signal
plotCDexperimentsHT <- function(cdAnalyzer,
                              plot_width,plot_height,
                              plot_type,axis_size,
                              legends,colorPalette,sels) {
  
  fig          <- plot_ly()
  
  signalsAll   <- cdAnalyzer$get_experiment_properties_modif('signalHT')
  wlsAll       <- cdAnalyzer$get_experiment_properties_modif('wavelength')
  
  counter <- 0
  i       <- 0
  
  for (expName in cdAnalyzer$experimentNames) {
    
    i <- i+1
    
    signals     <- signalsAll[[i]]
    wavelength  <- wlsAll[[i]]
    
    for (ii in 1:ncol(signals)) {
      counter <- counter + 1
      signal  <- signals[,ii]
      
      if (length(signal) != length(wavelength)) return(NULL)
      
      df      <- data.frame('wavelength'=wavelength,signal)
      
      # Skip trace if constant value only
      if (length(unique(df[,'signal'])) == 1) {
        next
      }
      
      if (sels[counter]) {
        fig <- fig %>% add_trace(data=df,color=I(colorPalette[counter]),x=~wavelength,y=~signal,
                                 type = 'scatter', mode = 'markers+lines',
                                 marker = list(size=4),
                                 line = list(width = 1.8),
                                 name = legends[counter])
      }
    }
  }
  
  minWL <- min(sapply(wlsAll, min)) - 5
  maxWL <- max(sapply(wlsAll, max)) + 5

  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  y <- list(title = 'High Tension (HT) Voltage',
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),showgrid = F)
  
  fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))
  
  fig <- configFig(fig,paste0("CDspectra_voltage_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  return(fig)
}

# Plot the CD signal and the voltage curves using a double y-axis
plot_cd_and_voltage <- function(cdAnalyzer,workingUnits,
                              plot_width,plot_height,plot_type,axis_size,
                              legends,colorPalette,sels) {
  
  fig         <- plot_ly()
  
  signalsHTAll    <- cdAnalyzer$get_experiment_properties_modif('signalHT')
  signalsAll      <- cdAnalyzer$get_experiment_properties_modif('signalDesiredUnit')
  wlsAll          <- cdAnalyzer$get_experiment_properties_modif('wavelength')
  spectraNamesAll <- cdAnalyzer$get_experiment_properties_modif('spectraNames')
  
  isFake                 <- unlist(cdAnalyzer$get_experiment_properties_modif('isFakeExperiment'))
  fakeExperimentSignal   <- unlist(cdAnalyzer$get_experiment_properties_modif('fakeExperimentSignal'))
  
  counter <- 0
  i       <- 0
  
  for (expName in cdAnalyzer$experimentNames) {
    
    i <- i+1
    
    signals     <- signalsAll[[i]]
    signalsHT   <- signalsHTAll[[i]]

    # Check that 1 ) the signal in the desired unit was computed. 
    # E.g, setting concentration to zero will give NA for the molar ellipticity
    
    # and 2) if the experiment is of type 'fake', the 'fakeExperimentSignal'
    # matches the desired unit to plot
    # Otherwise, skip to next iteration
    
    c1 <- all(is.na(signals))
    c2 <- isFake[i] & (fakeExperimentSignal[i] != workingUnits)
    
    if (c1 | c2) {
      counter <- counter + length(spectraNamesAll[[i]])
      next
    }
    
    wavelength  <- wlsAll[[i]]
    
    for (ii in 1:ncol(signals)) {
      
      counter   <- counter + 1
      signal    <- signals[,ii]
      signalHT  <- signalsHT[,ii]
      
      if (sels[counter]) {
        
        df      <- data.frame('wavelength'=wavelength,signal)
        dfHT    <- data.frame('wavelength'=wavelength,signalHT)
        
        fig <- fig %>% add_trace(data=df,color=I(colorPalette[counter]),x=~wavelength,y=~signal,
                                 type = 'scatter', mode = 'lines',
                                 name = legends[counter],list(width = 2))
        
        # Plot voltage data if it is not constant
        if (length(unique(dfHT[,'signalHT'])) > 1) {
          fig <- fig %>% add_trace(data=dfHT,color=I(colorPalette[counter]),x=~wavelength,y=~signalHT,
                                   type = 'scatter', mode = 'markers+lines',
                                   line = list(width = 1.5,opacity = 0.6),
                                   marker = list(size=4),
                                   yaxis = "y2", name = legends[counter],
                                   showlegend = F)
        }
        
      }
      
    }
    
  }
  
  minWL <- min(sapply(wlsAll, min)) - 5
  maxWL <- max(sapply(wlsAll, max)) + 5
  
  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  y <- list(title = workingUnits2ProperLabel(workingUnits),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),showgrid = F)
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend=TRUE,
                        legend = list(x = 1.05, y = 1,font = list(size = axis_size-1)),
                        yaxis2 = list(overlaying = "y", side = "right",
                                      title = "High Tension Voltage",
                                      titlefont = list(size = axis_size), 
                                      tickfont = list(size = axis_size)))
  
  fig <- configFig(fig,paste0("CDspectra_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  return(fig)
}

# Requires a dataframe with 4 columns:
#   'wavelength' / 'k', 'temperature / chem_con', 'value' (CD signal), 'legend'

# Returns the plot of the dataframe
# We plot the temperature versus signal, one subplot per legend, 
# and as many lines as wavelengths 

plot_unfolding_exp <- function(unfolding_exp_data,workingUnits,
                               plot_width=12,plot_height=8,
                               plot_type='svg', axis_size=12,
                               spectra_decomposition_method='fixedWL',
                               useLogAxis=F,
                               xLegend=1,yLegend=1,showTitle=T) {
  
  # Return null if there is no data
  if (is.null(unfolding_exp_data)) return(NULL)

  xaxisLabel <- get_x_axis_label(unfolding_exp_data)
  
  xaxis <- list(title = xaxisLabel,titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),showgrid = F)
  
  yAxisTitle <- workingUnits2ProperLabel(workingUnits)
  svd_or_pca_based <- spectra_decomposition_method %in% c('PCA','SVD')
  
  if (svd_or_pca_based) {
    yAxisTitle <- paste0(spectra_decomposition_method,' coefficient')
  }
  
  yaxis <- list(title = yAxisTitle,titlefont = list(size = axis_size),
                tickfont = list(size = axis_size),showgrid = F)
  
  plot_list   <- list()
  
  unfolding_exp_data <- add_measurement_factor_column(unfolding_exp_data)
  
  if (useLogAxis) {
    
    xaxis[['type']] <- 'log'
    xaxis[['exponentformat']] <- 'power'
    
  } 
  
  tot_cond    <- length(unique(unfolding_exp_data$legend))
  
  i <- 0
  for (leg in unique(unfolding_exp_data$legend)) {
    
    fig <- plot_ly()
    
    i   <- i + 1
    
    df_temp <- unfolding_exp_data[unfolding_exp_data$legend == leg,]

    # Assuming df_temp has a 'group' column indicating different groups
    if (svd_or_pca_based) {
      df_temp$group_var <- df_temp$k
    } else {
      df_temp$group_var <- df_temp$wavelength
    }
      
    unique_groups <- unique(df_temp$group_var)
    
    for (group_value in unique_groups) {
      
      subset_df <- df_temp[df_temp$group_var == group_value, ]
      
      fig <- add_trace(
        fig,
        x = subset_df$measurement_factor,
        y = subset_df$value,
        type = "scatter",
        mode = "markers",
        name = group_value
      )
    }
    
    fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,tot_cond,axis_size)
    
    legend_title <- ifelse(svd_or_pca_based,"Basis spectra","Wavelength (nm)")
    
    # Remove legend title to handle signal of unknown origin
    if (workingUnits == 'yOligomer') legend_title <- 'N-mer concentration'
    if (workingUnits == 'yMonomer')  legend_title <- 'Curve'

    fig <- fig %>% layout(legend = list(title = list(text = legend_title,
                                                     font = list(size = axis_size-1)),
                                        x = xLegend,y = yLegend))
    
    plot_list[[i]] <- fig
    
  }
  
  titleText <- unique(unfolding_exp_data$legend)
  if (!showTitle) titleText <- ''
  
  fig <- plot_list_to_fig(
    paste0('unfolding_curves_',Sys.Date()),plot_list,
    titleText,
    axis_size,plot_type,plot_width,plot_height)
  
  return( fig )
}

# Requires two dataframes with 4 columns (the original signal and the fitted one):
#   'wavelength', 'value' (CD signal), 'temperature', 'legend'
# Returns the plot of the fitted dataframe
# We plot the temperature versus signal, one subplot per legend, 
# and as many lines as wavelengths 

plot_unfolding_fitting <- function(
    unfolding_exp_data, unfolding_fitted_data,
    workingUnits,
    plot_width=12, plot_height=8,
    plot_type='svg', axis_size=12,
    fitted_coefficients_method='fixedWL',useLogAxis=F,
    xLegend=1,yLegend=1,showTitle=T) {
  
  # Return null if there is no data
  if (is.null(unfolding_exp_data)) return(NULL)
  
  xaxisLabel <- get_x_axis_label(unfolding_exp_data)
  
  xaxis <- list(title = xaxisLabel,titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),showgrid = F)
  
  if (useLogAxis) {
    
    xaxis[['type']] <- 'log'
    xaxis[['exponentformat']] <- 'power'
    
  } 
  
  yAxisTitle       <- workingUnits2ProperLabel(workingUnits)
  svd_or_pca_based <- fitted_coefficients_method %in% c('PCA','SVD')
    
  if (svd_or_pca_based) {
    yAxisTitle <- paste0(fitted_coefficients_method,' coefficient')
  }
  
  yaxis <- list(title = yAxisTitle,titlefont = list(size = axis_size),
                tickfont = list(size = axis_size),showgrid = F)
  
  plot_list   <- list()
  
  unfolding_exp_data    <- add_measurement_factor_column(unfolding_exp_data)
  unfolding_fitted_data <- add_measurement_factor_column(unfolding_fitted_data)
  
  tot_cond    <- length(unique(unfolding_exp_data$legend))
  
  i <- 0
  for (leg in unique(unfolding_exp_data$legend)) {
    
    fig <- plot_ly()
    
    i   <- i + 1
    
    df_temp     <- unfolding_exp_data[unfolding_exp_data$legend == leg,]
    df_temp_fit <- unfolding_fitted_data[unfolding_fitted_data$legend == leg,]
    
    # Assuming df_temp has a 'group' column indicating different groups
    unique_groups <- unique(df_temp$wavelength)
    
    for (group_value in unique_groups) {
      
      subset_df     <- df_temp[df_temp$wavelength == group_value, ]
      subset_df_fit <- df_temp_fit[df_temp_fit$wavelength == group_value, ]
      
      fig <- add_trace(
        fig,
        x = subset_df_fit$measurement_factor,
        y = subset_df_fit$value,
        type = "scatter",
        mode = "lines",
        line = list(color = "black",opacity = 0.6),
        showlegend = FALSE 
      )
      
      name <- ifelse(svd_or_pca_based,paste0("k = ", group_value),paste0(group_value))
      
      fig <- add_trace(
        fig,
        x = subset_df$measurement_factor,
        y = subset_df$value,
        type = "scatter",
        mode = "markers",
        name = name,
        showlegend = TRUE  
      )
      
    }
    
    fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,tot_cond,axis_size)
    
    legend_title <- ifelse(svd_or_pca_based,"Basis spectra","Wavelength (nm)")

    # Remove legend title to handle signal of unknown origin
    if (workingUnits == 'yOligomer') legend_title <- 'N-mer concentration'
    if (workingUnits == 'yMonomer')  legend_title <- 'Curve'

    fig <- fig %>% layout(legend = list(title = list(text = legend_title,
                                                     font = list(size = axis_size-1)),
                                        x = xLegend,y = yLegend))
    
    plot_list[[i]] <- fig
    
  }
  
  titleText <- unique(unfolding_exp_data$legend)
  if (!showTitle) titleText <- ''
  
  fig <- plot_list_to_fig(
    paste0('unfolding_fitted_curves_',Sys.Date()),
    plot_list,titleText,
    axis_size,plot_type,plot_width,plot_height)
  
  return( fig )
}

# Requires 
#  - unfolding_exp_data: a df with
#    four columns: 'wavelength', 'value' (CD signal), 'temperature / chem_con', 'legend'
# Optional:
#   - unfolding_fitted_data: a df with the same four named columns
#     that contains the fitted data

# Returns the plot of the dataframe
# We plot the wavelength versus signal, coloured by temperature, one subplot per legend, 

plot_unfolding_exp_spectra <- function(
    unfolding_exp_data, workingUnits,
    plot_width=12, 
    plot_height=8, plot_type='svg', axis_size=12,
    unfolding_fitted_data = NULL,
    plot_mode='markers',
    xLegend=1,yLegend=1,
    showTitle=T) {
  
  # Return null if there is no data
  if (is.null(unfolding_exp_data)) return(NULL)
  
  weHaveFittedData <- !is.null(unfolding_fitted_data)
  
  minWL <- min(unfolding_exp_data$wavelength) - 5
  maxWL <- max(unfolding_exp_data$wavelength) + 5
  
  xaxis <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  yAxisTitle <- workingUnits2ProperLabel(workingUnits)
  
  yaxis <- list(title = yAxisTitle,titlefont = list(size = axis_size),
                tickfont = list(size = axis_size),showgrid = F)
  
  plot_list   <- list()
  
  unfolding_exp_data <- add_measurement_factor_column(unfolding_exp_data)
  
  tot_cond    <- length(unique(unfolding_exp_data$legend))
  
  unfolding_exp_data$color <- get_colors_from_numeric_values(unfolding_exp_data$measurement_factor)
  
  temperatureBased <- any( grepl('temperature',colnames(unfolding_exp_data),  ignore.case = T) )
  chemBased        <- any( grepl('chem_conc',  colnames(unfolding_exp_data),  ignore.case = T) )
  
  i <- 0
  for (leg in unique(unfolding_exp_data$legend)) {
    
    fig <- plot_ly()
    
    i   <- i + 1
    
    df_temp <- unfolding_exp_data[unfolding_exp_data$legend == leg,]
    
    if (weHaveFittedData) {
      df_temp_fit <- unfolding_fitted_data[unfolding_exp_data$legend == leg,]
    }
    
    # Assuming df_temp has a 'group' column indicating different groups
    unique_groups <- unique(df_temp$measurement_factor)
    
    for (group_value in unique_groups) {
      
      subset_df <- df_temp[df_temp$measurement_factor == group_value, ]
      subset_df <- subset_df %>% arrange(wavelength)
      
      if (temperatureBased) {
        current_factor <- round(unique(subset_df$measurement_factor),1)
      } else {
        current_factor <- signif(unique(subset_df$measurement_factor),2)
      }
      name           <- paste0(current_factor)

      fig <- add_trace(
        fig,
        x = subset_df$wavelength,
        y = subset_df$value,
        color = I(unique(subset_df$color)),
        type = "scatter",
        mode = plot_mode,
        name = name,
        showlegend = T
      )
      
      if (weHaveFittedData) {
        subset_df_fit <- df_temp_fit[df_temp$measurement_factor == group_value, ]
        subset_df_fit <- subset_df_fit %>% arrange(wavelength)
        
        # Let's plot one line per replicate
        id_df <- subset_df_fit %>%
          group_by(wavelength) %>%
          mutate(duplicate_id = row_number())
        
        subset_df_fit <- inner_join(subset_df_fit,id_df,)
        
        for (id in unique(subset_df_fit$duplicate_id)) {
          
          subset_df_fit2 <- subset_df_fit[subset_df_fit$duplicate_id == id,]
          
          fig <- add_trace(
            fig,
            x = subset_df_fit2$wavelength,
            y = subset_df_fit2$value,
            type = "scatter",
            mode = "lines",
            line=list(color='black'),
            showlegend = FALSE
          )
          
        }
      }
    }
    
    fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,tot_cond,axis_size)
    
    legend_title <- ifelse(chemBased,"[Denaturant agent] (M)","Temperature (°C)")
    
    fig <- fig %>% layout(legend = list(title = list(text = legend_title,
                                                     font = list(size = axis_size-1)),
                                        x = xLegend,y = yLegend))
    
    plot_list[[i]] <- fig
    
  }
  
  plot_name <- ifelse(weHaveFittedData,
                      paste0('unfolding_fitted_spectra_',Sys.Date()),
                      paste0('unfolding_spectra_',Sys.Date()))
  
  titleTxt <- unique(unfolding_exp_data$legend)
  if (!showTitle) titleTxt <- ''
  
  fig <- plot_list_to_fig(
    plot_name,
    plot_list,titleTxt,
    axis_size,plot_type,plot_width,plot_height)
  
  return( fig )
}

# Requires a dataframe called 'res_df' with the following columns:
# 'wavelength', 'legend', 'residuals' and 'temperature'
plot_residuals <- function(res_df,
                           axis_size=12,
                           svd_or_pca_based=FALSE,
                           xlab="Temperature (ºC)",
                           use_log_axis=FALSE) {
  
  pre_str  <- ifelse(svd_or_pca_based,'k: ','WL: ')
  post_str <- ifelse(svd_or_pca_based,'',' nm')
  
  if (length(unique(res_df$legend))==1) {
    res_df$condition <- paste0(pre_str,res_df$wavelength,post_str)
  } else {
    res_df$condition <- paste0(pre_str,res_df$wavelength,post_str,' Dataset ',res_df$legend)
  }
  
  res_df <- add_measurement_factor_column(res_df)
  
  p <- ggplot(res_df,aes(x=measurement_factor,y=residuals))+
    geom_point(size=2)+
    theme_bw(base_size = axis_size)+
    xlab(xlab)+
    ylab("Residuals")+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  
  if (use_log_axis) {
    p <- p + scale_x_log10()
  }
  
  req_rows <- ceiling(length(unique(res_df$condition))/4)
  
  p2 <- p + facet_wrap(~ condition,
                       ncol = 4,nrow = req_rows)
  
  return(p2)
}

plot_basis_spectra <- function(basis_spectra_df,workingUnits,plot_width=12, 
                               plot_height=8, plot_type='svg', axis_size=12,
                               xLegend=1,yLegend=1,showTitle=T) {
  
  minWL <- min(basis_spectra_df$wavelength) - 5
  maxWL <- max(basis_spectra_df$wavelength) + 5
  
  xaxis <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  yAxisTitle <- workingUnits2ProperLabel(workingUnits)
  
  yaxis <- list(title = yAxisTitle,titlefont = list(size = axis_size),
                tickfont = list(size = axis_size),showgrid = F)
  
  plot_list   <- list()
  
  tot_cond    <- length(unique(basis_spectra_df$legend))
  
  i <- 0
  for (leg in unique(basis_spectra_df$legend)) {
    
    fig <- plot_ly()
    
    i   <- i + 1
    
    df_temp <- basis_spectra_df[basis_spectra_df$legend == leg,]
    df_temp <- df_temp %>% arrange(wavelength)
    
    fig <- add_trace(
      fig,
      x = df_temp$wavelength,
      y = df_temp$value,
      color = df_temp$k,
      mode = 'lines', 
      type = "scatter"
    )

    fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,tot_cond,axis_size)
    
    fig <- fig %>% layout(legend = list(title = list(text = "Basis spectra",
                                                     font = list(size = axis_size-1)),
                                        x = xLegend,y = yLegend))
    
    plot_list[[i]] <- fig
    
  }
  
  titleTxt <- unique(basis_spectra_df$legend)
  if(!showTitle) titleTxt <- '' 
    
  fig <- plot_list_to_fig(
    paste0('basis_spectra_',Sys.Date()),
    plot_list,titleTxt,
    axis_size,plot_type,plot_width,plot_height)
  
  return( fig )
  
}

# Plot the cumulative explained variance dataframe
# Requires a dataframe with 3 columns:
#   'k', 'explained_variance', 'legend'

# Returns the plot of the dataframe, one subplot per unique legend

plot_explained_variance <- function(df,plot_width=12,plot_height=8,
                                    plot_type='svg',axis_size=12,
                                    xLegend=1,yLegend=1,showTitle=T) {
  
  xaxis <- list(title = 'k',titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),showgrid = F)
  
  yaxis <- list(title = 'Cumulative explained variance (%)',
                titlefont = list(size = axis_size),
                tickfont = list(size = axis_size),
                range = c(0, 102),showgrid = F)
  
  plot_list   <- list()
  
  tot_cond    <- length(unique(df$legend))
  
  i <- 0
  for (leg in unique(df$legend)) {
    
    fig <- plot_ly()
    
    i   <- i + 1
    
    df_temp       <- df[df$legend == leg,]

    fig <- add_trace(
      fig,
      x = df_temp$k,
      y = df_temp$explained_variance,
      type = "scatter",
      mode = "markers",
      name = leg
    )
    
    fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,tot_cond,axis_size)

    legend_title <- "Dataset"
    
    fig <- fig %>% layout(legend = list(title = list(text = legend_title,
                                                     font = list(size = axis_size-1)),
                                        x = xLegend,y = yLegend))
    
    plot_list[[i]] <- fig
    
  }
  
  titleTxt <- unique(df$legend)
  if(!showTitle) titleTxt <- '' 
  
  fig <- plot_list_to_fig(
    paste0('explained_variance_',Sys.Date()),
    plot_list,titleTxt,
    axis_size,plot_type,plot_width,plot_height)
  
  return( fig )
}

# Input: List of lists. 
plot_fitted_spectra_sec_str <- function(list_of_signals,list_of_fitted_signals,
                                        max_wl_ref=240,step_wl_ref=1,
                                        plot_width=12  , plot_height=8,
                                        plot_type='svg', axis_size=18
                                        ) {
  
  # Initialize the plot 
  fig         <- plot_ly()
  minWL       <- 190
  
  for (i in 1:length(list_of_signals)) {
    
    exp_signals  <- list_of_signals[[c(i)]]
    exp_fittings <- list_of_fitted_signals[[c(i)]]
    
    for (ii in 1:length(exp_signals)) {
      
      exp_signal  <- exp_signals[[c(ii)]]
      exp_fitting <- exp_fittings[[c(ii)]]
      name        <- names(exp_signals)[c(ii)]
      
      wavelength <- max_wl_ref - (seq(0,length(exp_fitting)-1) * step_wl_ref)
      
      df  <- data.frame('wavelength'=wavelength,
                        'query'=exp_signal,'fitting'=exp_fitting)
      
      minWL <- min(minWL,df$wavelength)
      
      #color=I(colorPalette[counter])
      
      fig <- fig %>% add_trace(data = df,x=~wavelength,y=~fitting,
                               type = 'scatter', mode = 'lines',
                               line = list(color = "black",opacity = 0.6),
                               showlegend = FALSE)
      
      fig <- fig %>% add_trace(data = df,x=~wavelength,y=~query,
                               type = 'scatter', mode = 'markers',
                               name = name)
      
    }
    
  }
  
  minWL <- minWL      - 5
  maxWL <- max_wl_ref + 5
  
  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  y <- list(title = workingUnits2ProperLabel('meanUnitMolarExtinction'),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),showgrid = F)
  
  fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))
  
  fig <- configFig(fig,paste0("CD_spectra_fitted_delta_epsilon_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  
  return(fig)
}

# Requires 
#  - unfolding_exp_data: a df with
#    four columns: 'temperature / chem_conc', 'variable',value', 'legend'
#  - the x-axis label

# Returns the plot of the dataframe
# We plot the value versus temperature (or chemical agent concentration), 
# coloured by variable, one subplot per legend, 

plot_unfolding_fractions <- function(unfolding_fractions,
                               plot_width=12,plot_height=8,
                               plot_type='svg', axis_size=12,
                               xAxisLabel='Set',
                               xLegend=1,yLegend=1,showTitle=T) {
  
  # Return null if there is no data
  if (is.null(unfolding_fractions)) return(NULL)
  
  minX  <- min(unfolding_fractions[,1])
  maxX  <- max(unfolding_fractions[,1])
  delta <- maxX - minX
    
  xaxis <- list(title = xAxisLabel,titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),showgrid = F,
                range = c(minX-delta*0.1, maxX+delta*0.1))
  
  yaxis <- list(title = 'Fraction',titlefont = list(size = axis_size),
                tickfont = list(size = axis_size),showgrid = F)
  
  plot_list   <- list()
  
  tot_cond    <- length(unique(unfolding_fractions$legend))
  
  i <- 0
  for (leg in unique(unfolding_fractions$legend)) {
    
    fig <- plot_ly()
    
    i   <- i + 1
    
    df_temp <- unfolding_fractions[unfolding_fractions$legend == leg,]
    
    unique_groups <- unique(df_temp$variable)
    
    for (group_value in unique_groups) {
      
      subset_df <- df_temp[df_temp$variable == group_value, ]
      
      subset_df <- subset_df[order(subset_df[,1]), ]
      
      fig <- add_trace(
        fig,
        x = subset_df[,1],
        y = subset_df$value,
        type = "scatter",
        mode = "lines",
        name = group_value
      )
    }
    
    fig <- add_layout_to_subplot(fig,xaxis,yaxis,leg,tot_cond,axis_size)
    
    fig <- fig %>% layout(legend = list(
      title = list(text = "",font = list(size = axis_size-1)),
      x = xLegend,y = yLegend))
    
    plot_list[[i]] <- fig
    
  }
  
  titleTxt <- unique(unfolding_fractions$legend)
  if (!showTitle) titleTxt <- ''
  
  fig <- plot_list_to_fig(
    paste0('unfolding_fractions_',Sys.Date()),plot_list,
    titleTxt,
    axis_size,plot_type,plot_width,plot_height)
  
  return( fig )
}

plot_helicity <- function(x,y,color) {
  
  fig <- plot_ly(x = ~x, y = ~y, color = ~color,
                 type = 'scatter', mode = 'markers') %>%
    layout(
      xaxis = list(
        title = 'Temperature °C',
        tickfont = list(size  = 18),
        titlefont = list(size = 18)  # Set font size for x-axis label
      ),
      yaxis = list(
        title = 'Helicity fraction',
        range = c(0, 1),
        tickfont = list(size  = 18),
        titlefont = list(size = 18)  # Set font size for y-axis label
      )
    )
  
  return(fig)
  
}

plot_heatmap_bayes_posterior <- function(df,i=1,j=2,plot_width=12,plot_height=12,plot_type='svg',plot_axis_size=18) {

    if (i+1 > ncol(df) || j+1 > ncol(df)) {
      return(NULL)
    }

    df <- df[,c(i,j,ncol(df))]

    ss_classes <- colnames(df)

    df_grouped <- df %>%
      group_by(across(1:2)) %>%
      summarise(value = sum(across(everything()), na.rm = TRUE))

    # Normalise bins to probabilities
    df_grouped[,3] <- df_grouped[,3] / sum(df_grouped[,3]) * 100

    colnames(df_grouped) <- c('x','y','Percentage')

    bu_gn_palette <- c("#F7FCFD", "#E5F5F9", "#CCECE6", "#99D8C9", "#66C2A4", "#41AE76", "#238B45", "#006D2C", "#00441B")
    color_ramp    <- colorRamp(bu_gn_palette)

    df_grouped <- df_grouped[df_grouped$Percentage > 0.5,]

    df_grouped <- na.omit(df_grouped)

    fig <- plot_ly(
      data = df_grouped,
      x = ~x,
      y = ~y,
      z = ~Percentage,
      type = "heatmap",
      colors = color_ramp
    ) %>%
      layout(
        font="Roboto",
        xaxis = list(title = ss_classes[1], automargin = TRUE, showgrid = FALSE, range = c(0, 1),
        tickfont = list(size = plot_axis_size), titlefont = list(size = plot_axis_size)),
        yaxis = list(title = ss_classes[2], automargin = TRUE, showgrid = FALSE, range = c(0, 1),
        tickfont = list(size = plot_axis_size), titlefont = list(size = plot_axis_size)),
        legend = list(font = list(size = plot_axis_size-1)),
        margin = list(l = 60, r = 50, b = 60, t = 50)
      )

      #fig <- configFig(fig,paste0("sesca_posterior_probability_proj_",strsplit(as.character(Sys.time())," ")[[1]][1]),
     #              plot_type,plot_width,plot_height)

    return(fig)

}

# Input: 
# 'wavelength'    : vector    (float),     length n
# 'signal_matrix' : 2D matrix (float),     size n*m
# 'sample_names'  : vector    (character), size m
plot_gQuadruplexReferences <- function(wavelength,signal_matrix,sample_names,
                                       axis_size=18,plot_type='png',
                                       plot_width=12,plot_height=12) {
  
  fig          <- plot_ly()
  colorPalette <- getPalette(length(sample_names))
  
  for (i in 1:ncol(signal_matrix)) {
    
    signal  <- signal_matrix[,i]
    
    df      <- data.frame('wavelength'=wavelength,signal)
    
    fig <- fig %>% add_trace(data=df,color=I(colorPalette[i]),x=~wavelength,y=~signal,
                             type = 'scatter', mode = 'markers+lines',
                             marker = list(size=4),
                             line = list(width = 1.8),
                             name = sample_names[i])
    
  }
  
  minWL <- min(wavelength) - 5
  maxWL <- max(wavelength) + 5
  
  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(minWL, maxWL),showgrid = F)
  
  y <- list(title = workingUnits2ProperLabel('molarExtinction'),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),showgrid = F)
  
  fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))
  
  fig <- configFig(fig,paste0("CDspectra_G-Quadruplex_references_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  return(fig)  
  
}

#---------------------------------------------------------------------------------------------------------------------------------------
# Principal Cluster Analysis plotting function, provided by Rafael Del Villar-Guerra and minimally edited by O.B.
#---------------------------------------------------------------------------------------------------------------------------------------
plot_pca_analysis <- function( X,wavelength,spectraNames, text_graph_title = NULL, slider_graph_label_size = 5, slider_axis_label_size = 14,
                               slider_point_size = 2, slider_arrow_width = 1, slider_ellipse_confidence = 0.85, checkbox_show_ellipse = TRUE,
                               checkbox_show_legend = FALSE, slider_pca_num_clusters = 3, slider_pca_num_components = 8, PC1 = 1, PC2 = 2 ) {
  
  rownames(X) <- spectraNames
  colnames(X) <- wavelength
  
  axes = c( as.numeric( PC1 ), as.numeric( PC2 ) )
  if( checkbox_show_legend == FALSE ) { show_legend = "none" } else { show_legend = "right" }
  
  # Principal Component Analysis (PCA)
  result.pca  <- FactoMineR::PCA( X, scale.unit = TRUE, ncp = slider_pca_num_components, graph = FALSE )
  # Hierarchical Clustering on Principle Components (HCPC)
  result.hcpc <- FactoMineR::HCPC( result.pca, nb.clust = slider_pca_num_clusters, consol = TRUE, 
                                   min = 3, max = 10, graph = FALSE )
  
  nClust <- n_distinct(result.hcpc$data.clust$clust)

  wavelength_labels <- list()
  for( n in 1:nClust ) {
    pasted_string = paste0( "result.hcpc$desc.var$quanti$`", n, "`" )
    v.test_data <- eval( parse( text = pasted_string ) )
    max_value   <- abs( v.test_data[1,1] )
    min_value   <- abs( v.test_data[nrow(v.test_data),1] )
    
    v.test_data <- as.data.frame( v.test_data )

    if( max_value > min_value ) {
      wavelength_labels$list <- c( ( wavelength_labels$list ), rownames( v.test_data[1,] ) )
    } else {
      wavelength_labels$list <- c( ( wavelength_labels$list ), rownames( v.test_data[nrow(v.test_data),] ) )
    }
  }
  
  cluster = result.hcpc$data.clust[,ncol(result.hcpc$data.clust), drop = FALSE]
  
  colnames( cluster ) <- c( "Cluster" )
  df.cluster.bind     <- cbind.data.frame( cluster, X )
  
  result.pca <- FactoMineR::PCA( df.cluster.bind, scale.unit = TRUE, ncp = slider_pca_num_components, graph = FALSE, quali.sup = 1 )
  
  label.biplot.axis.x <- paste ("PC", axes[1], " (", round(result.pca$eig[PC1,2], digits = 1)," %", ")", sep="")
  label.biplot.axis.y <- paste ("PC", axes[2], " (", round(result.pca$eig[PC2,2], digits = 1)," %", ")", sep="")
  textElement_labels = element_text( face = "bold", color = "black", size = slider_axis_label_size )
  
  pca_plot <- factoextra::fviz_pca_biplot(
    result.pca,
    axes = axes,
    col.var = "blue",
    repel = TRUE,
    habillage = 1,
    addEllipses = checkbox_show_ellipse,
    labelsize = slider_graph_label_size,
    pointsize = slider_point_size,
    arrowsize = slider_arrow_width,
    ellipse.level = slider_ellipse_confidence,
    select.var =  list( name = wavelength_labels$list )
  ) +
    labs( title = text_graph_title, x = label.biplot.axis.x, y = label.biplot.axis.y ) +
    theme( legend.position = show_legend ) +
    theme( text = textElement_labels )
  
  return(pca_plot)
  
}

#---------------------------------------------------------------------------------------------------------------------------------------
# Principal Cluster Analysis plotting function, provided by Rafael Del Villar-Guerra and minimally edited by O.B.
#---------------------------------------------------------------------------------------------------------------------------------------
plot_cluster_analysis <- function( X, spectraNames,cluster_num_clusters = 3, cluster_graph_title = NULL, 
                                   cluster_text_size = 1, cluster_num_components = 8, axis_text_size = 14 ) {
  
  rownames(X) <- spectraNames
  
  result.pca  <- FactoMineR::PCA( X, scale.unit = TRUE, ncp = cluster_num_components, graph = FALSE )
  result.hcpc <- FactoMineR::HCPC( result.pca, nb.clust = cluster_num_clusters, consol = TRUE, min = 3, max = 10, graph = FALSE )
  
  tree_height <- max(result.hcpc$call$t$tree$height)
  upper_lim   <- ceiling(tree_height / 10) * 10
  breaks      <- seq(0,upper_lim,10)
  
  cluster_plot <- fviz_dend( result.hcpc,  k = cluster_num_clusters, cex = cluster_text_size, 
                              main = cluster_graph_title) +
    theme(axis.text.y  = element_text(size = axis_text_size),
          axis.title.y = element_text(size = axis_text_size+1))+
    scale_y_continuous(limits = c(-upper_lim/3,upper_lim),breaks = breaks)

  cluster_plot
}


