configPlot <- function(fig,plot_type,plot_width,plot_height,filename) {
  
  fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = filename,
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE,
    modeBarButtonsToRemove = list(
      'sendDataToCloud',
      'hoverClosestCartesian',
      'hoverCompareCartesian',
      'lasso2d',
      'select2d'))
  
}

## Plot fluorescence signal versus time (color by ligand concentration)
## fluo_m is a 3 column dataframe: time, fluo and conc. 
## which is the signal 
## Returns the plot

plot_fluo_signal <- function(fluo_m,cold_min,cold_max,hot_min,hot_max,
                             which,plot_width, plot_height, plot_type, axis_size) {
  
  x <- list(title = "Time (sec)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = which,titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  delta <- max(fluo_m$fluo) - min(fluo_m$fluo)
  min_y <- min(fluo_m$fluo) - 0.1*delta
  max_y <- max(fluo_m$fluo) + 0.2*delta
  
  fluo_m0 <- fluo_m %>% filter(conc == 0)
  fluo_m  <- fluo_m %>% filter(conc > 0)
  
  fig <- plot_ly()
  
  fig <- fig %>% 
    add_polygons(x=c(cold_min, cold_min, cold_max, cold_max),
                   y=c(min_y, max_y, max_y, min_y),
                   line=list(width=0),
                   fillcolor="#8bb8e8", inherit = FALSE,
                   name = 'Cold Region',opacity = 0.5) %>% 
    add_polygons(x=c(hot_min, hot_min, hot_max, hot_max),
                 y=c(min_y, max_y, max_y, min_y),
                 line=list(width=0),
                 fillcolor="#e40046", inherit = FALSE,
                 name = 'Hot Region',opacity = 0.3)
  
  if (nrow(fluo_m) == 0) { 
    
    fig <- fig %>% add_trace(data = fluo_m0 ,name = ~conc_, 
                             color = I("#389196"),y    = ~fluo   , x   = ~time) 
    
    # if we ONLY have positions with no ligand return this plot
    return(fig) 
    }
  
  for (position in unique(fluo_m0$conc_)) {
    fluo_m0_time <- fluo_m0 %>% filter(conc_ == position) 

    # set as black the positions with no ligand
    fig <- fig %>% add_trace(data  = fluo_m0_time,
                             color = I("black"),
                             y     = ~fluo, x=~time,name=position,
                             type = 'scatter', mode = 'lines', showlegend = FALSE) 
  }
  
  conc_vector <- unique(fluo_m$conc)
  n_colors    <- sapply(conc_vector,get_color_from_conc,min(fluo_m$conc),max(fluo_m$conc))
  
  for (position in unique(fluo_m$conc_)) {
    
    fluo_m_time <- fluo_m %>% filter(conc_== position) 
    idx         <- which(conc_vector == unique(fluo_m_time$conc))
    
    fig <- fig %>% add_trace(data  = fluo_m_time,
                             color = I(n_colors[idx]),
                             y = ~fluo, x=~time,name=position,
                             type = 'scatter', mode = 'lines', showlegend = FALSE) 
  }
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend = FALSE)
  
  fig <-  configPlot(fig,plot_type,plot_width,plot_height,
  paste0('mstTraces_',Sys.Date()))
  
  return(fig)
}


plot_lig_fluo <- function(f_cold,conc,expID_vector,
                          plot_width, plot_height, plot_type,
                          axis_text_size,legend_text_size) {
  
  df <- get_df_min_max_fluo(f_cold,conc,expID_vector)

  fig <- plot_ly(df, x = ~x, y = ~y, type = 'bar', color = ~z)
  fig <- fig %>% layout(title = "",
                        xaxis = list(title = "",tickfont = list(size = axis_text_size)),
                        yaxis = list(title = "",tickfont = list(size = axis_text_size)))
  
  fig <- fig %>%  layout(legend = list(font = list(size = legend_text_size)))
  
  fig <-  configPlot(fig,plot_type,plot_width,plot_height,
                     paste0('Fluorescence_at_the_max_and_min_ligand_concentration_',Sys.Date()))
  
  return(fig)
  
}

plot_fluo_vs_lig <- function(signal,conc,expID_vector,which,
                             plot_width, plot_height, plot_type, axis_size,
                             legend_text_size,
                             logScaleType,
                             showSD=T) {
  
  yAxisText <- "Fnorm"
  
  if (which == "cold")         yAxisText <-  "Fluorescence (AU)"
  if (which == "simulation")   yAxisText <-  "Signal (AU)"
  if (which == "simulationFO") yAxisText <-  "Fraction of occupied binding sites (%)"
  
  conc  <- conc * unlist(factorList[logScaleType])

  titleList <- list("defaultPlotly" = "Binding Partner Concentration (M)"   , 
                    "molar"         = "Binding Partner Concentration (M)"   ,
                    "milimolar"     = "Binding Partner Concentration (mM)"  , 
                    "micromolar"    = "Binding Partner Concentration (μM)")
  
  xTitle <- unlist(titleList[logScaleType])
  
  if (logScaleType == "defaultPlotly"){
    exponentformat <- "B"
  } else {
    exponentformat <- 'power'
  }
  
  x <- list(title = xTitle,titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat,
            showgrid = F,linecolor='black')
  
  y <- list(title = yAxisText,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),showgrid = F,linecolor='black')
  
  df <- data.frame(signal,conc,expID_vector)
  
  if (showSD) {
    
    df <- df %>% group_by(conc,expID_vector) %>% 
      summarise(err=sd(signal),signal = mean(signal))
  }

  showlegend <- F
  if (length(unique(expID_vector)) > 1) showlegend <- T
  
  fig <- plot_ly(data = df, x = ~conc, y = ~signal,
                 color = ~expID_vector,
                 marker = list(size = 9,
                               line = list(width = 0)),
                 showlegend = showlegend,type = 'scatter')
  
  # add error bars
  if (showSD) {
    
    # Remove conditions where it was not possible to estimate the sd
    temp_df <- df[!(is.na(df$err)),]
    
    fig <- fig %>% add_trace(data = temp_df,x = ~conc, y = ~signal,
                             error_y =  list(array = ~err),
                             opacity=0.75,showlegend = F)
    } 
  
  m   <- list(l = 10,r = 10,b = 20,t = 40,pad = 4)
  fig <- fig %>%  layout(xaxis = x, yaxis = y,margin=m)
  fig <- fig %>%  layout(legend = list(font = list(size = legend_text_size)))
  
  fig <-  configPlot(fig,plot_type,plot_width,plot_height,
                     paste0('mstSignalVsLigand_',Sys.Date()))
  
  return(fig)
  
}

formatAxis <- function(title,axis_size,type="-",exponentformat="B"){
  
  axis <- list(title = title,titlefont = list(size = axis_size), 
       tickfont = list(size = axis_size),type = type,
       exponentformat=exponentformat,showgrid = F,linecolor='black')
  
  return(axis)
}

# Force the signal to start at 0. Whenever possible use the value of the parameter RF1
# If not, use the mean signal of the first 3 points
# df is a data frame of 4 columns: signal, conc (ligand concentration), expID_vector (experiment ID), protConc (protein concetration)
# tidy_fit_list is a list with the fitted parameters
# uniqueExpIDs contains the experiment ids only of the fitted experiments (not selected as control)
set_initial_signal_to_zero <- function(df,tidy_fit_list,uniqueExpIDs) {
  
  dfs <- list()
  ids <- unique(df$expID_vector)
  
  c <- 0
  for (i in 1:length(ids)) {
    dfTemp <- df[df$expID_vector == ids[i],]
    
   
    proteinConcentration <- dfTemp$protConc
    
    if (ids[i] %in% uniqueExpIDs) {
      
      c <- c + 1

      fittedRF1        <- as.numeric(tidy_fit_list[[c]][1,2])
      valuesToSubtract <- fittedRF1 * proteinConcentration * 1e6

      dfTemp$signal   <- dfTemp$signal - valuesToSubtract  
      
    } else {
      
      dfTemp <- dfTemp %>% arrange(conc)
      dfTemp$signal   <- dfTemp$signal - mean(dfTemp$signal[1:min(3,nrow(dfTemp))])  
      
    }
    dfs[[i]] <- dfTemp
  }
  
  df <- do.call(rbind,dfs)
  
  return( df )
  
}

plot_fluo_vs_ligFit <- function(signal,conc,expID_vector,protConc,
                                plot_width, plot_height, plot_type,
                                legend_text_size,logScaleType,
                                tidy_fit_list,uniqueExpIDs,
                                legends,colorPalette,sels,ids,
                                showSD=T,scaleToZero=F) {
  
  conc  <- log10(conc * unlist(factorList[logScaleType]))
  df    <- data.frame(signal,conc,expID_vector,protConc)
  
  if (showSD) {
    
    df <- df %>% group_by(conc,expID_vector,protConc) %>% 
      summarise(err=sd(signal),signal = mean(signal))
  }
  
  if (scaleToZero) {
    df <- set_initial_signal_to_zero(df,tidy_fit_list,uniqueExpIDs) 
  }
  
  showlegend <- length(unique(expID_vector)) > 1
  
  fig <- plot_ly()
  
  for (var in unique(df$expID_vector)) {
    
    counterColorPalette <- which(ids == var)
    
    if (sels[counterColorPalette]) {
      hexColor            <- colorPalette[counterColorPalette]
      
      tempDf <- df[df$expID_vector == var,]
      fig    <- fig %>% add_trace(data=tempDf,x=~conc,y=~signal,
                                  color=I(hexColor),
                                  type = 'scatter', 
                                  name=legends[counterColorPalette],
                                  showlegend = showlegend,
                                  marker = list(size = 9,
                                                line = list(width = 0)))
      
      # add error bars
      if (showSD) {
        
        # Remove conditions where it was not possible to estimate the sd
        tempDf <- tempDf[!(is.na(tempDf$err)),]
        
        fig <- fig %>% add_trace(data = tempDf,x = ~conc, y = ~signal,
                                 error_y =  list(array = ~err),
                                 color=I(hexColor),
                                 opacity=0.75,showlegend = F)
      } 
      
    }
    
  }
  
  m   <- list(l = 10,r = 10,b = 20,t = 40,pad = 4)
  
  fig <- fig %>%  layout(legend = list(font = list(size = legend_text_size),
                                       margin=m))
  
  fig <-  configPlot(fig,plot_type,plot_width,plot_height,
                     paste0('mstSignalVsLigand_Fitting_',Sys.Date()))
  
  return(fig)
  
}

add_axis_to_plot <- function(fig,axis_size,
                             logScaleType,xAxisText,yAxisText,
                             xAxis_nticks,yAxis_nticks,
                             xAxis_lowerLimit,xAxis_upperLimit,
                             yAxis_lowerLimit,yAxis_upperLimit) {
  
  titleList <- list("molar"         = "(M)"   ,
                    "milimolar"     = "(mM)"  , 
                    "micromolar"    = "(μM)")
  
  xTitle         <- paste(xAxisText,unlist(titleList[logScaleType]))
  
  #x <- formatAxis(xTitle,axis_size,"log",'power')
  x <- formatAxis(xTitle,axis_size)
  y <- formatAxis(yAxisText,axis_size)
  
  x[['range']] <- list(log10(xAxis_lowerLimit),log10(xAxis_upperLimit))
  y[['range']] <- as.list(extendrange(c(yAxis_lowerLimit,yAxis_upperLimit)),0.02)
  
  # Select how many decimals to show in the plot
  diff <- abs(yAxis_upperLimit-yAxis_lowerLimit)
  if (diff > 1) {
    signifPlaces <- 0
  } else if (diff > 0.01) {
    signifPlaces <- 2
  } else {
    signifPlaces <- 4
  }
  
  y[['dtick']] <- signif(abs(yAxis_upperLimit-yAxis_lowerLimit) / (yAxis_nticks-1),signifPlaces)
  y[['tick0']] <- signif(yAxis_lowerLimit,signifPlaces)
  y[['zeroline']]    <- FALSE
  
  # max and min value of the x-axis (in logarithm base 10 scale)
  minX <- ceiling(log10(xAxis_lowerLimit))
  maxX <- floor(log10(xAxis_upperLimit))
  
  x[['tickmode']]    <- "array"
  
  tickvals           <- seq(minX,maxX,abs(maxX-minX) / (xAxis_nticks-1))
  x[['tickvals']]    <- as.list(tickvals)
  
  x[['ticktext']]    <- as.list(paste0("10","<sup>",signif(tickvals,2),"</sup>"))
  x[['zeroline']]    <- FALSE
  
  fig <- fig %>%  layout(xaxis = x, yaxis = y)
  
  return(fig)
  
}

plot_add_df <- function(fig,signal,conc,color,name) {
  
  df <- data.frame(signal,conc)
  fig <- fig %>% add_trace(data=df,y = ~signal,x=~conc,
                           name = name,mode = 'lines',type = "scatter",
                           inherit=FALSE,color=I(color))
  
  return(fig)
}

format_kd_for_plot_title <- function(fit_obj_list,tidy_fit_list,ci95,uniqueExpIDs,model_has_2_kds) {
  
  # Create title with fitted Kds
  titleXs <- ""
  for (iter in 1:length(tidy_fit_list)) {
    
    tidy_fit <- format_fitting_table_kd(tidy_fit_list[[iter]],fit_obj_list[[iter]])
    titleX   <- paste0(tidy_fit$Term[3]," : ",specify_decimal(tidy_fit$Estimate[3],1))
    
    if (!(model_has_2_kds)) ( titleX    <- paste0(titleX," Asymmetric CI95 : ",ci95[[iter]]))
    
    if (model_has_2_kds) {
      titleX2    <- paste0(tidy_fit$Term[4]," : ",specify_decimal(tidy_fit$Estimate[4],1))
      titleX     <- paste0(titleX,"  ",titleX2)
    }
    
    titleX  <- paste0(uniqueExpIDs[[iter]]," ",titleX)
    
    if (iter > 1) {
      titleXs <- paste0(titleXs,"\n",titleX)
    } else {
      titleXs <- paste0(titleXs,titleX)
    }
  }
  
  return(titleXs)
}

# Add the fitted line to the plot with the experimental data
plot_add_fitting <- function(fig,fit_obj_list,tidy_fit_list,uniqueExpIDs,
                             pconcs_list,
                             title,ci95,allExpIDs,
                             model_has_2_kds,axis_size,logScaleType,
                             legends,colorPalette,sels,ids,
                             kdInTitle,scaleToZero) {
  
  df <- get_fitting_plot(uniqueExpIDs,fit_obj_list,pconcs_list)

  df$Conc     <- df$Conc / 1e6                            # Back to molar from micromolar (used in the fitting process)

  df$Conc  <- log10(df$Conc  * unlist(factorList[logScaleType]))
  
  df <- df %>% rename(conc=Conc,expID_vector=expID)
  
  if (scaleToZero) {
    #TO DO - add protein concentration info
    df <- set_initial_signal_to_zero(df,tidy_fit_list,uniqueExpIDs) 
  }
  
  for (var in unique(allExpIDs)) {
    
    counterColorPalette <- which(ids == var)
    
    # Check we want to plot experiment and that it was fitted
    if ( sels[counterColorPalette] & (var %in% uniqueExpIDs) ) {
      hexColor            <- colorPalette[counterColorPalette]
      
      tempDf <- df[df$expID_vector == var,]
      fig    <- fig %>% add_trace(data=tempDf,x=~conc,y=~signal,
                                  color=I(hexColor),
                                  type = 'scatter', 
                                  mode = 'lines',
                                  name=legends[counterColorPalette],
                                  showlegend = FALSE,
                                  inherit=FALSE,
                                  line = list(width = 2))
    }
  }
  
  titleXs  <- ""
  if (kdInTitle) titleXs <- format_kd_for_plot_title(fit_obj_list,tidy_fit_list,ci95,uniqueExpIDs,model_has_2_kds)

  title <- ifelse(title != "", paste0(title,"\n",titleXs), titleXs)
  
  t   <- list(family = "Roboto",size = axis_size)
  fig <- fig %>% layout(title=list(text = title,font = t))
  
  return(fig)
}

plot_kd_advanced_fit <- function(fluo_fit_data_advance_nice,times,
                                 plot_width, plot_height, plot_type, axis_size,legend_text_size) {
  
  idx       <- grepl("Kd",fluo_fit_data_advance_nice$Term)
  
  kds       <- fluo_fit_data_advance_nice[idx,2]
  std_error <- fluo_fit_data_advance_nice[idx,3]
  low_CI95  <- fluo_fit_data_advance_nice[idx,4]
  
  expIds    <- fluo_fit_data_advance_nice[idx,7]
  
  acceptable_error <- kds > (std_error*2) & low_CI95 > 0
  
  kds       <- kds[acceptable_error]
  std_error <- std_error[acceptable_error]
  low_CI95  <- low_CI95[acceptable_error]
  upp_CI95  <- fluo_fit_data_advance_nice[idx,5][acceptable_error]
  
  nExps <- length(unique(expIds))
  if (nExps > 1) {
    for (i in 2:nExps) {
      times     <- c(times,times+0.1)
    }
  }

  times     <- times[acceptable_error]
  expIds    <- expIds[acceptable_error]
  
  # All to micromolar units
  factor_nM <- grepl("nM",fluo_fit_data_advance_nice$Term[idx][acceptable_error])
  factor_mM <- grepl("mM",fluo_fit_data_advance_nice$Term[idx][acceptable_error])
  factor_uM <- grepl("µM",fluo_fit_data_advance_nice$Term[idx][acceptable_error])
  
  kds       <-  kds * (factor_nM/1e3 + factor_mM*1e3 + factor_uM)
  std_error <-  std_error * (factor_nM/1e3 + factor_mM*1e3 + factor_uM)
  
  kd_units <- "µM"
  if (min(kds) < 1) {
    kds <- kds*1e3
    std_error <- std_error*1e3
    kd_units <- "nM"
  } 
  
  if (min(kds) > 1000) {
    kds <- kds / 1e3
    std_error <- std_error*1e3
    kd_units <- "mM"
  } 
  
  df <- data.frame(kds,times,std_error,expIds)

  y_title <- paste0("Kd ","(",kd_units,")")
  
  x <- list(title = "Time (sec)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = y_title,     titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fig <- plot_ly(data = df, x = ~times, y = ~kds,color=~expIds,
                 type = 'scatter', mode = 'markers',
                 marker = list(symbol = "square",size=8),
                 error_y = ~list(array = std_error,
                                 color = '#000000',opacity=0.7))
  
  fig <- fig %>%  layout(xaxis = x, yaxis = y,showlegend = TRUE)
  fig <- fig %>%  layout(legend = list(font = list(size = legend_text_size)))
  
  fig <-  configPlot(fig,plot_type,plot_width,plot_height,
                     paste0('mstKdEstimation_vs_FhotSelection_',Sys.Date()))
  
  return(fig)
  
}


