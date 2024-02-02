configFig <- function(fig,plot_width, plot_height, plot_type,plot_title_for_download) {
  
  fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = plot_title_for_download,
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

get_log_axis <- function(logScaleType,axis_size) {
  
  titleList <- list("defaultPlotly"   = "Time (s)"   , 
                    "seconds"         = "Time (s)"   ,
                    "miliseconds"     = "Time (ms)"  , 
                    "microseconds"    = "Time (μs)")
  
  xTitle <- unlist(titleList[logScaleType])
  
  if (logScaleType == "defaultPlotly"){
    exponentformat <- "B"
  } else {
    exponentformat <- 'power'
  }
  
  x <- list(title = xTitle,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat,
            showgrid = F,linecolor='black')
  
  y <- list(title = "Autocorrelation",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),showgrid = F,linecolor='black')
  
  return(list("xaxis"=x,"yaxis"=y))
  
}

plotAutocorrelation <- function(autocorrelationData,plot_width, plot_height, plot_type, axis_size,
                                logScaleType,splitFactor,expNames) {
  
  # autocorrelationData has at least 3 columns - 'variable', 'value', 'time'
  # 'value' versus 'time' represents the autocorrelation curve
  # 'variable' allows using different colours 
  
  # Check we have data to plot
  if (nrow(autocorrelationData) == 0) return(NULL)

  fig <- plot_ly()

  factorList <- list("defaultPlotly"=1,"seconds"=1,"miliseconds"=1e3,"microseconds"=1e6)
  
  autocorrelationData$time  <- autocorrelationData$time * unlist(factorList[logScaleType])
  
  autocorrelationData <- add_factor_column(autocorrelationData,splitFactor,expNames)
  
  factors  <- unique(autocorrelationData$factor)
  nb.cols  <- length(factors)
  mycolors <- getPalette(nb.cols)
  
  for (f in factors) {
    
    tempDF <- autocorrelationData[autocorrelationData$factor == f,]
    color  <- mycolors[which(factors == f)]
      
    fig <- fig %>% add_trace(data = tempDF ,name = ~variable, 
                             color = I(color),y    = ~value   , x   = ~time,
                             type = 'scatter', mode = 'markers',name=f)
    
  }
  
  axis <- get_log_axis(logScaleType,axis_size)
  
  fig  <- fig %>% layout(xaxis = axis$xaxis, yaxis = axis$yaxis,showlegend = TRUE,
                        font="Roboto",legend = list(font = list(size = axis_size)))
  
  fig <- configFig(fig,plot_width, plot_height, plot_type,paste0("AutocorrelationPlot_",Sys.Date()))
  
  return(fig)
}

plotFittedAutocorrelation1 <- function(autocorrelationData,predictedAutocorrelation,
                                      plot_width, plot_height, 
                                      plot_type, axis_size,
                                      logScaleType,splitFactor,expNames) {
  
  # autocorrelationData and predictedAutocorrelation have 3 columns - 'variable', 'value', 'time'
  # 'value' versus 'time' represents the autocorrelation curve
  # 'variable' allows using different colours 
  # predictedAutocorrelation is a named list with key-value pairs:
  #                                               fitted curve id - fitted data points
  
  fig <- plot_ly()
  
  factorList <- list("defaultPlotly"=1,"seconds"=1,"miliseconds"=1e3,"microseconds"=1e6)
  
  autocorrelationData$time       <- autocorrelationData$time * unlist(factorList[logScaleType])
  predictedAutocorrelation$time  <- predictedAutocorrelation$time * unlist(factorList[logScaleType])
  
  fig <- fig %>% add_trace(data = predictedAutocorrelation ,name = ~variable, 
                           color = I('black'),y    = ~value   , x   = ~time,
                           type = 'scatter', mode = 'lines',
                           showlegend = F)
  
  autocorrelationData <- add_factor_column(autocorrelationData,splitFactor,expNames)
  
  factors  <- unique(autocorrelationData$factor)
  nb.cols  <- length(factors)
  mycolors <- getPalette(nb.cols)
  
  for (f in factors) {
    
    tempDF <- autocorrelationData[autocorrelationData$factor == f,]
    color  <- mycolors[which(factors == f)]
    
    fig <- fig %>% add_trace(data = tempDF ,name = ~variable, 
                             color = I(color),y    = ~value   , x   = ~time,
                             type = 'scatter', mode = 'markers',name=f)
    
  }

  axis <- get_log_axis(logScaleType,axis_size)
  
  fig  <- fig %>% layout(xaxis = axis$xaxis, yaxis = axis$yaxis,
                         font="Roboto",legend = list(font = list(size = axis_size)))
  
  fig <- configFig(fig,plot_width, plot_height, plot_type,paste0("AutocorrelationPlotFitted_",Sys.Date()))
  
  return(fig)
}

formatHrDistributionPlot <- function(fig,axis_size,useDiffCoeff=F) {
  
  xLabel <- ifelse(useDiffCoeff,"Diffusion Coefficient (m^2/s)",
                   'Hydrodynamic radius (nm)')
  
  fig <- fig  +
    theme_bw(base_size = axis_size) +
    theme(legend.position = "none") +
    xlab(xLabel)
  
  return(fig + ylab("Intensity weighted contribution (%)"))
  
}

plothydrodynamicRadiusDistributionCollapsed <- function(
    estimatedContributions,plot_width, 
    plot_height,plot_type,axis_size,
    peakFilter=T) {
  
  estimatedContributions$value <- estimatedContributions$value*100
  
  if (peakFilter) {
    estimatedContributions$value[estimatedContributions$value < 1] <- 0
  }
  
  # Detect if we should use the diff coefficients
  useDiffCoeff <- 'diff' %in% colnames(estimatedContributions)
  if (useDiffCoeff) {
    estimatedContributions$hr <- estimatedContributions$diff
  }

  ids <- unique(estimatedContributions$variable)
  
  nb.cols  <- length(ids)
  mycolors <- getPalette(nb.cols)
  
  fig <- plot_ly(type = 'scatter', mode = 'lines')
  
  for (id in ids) {
    color <- mycolors[which(id == ids)]
    tempDF <- estimatedContributions[estimatedContributions$variable == id,]
    fig <- fig %>% add_trace(data=tempDF,x = ~hr, y = ~value, 
                             name = id,color=I(color))
  }

  xLabel <- ifelse(useDiffCoeff,"Diffusion Coefficient (m^2/s)",
                   'Hydrodynamic radius (nm)')
  
  x <- list(title = xLabel,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),type = "log",exponentformat='power',
            showgrid = F,linecolor='black')
  
  yAxisTitle <- "Intensity weighted contribution (%)"
  
  y <- list(title = yAxisTitle,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),showgrid = F,linecolor='black')
  
  fig  <- fig %>% layout(xaxis = x, yaxis = y,
                         font="Roboto",legend = list(font = list(size = axis_size)))
  
  title <- "HrDistributionIntensity_"
  
  fig  <- configFig(fig,plot_width, plot_height, plot_type,paste0(title,Sys.Date()))
  
  return(fig)

}

plothydrodynamicRadiusDistribution <- function(
    estimatedContributions,overlapFactor) {
  
  ggplot(estimatedContributions, aes(hr, variable, height = value, group = variable,
                                     color=variable)) + 
    geom_ridgeline(scale = overlapFactor,alpha=0)+
    scale_color_viridis_d()

}

add_peak_areas_to_plot <- function(plot,leftBounds,rightBounds) {
  
  # We will detect peak present in different regions 
  # defined in the vectors leftBounds and rightBounds
  
  palette <- brewer.pal(n=8,"Dark2")
  # Repeat palette if we have more than 12 peak areas to plot
  palette <- rep(palette,ceiling(length(leftBounds)/8)) 
  
  for (i in 1:length(leftBounds)) {
    plot <- plot + geom_vline(xintercept = rightBounds[i],color=palette[i])
    plot <- plot + geom_vline(xintercept = leftBounds[i],color=palette[i])
  }
  
  return(plot)
}

plothydrodynamicRadiusDistributionHistogram <- function(
    estimatedContributions,overlapFactor) {
  
  # Detect if the 'hr' column are the diffusion coefficients
  useDiffCoeff <- min(estimatedContributions$hr) < 1e-10
  
  # Plot the estimated hydrodynamic radius using histograms
  # I couldn't find a way to give the exact bin values for ggridges
  # so we divide the   
  logScale <- estimatedContributions$hr[2] / estimatedContributions$hr[1] 
  
  limit1LogHr <- log(min(estimatedContributions$hr),base=logScale) - 0.5
  limit2LogHr <- log(max(estimatedContributions$hr),base=logScale) + 0.5
  
  estimatedContributions <- estimatedContributions %>% filter(value  > 0.01)
  
  # Used to build the rectangles of the custom-built histogram, 
  # one rectangle for each contribution > 0.001
  estimatedContributionsPre1  <- estimatedContributions # Point (x1,0) x-axis is the Hr
  estimatedContributionsPre2  <- estimatedContributions # Point (x1,y) y-axis is the % contribution
  estimatedContributionsPost1 <- estimatedContributions # Point (x2,y)
  estimatedContributionsPost2 <- estimatedContributions # Point (x2,0)
  
  nrows <- 1:nrow(estimatedContributions)
  
  estimatedContributionsPre1$intercalateNumber  <- nrows
  estimatedContributionsPre2$intercalateNumber  <- nrows
  estimatedContributionsPost1$intercalateNumber <- nrows
  estimatedContributionsPost2$intercalateNumber <- nrows
  
  estimatedContributions$logHr <- log(estimatedContributions$hr,base = logScale)
  
  # Detect if we used the hr (and not the diffusion coefficients!)
  if (useDiffCoeff) {
    # Get the x values of the rectangles (hr) - the hr is evenly spaced in log space
    estimatedContributionsPre1$logHr  <- estimatedContributions$logHr + 0.5 
    estimatedContributionsPre2$logHr  <- estimatedContributions$logHr + 0.5 
    estimatedContributionsPost1$logHr <- estimatedContributions$logHr - 0.5 
    estimatedContributionsPost2$logHr <- estimatedContributions$logHr - 0.5 
  } else {
    # Get the x values of the rectangles (hr) - the hr is evenly spaced in log space
    estimatedContributionsPre1$logHr  <- estimatedContributions$logHr - 0.5 
    estimatedContributionsPre2$logHr  <- estimatedContributions$logHr - 0.5 
    estimatedContributionsPost1$logHr <- estimatedContributions$logHr + 0.5 
    estimatedContributionsPost2$logHr <- estimatedContributions$logHr + 0.5 
  }
  
  # Get the y values of the rectangles (contribution)
  estimatedContributionsPre1$value  <- 0 
  estimatedContributionsPre2$value  <- estimatedContributions$value
  estimatedContributionsPost1$value <- estimatedContributions$value
  estimatedContributionsPost2$value <- 0
  
  # Merge dataframes
  l <- list(estimatedContributionsPre1,estimatedContributionsPre2,
            estimatedContributionsPost1,estimatedContributionsPost2)
  
  estimatedContributions2 <- do.call(rbind,l) %>% 
    arrange(hr) %>% 
    arrange(intercalateNumber) %>% select(-hr)
  
  plot <- ggplot()
  
  factors  <- unique(estimatedContributions2$variable)
  nb.cols  <- length(factors)
  mycolors <- getPalette(nb.cols)
  
  # Split them by variable, useful for plotting
  estContrSplit <- split(estimatedContributions2, estimatedContributions2$variable)
  maxRows       <- max(sapply(estContrSplit, nrow))
  
  # geom_line doesn't follow row order to plot
  # Therefore, we can plot at most, 4 rows of each dataframe
  for (i in 1:(maxRows / 4)) {
    
    begin <- (i-1)*4+1
    end   <- begin + 3
    idx   <- begin:end
    
    # Select only the conditions with enough rows
    selectedDFs <- sapply(1:length(estContrSplit),function(x) nrow(estContrSplit[[x]]) >= end)
    
    # Merge dataframes
    estContrSplit2 <- estContrSplit[selectedDFs]
    estContrSplit3 <- lapply(estContrSplit2, function(x) x[idx,])
    tempDF         <- do.call(rbind,estContrSplit3)
    
    plot <- plot +
      geom_ridgeline(data=tempDF,
                     aes(logHr, variable, height = value,
                         fill=variable),color='black',scale=overlapFactor) +
      scale_fill_manual(breaks = factors, 
                        values = mycolors)
      
 
    
  }
  
  # Create the dataframe to add a baseline 
  vars    <- unique(estimatedContributions2$variable)
  lastDF1 <- data.frame('var'=vars,'logHr'=limit1LogHr)
  lastDF2 <- data.frame('var'=vars,'logHr'=limit2LogHr)
  lastDF  <- rbind(lastDF1,lastDF2)
  lastDF$value <- 0
  
  plot <- plot +
    geom_ridgeline(data=lastDF,
                   aes(logHr, var, height = value),
                   color='black',scale=overlapFactor) 

  # Detect if we used the hr (and not the diffusion coefficients!)
  if (!useDiffCoeff) {
    tickHr <- c(0.1,1,10,1e2,1e3,1e4,1e5,1e6)
    breaks <- sapply(tickHr, function(x) log(x,base=logScale))
    labels <- c(expression(paste("10"^"-1")),expression(paste("10"^"0")),
                expression(paste("10"^"1")),
                expression(paste("10"^"2")),expression(paste("10"^"3")),
                expression(paste("10"^"4")),expression(paste("10"^"5")),
                expression(paste("10"^"6")))
    
    plot <- plot + scale_x_continuous(breaks=breaks, label=labels)
  } else {
    tickHr <- c(1e-10,1e-12,1e-14)
    breaks <- sapply(tickHr, function(x) log(x,base=logScale))
    labels <- c(expression(paste("10"^"-10")),
                expression(paste("10"^"-12")),
                expression(paste("10"^"-14"))
                )
    plot <- plot + scale_x_continuous(breaks=breaks, label=labels,trans = "reverse")
  }
  
  return(plot) 
}

plothydrodynamicRadiusDistributionGray <- function(estimatedContributions,
                                                   overlapFactor=NULL) {
  
  # overlapFactor is not used, but we need to leave this argument
  # because of the way we call this function
  
  estimatedContributions$height <- 1
  estimatedContributions$value  <- as.numeric(estimatedContributions$value)

  ggplot(estimatedContributions, aes(hr, variable,
                                     height = height, group = variable,fill = factor(value))) + 
    geom_ridgeline_gradient(scale=1)   +
    scale_fill_grey(start = 1,end = 0) 
  
}

plothydrodynamicRadiusDistributionDual <- function(
    estimatedContributions,overlapFactor) {
  
  ggplot(estimatedContributions, aes(hr, variable, 
                                              height = value, group = variable,
                                     fill=variable)) + 
    geom_ridgeline(scale=overlapFactor,alpha=0.7)+
    scale_fill_cyclical(values = c("#FC766AFF", "#5B84B1FF"))
  
}

plotResiduals <- function(autocorrelationData,residuals,
                          axis_size,
                          residualsPlotSelection) {
  
  # autocorrelationData has 4 columns - 'variable', 'value', 'time'
  # 'value' versus 'time' represents the autocorrelation curve
  # 'variable' allows using different colours 
  # predictedAutocorrelation is a named list with key-value pairs:
  #                                               fitted curve id - fitted data points
  
  selection          <- as.numeric(rev(strsplit(residualsPlotSelection,'-')[[1]])[1])
  selectedConditions <- unique(autocorrelationData$variable)[(selection-24):selection]
  
  autocorrelationData$residuals <- residuals
    
  autocorrelationData <- autocorrelationData[
    autocorrelationData$variable %in% selectedConditions,]
  
  p <- ggplot(autocorrelationData,aes(x=time,y=residuals))+
    geom_point(size=0.4)+
    theme_bw(base_size = axis_size-3)+
    xlab("Time (s)")+
    ylab("Residuals")+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  p2 <- p + facet_wrap(~ variable,scales = "free",
                       ncol = 5,nrow = 5)
  
  return(p2)
  
}

plot_type_hr_distribution <- function(hr_plot_type,estimatedContributions,
                                      axis_size,overlapFactor,rowOrder='default',
                                      peakFilter=T,useAverage=F,mapping=NULL)  {
  
  # estimatedContributions has 3 columns - 'variable', 'value', 'hr' (or 'diff')
  # 'value' versus 'hr' (or 'diff') represents the hydrodynamic radii 
  # (or diffusion coefficients ) distribution
  # 'variable' allows using different colours 
  
  ## To avoid rewriting all functions and use the diffusion coefficients
  # we use this small hack 
  useDiffCoeff <- 'diff' %in% colnames(estimatedContributions)
  if (useDiffCoeff) {
    estimatedContributions$hr <- estimatedContributions$diff
  }
  
  estimatedContributions$value <- estimatedContributions$value*100 # Transform from 0 - 1 to 0 - 100
  
  # Average the samples, if required
  if (useAverage) {
    
    estimatedContributions$variable <- mapping$Group[match(estimatedContributions$variable, mapping$Condition)]
    
    # Define the order vector
    order_vector <- unique(estimatedContributions$variable)
    
    estimatedContributions <- estimatedContributions %>% 
      group_by(hr,variable) %>% 
      summarise(value = mean(value))
    
    # Sort by diff
    if (useDiffCoeff) {
      estimatedContributions <- estimatedContributions %>% arrange(-hr)
    }
    
    # Order the DataFrame using the order vector
    estimatedContributions <- estimatedContributions[order(match(estimatedContributions$variable, order_vector)), ]
    
  }
  
  samplesN <- length(unique(estimatedContributions$variable))
  
  if (peakFilter) {
    estimatedContributions$value[estimatedContributions$value < 1] <- 0
  }
  
  functionsList  <-list("rainbow"="plothydrodynamicRadiusDistribution",
                        "gray"="plothydrodynamicRadiusDistributionGray",
                        "dual"="plothydrodynamicRadiusDistributionDual",
                        'histogram'='plothydrodynamicRadiusDistributionHistogram')
  
  # Call the function using a string
  plottingFunction <- get(functionsList[[hr_plot_type]])
  
  factors <- unique(estimatedContributions$variable)
  
  if (rowOrder == "alphanumeric") {
    # Match row order and plot order
    estimatedContributions$variable <- factor(
      estimatedContributions$variable,
      str_sort(factors, numeric = TRUE), ordered = TRUE)
  } else {
    # Match row order and plot order
    estimatedContributions$variable <- factor(
      estimatedContributions$variable,factors, ordered = TRUE)
  }
  
  fig <- plottingFunction(estimatedContributions,overlapFactor * 0.01)
  fig <- formatHrDistributionPlot(fig,axis_size,useDiffCoeff)
  
  if (hr_plot_type != "histogram") {
    fig <- fig +
      scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                    labels = trans_format("log10", math_format(10^.x)))
  }
  
  return( fig )
  
}

# Plot only one simulated distribution
plot_simulated_distribution <- function(df) {
  
  # df should have one column called 'hr'(in nanometers)  
  # and one called 'value' (contributions, which can be intensity, volume
  # or number weighted)
  
  logScale <- df$hr[2] / df$hr[1]
  df$logHr <- log(df$hr,base = logScale)
  df$value <- df$value*100 # to percentage
  
  plot <- ggplot(data=df, aes(x=logHr, y=value)) +
    geom_bar(stat="identity",fill='red') +
    theme_classic(base_size = 18)
  
  tickHr <- c(0.1,1,10,1e2,1e3,1e4,1e5,1e6)
  breaks <- sapply(tickHr, function(x) log(x,base=logScale))
  labels <- c(expression(paste("10"^"-1")),expression(paste("10"^"0")),
              expression(paste("10"^"1")),
              expression(paste("10"^"2")),expression(paste("10"^"3")),
              expression(paste("10"^"4")),expression(paste("10"^"5")),
              expression(paste("10"^"6")))
  
  plot <- plot + scale_x_continuous(breaks=breaks, label=labels)+
    xlab("Hydrodynamic radius (nm)")
  
  return(plot)
}

plotLcurve <- function(lCurveDf,dfSel,axis_size,
                       residualsPlotSelection) {
  
  # lCurveDf (and dfSel) have 3 columns - 'variable', 'penalty', 'fidelity'
  # 'penalty' is the norm of the penalization term norm (Mx) where M
  # is the second order derivative matrix
  # 'fidelity' is the norm of the residuals (of the first order autocorrelation)
  # dfSel contains only the selected corners
  
  selection          <- as.numeric(rev(strsplit(residualsPlotSelection,'-')[[1]])[1])
  selectedConditions <- unique(lCurveDf$variable)[(selection-24):selection]
  
  lCurveDf <- lCurveDf[lCurveDf$variable %in% selectedConditions,]
  dfSel    <- dfSel[dfSel$variable %in% selectedConditions,]
  
  p <- ggplot(lCurveDf,aes(x=fidelity,y=penalty))+
    geom_point(size=0.4)+
    geom_point(data=dfSel,shape=4,color='blue',size=4) +
    theme_bw(base_size = axis_size-3)+
    xlab('Fidelity term')+
    ylab('Penalty term')+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))+
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x)))
  
  p2 <- p + facet_wrap(~ variable,scales = "free",
                       ncol = 5,nrow = 5)
  
  return(p2)
  
}




