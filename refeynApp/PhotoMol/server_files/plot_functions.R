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

plotRefeynFitNormalized <- function(refeyn,baseline,plot_width, plot_height, plot_type, axis_size,
                                    legends,colorPalette,sels,addMassesToLegend){
  
  dfMass <- data.frame("mass"=refeyn$masses_kDa)
  dfMass <- dfMass %>% 
    filter(mass >= refeyn$hist_window[1]) %>% 
    filter(mass <= refeyn$hist_window[2])
  
  # refeynFit is a matrix with n columns whose 1st column is the x axis, 
  # the 2nd to (n-1) column are the predicted values using the individual gaussians
  # the n column is the predicted value from the gaussian sum
  
  refeynFit  <- as.data.frame(refeyn$fit)
  prevNames  <- colnames(refeynFit)[-ncol(refeynFit)]
  colnames(refeynFit) <- c("x",prevNames[-1],"y")
  
  fig       <- plot_ly()
  baseline <- baseline / nrow(dfMass)# normalize
  
  # Only plot sum if we have more than one gaussian
  counterColorPalette <- 0
  if (ncol(refeynFit) > 3) {
    
    counterColorPalette <- counterColorPalette+1
    if (sels[counterColorPalette]) {
      gaussianSum   <- refeynFit[,c(1,ncol(refeynFit))] %>% filter(y > baseline+0.05 | y == 0)
      gaussianSum$y <- gaussianSum$y / nrow(dfMass) # Normalize
      
      fig <- fig %>% add_trace(data=gaussianSum,color=I(colorPalette[1]),x=~x,y=~y,
                               type = 'scatter', mode = 'lines',
                               name=legends[1],list(width = 3))
      
      
    }
    legends <- legends[-1]
    
  }
  
  fitted_means <- refeyn$fit_table[,1]
  
  if (addMassesToLegend)   legends <- paste0(legends," : ",round(fitted_means))
  legends <- paste0(legends," kDa")
  
  gaussianInd <- refeynFit[,-ncol(refeynFit)]
  colnames(gaussianInd) <- c("x",legends)
  
  gaussianInd <- reshape2::melt(gaussianInd,id.vars="x") %>% filter(value > baseline+0.05)
  gaussianInd$variable <- as.factor(gaussianInd$variable)
  
  for (var in unique(gaussianInd$variable)) {
    
    counterColorPalette <- counterColorPalette + 1
    
    if (sels[counterColorPalette]) {
      hexColor            <- colorPalette[counterColorPalette]
      
      tempDf <- gaussianInd[gaussianInd$variable == var,]
      
      tempDf$value <- tempDf$value / nrow(dfMass)# normalize
      
      fig    <- fig %>% add_trace(data=tempDf,x=~x,y=~value,
                                  color=I(hexColor),
                                  type = 'scatter', mode = 'lines',
                                  name=var,
                                  line = list(width = 3))
    }
    
  }
  
  start <- refeyn$hist_window[1]
  end   <- refeyn$hist_window[2]
  size  <- refeyn$bin_width
  
  fig <- fig %>% add_histogram(data=dfMass,x=~mass, color = I("#389196"),
                               alpha = 0.4, xbins=list(start=start,
                                                       end=end,
                                                       size=size),
                               histnorm = "probability",
                               showlegend=FALSE)
  
  y <- list(title = "Normalized counts",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size))
  
  xtitle <- "Mass (kDa)"

  x <- list(title = xtitle,titlefont = list(size = axis_size),
            tickfont = list(size = axis_size))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,showlegend = TRUE,
                        font="Roboto",legend = list(font = list(size = axis_size)))
  

  fig <- configFig(fig,plot_width, plot_height, plot_type,paste0("refeynPlotNormalized-",Sys.Date()))
  
  return(fig)
  
}

plotRefeynFit <- function(refeyn,baseline,plot_width, plot_height, plot_type, axis_size,
                          legends,colorPalette,sels,addMassesToLegend,contrasts=FALSE) {

  dfMass <- data.frame("mass"=refeyn$contrasts)
  if (!contrasts) {
    dfMass <- data.frame("mass"=refeyn$masses_kDa) 
  } 
 
  dfMass <- dfMass %>% 
    filter(mass >= refeyn$hist_window[1]) %>% 
    filter(mass <= refeyn$hist_window[2])
  
  # Replace mass with contrasts
  if (contrasts) {
    dfMass$mass <- dfMass$mass*factorForContrast 
  } 
  
  # refeynFit is a matrix with n columns whose 1st column is the x axis, 
  # the 2nd to (n-1) column are the predicted values using the individual gaussians
  # the n column is the predicted value from the gaussian sum
  
  refeynFit  <- as.data.frame(refeyn$fit)
  prevNames  <- colnames(refeynFit)[-ncol(refeynFit)]
  colnames(refeynFit) <- c("x",prevNames[-1],"y")
  
  if (contrasts) refeynFit$x <- refeynFit$x * factorForContrast
  
  fig <- plot_ly()
  # Only plot sum if we have more than one gaussian
  counterColorPalette <- 0
  if (ncol(refeynFit) > 3) {
    
    counterColorPalette <- counterColorPalette+1
    if (sels[counterColorPalette]) {
      gaussianSum <- refeynFit[,c(1,ncol(refeynFit))] %>% filter(y > baseline+0.05 | y == 0)
      
      fig <- fig %>% add_trace(data=gaussianSum,color=I(colorPalette[1]),x=~x,y=~y,
                               type = 'scatter', mode = 'lines',
                               name=legends[1],list(width = 3))
      
      
    }
    legends <- legends[-1]
    
  }
  
  fitted_means <- refeyn$fit_table[,1]
  
  if (contrasts) {
    fitted_means <- fitted_means*factorForContrast
  }
  
  if (addMassesToLegend)   legends <- paste0(legends," : ",round(fitted_means))
  
  if (!contrasts) legends <- paste0(legends," kDa")
  
  gaussianInd <- refeynFit[,-ncol(refeynFit)]
  colnames(gaussianInd) <- c("x",legends)
  
  gaussianInd <- reshape2::melt(gaussianInd,id.vars="x") %>% filter(value > baseline+0.05)
  gaussianInd$variable <- as.factor(gaussianInd$variable)
  
  for (var in unique(gaussianInd$variable)) {
    
    counterColorPalette <- counterColorPalette + 1
    
    if (sels[counterColorPalette]) {
      hexColor            <- colorPalette[counterColorPalette]
      
      tempDf <- gaussianInd[gaussianInd$variable == var,]
      fig    <- fig %>% add_trace(data=tempDf,x=~x,y=~value,
                                  color=I(hexColor),
                                  type = 'scatter', mode = 'lines',
                                  name=var,
                                  line = list(width = 3))
    }
    
  }
    
  
  
  start <- refeyn$hist_window[1]
  end   <- refeyn$hist_window[2]
  size  <- refeyn$bin_width
  
  if (contrasts) {
    start <- start*factorForContrast
    end   <- end*factorForContrast
    size  <- size*factorForContrast
  }
  
  fig <- fig %>% add_histogram(data=dfMass,x=~mass, color = I("#389196"),
                               alpha = 0.4, xbins=list(start=start,
                                                       end=end,
                                                       size=size),
                               showlegend=FALSE)
  
  y <- list(title = "Counts",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size))
  
  xtitle <- "Mass (kDa)"
  if (contrasts) xtitle <- "Ratiometric contrast * 1e3"
  
  x <- list(title = xtitle,titlefont = list(size = axis_size),
            tickfont = list(size = axis_size))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,showlegend = TRUE,
                        font="Roboto",legend = list(font = list(size = axis_size)))
  
  
  fig <- configFig(fig,plot_width, plot_height, plot_type,paste0("refeynPlot-",Sys.Date()))
  
  return(fig)
  
}

truncatedGaussian <- function(x,mean,std,amplitude,leftBound) {
  
  return( (x >= leftBound) * (amplitude * exp( - (x-mean)**2 / (2 * (std**2) ))))
}

addSimulation2plotRefeynFit <- function(fig,mean,std,amplitude,leftBound) {
  
  xSeq <- seq(max(leftBound,mean-std*2.5),mean+std*2.5)
  y    <- truncatedGaussian(xSeq,mean,std,amplitude,leftBound)
  
  df   <- data.frame(xSeq,y)
  
  fig %>% add_trace(data=df,color=I("black"),x=~xSeq,y=~y,
                           type = 'scatter', mode = 'lines',
                           name="Simulation",list(width = 3))
  
}

addLabels2plotRefeynFit <- function(fig,means,amplitudes,selected,text_size,
                                    contrasts=FALSE) {
  
  # means gives the x position of the labels and also the text
  # amplitudes gives the y position of the labels
  # selected is a boolean vector to select which labels to add
  # text size controls the labels size

  t <- list(
    family = "Roboto",
    size = text_size-2,
    color = toRGB("black"))
  
  if (sum(selected) == 0) {return(fig)} # No labels
  
  if (length(selected) > 1) {
    
    if (sum(selected) == 1 && selected[1]) {return(fig)} # No labels
    
    selected <- selected[-1] # Remove gaussian sum from the selection
  } 
  
  yShift <- max(amplitudes)*0.02

  if (contrasts) means <- means * factorForContrast
  
  dfLabel <- data.frame("means"=means,"amplitudes"=amplitudes+yShift,
                        "labels"=paste0("   ",round(means)))

  dfLabel <- dfLabel[selected,]

  if (!contrasts) dfLabel$labels <- paste0(dfLabel$labels," kDa")
    
  fig <- fig %>% add_text(data = dfLabel, x = ~means, y = ~amplitudes, 
                          text = ~labels,
                          textfont = t, textposition = "top right",
                          showlegend=FALSE)
  
}

plotMass_vs_contrast <- function(mass,contrast,slope,intercept,
                                 plot_width, plot_height, plot_type,axis_size) {
  
  fig <- plot_ly()
  
  df <- data.frame(mass,contrast)

  fig    <- fig %>% add_trace(data=df,x=~mass,y=~contrast,
                              color=I("#377EB8"),
                              type = 'scatter', mode = 'markers'
                              )
  
  x1 <- min(mass) 
  y1 <- x1*slope + intercept
  x2 <- max(mass) 
  y2 <- x2*slope + intercept
  
  dfPred <- data.frame(x=c(x1,x2),y=c(y1,y2))

  fig    <- fig %>% add_trace(data=dfPred,x=~x,y=~y,
                              color=I("#ffa500"),
                              type = 'scatter', mode = 'lines',
                              list(width = 3)
  )
  
  y <- list(title = "Ratiometric contrast",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size))
  
  x <- list(title = "Mass (kDa)",titlefont = list(size = axis_size),
            tickfont = list(size = axis_size))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,showlegend = FALSE,
                        font="Roboto",legend = list(font = list(size = axis_size)))
  
  fig <- configFig(fig,plot_width, plot_height, plot_type,paste0("calibrationMassVsContrast-",Sys.Date()))
  return(fig)
  
}

plotRefeynMassHist <- function(refeyn,plot_width, plot_height, plot_type,axis_size) {
  
  dfMass <- data.frame("mass"=refeyn$masses_kDa) 
  fig     <- plot_ly()

  fig <- fig %>% add_histogram(data=dfMass,x=~mass, color = I("#389196"),
                               alpha = 0.4, xbins=list(size=refeyn$bin_width),
                               showlegend=FALSE)
  
  y <- list(title = "Counts",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size))
  
  xtitle <- "Mass (kDa)"

  x <- list(title = xtitle,titlefont = list(size = axis_size),
            tickfont = list(size = axis_size))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,showlegend = TRUE,
                        font="Roboto",legend = list(font = list(size = axis_size)))
  
  
  fig <- configFig(fig,plot_width, plot_height, plot_type,paste0("refeynPlotBindingAndUnbindingEvents-",Sys.Date()))
  
  return(fig)
  
}

