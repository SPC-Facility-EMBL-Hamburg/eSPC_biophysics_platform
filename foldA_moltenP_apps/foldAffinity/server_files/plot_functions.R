## Plot fluorescence signal versus temperature (color by ligand concentration)
## fluo_m is a 3 column dataframe: temp, fluo and conc. 
## which is the signal (350,330 or Ratio)
## Returns the plot

plot_fluo_signal <- function(fluo_m,which="Ratio",
                             plot_width, plot_height, plot_type, axis_size) {
  
  x <- list(title = "Temperature (°C)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = which,titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fluo_m0 <- fluo_m %>% filter(conc == 0)
  fluo_m  <- fluo_m %>% filter(conc > 0)
  
  fig <- plot_ly(type = 'scatter', mode = 'lines')
  
  if (nrow(fluo_m) == 0) { 
    
    fig <- fig %>% add_trace(data = fluo_m0 ,name = ~conc_, 
                             color = I("#389196"),y    = ~fluo   , x   = ~temp) 
    
    # if we ONLY have positions with no ligand return this plot
    return(fig) 
    }
  
  for (position in unique(fluo_m0$conc_)) {
    fluo_m0_temp <- fluo_m0 %>% filter(conc_== position) 

    # set as black the positions with no ligand
    fig <- fig %>% add_trace(data  = fluo_m0_temp,
                             color = I("black"),
                             y     = ~fluo, x=~temp,name=position) 
  }
  
  conc_vector <- unique(fluo_m$conc)
  n_colors   <- sapply(conc_vector,get_color_from_conc,min(fluo_m$conc),max(fluo_m$conc))
  
  for (position in unique(fluo_m$conc_)) {
    
    fluo_m_temp <- fluo_m %>% filter(conc_== position) 
    idx        <- which(conc_vector == unique(fluo_m_temp$conc))
    
    fig <- fig %>% add_trace(data  = fluo_m_temp,
                             color = I(n_colors[idx]),
                             y = ~fluo, x=~temp,name=position) 
  }
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend = FALSE)
  
  fig <- fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE)
  
  return(fig)
}


# Plot maximum of derivative
generate_der_plot <- function(tms,concentrations,
                              plot_width, plot_height, plot_type, axis_size,
                              logScaleType){
  
  #logScaleType can be defaultPlotly, micromolar, milimolar, molar
  
  df2 <- generate_der_df(tms,concentrations)
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  
  df2$concentration <- df2$concentration * unlist(factorList[logScaleType])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[logScaleType])
  
  if (logScaleType == "defaultPlotly"){
    exponentformat <- "B"
  } else {
    exponentformat <- 'power'
  }
  
  x <- list(title = xTitle,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat)
  
  y <- list(title = "Tm (°C) (using the 1st derivative)",titlefont = list(size = axis_size),
            tickfont = list(size = axis_size))
  
  fig <- plot_ly(df2, y = ~Tm, x = ~concentration, type = 'scatter',
                 mode = "markers", marker = list(size = 10,color = "rgba(255, 182, 193, .9)"))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend = FALSE)
  
  fig <- fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE)
  
  return(fig)
}

## Plot initial fluorescence signal versus ligand concentration
## fluo_m is a 3 column dataframe: temp, fluo and conc. 
## which is the signal (350,330 or Ratio)
## Returns the plot

plot_initialFluo_signal <- function(fluo_m,which="Ratio",
                             plot_width, plot_height, plot_type, axis_size) {
  
  x <- list(title = "Ligand concentration (M)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),type = "log")
  y <- list(title = which,titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fluo_m  <- fluo_m %>% filter(conc > 0)
  
  fig <- plot_ly(type = 'scatter', mode = 'markers')
  
  if (nrow(fluo_m) == 0) { return(NULL) }
  
  fluo_m <- fluo_m[fluo_m$temp==min(fluo_m$temp),]
  
  # set as black the positions with no ligand
  fig <- fig %>% add_trace(data    = fluo_m,
                           color = I("black"),
                           y     = ~fluo, x=~conc) 
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto",showlegend = FALSE)
  
  fig <- fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE)
  
  return(fig)
}

## Plot the fluorescence fit.
## Requires the list of dataframes from both the experimental and fitted fluorescence
## Plots the selected element from that list

plot_fluorescence_fit <- function(fluo_m_real_list,fluo_m_pred_list,selected) {
  
  fluo_m <- fluo_m_real_list
  
  fluo_m <- fluo_m[[selected]]
  
  if (max(fluo_m$conc) > 1e-3) {
    fluo_m$conc <- fluo_m$conc*1e3
    scale_f <- "mM"
  } else {
    fluo_m$conc <- fluo_m$conc*1e6
    scale_f <- "µM"
  }
  
  fluo_m$legend <- paste0("Rep ",fluo_m$rep, ". ",signif(fluo_m$conc,2)," ",scale_f)
  neworder <- unique(fluo_m$legend)
  fluo_m <- arrange(transform(fluo_m,legend=factor(legend,levels=neworder)),legend)
  
  p <- ggplot(fluo_m,aes(x=temp,y=fluo))+
    geom_point(size=0.4)+
    theme_bw(base_size = 14)+
    xlab("Temperature (ºC)")+
    ylab("Fluorescence (AU)")+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  
  fluo_m <- fluo_m_pred_list[[selected]]
  
  if (max(fluo_m$conc) > 1e-3) {
    fluo_m$conc <- fluo_m$conc*1e3
    scale_f <- "mM"
  } else {
    fluo_m$conc <- fluo_m$conc*1e6
    scale_f <- "µM"
  }
  
  fluo_m$legend <- paste0("Rep ",fluo_m$rep, ". ",signif(fluo_m$conc,2)," ",scale_f)
  neworder <- unique(fluo_m$legend)
  fluo_m <- arrange(transform(fluo_m,legend=factor(legend,levels=neworder)),legend)
  
  p2 <- p +
    geom_line(data=fluo_m,color="red",linewidth=0.3,
              aes(x=temp,y=fluo),inherit.aes = FALSE)
  
  p3 <- p2 + facet_wrap(~ legend,scales = "free",
                        ncol = global_plot_columns,nrow = global_plot_rows,drop = FALSE)
  
  return(p3)
}

## Plot selected fluorescence fit if possible

plot_sel_fluo_fit <- function(real_data,model_data,selected){
  if (length(real_data) < selected) {return(NULL)}
  
  p <- plot_fluorescence_fit(real_data,model_data,selected)
  return(p)
}

## Plot the isothermal exp data
## iso_real_data has two columns: conc and fu (unfolded_fraction, "experimental")
## Returns the plot

plot_isothermal_fitting_exp <- function(iso_real_data,plot_width, plot_height, plot_type, axis_size,
                                        logScaleType) {
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  
  iso_real_data$conc <- iso_real_data$conc * unlist(factorList[logScaleType])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[logScaleType])
  
  if (logScaleType == "defaultPlotly"){
    exponentformat <- "B"
  } else {
    exponentformat <- 'power'
  }
  
  xaxis <- list(title = xTitle,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat)
  
  yaxis <- list(title = "Unfolded Fraction",titlefont = list(size = axis_size),
            tickfont = list(size = axis_size))
  
  plot_list <- list()
  x_pos_annot <- min(log10(iso_real_data$conc)) + (max(log10(iso_real_data$conc)) - min(log10(iso_real_data$conc))) / 2
  
  tot_cond <- length(unique(iso_real_data$legend))
  
  i <- 0
  for (leg in unique(iso_real_data$legend)) {
    i <- i + 1
    
    df_temp <- iso_real_data %>% filter(legend == leg)
    
    fig <- plot_ly(df_temp, x = ~conc, y = ~fu, color = I("#00AFBB"), type = "scatter") %>% 
      layout(xaxis = xaxis,yaxis=yaxis,
             annotations = list(x = x_pos_annot , y = 1.12, 
                                text = ifelse(tot_cond > 1,leg,""), showarrow = F, 
                                xref='x', yref='paper',font = list(size = axis_size)),
             showlegend=FALSE)
      
    plot_list[[i]] <- fig
    
  }
  
  if (i > 1) {
    
    return( subplot(plot_list,nrows = 2,margin = c(0.03,0.03,0.1,0.1),
                    shareY = TRUE, shareX = TRUE) %>%  config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  } else {
    return( plot_list[[1]] %>%  layout(
      title=list(text=unique(iso_real_data$legend),font = list(size = axis_size))) %>% config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  }

}

## Plot the isothermal fitting
## iso_real_data has two columns: conc and fu (unfolded_fraction, "experimental")
## iso_model_data has four columns: conc fu (predicted) fu_low fu_upp (lower and upper bound of the error)
## Returns the plot

plot_isothermal_fitting <- function(iso_real_data,iso_model_data,plot_width, plot_height, plot_type,
                                    axis_size,logScaleType) {
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  
  iso_real_data$conc  <- iso_real_data$conc * unlist(factorList[logScaleType])
  iso_model_data$conc <- iso_model_data$conc * unlist(factorList[logScaleType])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[logScaleType])
  
  if (logScaleType == "defaultPlotly"){
    exponentformat <- "B"
  } else {
    exponentformat <- 'power'
  }
  
  xaxis <- list(title = xTitle,titlefont = list(size = axis_size), 
                tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat)
  
  yaxis <- list(title = "Unfolded Fraction",titlefont = list(size = axis_size),
                tickfont = list(size = axis_size))
  
  plot_list <- list()
  
  x_pos_annot <- min(log10(iso_real_data$conc)) + (max(log10(iso_real_data$conc)) - min(log10(iso_real_data$conc))) / 2
  
  tot_cond <- length(unique(iso_real_data$legend))
  
  i <- 0
  for (leg in unique(iso_model_data$legend)) {
    i <- i + 1
    
    df_temp <- iso_model_data %>% filter(legend == leg) %>% arrange(conc)
    
    fig <- plot_ly(type = 'scatter',mode = 'lines+markers')
    
    # Uncomment this lines when fu_upp and fu_low are calculated properly!!!! 
    # Real confidence bands from non linear regression
    
    #fig <- fig %>% add_trace(data=df_temp, x = ~conc, y = ~fu_upp, type = 'scatter',mode = 'lines',
    #                         line = list(color = 'rgba(255,163,0,0.2)'))
    
    #fig <- fig %>% add_trace(data=df_temp, x = ~conc,y = ~fu_low, type = 'scatter', mode = 'lines',
    #                         fill = 'tonexty', fillcolor='rgba(255,163,0,0.2)', 
    #                        line = list(color = 'rgba(255,163,0,0.2)'))
    
    fig <- fig %>% add_trace(data=df_temp, x = ~conc,y = ~fu, type = 'scatter', mode = 'lines',
                             line = list(color = 'rgba(255,163,0,0.5)'))
    
    df_temp2 <- iso_real_data %>% filter(legend == leg) %>% arrange(conc)
    
    fig <- fig %>% add_trace(data=df_temp2, x = ~conc, y = ~fu, color = I("#00AFBB"),
                             type = "scatter",mode = "markers")
    
    m   <- list(l = 10,r = 10,b = 20,t = 40,pad = 4)
    fig <- fig %>% layout(xaxis = xaxis,yaxis=yaxis,
                          annotations = list(x = x_pos_annot , y = 1.13,
                                             text = ifelse(tot_cond > 1,leg,""), showarrow = F,
                                             xref='x', yref='paper',font = list(size = axis_size)),
                          showlegend=FALSE)
    
    plot_list[[i]] <- fig
    
  }
  
  if (i > 1) {
    return( subplot(plot_list,nrows = 2,margin = c(0.03,0.03,0.1,0.1),
                    shareY = TRUE, shareX = TRUE) %>%  config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  } else {
    return( plot_list[[1]] %>%  layout(title=list(text=unique(iso_real_data$legend),
                                                  font = list(size = axis_size)),
                                       margin=m) %>%  config(
      toImageButtonOptions = list(
        format = plot_type,
        filename = "myplot",
        width = plot_width * 50,
        height = plot_height * 50
      ), displaylogo = FALSE,
      modeBarButtonsToRemove = list('sendDataToCloud',
                                    'hoverClosestCartesian','hoverCompareCartesian',
                                    'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                    'zoom2d','pan2d','autoScale2d','resetScale2d')))
  }
  
}

plot_tm_shift <- function(conc,tms,plot_width, plot_height, plot_type, axis_size,
                          logScaleType) {
  
  df_model <- get_tm_df(conc,tms)
  
  df_model$tms <- df_model$tms - 273.15 # To centigrade
  
  factorList <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  df_model$l_conc <- df_model$l_conc * unlist(factorList[logScaleType])
  
  titleList <- list("defaultPlotly" = "[Ligand] (M)"   , "molar"      = "[Ligand] (M)",
                    "milimolar"     = "[Ligand] (mM)"  , "micromolar" = "[Ligand] (μM)")
  
  xTitle <- unlist(titleList[logScaleType])
  
  if (logScaleType == "defaultPlotly"){
    exponentformat <- "B"
  } else {
    exponentformat <- 'power'
  }
  
  xaxis <- list(title = xTitle,titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),type = "log",exponentformat=exponentformat)
  
  yaxis <- list(title = "TmObs (°C)",titlefont = list(size = axis_size),
                tickfont = list(size = axis_size))
  
  fig <- plot_ly(df_model, x = ~l_conc, y = ~tms, color = I("#00AFBB"), type = "scatter") %>% 
    layout(xaxis = xaxis,yaxis=yaxis,showlegend=FALSE) %>%  
    config(toImageButtonOptions = list(
      format = plot_type,
      filename = "myplot",
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE,
    modeBarButtonsToRemove = list('sendDataToCloud',
                                  'hoverClosestCartesian','hoverCompareCartesian',
                                  'lasso2d','select2d','zoomIn2d','zoomOut2d',
                                  'zoom2d','pan2d','autoScale2d','resetScale2d'))
  
  return(fig)
}

plot_tm_shift_fit <- function(conc,tms,df_pred,fit_info,asymmetricCI95,
                              plot_width, plot_height, plot_type, axis_size,
                              logScaleType) {
  
  legend <- paste0("Kd1 (M) = ",
                   signif(fit_info$estimate[2],2)," ± ",
                   round(fit_info$std.error[2]/fit_info$estimate[2]*100,0),"%"
  )
  
  if (!(is.null(asymmetricCI95))) {
    
    legend <- paste0("Kd (M) : ",
                     signif(fit_info$estimate[2],2)," Asymmetric CI95 : [",
                     signif(asymmetricCI95$kd_min95,2)," ; ",
                     signif(asymmetricCI95$kd_max95,2),"]"
                     
    )
    
  }
  
  if (nrow(fit_info) == 3){
    legend <- paste0(legend,
                     " Kd2 (M) = ",
                     signif(fit_info$estimate[3],2)," ± ",
                     round(fit_info$std.error[3]/fit_info$estimate[3]*100,0),"%"
    )
  } 
  
  df_pred$tms <- df_pred$tms - 273.15 # To centigrade
  
  df_pred <- df_pred %>% arrange(l_conc)
  
  factorList     <- list("defaultPlotly"=1,"molar"=1,"milimolar"=1e3,"micromolar"=1e6)
  df_pred$l_conc <- df_pred$l_conc * unlist(factorList[logScaleType])
  
  m   <- list(l = 10,r = 10,b = 20,t = 40,pad = 4)
  
  fig <- plot_tm_shift(conc,tms,plot_width, plot_height, plot_type, axis_size,
                       logScaleType) %>% 
    add_trace(data=df_pred, x = ~l_conc,y = ~tms, type = 'scatter', mode = 'lines',
              line = list(color = 'rgba(255,163,0,0.5)')) %>% 
    layout(title=legend,margin=m)
    
  return(fig)
}

