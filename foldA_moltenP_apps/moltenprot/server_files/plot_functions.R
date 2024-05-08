# change width, height and plot type 
config_figure <- function(fig,plot_type,plot_width,plot_height) {
  
  fig %>%  config(toImageButtonOptions = list(
      format = plot_type,filename = "myplot",
      width = plot_width * 50,height = plot_height * 50
    ), displaylogo = FALSE)
  
}

## Plot fluorescence signal versus temperature (color by condition)
## fluo_m is a 3 column dataframe: temp, fluo and cond 

plot_fluo_signal <- function(fluo_m,y_label, plot_width, plot_height, plot_type, 
                             legend_text_size,axis_size) {
  
  # Avoid duplicate names
  fluo_m <- avoid_positions_with_the_same_name_in_df(fluo_m)
  fluo_m$temp <- fluo_m$temp - 273.15 # To celsius 
  
  x <- list(title = "Temperature (°C)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = y_label,titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fig <- plot_ly(data = fluo_m, x = ~temp, y = ~fluo,color=~Condition,
                 type = 'scatter', mode = 'lines')
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  fig <- fig %>% layout(showlegend = TRUE, legend = list(font = list(size = legend_text_size)))
                    
  return(fig ) 
}

#plot_spectra <- function(signal_data_dictionary,temp_data_dictionary,
#                         min_temp,max_temp,conditions) {
  
  #dfs <- lapply(signal_data_dictionary, function(l) {
    
    #  })
  
#}

# Plot the whole spectral data
# Requires: a dataframe with columns named 'wl', 'value', 'temps' and 'name'

# 'wl'    : wavelength in nanometers,       for the x-axis
# 'value' : the measured fluorescence,      for the y-axis
# 'temps' : temperature in degree Celsius,  for coloring
# 'name'  : name of the condition,          for dividing the plot

plot_whole_spectra <- function(df,font_size=18) {
  
  ggplot(df,aes(x=wl,y=value,color=temps,group=temps))+
    geom_line()+
    facet_wrap(~name)+
    theme_classic(base_size = font_size)+
    ylab('Fluorescence (AU)')+
    xlab('Wavelength (nm)')+
    scale_color_viridis_c(name='Temperature (°C)')
  
}

# Plot maximum of derivative
generate_max_der_plot <- function(tms,conditions,plot_width, plot_height, plot_type, 
                                  legend_text_size,axis_size){

  df2 <- generate_max_der_df(tms,conditions)

  df2$condition <- factor(df2$condition, levels = unique(df2$condition))
  
  x <- list(title = "",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = "Tm (°C) (using the 1st derivative)",
            titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fig <- plot_ly(df2, y = ~Tm, x = ~condition, type = 'scatter',
                 mode = "markers", marker = list(size = 10,color = "rgba(255, 182, 193, .9)"))
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  fig <- fig %>% layout(showlegend = FALSE, legend = list(font = list(size = legend_text_size)))
  
  return(fig)
}

## Plot the fluorescence fit.
## Requires the list of dataframes from both the experimental and fitted fluorescence
## Plots the selected element from that list

plot_fluorescence_fit <- function(fluo_m_real_list,fluo_m_pred_list,selected,is_report_plot=FALSE) {
  
  fluo_m      <- fluo_m_real_list[[selected]]
  fluo_m_pred <- fluo_m_pred_list[[selected]]
  
  # Avoid duplicate names
  fluo_m        <- avoid_positions_with_the_same_name_in_df(fluo_m)
  
  fluo_m$legend <- paste0(fluo_m$Condition)
  fluo_m$legend <- factor(fluo_m$legend , levels = unique(fluo_m$legend ))
  fluo_m$temp   <- fluo_m$temp - 273.15 # To celsius 
  
  # Avoid duplicate names (predicted data)
  fluo_m_pred <- avoid_positions_with_the_same_name_in_df(fluo_m_pred)
  
  fluo_m_pred$legend <- paste0(fluo_m_pred$Condition)
  fluo_m_pred$legend <- factor(fluo_m_pred$legend , levels = unique(fluo_m_pred$legend ))
  fluo_m_pred$temp   <- fluo_m_pred$temp - 273.15 # To celsius
  
  # Return plot for shiny app (facet wrap is fine here)
  if (!(is_report_plot)) {
    
    p <- ggplot(fluo_m,aes(x=temp,y=fluo))+
      geom_point(size=0.4)+
      theme_bw(base_size = 14)+
      xlab("Temperature (ºC)")+
      ylab("Fluorescence (AU)")+
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"))
    
    p2 <- p +
      geom_line(data=fluo_m_pred,color="red",linewidth=0.4,
                aes(x=temp,y=fluo),inherit.aes = FALSE)
    
    p3 <- p2 + facet_wrap(~ legend,scales = "free",
                          ncol = global_plot_columns,nrow = global_plot_rows,drop = FALSE)
   
    return(p3)         
  }
  
  # In the report we want all the plots with the same size so we change the plotting routine
  
  plots_list    <- list()
  total_plots   <- length(levels(fluo_m$legend))
  missing_plots <- global_chunck_n - total_plots
  
  for (i in 1:total_plots) {
    
    selected_legend <- levels(fluo_m$legend)[i]
    
    xlab <- ylab <- ""
    
    if (i %in% c(1,5,9,13))                 { ylab <- "Fluorescence (AU)" }
    if (i %in% (total_plots-4):total_plots) { xlab <- "Temperature (ºC) " }
    
    p <- ggplot(fluo_m %>% filter(legend == selected_legend),aes(x=temp,y=fluo))+
      geom_point(size=0.4)+
      theme_bw(base_size = 7)+
      xlab(xlab)+
      ylab(ylab)+
      theme(panel.grid = element_blank(),
            strip.background = element_blank(),
            panel.border = element_rect(colour = "black"))+
      ggtitle(substr(selected_legend, 1, 19))
    
    p2 <- p + geom_line(data=fluo_m_pred %>% filter(legend == selected_legend),
                        aes(x=temp,y=fluo),inherit.aes = FALSE,
                        color="red",linewidth=0.4)

    plots_list[[i]] <- p2
  }
  
  # Add dummy plots so every plot has the same size
  if (missing_plots > 0) {
    dummy_plot <- ggplot() + theme_void()
    plots_list[((global_chunck_n-missing_plots)+1):global_chunck_n] <- rep(list(dummy_plot), missing_plots, simplify=FALSE)
  }

  grid.newpage()

  return(egg::ggarrange(plots=plots_list))
}

## Plot the fluorescence fit standarized residuals
## Requires the list of dataframes from both the experimental and fitted fluorescence. 
## And the standard error of the fitting and fitted conditions
## Plots the selected element from that list

plot_fluorescence_residuals <- function(fluo_m_real_list,fluo_m_pred_list,selected,standard_error,fitted_conditions) {
  
  fluo_m        <- fluo_m_real_list[[selected]]
  
  fluo_m$legend <- paste0(fluo_m$Condition)
  neworder      <- unique(fluo_m$legend)
  fluo_m        <- dplyr::arrange(transform(fluo_m,legend=factor(legend,levels=neworder)),legend)
  fluo_m$temp   <- fluo_m$temp - 273.15 # To celsius 
  
  fluo_m2        <- fluo_m_pred_list[[selected]]
  fluo_m2$legend <- paste0(fluo_m2$Condition)
  neworder       <- unique(fluo_m2$legend)
  fluo_m2        <- dplyr::arrange(transform(fluo_m2,legend=factor(legend,levels=neworder)),legend)
  fluo_m2$temp   <- fluo_m2$temp - 273.15 # To celsius 
  
  std_error_df     <- data.frame(standard_error=unlist(standard_error),
                                 cond_=c(paste0(1:length(fitted_conditions),"_",fitted_conditions)))
  
  fluo_residual    <- fluo_m$fluo - fluo_m2$fluo
  fluo_residual_df <- data.frame(temp=fluo_m$temp,residual=fluo_residual,
                                 Condition=fluo_m$Condition,legend=fluo_m$legend,cond_=fluo_m$cond_)
  
  fluo_residual_df <- inner_join(fluo_residual_df,std_error_df)
  fluo_residual_df$standarized_residual <- fluo_residual_df$residual / fluo_residual_df$standard_error
  
  fluo_residual_df <- avoid_positions_with_the_same_name_in_df(fluo_residual_df)
  
  p <- ggplot(fluo_residual_df,aes(x=temp,y=standarized_residual))+
    geom_point(size=0.4)+
    theme_bw(base_size = 14)+
    xlab("Temperature (ºC)")+
    ylab("Standarized Residual")+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))
  
  p2 <- p + facet_wrap(~ Condition,scales = "free",
                       ncol = global_plot_columns,nrow = global_plot_rows)
  
  return(p2)
}

# Generate the fraction unfolded versus temperature plot
generate_fractions_plot <- function(dhs,tms,conditions,t_onset,
                                    plot_width, plot_height, plot_type, 
                                    legend_text_size,axis_size,color_by_capillary) {
  
  dfs <- list()
  
  i <- 0
  for (temp in seq(25,95,1)) {
    
    i <- i+1
    fus <- mapply(get_fraction_unfolded_EquilTwoState, dhs, tms, temp+273.15)
    df_unfolded <- data.frame("fu"=fus,"Condition"=conditions)
    df_unfolded$temp <- temp
    
    df_unfolded$Condition <- sapply(df_unfolded$Condition,trimws,"both")
    df_unfolded$cond_     <- c(paste0(1:nrow(df_unfolded),"_",df_unfolded$Condition))
    
    dfs[[i]] <- df_unfolded
    
  }
  
  df_tog <- do.call(rbind,dfs)
  
  # Avoid duplicate names
  df_tog <- avoid_positions_with_the_same_name_in_df(df_tog)
  
  # Generate t_onset versus t_m dataframe. Only useful to set the colors of the plot
  df     <- generate_tm_tonset_df(conditions,tms,t_onset)
  # Set the new order of the dataframe
  neworder       <- unique(df$condition)
  df_tog         <- dplyr::arrange(
    transform(df_tog,original_condition=factor(original_condition,levels=neworder)),original_condition)
  
  x <- list(title = "Temperature (°C)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = "Unfolded fraction ",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  # if color_by_capillary, we will have one color by capillary / position
  # otherwise, we have one color by condition name
  if (color_by_capillary) {
    fig <- plot_ly(data = df_tog, x = ~temp, y = ~fu,color=~Condition,
                   type = 'scatter', mode = 'markers',marker = list(size = 8))
    fig <- fig %>% layout(showlegend = TRUE, legend = list(font = list(size = legend_text_size)))
  } else {
    fig <- plot_ly(type = 'scatter', mode = 'markers',marker = list(size = 8))
    
    unique_cond  <- levels(df_tog$original_condition)
    total_colors <- length(unique_cond)
    
    if (total_colors <= 8) {
      palette      <- RColorBrewer::brewer.pal(total_colors, "Set2")
    } else {
      palette      <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(total_colors)
    }
    
    for (cond in unique_cond) {
      
      df2_temp <- df_tog %>% filter(original_condition == cond) 
      
      current_original_condition <- unique(df2_temp$original_condition)
      idx        <- which(current_original_condition == unique_cond)
      
      fig <- fig %>% add_trace(data  = df2_temp,
                               color = I(palette[idx]),
                               y = ~fu, x= ~temp, showlegend = T,
                               name=current_original_condition) 
    }
  }
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  fig <- fig %>% layout(showlegend = TRUE, legend = list(font = list(size = legend_text_size)))
  
  return(fig)
  
}

generate_tm_plot <- function(tms,conditions,
                             plot_width, plot_height, plot_type, 
                             axis_size){
  
  df <- data.frame("condition"=conditions,"Tm"=tms)
  df$Tm <-  df$Tm - 273.15
  df <- df %>% dplyr::arrange(-Tm)
  
  df2 <- df[1:25,]
  df2$condition  <- as.character(df2$condition)
  neworder       <- unique(df2$condition)
  df2            <- dplyr::arrange(transform(df2,condition=factor(condition,levels=neworder)),condition)
  
  x <- list(title = "",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = "Tm (°C)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fig <- plot_ly(df2, y = ~Tm, x = ~condition, type = 'scatter',
                 mode = "markers", marker = list(size = 10,color = "rgba(255, 182, 193, .9)"))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  fig <- fig %>% layout(showlegend = FALSE)
  
  return(fig)
}

generate_tm_tonset_plot <- function(tms,conditions,t_onset,
                                    plot_width, plot_height, plot_type, 
                                    legend_text_size,axis_size){
  
  df <- generate_tm_tonset_df(conditions,tms,t_onset)
  
  df2 <- df[1:25,]
  df2$condition  <- as.character(df2$condition)
  neworder       <- unique(df2$condition)
  df2            <- dplyr::arrange(transform(df2,condition=factor(condition,levels=neworder)),condition)
  
  x <- list(title = "T onset (1 % unfolded)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = "Tm (°C)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fig <- plot_ly(type = 'scatter', mode = 'markers',marker = list(size = 10))
  unique_cond  <- levels(df2$condition)
  total_colors <- length(unique_cond)
  
  if (total_colors <= 8) {
    palette      <- RColorBrewer::brewer.pal(total_colors, "Set2")
  } else {
    palette      <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(total_colors)
  }
  
  for (cond in unique_cond) {
    
    df2_temp <- df2 %>% filter(condition == cond) 
    
    current_original_condition <- unique(df2_temp$condition)
    idx        <- which(current_original_condition == unique_cond)
    
    fig <- fig %>% add_trace(data  = df2_temp,
                             color = I(palette[idx]),
                             y = ~Tm, x= ~t_onset, showlegend = T,
                             name=current_original_condition) 
  }
 
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  fig <- fig %>% layout(showlegend = TRUE, legend = list(font = list(size = legend_text_size)))
  
  return(fig)
}

generate_2t_plot <- function(t1,t2,conditions,score,
                             plot_width, plot_height, plot_type, 
                             legend_text_size,axis_size){
  
  df    <-  data.frame(condition=conditions,T1=t1,T2=t2,score=score)
  df$T1 <-  df$T1 - 273.15
  df$T2 <-  df$T2 - 273.15
  df <- df %>% dplyr::arrange(-score)
  
  df2 <- df[1:25,]
  df2$condition  <- as.character(df2$condition)
  neworder       <- unique(df2$condition)
  df2            <- dplyr::arrange(transform(df2,condition=factor(condition,levels=neworder)),condition)
  
  x <- list(title = "",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = "T1 and T2 (°C)",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  fig <- plot_ly(df2, y = ~T1, x = ~condition, type = 'scatter',showlegend = FALSE,
                 mode = "markers", marker = list(size = 10,color = "rgba(255, 182, 193, .9)"))
  
  fig <- fig %>% add_trace(y = ~T2, x = ~condition, type = 'scatter',showlegend = FALSE,
                           mode = "markers", marker = list(color = "rgba(255, 182, 193, .9)"))
  
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  
  return(fig)
}

# score_table is a data frame where the first column has the 'score'
# and the last column has the 'condition name'
generate_score_plot <- function(score_table,
                                plot_width,plot_height, 
                                plot_type,axis_size){
  
  x <- list(title = "",titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  y <- list(title = colnames(score_table)[1],titlefont = list(size = axis_size), tickfont = list(size = axis_size))
  
  score_table           <- score_table[,c(1,ncol(score_table))]
  colnames(score_table) <- c("score","condition")
  
  score_table$condition  <- as.character(score_table$condition)
  neworder       <- unique(score_table$condition)
  score_table    <- dplyr::arrange(transform(score_table,condition=factor(condition,levels=neworder)),condition)
  
  fig <- plot_ly(score_table, y = ~score, x = ~condition, type = 'scatter',
                 mode = "markers", marker = list(size = 10,color = "rgba(255, 182, 193, .9)"))
  fig <- fig %>% layout(xaxis = x, yaxis = y,font="Roboto")
  
  fig <- config_figure(fig,plot_type,plot_width,plot_height)
  fig <- fig %>% layout(showlegend = FALSE)
  
  return(fig)
}
