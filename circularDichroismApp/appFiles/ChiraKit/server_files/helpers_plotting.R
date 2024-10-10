# Helpers functions called in the script plotFunctions.R

# Based on the input units, e.g. millidegrees, decide the text
# that will be shown in the y-axis of the plots
workingUnits2ProperLabel <- function(choice) {
  
  if (choice == 'millidegrees')   return("Millidegrees (m°)" )
  if (choice == 'degrees')       return("Degrees (°)"      )
  
  if (choice == 'molarEllipticity')  {
    return("Molar ellipticity (deg\u2219cm<sup>2</sup>\u2219dmol<sup>-1</sup>)")
  }         

  if (choice == 'meanUnitMolarEllipticity')  {
    return("MRE (deg\u2219cm<sup>2</sup>\u2219dmol<sup>-1</sup>\u2219chromophores<sup>-1</sup>)")
  }   
  
  if (choice == 'absorbance')                 return("Differential absorbance (\u0394A)"     )
  if (choice == 'milliabsorbance')            return("Differential milliabsorbance (\u0394mA)")
  
  if (choice == 'molarExtinction')  {
    return("Molar extinction (L\u2219mol<sup>-1</sup>\u2219cm<sup>-1</sup>)")
  }         
  
  if (choice == 'meanUnitMolarExtinction')  {
    return("\u0394\u03B5 (L\u2219mol<sup>-1</sup>\u2219cm<sup>-1</sup>\u2219chromophores<sup>-1</sup>)")
  }
  
  return('Signal')
}

# - remove unnecessary options shown in the plotly figures
# - set the figure size and name
configFig <- function(fig,nameForDownload,plot_type,plot_width,plot_height) {
  fig %>%  config(
    toImageButtonOptions = list(
      format = plot_type,
      filename = nameForDownload,
      width = plot_width * 50,
      height = plot_height * 50
    ), displaylogo = FALSE,
    modeBarButtonsToRemove = list('sendDataToCloud','hoverClosestCartesian',
                                  'hoverCompareCartesian','lasso2d','select2d'))
}

# Combine a list of plotly plots into one figure
plot_list_to_fig <- function(plot_name,plot_list,titleText,axis_size,plot_type,plot_width,plot_height,
                             nrows=2,shareAxis = TRUE) {
  
  i <- length(plot_list)
  
  if (i > 1) {
    
    if (shareAxis) {
      
      fig <- subplot(plot_list,nrows = nrows,margin = c(0.03,0.03,0.1,0.2),
                     shareY = TRUE, shareX = TRUE)
      
    } else {
      
      fig <- subplot(plot_list,nrows = nrows,margin = c(0.03,0.03,0.1,0.2),
                     shareY = FALSE, shareX = FALSE,titleX = TRUE,titleY = TRUE)
      
    }
    

  } else {
    
    fig <- plot_list[[1]] %>%  layout(
      title=list(text=titleText,font = list(size = axis_size)))
    
  }
  
  fig <- configFig(fig,plot_name,plot_type,plot_width,plot_height)
  return(fig)
}

# This function will:
#  add the x axis title to the plot     argument 'xaxis'
#  add the y axis title to the plot     argument 'yaxis'
#  add the legend                       argument 'leg'
#  at a certain position                argument 'x_pos_annot'
#  if there is more than one condition  argument 'tot_cond'
# set the size font                     argument 'axis_size'
add_layout_to_subplot <- function(fig,xaxis,yaxis,leg,tot_cond,axis_size,
                                  showlegend=TRUE) {
  
  fig <- fig %>% 
    layout(xaxis = xaxis,yaxis=yaxis,
           legend = list(font = list(size = axis_size-1)),
           annotations = list(
             x = 0.5 , y = 1.12, 
             text = ifelse(tot_cond > 1,leg,""), showarrow = F, 
             xref='paper', yref='paper',font = list(size = axis_size)),
           showlegend=showlegend)
  
  return(fig)
}

# obtain the name of the x-axis
get_x_axis_label <- function(df) {
  
  temperatureBased <- any(grepl('temperature',colnames(df),ignore.case = T))
  chemBased        <- any(grepl('chem_conc',colnames(df),ignore.case = T))
  customBased      <- !chemBased & !temperatureBased
  
  xaxisLabel <- 'Temperature (°C)'
  
  if (chemBased)   xaxisLabel <- '[Denaturant agent] (M)'
  if (customBased) xaxisLabel <- colnames(df)[2]
  
  return(xaxisLabel)
}

# Add a column to a dataframe named measurement_factor
# This column will be the same as the one named temperature, chem_conc, or the second column
add_measurement_factor_column <- function(df) {
  
  df_colnames <- colnames(df)
  
  if ('temperature' %in% df_colnames) {
    df$measurement_factor <- df$temperature
  } else if ('chem_con' %in% df_colnames) {
    df$measurement_factor <- df$chem_con
  } else {
    # Second column as default
    df$measurement_factor <- df[,2] 
  }
  
  return(df)
}

find_y_plotting_range <- function(y_values_matrix,y_err_matrix) {
  
  yMax <- max(y_values_matrix + y_err_matrix)
  yMin <- min(y_values_matrix - y_err_matrix)
  
  totalRange <- yMax - yMin
  
  factor <- 0.1
  
  return(list('yMin'=yMin - totalRange*factor,'yMax'=yMax + totalRange*factor))
}

# Find if the label contains the selected reference
reference_is_in_label <- function(reference,label,separator=' - ') {
  
  cond0 <- reference == "No reference"
  
  # Verify if the reference is in the selected label 
  cond1 <- grepl(paste0(separator,reference),label) 
  cond2 <- grepl(paste0(reference,separator),label) 
  # To avoid problems such as selecting '1' as reference, but also '10' is in the labels
  
  vec   <- strsplit(label,paste0(separator,reference))[[1]]
  vec   <- vec[vec != ""]
  cond3 <- length(vec) == 1   
  
  vec   <- strsplit(label,paste0(reference,separator))[[1]]
  vec   <- vec[vec != ""]
  cond4 <- length(vec) == 1 

  return(cond0 | (cond1 & cond3) | (cond2 & cond4))
  
}





