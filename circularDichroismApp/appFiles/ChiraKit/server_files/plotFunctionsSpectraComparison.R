# Scripts to generate the figures of the comparison panel

# Input:
# dataframe with columns 'x', 'y' and 'err'
# 'x' is the wavelength, 'y' the average CD signal, and 'err' the standard deviation
add_average_spectrum_to_fig <- function(df,hex_color,lbl,fig) {
  
  fig    <- fig %>% add_trace(data=df,x=~x,y=~y,
                              type = 'scatter', 
                              name=lbl,
                              color=I(hex_color),
                              showlegend = TRUE,
                              marker = list(size = 9,
                                            line = list(width = 0)))
  
  df_temp <- df[!is.na(df$err),]
  
  if (nrow(df_temp) > 0) {
    
    fig <- fig %>% add_trace(data = df_temp,x = ~x, y = ~y,
                             error_y =  list(array = ~err),
                             color=I(hex_color),
                             opacity=0.75,showlegend = F)
    
  }
  
  return(fig)
  
}

# Input:
# means - matrix with the averages
# sds   - matrix with the standard deviation
# lbls  - labels of the previous two matrices
# wl    - wavelength vector, determines the number of rows in 'means' and 'sds'
# selectedUnits - the working units used to create the comparison dataset
plot_average_spectra <- function(
    means,sds,lbls,wl,selectedUnits,
    plot_width=12,plot_height=8,plot_type='svg',axis_size=16) {
  
  colors <- getPalette(length(lbls)) 
  
  fig <- plot_ly()
  
  for (i in 1:length(lbls)) {
    
    tempDf <- data.frame(x=wl,y=means[,i],err=sds[,i])
    fig    <- add_average_spectrum_to_fig(tempDf,colors[i],lbls[i],fig)
    
  }
  
  yRange <- find_y_plotting_range(means,sds)
  
  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(minWL <- min(wl) - 5, max(wl) + 5),
            showgrid = F)
  
  y <- list(title = workingUnits2ProperLabel(selectedUnits),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),
            showgrid = F, range = c(yRange$yMin,yRange$yMax))
  
  fig <- fig %>% layout(showlegend = TRUE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))
  
  fig <- configFig(fig,paste0("CDspectraAvg_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  
  return(fig)
}

# Input:
# dS    - matrix with the difference spectra
# dsE   - matrix with the associated errors (calculated using error propagation rules)
# dsLbl - labels of the comparison. One label per column in 'dS'
# wl    - wavelength vector, determines the number of rows in 'dS' and 'dsE' 
plot_difference_spectra <- function(
    dS,dSe,dsLbl,wl,wU,
    plot_width=12,plot_height=8,plot_type='svg',axis_size=16,
    reference="No reference") {
  
  plot_list   <- list()
  colors      <- getPalette(length(dsLbl)) 

  x <- list(title = "Wavelength (nm)",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),range = c(min(wl) - 5, max(wl) + 5),
            showgrid = F)

  yRange <- find_y_plotting_range(dS,dSe)
  
  y <- list(title = workingUnits2ProperLabel(wU),
            titlefont = list(size = axis_size), tickfont = list(size = axis_size),
            range = c(yRange$yMin,yRange$yMax),
            showgrid = F)
  
  # counter i and j
  # 'i' iterates over the mean matrix, 'j' is used to fill the list 'plot_list'
  i     <- 0
  j     <- 0
  
  cond0 <- reference == "No reference"

  for (cat in dsLbl) {
    
    i      <- i + 1
  
    if (reference_is_in_label(reference,cat)) {
      
      tempDf <- data.frame(x=wl,y=dS[,i],err=dSe[,i])
      fig    <- plot_ly()
      fig    <- add_average_spectrum_to_fig(tempDf,colors[i],cat,fig)
      fig    <- add_layout_to_subplot(fig,x,y,cat,length(dsLbl),axis_size,F)
      
      j      <- j + 1
      plot_list[[j]] <- fig
    } 
    
  }
  
  fig <- plot_list_to_fig(
    paste0('CD_difference_spectra_',Sys.Date()),
    plot_list,dsLbl,
    axis_size,plot_type,plot_width,plot_height)
  
  return(fig)
  
}

# Input:
# means  - matrix with the average spectra
# labels - labels, one per column in 'means'
# wU - the working units used to create the comparison dataset
plot_similarity_spectra <- function(
    means,labels,wU,
    plot_width=12,plot_height=8,plot_type='svg',axis_size=16,
    reference='No reference') {
  
  plot_list   <- list()
  nMeans      <- length(labels)
  cond0       <- reference == 'No reference'
  
  c <- 0
  for (i in 1:(nMeans-1)) { 
    
    x <- list(title = paste0(labels[i],'; ',workingUnits2ProperLabel(wU)),titlefont = list(size = axis_size), 
              tickfont = list(size = axis_size),
              showgrid = F)
    
    cond1 <- reference == labels[i] 
    
    for (ii in (i+1):nMeans) {
      
      y <- list(title = paste0(labels[ii],'; ',workingUnits2ProperLabel(wU)),
                titlefont = list(size = axis_size), tickfont = list(size = axis_size),
                showgrid = F)
      
      cond2 <- reference == labels[ii] 
      
      if (cond0 | cond1 | cond2) {
        
        c   <- c + 1
        fig <- plot_ly()
        tempDf <- data.frame(x=means[,i],y=means[,ii])
        
        model       <- lm(y ~ x , data=tempDf)
        adjRsquared <- signif(summary(model)$adj.r.squared,4)
        
        fig    <- fig %>% add_trace(data=tempDf,x=~x,y=~y,
                                    type = 'scatter', 
                                    marker = list(size = 9,
                                                  line = list(width = 0)))
        
        fig    <- add_layout_to_subplot(fig,x,y,paste("Adj. R-squared: ", adjRsquared),nMeans,axis_size,F)
        
        plot_list[[c]] <- fig
        
      }
    }
    
  }
  
  nrows <- min(length(plot_list),3)
  
  fig <- plot_list_to_fig(
    paste0('CD_similarity_spectra_',Sys.Date()),
    plot_list,paste("Adj. R-squared: ", adjRsquared),
    axis_size,plot_type,plot_width,plot_height,nrows,F)
  
  return(fig)
  
}

# Input:
# comparison_labels - vector of strings, has the name of the comparison that was done
# distances_lst     - list of vectors containing the (normalised) euclidean distances of the comparisons
plot_distances <- function(comparison_labels,distances_lst,
                           plot_style = 'boxplot',
                           plot_width=12,plot_height=8,plot_type='svg',axis_size=16,
                           reference="No reference") {
  
  fig <- plot_ly(orientation='h')
  
  i <- 1
  for (cat in comparison_labels) {
    
    cond1 <- reference_is_in_label(reference,cat,' versus ')
    vec   <- strsplit(cat,' versus ')[[1]] 
    cond2 <- vec[1] == vec[2] # show always intra-group distances
    
    if (cond1  | cond2 ) {
      
      dists <- distances_lst[[i]]
      
      if (plot_style == 'boxplot') {
        
        fig   <- fig %>% add_trace(x=dists,y=rep(cat,length(dists)),
                                   boxpoints = "outliers",
                                   type = "box")
        
      } else {
        
        fig   <- fig %>% add_trace(x=dists,y=rep(cat,length(dists)),
                                   type = "scatter",
                                   marker = list(size = 9))
        
      }
    }

    i <- i + 1
    
  }
  
  x <- list(title = "Euclidean distance",titlefont = list(size = axis_size), 
            tickfont = list(size = axis_size),showgrid = F)
  
  y <- list(tickfont = list(size = axis_size),
            showgrid = F)
  
  fig <- fig %>% layout(showlegend = FALSE,xaxis = x, yaxis = y,font="Roboto",
                        legend = list(font = list(size = axis_size-1)))
  
  fig <- configFig(fig,paste0("CDspectraDistances_",strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  
  return(fig)
  
}

# Requires:
# 'distance_matrix': the distance matrix, with at least filled values in the upper right triangle.
# 'groups'         : vector with the spectra groups. e.g., c('WT','WT','Mut','Mut')    
plot_dendogram <- function(distance_matrix,groups,
                           plot_width=12,plot_height=8,
                           plot_type='svg',axis_size=16) {
  
  m              <- as.matrix(distance_matrix)
  upper_triangle <- m[upper.tri(m)]
  full_matrix    <- matrix(0, nrow(m), nrow(m))
  
  # Fill the upper triangle of the full matrix with the given values
  full_matrix[upper.tri(full_matrix, diag = FALSE)] <- upper_triangle
  
  # Make the matrix symmetric by copying the upper triangle to the lower triangle
  full_matrix <- full_matrix + base::t(full_matrix)

  model       <- hclust(as.dist(full_matrix), method = "ward.D2")
  
  dhc       <- as.dendrogram(model)
  dend_data <- dendro_data(dhc, type = "rectangle")
  
  segments <- dend_data$segments
  labels   <- dend_data$labels
  labels$y <- -max(abs(dend_data$segments$y))*0.2
  
  labels$label <- groups[model$order]
  
  maxChars <- max(sapply(labels$label, nchar))
  
  # Left-align labels
  labels$label <- sapply(labels$label, function (x) {
    
    chars    <- nchar(x)
    reqChars <- maxChars - chars
    strToAdd <- paste0(rep(' ',1+reqChars),collapse = '')
  
    return(paste0(x,strToAdd))
      
  })
  
  p <- ggplot() + 
    geom_segment(data = segments, aes(x = x, y = y, xend = xend, yend = yend)) + 
    geom_text(data = labels, aes(x = x, y = y, label = label),
              size = axis_size / 3.4) +
    coord_flip() + 
    scale_y_reverse(limits = c(NA,labels$y*2.2)) +
    theme_classic(base_size = axis_size) +
    labs(title = "Hierarchical Clustering Dendrogram (Ward.D2)", 
         x = "Spectra", y = "Height")+
    theme(axis.ticks.y = element_blank(),
          axis.text.y  = element_blank()) 

  
  fig <- ggplotly(p)
  
  fig <- configFig(fig,
                   paste0("CDspectraTree_",
                          strsplit(as.character(Sys.time())," ")[[1]][1]),
                   plot_type,plot_width,plot_height)
  
  return(fig)
  
}





