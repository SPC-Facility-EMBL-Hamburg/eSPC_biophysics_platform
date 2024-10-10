observeEvent(input$launchReferenceGquad,{
  
  reactives$GQ_ref_load <- FALSE
  
  if (!input$useDefaultReferenceSetGQuad) {
    
    file1   <- input$cd_ref_spectra_GQ
    file2   <- input$cd_ref_sec_params_GQ
    file3   <- input$cd_ref_ter_params_GQ
    
    if (is.null(file1)) {
      shinyalert(
      text = "Please load the reference spectra: A CSV file with the wavelength 
      data in the first column and the 
      spectral data in the subsequent columns (in molar extinction units). 
      The files should have a header and the spectra names should match the column names
      of the reference spectra file.",
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
      return(NULL) 
    }
    
    if (is.null(file2) & is.null(file3)) {
      shinyalert(
        text = "Please load the secondary and/or tertiary parameters: CSV files with the
        spectra names in the first column and the fractions data in the subsequent columns. 
        The files should have a header and the spectra names should match the column names
        of the reference spectra file.",
        type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
      return(NULL) 
    }
    
    # In molar extinction units!
    gQuadRefPyClass$load_data(file1$datapath)
    gQuadRefPyClass$spectraNames <- trimws(gQuadRefPyClass$spectraNames)
    
    if (!is.null(file2)) {
      reactives$secondary_parameters <- read.csv(file2$datapath,check.names = F)
    }
    
    if (!is.null(file3)) {
      reactives$tertiary_parameters <- read.csv(file3$datapath,check.names = F)
    }
    
  } else {
    
    # Load structure default parameters
    reactives$secondary_parameters <- read.csv('www/default_secondary_parameters.csv',check.names = F)
    reactives$tertiary_parameters  <- read.csv('www/default_tertiary_parameters.csv', check.names = F)
    
    # In molar extinction units!
    gQuadRefPyClass$load_data('www/default_spectra_GQuadruplex.csv')
    gQuadRefPyClass$spectraNames <- trimws(gQuadRefPyClass$spectraNames)
    
  }
  
  output$cdSpectraGQ <- renderPlotly({
    
    plot_gQuadruplexReferences(
      gQuadRefPyClass$wavelength,gQuadRefPyClass$signalInput,
      gQuadRefPyClass$spectraNames)
    
  })
  
  if (!is.null(reactives$secondary_parameters)) {
    output$secondary_params_GQ <- renderTable({reactives$secondary_parameters})
  }
  
  if (!is.null(reactives$tertiary_parameters)) {
    output$tertiary_params_GQ <- renderTable({reactives$tertiary_parameters})
  }
  
  reactives$GQ_ref_load <- TRUE
  
})

# Reactive expression to convert input to numeric
pca_num_clusters <- reactive({
  if (input$GQnClust == 'Auto') {
    return(-1)
  } else {
    return(as.numeric(input$GQnClust))
  }
})

reference_cd_plot <- reactive({
  
  req(reactives$GQ_ref_load)
  signal <- base::t(gQuadRefPyClass$signalInput)
  wl     <- gQuadRefPyClass$wavelength
  names  <- gQuadRefPyClass$spectraNames
  
  pca_num_clusters <- pca_num_clusters()
  
  reference_cd_plot <- plot_pca_analysis( signal,wl,names,PC1=input$PC_ha,PC2=input$PC_va,slider_pca_num_clusters = pca_num_clusters,
                                          text_graph_title = input$pca_graph_title,slider_graph_label_size = input$pca_graph_label,
                                          slider_axis_label_size = input$pca_axis_label,slider_point_size = input$pca_point_size,
                                          slider_arrow_width = input$pca_arrow_width,slider_ellipse_confidence = input$pca_ellipse_confidence / 100,
                                          checkbox_show_ellipse =   input$show_ellipse_pcaGQ, checkbox_show_legend = input$show_legend_pcaGQ)
  
  return(reference_cd_plot)
  
})

output$pca_results_GQ <- renderPlot({
  
  reference_cd_plot()
  
})

cluster_plot <- reactive({
  
  req(reactives$GQ_ref_load)
  signal <- base::t(gQuadRefPyClass$signalInput)
  names  <- gQuadRefPyClass$spectraNames
  
  pca_num_clusters <- pca_num_clusters()
  
  cluster_plot <- plot_cluster_analysis( signal,names,cluster_num_clusters = pca_num_clusters,
                                         cluster_graph_title = input$cluster_graph_title,
                                         cluster_text_size   = input$cluster_text_size / 10,
                                         axis_text_size      = input$cluster_text_size + 5)
  
  return(cluster_plot)
  
})

output$pca_clustering_GQ <- renderPlot({
  
  cluster_plot()
  
})

observeEvent(input$launchSamplesPCAGquad,{
  
  if (!reactives$GQ_ref_load) {
    shinyalert(text = "Please load the references first.",
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    return(NULL)
  }
  
  reactives$GQ_sample_load <- FALSE

  exps <- cdAnalyzer$experimentNames
  
  if(length(exps)==0){
    shinyalert(text = "Please import sample data using the module '1.Import data'",
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    return(NULL)
  } 
  
  if(input$workingUnits != 'molarExtinction'){
    shinyalert(
      text = "Please choose 'Molar extinction' as the working units (Module '1.Import data')",
      type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
    return(NULL)
  } 
  
  signals <- cdAnalyzer$get_experiment_properties_modif('signalDesiredUnit')  # In molar extinction units!
  
  # Find which experiments have data
  sel_exps_ids <- c()
  
  counter <- 1
  for (s in signals) {
    nas <- sum(is.na(s))
    if (nas == 0) {
      sel_exps_ids <- c(sel_exps_ids,counter)
    }
    counter <- counter + 1
  }
  
  if(length(sel_exps_ids)==0){
    shinyalert(
    text = "No data available in 'Molar extinction' units. 
    Please import data and revise the Table '2. Parameters for molar ellipticity / extinction'. 
    For example, set the input units to 'Molar extinction'",
    type = "warning",closeOnEsc = T,closeOnClickOutside = T,html=T)
    return(NULL)
  } 
  
  relevantSpectra <- unlist(cdAnalyzer$get_experiment_properties('internalID')[sel_exps_ids])
  
  merged  <- get_signal_dfs_from_selected_spectra(relevantSpectra,cdAnalyzer)
  signal  <- as.matrix(merged[,-1],drop = FALSE)
  
  gQuadSamplePyClass$wavelength   <- np_array(merged[,1])
  gQuadSamplePyClass$signalInput  <- np_array(signal)
  gQuadSamplePyClass$name         <- ''
  
  # Filter according to the reference wavelength range
  minWL <- min(gQuadRefPyClass$wavelength)
  maxWL <- max(gQuadRefPyClass$wavelength)
  
  gQuadSamplePyClass$signalInput <- filter_matrix_by_vector(gQuadSamplePyClass$signalInput,gQuadSamplePyClass$wavelength,minWL,maxWL)
  gQuadSamplePyClass$wavelength  <- filter_vector_by_values(gQuadSamplePyClass$wavelength,minWL,maxWL)
  
  gQuadSamplePyClass$spectraNames <- relevantSpectra

  output$cdSpectraGQ_samples <- renderPlotly({
    
    # We reuse the plotting function for the references
    plot_gQuadruplexReferences(
      gQuadSamplePyClass$wavelength,gQuadSamplePyClass$signalInput,
      gQuadSamplePyClass$spectraNames)
    
  })

  reactives$GQ_sample_load <- TRUE

})

# Reactive expression to convert input to numeric
pca_num_clusters_samples <- reactive({
  if (input$GQnClust_samples == 'Auto') {
    return(-1)
  } else {
    return(as.numeric(input$GQnClust_samples))
  }
})

samples_pca_plot <- reactive({
  
  req(reactives$GQ_sample_load)
  
  signal <- base::t(gQuadSamplePyClass$signalInput)
  wl     <- gQuadSamplePyClass$wavelength
  names  <- gQuadSamplePyClass$spectraNames
  
  pca_num_clusters <- pca_num_clusters_samples()
  
  samples_pca_plot <- plot_pca_analysis(
    signal,wl,names,PC1=input$PC_ha_samples,PC2=input$PC_va_samples,
    slider_pca_num_clusters   = pca_num_clusters,
    text_graph_title          = input$pca_graph_title_samples,
    slider_graph_label_size   = input$pca_graph_label_samples,
    slider_axis_label_size    = input$pca_axis_label_samples,
    slider_point_size         = input$pca_point_size_samples,
    slider_arrow_width        = input$pca_arrow_width_samples,
    slider_ellipse_confidence = input$pca_ellipse_confidence_samples / 100,
    checkbox_show_ellipse     = input$show_ellipse_pcaGQ_samples, 
    checkbox_show_legend      = input$show_legend_pcaGQ_samples)
  
  return(samples_pca_plot)
  
})

output$pca_results_GQ_samples <- renderPlot({
  
  samples_pca_plot()
  
})

combined_pca_plot <- reactive({
  
  req(reactives$GQ_sample_load)
  
  signal1 <- base::t(gQuadSamplePyClass$signalInput)
  wl1     <- gQuadSamplePyClass$wavelength
  names1  <- gQuadSamplePyClass$spectraNames
  
  signal2 <- base::t(gQuadRefPyClass$signalInput)
  wl2     <- gQuadRefPyClass$wavelength
  names2  <- gQuadRefPyClass$spectraNames
  
  wl      <- intersect(wl1,wl2)

  signal1 <- signal1[,wl1 %in% wl2]
  signal2 <- signal2[,wl2 %in% wl1]

  signal  <- rbind(signal1,signal2)
  names   <- c(names1,names2)

  pca_num_clusters <- pca_num_clusters_samples()
  
  combined_pca_plot <- plot_pca_analysis(
    signal,wl,names,PC1=input$PC_ha_samples,PC2=input$PC_va_samples,
    slider_pca_num_clusters   = pca_num_clusters,
    text_graph_title          = input$pca_graph_title_samples,
    slider_graph_label_size   = input$pca_graph_label_samples,
    slider_axis_label_size    = input$pca_axis_label_samples,
    slider_point_size         = input$pca_point_size_samples,
    slider_arrow_width        = input$pca_arrow_width_samples,
    slider_ellipse_confidence = input$pca_ellipse_confidence_samples / 100,
    checkbox_show_ellipse     = input$show_ellipse_pcaGQ_samples, 
    checkbox_show_legend      = input$show_legend_pcaGQ_samples)
  
  return(combined_pca_plot)
  
})


output$pca_results_GQ_combined <- renderPlot({
  
  combined_pca_plot()
  
})

cluster_plot_samples <- reactive({
  
  req(reactives$GQ_sample_load)
  signal <- base::t(gQuadSamplePyClass$signalInput)
  names  <- gQuadSamplePyClass$spectraNames
  
  pca_num_clusters <- pca_num_clusters_samples()
  
  cluster_plot_samples <- plot_cluster_analysis( signal,names,cluster_num_clusters = pca_num_clusters,
                                         cluster_graph_title = input$cluster_graph_title_samples,
                                         cluster_text_size   = input$cluster_text_size_samples / 10,
                                         axis_text_size      = input$cluster_text_size_samples + 5)
  return(cluster_plot_samples)
  
})

output$pca_clustering_GQ_samples <- renderPlot({
  
  cluster_plot_samples()
  
})

cluster_plot_combined <- reactive({
  
  req(reactives$GQ_sample_load)

  signal1 <- base::t(gQuadSamplePyClass$signalInput)
  wl1     <- gQuadSamplePyClass$wavelength
  names1  <- gQuadSamplePyClass$spectraNames
  
  signal2 <- base::t(gQuadRefPyClass$signalInput)
  wl2     <- gQuadRefPyClass$wavelength
  names2  <- gQuadRefPyClass$spectraNames
  
  signal1 <- signal1[,wl1 %in% wl2]
  signal2 <- signal2[,wl2 %in% wl1]
  
  signal  <- rbind(signal1,signal2)
  names   <- c(names1,names2)
  
  pca_num_clusters <- pca_num_clusters_samples()
  
  cluster_plot_combined <- plot_cluster_analysis( signal,names,cluster_num_clusters = pca_num_clusters,
                                                 cluster_graph_title = input$cluster_graph_title_samples,
                                                 cluster_text_size   = input$cluster_text_size_samples / 10,
                                                 axis_text_size      = input$cluster_text_size_samples + 5)
  return(cluster_plot_combined)
  
})

output$pca_clustering_GQ_combined <- renderPlot({
  
  cluster_plot_combined()
  
})

observeEvent(input$launchSVD_GQ,{
  
  req(reactives$GQ_sample_load)
  
  cd_reference           <- base::t(gQuadRefPyClass$signalInput)
  rownames(cd_reference) <- gQuadRefPyClass$spectraNames
  colnames(cd_reference) <- gQuadRefPyClass$wavelength
  
  do_sec_params <- !is.null(reactives$secondary_parameters)
  do_ter_params <- !is.null(reactives$tertiary_parameters)
  
  if (do_sec_params) {
    parameters        <- reactives$secondary_parameters
    parameters_sorted <- parameters[order(match(parameters[,1], gQuadRefPyClass$spectraNames)),][,-1]
  }
  
  if (do_ter_params) {
    parameters_ter        <- reactives$tertiary_parameters
    parameters_sorted_ter <- parameters_ter[order(match(parameters_ter[,1], gQuadRefPyClass$spectraNames)),][,-1]
  }

  wlSample <- gQuadSamplePyClass$wavelength
  
  if (!identical(wlSample,gQuadRefPyClass$wavelength)) {
    shinyalert(text = "The estimation was cancelled.
               Please use the same wavelength data for the reference and the sample datasets.",
               type = "warning",closeOnEsc = T,closeOnClickOutside = T,
               html=T)
    return(NULL)
  }
  
  cd_sample     <- base::t(gQuadSamplePyClass$signalInput)
  names_sample  <- gQuadSamplePyClass$spectraNames
  
  est_results_all     <- list()
  est_results_all_ter <- list()
  
  for (i in 1:nrow(cd_sample)) {
    name                 <- names_sample[i]
    test_spectra         <- matrix(cd_sample[i,],nrow = 1)
    
    if (do_sec_params) {
      res                  <- svd_spectra_analysis(cd_reference,parameters_sorted,test_spectra,name,tertiary = F)
      est_results_all[[i]] <- res
    }

    if (do_ter_params) {
      res                      <- svd_spectra_analysis(cd_reference,parameters_sorted_ter,test_spectra,name,tertiary = T)
      est_results_all_ter[[i]] <- res
    }
    
  }
  
  if (do_sec_params) {
    sec_est_all <- lapply(est_results_all, function(x) x$predicted)
    sec_est_all <- data.frame(do.call(rbind,sec_est_all),check.names = F)
    sec_est_all$name <- names_sample
    output$fitted_secondary_str_GQ <- renderTable({sec_est_all})
  }
  
  if (do_ter_params) {
    sec_est_all_ter      <- lapply(est_results_all_ter, function(x) x$predicted)
    sec_est_all_ter      <- data.frame(do.call(rbind,sec_est_all_ter),check.names = F)
    sec_est_all_ter$name <- names_sample  
    output$fitted_tertiary_str_GQ  <- renderTable({sec_est_all_ter})
  }

})


