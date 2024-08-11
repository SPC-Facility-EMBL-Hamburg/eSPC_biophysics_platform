output$download_raw_signal_plot <- downloadHandler(filename = function() {
  paste0('raw_signal_',reactives$exportName,'.csv') 
  },content = function(file) {
    
    fluo_m <- make_df4plot(mst$signal,mst$concs,mst$times)
    colnames(fluo_m) <- c("Time","Replicate_Concentration","Signal","Concentration")
    
    write.csv(fluo_m,file,row.names = F,quote = F)
  })

output$download_lig_signal_plot <- downloadHandler(filename = function() {
  paste0('ligand_signal_',reactives$exportName,'.csv') 
  },content = function(file) {
    
    df <- get_df_min_max_fluo(mst$F_cold,mst$concs,mst$expID_vector)
    colnames(df) <- c("concentration","signal","experimentID")
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_init_fluo_plot <- downloadHandler(filename = function() {
  paste0('fluo_init_vs_ligand_conc_',reactives$exportName,'.csv') 
  },content = function(file) {
    
    df <- data.frame(mst$concs,mst$F_cold,mst$expID_vector)
    colnames(df) <- c("Concentration","InitialFluorescence","ExperimentID")
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_fnorm_plot <- downloadHandler(filename = function() {
  paste0('fnorm_vs_ligand_conc_',reactives$exportName,'.csv') 
  },content = function(file) {
    
    df <- data.frame(mst$concs,mst$F_norm,mst$expID_vector)
    colnames(df) <- c("Concentration","Fnorm","ExperimentID")
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_fitting_plot <- downloadHandler(filename = function() {
  paste0('fitted_curves_',reactives$exportNameFittings,'_',Sys.Date(),'.csv')
  },content = function(file) {
    
    fluo_fit_data <- fluo_fit_data()

    df <- get_fitting_plot(fluo_fit_data$uniqueExpIDs,fluo_fit_data$fit_obj,
                           fluo_fit_data$pconcs)[,-4]
    
    colnames(df) <- c('ligand_conc','measured_signal','fitted_signal','expID','protein_conc')
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_params_table <- downloadHandler(filename = function() {
  paste0('fitted_params_',reactives$exportNameFittings,'_',Sys.Date(),'.csv')
  },content = function(file) {
    
    fluo_fit_data <- fluo_fit_data()
    df <- get_parameters_table(fluo_fit_data,input$model_selected) 

    write.csv(df,file,row.names = F,quote = F)
  })
