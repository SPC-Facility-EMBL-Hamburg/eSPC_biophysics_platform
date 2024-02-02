output$download_raw_signal_plot <- downloadHandler(filename = function() {
  "raw_Signal.csv"},content = function(file) {
    
    fluo_m <- make_df4plot(mst$signal,mst$concs,mst$times)
    colnames(fluo_m) <- c("Time","Replicate_Concentration","Signal","Concentration")
    
    write.csv(fluo_m,file,row.names = F,quote = F)
  })

output$download_lig_signal_plot <- downloadHandler(filename = function() {
  "ligand_Signal.csv"},content = function(file) {
    
    df <- get_df_min_max_fluo(mst$F_cold,mst$concs,mst$expID_vector)
    colnames(df) <- c("concentration","signal","experimentID")
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_init_fluo_plot <- downloadHandler(filename = function() {
  "initial_Fluorescence_versus_[Ligand].csv"},content = function(file) {
    
    df <- data.frame(mst$concs,mst$F_cold,mst$expID_vector)
    colnames(df) <- c("Concentration","InitialFluorescence","ExperimentID")
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_fnorm_plot <- downloadHandler(filename = function() {
  "Fnorm_versus_[Ligand].csv"},content = function(file) {
    
    df <- data.frame(mst$concs,mst$F_norm,mst$expID_vector)
    colnames(df) <- c("Concentration","Fnorm","ExperimentID")
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_fitting_plot <- downloadHandler(filename = function() {
  "fitted_Curve.csv"},content = function(file) {
    
    fluo_fit_data <- fluo_fit_data()
    
    df <- get_fitting_plot(fluo_fit_data$uniqueExpIDs,fluo_fit_data$fit_obj,
                           fluo_fit_data$pconcs)[,c(1,3,6)]
    
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_params_table <- downloadHandler(filename = function() {
  "Fitted_Parameters.csv"},content = function(file) {
    
    fluo_fit_data <- fluo_fit_data()
    df <- get_parameters_table(fluo_fit_data,input$model_selected) 

    write.csv(df,file,row.names = F,quote = F)
  })


