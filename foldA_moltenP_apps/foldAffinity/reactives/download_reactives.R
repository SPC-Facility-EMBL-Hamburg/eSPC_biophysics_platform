
## Download section ##

output$download_params_table        <-   downloadHandler(
  filename = function() {paste0("Fitted_Parameters_Fluorescence_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(fluo_fit_data())
    write.csv(get_params_df(fluo_fit_data()$fit_params,dsf$concentrations),
              file,row.names = F,quote = F)
})
  
output$download_params_errors_table <-   downloadHandler(
  filename = function() {paste0("Fitted_Parameters_Relative_Errors_Fluorescence_",Sys.Date(),".csv")},
  content  = function(file) {
    
    write.csv(get_params_errors_df(fluo_fit_data()$fit_params,dsf$concentrations,fluo_fit_data()$fit_errors),
              file,row.names = F,quote = F)
})
  
output$download_params_table_fu <-   downloadHandler(
  
  filename = function() {paste0("Fitted_Parameters_UnfoldedFraction_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(ist_fit_data_pred())
    write.csv(bind_params2df(dsf$bind_params,dsf$isothermal_ts),file,row.names = F,quote = F)
    
  })
  
output$download_params_errors_table_fu <-   downloadHandler(
  
  filename = function() {paste0("Fitted_Parameters_Relative_Errors_UnfoldedFraction_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(ist_fit_data_pred())
    
    df <- bind_err_params2df(dsf$bind_params,dsf$bind_errors,dsf$isothermal_ts,
                             dsf$bind_ci95_asymmetric_low, dsf$bind_ci95_asymmetric_up)
    
    colnames(df) <- c("Ku relative error","Kd relative error",
                      "Temperature","AsymmetricCI95_lower","AsymmetricCI95_upper")
    
    write.csv(df,file,row.names = F,quote = F)
    
  })

output$download_signal_plot   <-   downloadHandler(
  
  filename = function() { paste0("fluo_signal_data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(fluo_fit_data())  
    df <- make_df4plot(dsf$fluo,dsf$concentrations,dsf$temps)
    colnames(df) <- c("Temperature","ConcentrationID","Signal","Concentration (M)")
    write.csv(df,file,row.names = F,quote = F)
    
  })

output$download_signalDerivative_plot   <-   downloadHandler(
  filename = function() { paste0("fluo_signal_derivative_data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(fluo_fit_data()) 
    df <- make_df4plot(dsf$derivative,dsf$concentrations,dsf$temps)
    colnames(df) <- c("Temperature","ConcentrationID","Signal_FirstDerivative","Concentration (M)")
    write.csv(df,file,row.names = F,quote = F)
    
    })
  
output$download_initialSignal_plot   <- downloadHandler(
  
  filename = function() { paste0("initialSignal_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(fluo_fit_data()) 
    write.table(getInitialSignalDF(dsf$fluo,dsf$concentrations,dsf$temps),file,row.names = F,quote = F,col.names=F)
    
    })
  
output$download_fluo_fit_plot <- downloadHandler(
  
  filename = function() { paste0("fluo_signal_fitted_data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(fluo_fit_data())  
    df <- make_df4plot(dsf$fit_fluo_pred,dsf$concentrations,dsf$temps)
    colnames(df) <- c("Temperature","ConcentrationID","Signal","Concentration (M)")
    write.csv(df,file,row.names = F,quote = F)
    
  })
  

output$download_fraction_unfolded_plot     <-   downloadHandler(
  filename = function() {paste0("UnfoldedFraction_data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(ist_fit_data_pred())
    df           <- format_ist_data_exp(dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts)
    colnames(df) <- c("Temperature","Concentration (M)","Unfolded fraction","Legend")
    write.csv(df,file,row.names = F,quote = F)
    
    })
  
output$download_fraction_unfolded_fit_plot <-   downloadHandler(
  filename = function() {paste0("UnfoldedFraction_fitting_data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    req(ist_fit_data_pred())
    df <- format_ist_data_exp_and_pred(
      dsf$isothermal_data,dsf$concentrations,dsf$isothermal_ts,
      dsf$kd_models,dsf$kd_models_lower,dsf$kd_models_upper,dsf$kd_model_conc,
      dsf$bind_params,dsf$bind_errors,0,0)$data_pred[,c(1,2,3)]
    
    colnames(df) <- c("Temperature","Concentration (M)","Unfolded fraction")
    write.csv(df,file,row.names = F,quote = F)
    
    })
  
# Download the plots
output$download_fit_plots = downloadHandler(
  filename = 'fitting_plots.zip',
  content  = function(file){
    
    # Set temporary working directory
    owd <- setwd( tempdir())
    on.exit( setwd( owd))
    
    withBusyIndicatorServer("download_fit_plots",{
      total_plots <- ceiling(ncol(dsf$fluo)/global_chunck_n)
      
      iter <- 0
      fns <- c()
      for (selected in 1:total_plots) {
        iter <- iter + 1
        real_data  <- fluo_fit_data()$fluo_fit_real
        model_data <- fluo_fit_data()$fluo_fit_pred
        fname <- paste0("fitting_plot_",iter,".png")
        fns   <- c(fns,fname)
        plot_fluorescence_fit(real_data,model_data,selected) %>% 
          ggpubr::ggexport(filename = fname,res=240,width = 1800,height = 1800)
      }
    })
    
    # Zip them up
    zip( file, fns)
  }
)

output$download_session   <-   downloadHandler(
  filename = function() { paste0("fAffinitySession_",Sys.Date(),".json")},
  content  = function(file) {
    
    dsf$export_JSON_file(file)
    
  })

