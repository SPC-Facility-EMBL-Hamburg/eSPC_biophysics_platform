# Download the plots
output$download_plots = downloadHandler(
  filename = paste0('dls_plots_',Sys.Date(),'.zip'),
  content = function(file){
    
    # Set temporary working directory
    owd <- setwd( tempdir())
    on.exit( setwd( owd))
    
    withBusyIndicatorServer("download_plots",{

      fileExt <- paste0(".",input$plot_type)
      
      fns <- c()

      req(dlsDataUpdated())
      data <- estimatedContributions()
      
      # TODO - add volume distribution using MIE theory
      for (plotType in c('rainbow','gray','dual','histogram')) {
        
        fn  <- paste0('Rh_distribution_intensity_',plotType,fileExt)
        fns <- c(fns,fn)
       
        plot_type_hr_distribution(plotType,data,input$plot_axis_size,
                                  input$overlapFactor,
                                  input$rowOrder,
                                  input$peak_filter) %>%
          ggsave(file=fn,width = 10,height = 10,dpi = 240)
        
        # Detect the mapping for the average
        mapping <- tryCatch(
          hot_to_r(input$tableForAvg),
          error=function(e) return(NULL)
        )
        
        # Plot the average, if possible
        if (!is.null(mapping)) {
          fn  <- paste0('Average_Rh_distribution_intensity_',plotType,fileExt)
          fns <- c(fns,fn)
          plot_type_hr_distribution(plotType,data,input$plot_axis_size,
                                    input$overlapFactor,
                                    input$rowOrder,
                                    input$peak_filter,
                                    TRUE,mapping) %>%
            ggsave(file=fn,width = 10,height = 10,dpi = 240)
        }
      }
      
      # Zip them up
      zip( file, fns)
    })

  }
)

output$download_fitted_hr_rowWise <- downloadHandler(filename = function() {
  
  alphaTerm <- ifelse(input$regMethod == 'fixedAlpha',
                      paste0('alpha',input$alphaRegTerm),
                      paste0('Lcurve_start',input$alphaRegTermLcurve_start,
                             '_stop',input$alphaRegTermLcurve_stop,
                             '_step',input$alphaRegTermLcurve_step))
  
  paste0("Estimated_Rh_Raynals_RowWise_",alphaTerm,'_',Sys.Date(),".csv")
},content = function(file) {
  
  req(estimatedContributions())
  estimatedContributions       <- estimatedContributions()
  estimatedContributions$value <- estimatedContributions$value*100
  
  colnames(estimatedContributions) <- c('Rh','variable','intensity')
  write.csv(estimatedContributions,file,row.names = F,quote = F)
})

output$download_fitted_hr_colWise <- downloadHandler(filename = function() {
  alphaTerm <- ifelse(input$regMethod == 'fixedAlpha',
                      paste0('alpha',input$alphaRegTerm),
                      paste0('Lcurve_start',input$alphaRegTermLcurve_start,
                             '_stop',input$alphaRegTermLcurve_stop,
                             '_step',input$alphaRegTermLcurve_step))
  
  paste0("Estimated_Rh_Raynals_ColWise_",alphaTerm,'_',Sys.Date(),".csv")
},content = function(file) {
  
  req(estimatedContributions())
  estimatedContributions       <- estimatedContributions()
  estimatedContributions$value <- estimatedContributions$value*100
  
  colnames(estimatedContributions) <- c('Rh','variable','intensity')
  
  write.csv(dcast(estimatedContributions,Rh~variable),file,row.names = F,quote = F)
})

output$download_fitted_autocorrelation_rowWise <- downloadHandler(filename = function() {

  alphaTerm <- ifelse(input$regMethod == 'fixedAlpha',
                      paste0('alpha',input$alphaRegTerm),
                      paste0('Lcurve_start',input$alphaRegTermLcurve_start,
                             '_stop',input$alphaRegTermLcurve_stop,
                             '_step',input$alphaRegTermLcurve_step))
  
  paste0("Fitted_Autocorrelation_Raynals_RowWise_",alphaTerm,'_',Sys.Date(),".csv")
  
},content = function(file) {
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  data <- dlsDataUpdated()
  
  dlsPredData     <- data$predictedData

  dlsPredData[,1] <- dlsPredData[,1] * 1e6 # Back to microseconds
  colnames(dlsPredData) <- c('time (microseconds)','variable','FittedAutocorrelation')
  
  write.csv(dlsPredData,file,row.names = F,quote = F)
})

output$download_fitted_autocorrelation_colWise <- downloadHandler(filename = function() {
  alphaTerm <- ifelse(input$regMethod == 'fixedAlpha',
                      paste0('alpha',input$alphaRegTerm),
                      paste0('Lcurve_start',input$alphaRegTermLcurve_start,
                             '_stop',input$alphaRegTermLcurve_stop,
                             '_step',input$alphaRegTermLcurve_step))
  
  paste0("Fitted_Autocorrelation_Raynals_ColWise_",alphaTerm,'_',Sys.Date(),".csv")},content = function(file) {
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  data <- dlsDataUpdated()
  
  dlsPredData  <- data$predictedData
  
  dlsPredData[,1] <- dlsPredData[,1] * 1e6 # Back to microseconds
  colnames(dlsPredData) <- c('time_microseconds','variable','FittedAutocorrelation')
  
  write.csv(dcast(dlsPredData,time_microseconds~variable, value.var = "FittedAutocorrelation"),file,row.names = F,quote = F)
})

# Downloadable csv of selected dataset ----
output$exportSimData <- downloadHandler(
  
  filename = function() {
    paste0("SimulatedDataDLS_",Sys.Date(),".zip")
  },
  content = function(zipFile) {
    
    df <- dcast(simulatedAutocorrelation(), time ~ variable)
    df$time <- df$time * 1e6 # Save data in microseconds
    
    columnNamesInOrder <- unique(simulatedAutocorrelation()$variable)
    if (input$rowOrder == 'alphanumeric') {
      columnNamesInOrder <- str_sort(columnNamesInOrder, numeric = TRUE)
    }  
  
    df <- df[,c("time",columnNamesInOrder)]
    
    f1 <- paste0("autocorrelation_",Sys.Date(),".csv")
    f2 <- paste0("metadata_",Sys.Date(),".csv")
    f3 <- paste0("relativeContributions_",Sys.Date(),".csv")
    
    write.csv(df, f1, row.names = FALSE)
    write.csv(reactives$simulatedMetaData, f2, row.names = FALSE)
    
    estimatedContributions       <- simulatedDistributions()$allDistributions
    write.csv(estimatedContributions, f3, row.names = FALSE)
    
    zip(zipfile=zipFile, files=c(f1,f2,f3))    

  },
  contentType = "application/zip"
)

output$download_fidelityPenaltyNorms <- downloadHandler(filename = function() {
  paste0("Norms_Raynals_ColWise_",Sys.Date(),".csv")
},content = function(file) {
  
  req(reactives$data_loaded)
  req(dlsDataUpdated())
  data <- dlsDataUpdated()
  
  lCurveData <- dlsDataUpdated()$lCurveData
  
  if (is.null(lCurveData)) return(NULL)
  
  lCurveData$lCurveDf$logFidelity <- log(lCurveData$lCurveDf$fidelity)
  lCurveData$lCurveDf$logPenalty  <- log(lCurveData$lCurveDf$penalty)
  
  write.csv(lCurveData$lCurveDf,file,row.names = F,quote = F)
})

output$download_fitted_dfs_rowWise <- downloadHandler(filename = function() {
  
  alphaTerm <- ifelse(input$regMethod == 'fixedAlpha',
                      paste0('alpha',input$alphaRegTerm),
                      paste0('Lcurve_start',input$alphaRegTermLcurve_start,
                             '_stop',input$alphaRegTermLcurve_stop,
                             '_step',input$alphaRegTermLcurve_step))
  
  paste0("Estimated_Diff_Raynals_RowWise_",alphaTerm,'_',Sys.Date(),".csv")
},content = function(file) {
  
  req(estimatedContributions())
  estimatedContributions       <- estimatedContributionsDiff()
  estimatedContributions$value <- estimatedContributions$value*100
  
  colnames(estimatedContributions) <- c('diffusionCoefficient','variable','intensity')
  write.csv(estimatedContributions,file,row.names = F,quote = F)
})

output$download_fitted_dfs_colWise <- downloadHandler(filename = function() {
  alphaTerm <- ifelse(input$regMethod == 'fixedAlpha',
                      paste0('alpha',input$alphaRegTerm),
                      paste0('Lcurve_start',input$alphaRegTermLcurve_start,
                             '_stop',input$alphaRegTermLcurve_stop,
                             '_step',input$alphaRegTermLcurve_step))
  
  paste0("Estimated_Diff_Raynals_ColWise_",alphaTerm,'_',Sys.Date(),".csv")
},content = function(file) {
  
  req(estimatedContributions())
  estimatedContributions       <- estimatedContributionsDiff()
  estimatedContributions$value <- estimatedContributions$value*100
  
  colnames(estimatedContributions) <- c('diffusionCoefficient','variable','intensity')
  
  write.csv(dcast(estimatedContributions,diffusionCoefficient~variable),file,row.names = F,quote = F)
})

