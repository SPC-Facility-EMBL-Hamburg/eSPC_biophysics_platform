output$download_cd_data_row_wise        <-   downloadHandler(
  filename = function() {paste0("CD_spectra_",Sys.Date(),".csv")},
  content  = function(file) {
    
    # Filter experiments where the signalDesiredUnit matrix units does not match the working units
    id_to_keep <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
    
    signalsAll   <- cdAnalyzer$getExperimentPropertiesModif('signalDesiredUnit')[id_to_keep]
    wlsAll       <- cdAnalyzer$getExperimentPropertiesModif('wavelength')[id_to_keep]
    spectraNames <- cdAnalyzer$getExperimentProperties('spectraNames')[id_to_keep]

    concentration  <- unlist(cdAnalyzer$getExperimentProperties('concentration')[id_to_keep])
    pathLength     <- unlist(cdAnalyzer$getExperimentProperties('pathLength')[id_to_keep])
    molWeight      <- unlist(cdAnalyzer$getExperimentProperties('molecularWeight')[id_to_keep])
    nResidues      <- unlist(cdAnalyzer$getExperimentProperties('numberOfResidues')[id_to_keep])
    
    dfs <- list()
    for (i in 1:length(wlsAll)) {
      
      spectraData           <- as.data.frame(cbind(wlsAll[[i]],signalsAll[[i]]))
      colnames(spectraData) <- c('wavelength',spectraNames[[i]])
      
      spectraDataMelt <-  melt(spectraData,id.vars = c('wavelength'))
      spectraDataMelt <- spectraDataMelt[order(spectraDataMelt$wavelength),]
      dfs[[i]]        <- spectraDataMelt
        
    }
    
    df <- do.call(rbind,dfs)
    colnames(df) <- c('wavelength','Sample_name',input$workingUnits)
    
    header_lines <- df_usual_comments(input$workingUnits,concentration,
                                      pathLength,molWeight,nResidues)
    
    data_lines   <- df_to_lines(df)
    
    lines        <- c(header_lines,data_lines)
    
    cat(lines, file = file, sep = "\n")
  })

output$download_cd_data_col_wise        <-   downloadHandler(
  filename = function() {paste0("CD_spectra_",Sys.Date(),".csv")},
  content  = function(file) {
    
    # Filter experiments where the signalDesiredUnit matrix units does not match the working units
    id_to_keep <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
    
    signalsAll   <- cdAnalyzer$getExperimentPropertiesModif('signalDesiredUnit')[id_to_keep]
    wlsAll       <- cdAnalyzer$getExperimentPropertiesModif('wavelength')[id_to_keep]
    spectraNames <- cdAnalyzer$getExperimentProperties('spectraNames')[id_to_keep]

    concentration  <- unlist(cdAnalyzer$getExperimentProperties('concentration')[id_to_keep])
    pathLength     <- unlist(cdAnalyzer$getExperimentProperties('pathLength')[id_to_keep])
    molWeight      <- unlist(cdAnalyzer$getExperimentProperties('molecularWeight')[id_to_keep])
    nResidues      <- unlist(cdAnalyzer$getExperimentProperties('numberOfResidues')[id_to_keep])
    
    dfs <- list()
    for (i in 1:length(wlsAll)) {
      
      spectraData           <- as.data.frame(cbind(wlsAll[[i]],signalsAll[[i]]))
      colnames(spectraData) <- c('wavelength',spectraNames[[i]])
      dfs[[i]]              <- spectraData
      
    }
    
    df <- full_join_df_lst(dfs)
    df <- df[order(df$wavelength),]
    
    header_lines <- df_usual_comments(input$workingUnits,concentration,
                                      pathLength,molWeight,nResidues)
    
    data_lines   <- df_to_lines(df)
    lines        <- c(header_lines,data_lines)
    
    cat(lines, file = file, sep = "\n")
  })

output$download_generated_cd_data_row_wise        <-   downloadHandler(
  filename = function() {paste0("CD_ChiraKit_generated_spectra_",Sys.Date(),".csv")},
  content  = function(file) {
    
    # Filter experiments where the signalDesiredUnit matrix units does not match the working units
    id_to_keep   <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
    isGenerated  <- unlist(cdAnalyzer$getExperimentProperties('isGenerated')[id_to_keep])
    
    signalsAll   <- cdAnalyzer$getExperimentPropertiesModif('signalDesiredUnit')[id_to_keep][isGenerated]
    wlsAll       <- cdAnalyzer$getExperimentPropertiesModif('wavelength')[id_to_keep][isGenerated]
    spectraNames <- cdAnalyzer$getExperimentProperties('spectraNames')[id_to_keep][isGenerated]
    
    concentration  <- unlist(cdAnalyzer$getExperimentProperties('concentration')[id_to_keep][isGenerated])
    pathLength     <- unlist(cdAnalyzer$getExperimentProperties('pathLength')[id_to_keep][isGenerated])
    molWeight      <- unlist(cdAnalyzer$getExperimentProperties('molecularWeight')[id_to_keep][isGenerated])
    nResidues      <- unlist(cdAnalyzer$getExperimentProperties('numberOfResidues')[id_to_keep][isGenerated])
    
    dfs <- list()

    for (i in 1:length(wlsAll)) {
      
      spectraData           <- as.data.frame(cbind(wlsAll[[i]],signalsAll[[i]]))
      colnames(spectraData) <- c('wavelength',spectraNames[[i]])
      
      spectraDataMelt <-  melt(spectraData,id.vars = c('wavelength'))
      spectraDataMelt <- spectraDataMelt[order(spectraDataMelt$wavelength),]
      dfs[[i]]        <- spectraDataMelt
        
    }
    
    df <- do.call(rbind,dfs)
    colnames(df) <- c('wavelength','Sample_name',input$workingUnits)
    
    header_lines <- df_usual_comments(input$workingUnits,concentration,
                                      pathLength,molWeight,nResidues)
    
    data_lines   <- df_to_lines(df)
    lines        <- c(header_lines,data_lines)
    
    cat(lines, file = file, sep = "\n")
  })

output$download_generated_cd_data_col_wise        <-   downloadHandler(
  filename = function() {paste0("CD_ChiraKit_generated_spectra_",Sys.Date(),".csv")},
  content  = function(file) {
    
    # Filter experiments where the signalDesiredUnit matrix units does not match the working units
    id_to_keep   <- !find_non_matching_units_experiments(cdAnalyzer,input$workingUnits)
    isGenerated  <- unlist(cdAnalyzer$getExperimentProperties('isGenerated')[id_to_keep])
    
    signalsAll   <- cdAnalyzer$getExperimentPropertiesModif('signalDesiredUnit')[id_to_keep][isGenerated]
    wlsAll       <- cdAnalyzer$getExperimentPropertiesModif('wavelength')[id_to_keep][isGenerated]
    spectraNames <- cdAnalyzer$getExperimentProperties('spectraNames')[id_to_keep][isGenerated]

    concentration  <- unlist(cdAnalyzer$getExperimentProperties('concentration')[id_to_keep][isGenerated])
    pathLength     <- unlist(cdAnalyzer$getExperimentProperties('pathLength')[id_to_keep][isGenerated])
    molWeight      <- unlist(cdAnalyzer$getExperimentProperties('molecularWeight')[id_to_keep][isGenerated])
    nResidues      <- unlist(cdAnalyzer$getExperimentProperties('numberOfResidues')[id_to_keep][isGenerated])
    
    dfs <- list()
    for (i in 1:length(wlsAll)) {
      
      spectraData           <- as.data.frame(cbind(wlsAll[[i]],signalsAll[[i]]))
      colnames(spectraData) <- c('wavelength',spectraNames[[i]])
      dfs[[i]]              <- spectraData
      
    }
    
    df <- full_join_df_lst(dfs)
    df <- df[order(df$wavelength),]

    header_lines <- df_usual_comments(input$workingUnits,concentration,
                                      pathLength,molWeight,nResidues)
    
    data_lines   <- df_to_lines(df)
    lines        <- c(header_lines,data_lines)
    
    cat(lines, file = file, sep = "\n")
  })

output$download_selected_cd_exp        <-   downloadHandler(
  filename = function() {
    paste0("CD_spectra_",input$selected_cd_exp,".txt")},
  
  content  = function(file) {
    
    exp <- input$selected_cd_exp
    
    # Retrieve the experiment python object
    pyObject  <- cdAnalyzer$experimentsModif[[exp]]

    # Metadata
    metadata_info <- pyObject$metadata
    
    # Start the lines vector, to fill the output file
    lines <- c()
    
    if (!is.null(metadata_info)) {
      
      metadataFeature   <- names(metadata_info)
      metadataFeature   <- paste0('#',metadataFeature,' :')
      metadataFeature   <- identate(metadataFeature,48)
      metadataValue     <- unlist(metadata_info)
      metadata_df       <- data.frame(metadataFeature,metadataValue)
      
      for (i in 1:nrow(metadata_df)) {
        lines <- c(lines,paste0(metadata_df[i,],collapse = ''))
      }
      
    }
    
    # Find the units of the desiredSignal matrix
    
    # case 1, it is a fake experiment which can't be further processed 
    if (pyObject$isFakeExperiment) {
    
      workingUnits <- pyObject$fakeExperimentSignal
      
      # case 2, the signal matches the one selected by the user at the moment of
      # exporting this data
    } else {
     
      workingUnits <- input$workingUnits
      
    }
    
    # Convert to a nice format
    workingUnits <- workingUnits2ProperLabel(workingUnits)
    
    # Replace the html like <sup> tag with the character '^' 
    workingUnits <- clean_html_sup_tag(workingUnits)

    # Path length, molecular weight, concentration and number of residues
    parametersValue <- c(paste0(pyObject$temperature,collapse = ' '), 
                         pyObject$concentration, pyObject$numberOfResidues , 
                         pyObject$pathLength, pyObject$molecularWeight,
                         workingUnits)
    
    parametersInfo <- c('Temperature (°C)','Concentration (mg/ml)',
                        'Number of residues (for protein samples)',
                        'Path length (cm)','Molecular weight (Dalton)',
                        'Units of the CD signal')
    
    parametersInfo   <- paste0('#',parametersInfo,' :')
    parametersInfo   <- identate(parametersInfo,48) 
    
    parametersDF     <- data.frame(parametersInfo,parametersValue)
    
    for (i in 1:nrow(parametersDF)) {
      lines <- c(lines,paste0(parametersDF[i,],collapse = ''))
    }
    
    # Wavelength vector
    wl       <- pyObject$wavelength
    
    # Signal matrix
    signal <- pyObject$signalDesiredUnit

    # High tension voltage matrix
    ht     <- pyObject$signalHT
    
    df     <- data.frame(wl,signal)
    
    # Retrieve cd curves names
    cd_curves_names <- pyObject$spectraNames
    # Remove extra spaces
    cd_curves_names <- gsub("\\s+"," ",cd_curves_names)
    # Replace spaces with underscores
    cd_curves_names <- gsub(" ","_",cd_curves_names)
    
    if (sum(is.na(ht)) < 10) {
      
      dfHT <- data.frame(ht)
      df   <- cbind(df,dfHT)
      
      namesHT     <- paste0('HT_/_',cd_curves_names)
      
    } else {
      
      namesHT <- NULL
      
    }
    
    df <- signif(df,5)
    
    # order using the wavelength
    df <- df[rev(order(df$wl)),]
    
    # Convert all columns to character type
    df <- as.data.frame(lapply(df, as.character))
    
    # Add the first line (columns names)
    namesSignal     <- paste0('CD_/_',cd_curves_names)
    
    maxNchar <- max(c(sapply(namesSignal, nchar),15)) + 1
    
    # Find the spacing we need based on how long the cd curves names are
    custom_spacing <- function(x) sprintf(paste0('%-',maxNchar,'s'),x)
    
    df <- as.data.frame(lapply( df, custom_spacing ))
    
    line1 <- sapply(c('Wavelength_(nm)',namesSignal,namesHT), custom_spacing)
    line1 <- paste0(line1,collapse = '')
    # Create a character vector with custom text
    lines <- c(lines,line1)
    
    for (i in 1:nrow(df)) {
      lines <- c(lines,paste0(df[i,],collapse = ''))
    }
    
    # Specify the file and use the cat() function to create the text file
    cat(lines, file = file, sep = "\n")
    
  })

output$download_melting_data        <-   downloadHandler(
  filename = function() {paste0("CD_Melting_Data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df      <- generate_thermal_ramp_df(cdAnalyzer)
    df      <- reassign_unfolding_df_colnames(df,reactives$spectra_decomposition_method_thermal)
    
    write.csv(df,file,row.names = F)
  })

output$download_chemical_data        <-   downloadHandler(
  filename = function() {paste0("CD_ChemicalUnfolding_Data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df      <- generate_chemical_unfolding_df(cdAnalyzer)
    df      <- reassign_unfolding_df_colnames(df,reactives$spectra_decomposition_method_chemical)
    
    write.csv(df,file,row.names = F)
  })

output$download_custom_data        <-   downloadHandler(
  filename = function() {paste0("CD_Custom_Data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df      <- generate_custom_df(cdAnalyzer)
    df      <- reassign_unfolding_df_colnames(df,reactives$spectra_decomposition_method_custom)
    
    write.csv(df,file,row.names = F)
  })

output$download_melting_data_fit        <-   downloadHandler(
  filename = function() {paste0("CD_Melting_Fitted_Data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    dfFit   <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signal_predicted')
    dfFit   <- reassign_unfolding_df_colnames(dfFit,reactives$fitted_coefficients_method_thermal)
    
    colnames(dfFit)[3] <- paste0(colnames(dfFit)[3],'_fit')
    
    write.csv(dfFit,file,row.names = F)
  })

output$download_chemical_data_fit        <-   downloadHandler(
  filename = function() {paste0("CD_ChemicalUnfolding_Fitted_Data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    dfFit   <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='signal_predicted')
    dfFit   <- reassign_unfolding_df_colnames(dfFit,reactives$fitted_coefficients_method_chemical)
    
    colnames(dfFit)[3] <- paste0(colnames(dfFit)[3],'_fit')
    
    write.csv(dfFit,file,row.names = F)
  })

output$download_custom_data_fit        <-   downloadHandler(
  filename = function() {paste0("CD_Custom_Fitted_Data_",Sys.Date(),".csv")},
  content  = function(file) {
    
    dfFit   <- generate_custom_df(cdAnalyzer,signal_type='signal_predicted')
    dfFit   <- reassign_unfolding_df_colnames(dfFit,reactives$fitted_coefficients_method_custom)
    
    colnames(dfFit)[3] <- paste0(colnames(dfFit)[3],'_fit')
    
    write.csv(dfFit,file,row.names = F)
  })

output$download_melting_data_all_spectra        <-   downloadHandler(
  filename = function() {paste0("CD_Melting_Data_All_Wavelengths_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df      <- generate_thermal_ramp_df(cdAnalyzer,signal_type='signalDesiredUnit')
    df      <- reassign_unfolding_df_colnames(df)
    
    write.csv(df,file,row.names = F)
  })

output$download_chemical_data_all_spectra        <-   downloadHandler(
  filename = function() {paste0("CD_ChemicalUnfolding_Data_All_Wavelengths_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df      <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='signalDesiredUnit')
    df      <- reassign_unfolding_df_colnames(df)
    
    write.csv(df,file,row.names = F)
  })

output$download_custom_data_all_spectra        <-   downloadHandler(
  filename = function() {paste0("CD_Custom_Data_All_Wavelengths_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df      <- generate_custom_df(cdAnalyzer,signal_type='signalDesiredUnit')
    df      <- reassign_unfolding_df_colnames(df)
    
    write.csv(df,file,row.names = F)
  })

output$download_fitted_spectra_melting        <-   downloadHandler(
  filename = function() {paste0("CD_Melting_Data_All_Wavelengths_Reconstructed_",Sys.Date(),".csv")},
  content  = function(file) {
    
    dfFit      <- generate_thermal_ramp_df(cdAnalyzer,signal_type='fitted_spectra')
    dfFit      <- reassign_unfolding_df_colnames(dfFit)
    
    colnames(dfFit)[3] <- paste0(colnames(dfFit)[3],'_fit')
    
    write.csv(dfFit,file,row.names = F)
  })

output$download_fitted_spectra_chemical        <-   downloadHandler(
  filename = function() {paste0("CD_ChemicalUnfolding_Data_All_Wavelengths_Reconstructed_",Sys.Date(),".csv")},
  content  = function(file) {
    
    dfFit      <- generate_chemical_unfolding_df(cdAnalyzer,signal_type='fitted_spectra')
    dfFit      <- reassign_unfolding_df_colnames(dfFit)
    
    colnames(dfFit)[3] <- paste0(colnames(dfFit)[3],'_fit')
    
    write.csv(dfFit,file,row.names = F)
  })

output$download_fitted_spectra_custom        <-   downloadHandler(
  filename = function() {paste0("CD_Custom_Data_All_Wavelengths_Reconstructed_",Sys.Date(),".csv")},
  content  = function(file) {
    
    dfFit      <- generate_custom_df(cdAnalyzer,signal_type='fitted_spectra')
    dfFit      <- reassign_unfolding_df_colnames(dfFit)
    
    colnames(dfFit)[3] <- paste0(colnames(dfFit)[3],'_fit')
    
    write.csv(dfFit,file,row.names = F)
  })

output$download_melting_params    <-   downloadHandler(
  filename = function() {paste0("fitted_params_melting_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df <- get_fitted_params_unfolding(cdAnalyzer)
    
    write.csv(df,file,row.names = F)
  })

output$download_chemical_params    <-   downloadHandler(
  filename = function() {paste0("fitted_params_chemical_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df <- get_fitted_params_unfolding(cdAnalyzer,'Chemical')
    
    write.csv(df,file,row.names = F)
  })

output$download_custom_params    <-   downloadHandler(
  filename = function() {paste0("fitted_params_custom_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df <- get_fitted_params_unfolding(cdAnalyzer,'Custom')
    
    write.csv(df,file,row.names = F)
  })

output$download_melting_params_error  <-   downloadHandler(
  filename = function() {paste0(
    "fitted_params_relative_errors_percentage_melting_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df <- get_fitted_params_unfolding(cdAnalyzer,errors=TRUE)
    
    write.csv(df,file,row.names = F)
  })

output$download_chemical_params_error  <-   downloadHandler(
  filename = function() {paste0(
    "fitted_params_relative_errors_percentage_chemical_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df <- get_fitted_params_unfolding(cdAnalyzer,'Chemical',errors=TRUE)
    
    write.csv(df,file,row.names = F)
  })

output$download_custom_params_error  <-   downloadHandler(
  filename = function() {paste0(
    "fitted_params_relative_errors_percentage_custom_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df <- get_fitted_params_unfolding(cdAnalyzer,'Custom',errors=TRUE)
    
    write.csv(df,file,row.names = F)
  })

output$download_basis_spectra    <-   downloadHandler(
  filename = function() {
    paste0("basis_spectra_melting_",
           reactives$spectra_decomposition_method_thermal,"_",Sys.Date(),".csv")
    },
  content  = function(file) {
    
    df           <- get_basis_spectra_df(cdAnalyzer)
    colnames(df) <- c('wavelength_(nm)','k','CD_signal_value','legend')
    
    write.csv(df,file,row.names = F)
  })

output$download_basis_spectra_chemical    <-   downloadHandler(
  filename = function() {
    paste0("basis_spectra_chemical_",
           reactives$spectra_decomposition_method_chemical,"_",Sys.Date(),".csv")
  },
  content  = function(file) {
    
    df           <- get_basis_spectra_df(cdAnalyzer,'Chemical')
    colnames(df) <- c('wavelength_(nm)','k','CD_signal_value','legend')
    
    write.csv(df,file,row.names = F)
  })

output$download_basis_spectra_custom    <-   downloadHandler(
  filename = function() {
    paste0("basis_spectra_custom_",
           reactives$spectra_decomposition_method_custom,"_",Sys.Date(),".csv")
  },
  content  = function(file) {
    
    df           <- get_basis_spectra_df(cdAnalyzer,'Custom')
    colnames(df) <- c('wavelength_(nm)','k','CD_signal_value','legend')
    
    write.csv(df,file,row.names = F)
  })

output$download_explained_variance    <-   downloadHandler(
  filename = function() {
    paste0("explained_variance_melting_",
           reactives$spectra_decomposition_method_thermal,"_",Sys.Date(),".csv")
    },
  content  = function(file) {
    
    df  <- get_explained_variance_df(cdAnalyzer,'Thermal')
    
    write.csv(df,file,row.names = F)
  })

output$download_explained_variance_chemical    <-   downloadHandler(
  filename = function() {
    paste0("explained_variance_chemical_",
    reactives$spectra_decomposition_method_chemical,
    Sys.Date(),".csv")},
  content  = function(file) {
    
    df  <- get_explained_variance_df(cdAnalyzer,'Chemical')
    
    write.csv(df,file,row.names = F)
  })

output$download_explained_variance_custom    <-   downloadHandler(
  filename = function() {
    paste0("explained_variance_custom_",
           reactives$spectra_decomposition_method_custom,
           Sys.Date(),".csv")},
  content  = function(file) {
    
    df  <- get_explained_variance_df(cdAnalyzer,'Custom')
    
    write.csv(df,file,row.names = F)
  })

output$download_log_book <-   downloadHandler(
  filename = function() {
    paste0("logbook_ChiraKit_",
           Sys.Date(),".txt")},
  content  = function(file) {
    
    lines <- unlist(reactives$logbook)
    lines <- purge_logbook_lines(lines)
    
    cat(lines,file = file, sep = "\n")
  })

output$download_calculated_sec_str <-   downloadHandler(
  filename = function() {
    paste0("dssp_based_sec_str_",
           Sys.Date(),".txt")},
  content  = function(file) {
    
    sec_str_dfs <- list()
    
    cd_data_files   <- input$pdbFiles$datapath
    names           <- input$pdbFiles$name
    
    for (i in 1:length(cd_data_files)) {
      df <- run_dssp_workflow(cd_data_files[i])
      if (!is.null(df)) {
        
        df$Name     <- names[i]
        sec_str_dfs <- append(sec_str_dfs,list(df))
        
      }
    }
    
    all_dfs <- do.call(rbind,sec_str_dfs)
    
    write.csv(all_dfs,file,row.names = F,quote = F)
  })

output$download_estimated_sec_str <-   downloadHandler(
  filename = function() {
    paste0("spectra_derived_sec_str_",
           Sys.Date(),".txt")},
  content  = function(file) {
    
    sec_str_dfs <- list()
    
    expNames <- cdAnalyzer$experimentNames
    
    for (exp in expNames) {
      
      secondary_structure_content_lst <- cdAnalyzer$experimentsOri[[exp]]$secondary_structure_content
      
      if (length(secondary_structure_content_lst) > 0) {
        
        spectraNames <- cdAnalyzer$experimentsOri[[exp]]$spectraNames
        
        for (i in 1:length(spectraNames)) {
          
          sec_str_df      <- as.data.frame(secondary_structure_content_lst[i])
          sec_str_df$name <- spectraNames[i]
                                
          method <- gsub('Component_','',colnames(sec_str_df)[1])
      
          colnames(sec_str_df)[1] <- "Component"
          
          sec_str_df$method       <- method
          
          sec_str_dfs             <- append(sec_str_dfs,list(sec_str_df))
        }
      }
    }
    
    all_dfs <- do.call(rbind,sec_str_dfs)
    
    write.csv(all_dfs,file,row.names = F,quote = F)
  })

output$download_sec_str_spectra_fittings <-   downloadHandler(
  filename = function() {
    paste0("spectra_fittings_",
           Sys.Date(),".txt")},
  content  = function(file) {
    
    df <- lst_of_cd_spectra_to_df(reactives$list_of_signals,reactives$list_of_fittings)
      
    write.csv(df,file,row.names = F,quote = F)
  })

output$download_average_spectra        <-   downloadHandler(
  filename = function() {paste0("CD_average_spectra_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df           <- cbind(compareSpectraPyClass$wavelength,compareSpectraPyClass$means)
    colnames(df) <- c('wavelength_nm',compareSpectraPyClass$labels_unique)
      
    write.csv(df,file,row.names = F)
  })

output$download_standard_deviation        <-   downloadHandler(
  filename = function() {paste0("CD_standard_deviation_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df           <- cbind(compareSpectraPyClass$wavelength,compareSpectraPyClass$sds)
    colnames(df) <- c('wavelength_nm',paste0('std ',compareSpectraPyClass$labels_unique))
    
    write.csv(df,file,row.names = F)
  })

output$download_difference_spectra        <-   downloadHandler(
  filename = function() {paste0("CD_difference_spectra_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df           <- cbind(compareSpectraPyClass$wavelength,compareSpectraPyClass$difference_spectra)
    colnames(df) <- c('wavelength_nm',compareSpectraPyClass$difference_spectra_lbl)
  
    write.csv(df,file,row.names = F)
  })

output$download_difference_std        <-   downloadHandler(
  filename = function() {paste0("CD_difference_spectra_std_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df           <- cbind(compareSpectraPyClass$wavelength,compareSpectraPyClass$difference_spectra_sd)
    colnames(df) <- c('wavelength_nm',paste0('std ',compareSpectraPyClass$difference_spectra_lbl))

    write.csv(df,file,row.names = F)
  })

output$download_distance_data        <-   downloadHandler(
  filename = function() {paste0("CD_distance_matrix_normalised_euclidean_",Sys.Date(),".csv")},
  content  = function(file) {
    
    df           <- compareSpectraPyClass$distance_matrix
    colnames(df) <- compareSpectraPyClass$labels
    
    write.csv(df,file,row.names = F)
  })

