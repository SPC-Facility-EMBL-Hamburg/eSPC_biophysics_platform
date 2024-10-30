observeEvent(list(input$pdbFilesSESCA,input$sescaBasisSet,input$sescaReference),{
    output$sesca_plot <- NULL
})

find_ref_spec <- function (sescaReference) {

    # Find the reference spectrum in the CDAnalyzer object
    names          <- cdAnalyzer$experimentNames
    internalIDs    <- cdAnalyzer$get_experiment_properties('internalID')

    id                 <- which(unlist(internalIDs) == sescaReference)
    names(internalIDs) <- names
    id                 <- found_ids(internalIDs,c(id))

    exp <- names(id)[1]
    pos <- id[1]

    return(list(exp,pos))
}

# Obtain the wavelength and signal of the selected experiment and spectrum position
process_sesca_signal_wl <- function (exp,pos,scaling_factor=1) {

    if (file.exists("sesca_reference.dat")) file.remove("sesca_reference.dat")

    # Convert to delta epsilon
    cdAnalyzer$experimentsModif[[exp]]$experiment_from_abs_to_other_units('meanUnitMolarExtinction')

    if (all(is.na(cdAnalyzer$experimentsModif[[exp]]$signalDesiredUnit))) {
        cdAnalyzer$experimentsModif[[exp]]$experiment_from_abs_to_other_units(input$workingUnits)

        popUpWarning("The reference spectrum is not available in delta epsilon units.
        Please check that the require parameters are available (protein concentration, path length,
        number of peptide bonds and molecular weight).")
        return(NULL)
    }

    signal <- cdAnalyzer$experimentsModif[[exp]]$signalDesiredUnit[,pos] * scaling_factor

    wavelength_ori     <- cdAnalyzer$experimentsModif[[exp]]$wavelength

    # Interpolate the signal to every 1 nm, so SESCA can process the data
    wavelength <- seq(ceiling(min(wavelength_ori)),floor(max(wavelength_ori)),1)

    signal_interpolated <- approx(wavelength_ori, signal, xout = wavelength)$y

    signal_sesca_units  <- signal_interpolated * 3.298 # convert from delta epsilon units to kMRE units, used in the SESCA algorithm

    data_to_save <- data.frame(wavelength, signal_sesca_units)

    # Create a file called sesca_reference.dat with the interpolated data

    append_record_to_logbook(
    paste0("Creating the reference spectrum (",(input$sescaReference),") for the SESCA algorithm.The data was interpolated to be evenly spaced every 1 nm.")
    )

    append_record_to_logbook(paste0("The scaling factor applied to the reference spectrum is ",scaling_factor))

    write.table(data_to_save, file = "sesca_reference.dat", sep = " ", row.names = FALSE, col.names = FALSE,quote = FALSE)

    cdAnalyzer$experimentsModif[[exp]]$experiment_from_abs_to_other_units(input$workingUnits)
    Sys.sleep(0.5)

    return(list(wavelength,signal_interpolated))
}

observeEvent(input$runSESCA,{

    reactives$sesca_pred_was_run <- FALSE

    if (is.null(input$pdbFilesSESCA)) return(NULL)

    pdb_paths  <- input$pdbFilesSESCA$datapath
    pdb_names  <- input$pdbFilesSESCA$name

    # Check if the file is a PDB file
    if (length(grep(".pdb",pdb_names)) == 0) {
        popUpWarning("The file(s) must be in PDB format.")
        return(NULL)
    } else {

        basisSet <- input$sescaBasisSet

        if (input$sescaReference == 'None') {

            popUpInfo("The SESCA algorithm is running. This may take a few minutes.")
            sescaPyClass$predict(pdb_paths,pdb_names,basisSet)
            Sys.sleep(0.2)
            output$sesca_plot <- renderPlotly({sesca_plot(sescaPyClass,input$sesca_plot_width,
                                                          input$sesca_plot_height,input$sesca_plot_type,
                                                          input$sesca_axis_size,
                                                          average_ensemble = input$sescaAverageEnsemble)})
            output$sesca_comparison_stats <- NULL

        } else {

            # Find the reference spectrum in the CDAnalyzer object
            l <- find_ref_spec(input$sescaReference)

            exp <- l[[1]]
            pos <- l[[2]]

            scaling_factor <- as.numeric(input$sescaScalingFactor)

            l <- process_sesca_signal_wl(exp,pos,scaling_factor)

            if (is.null(l)) return(NULL) # Stop if no data in meanUnitMolarExtinction units

            wavelength_ref  <- l[[1]]
            signal_ref      <- l[[2]]

            popUpInfo("The SESCA algorithm is running. This may take a few minutes.")
            sescaPyClass$predict(pdb_paths,pdb_names,basisSet,"sesca_reference.dat")
            Sys.sleep(0.1)

            ref_name <- input$sescaReference
            if (scaling_factor != 1) {
                ref_name <- paste0(ref_name," (scaled by ",scaling_factor,")")
            }

            output$sesca_plot <- renderPlotly({sesca_plot(
              sescaPyClass,input$sesca_plot_width,input$sesca_plot_height,input$sesca_plot_type, input$sesca_axis_size,
              wavelength_ref,signal_ref,ref_name,input$sescaAverageEnsemble)})

            output$sesca_comparison_stats <- renderTable({sescaPyClass$comparison_stats},options = list(scrollX = TRUE))
        }

        output$sesca_sec_str <- renderTable({
            table <- sescaPyClass$SS_Comp
            colnames(table) <- fix_ss_labels(colnames(table),basisSet)
            return(table)
        })

        append_record_to_logbook("Predicting the CD spectra of the provided PDB file(s) using the SESCA algorithm.")

        append_record_to_logbook(paste0("Selected PDBs: ", paste(pdb_names, collapse = ", ")))
        append_record_to_logbook(paste0("Selected Basis Set: ", basisSet))

        reactives$sesca_pred_was_run <- TRUE

    }

})

observeEvent(input$runSESCA_est,{

    output$sesca_sec_str_est <- NULL

    reactives$sesca_est_was_run <- FALSE

    if (input$sescaReferenceEstimation == 'None') return(NULL)

        # Find the reference spectrum in the CDAnalyzer object
    l <- find_ref_spec(input$sescaReferenceEstimation)

    exp <- l[[1]]
    pos <- l[[2]]

    l <- process_sesca_signal_wl(exp,pos)
    if (is.null(l)) return(NULL) # Stop if no data in meanUnitMolarExtinction units

    # Check lower wavelength limit
    wavelength <- l[[1]]
    if (min(wavelength) > 190) {
        popUpWarning("The wavelength range of the current reference spectrum is too high.
        Please choose a reference spectrum that extends to at least 190 nm.")
        return(NULL)
    }

    basisSet <- input$sescaBasisSetEstimation

    # Check lower wavelength limit for the SESCA basis set
    # Basis set with 3 or 4 secondary structure elements
    basis_set_filt <- c('DS-dT','DSSP-T', 'DSSP-1','DS-dTSC3')

    if (min(wavelength) > 180 && !basisSet %in% basis_set_filt) {
        popUpWarning("The wavelength range of the current reference spectrum is too high.
        Please choose a reference spectrum that extends to at least 180 nm, or change the basis set
        to DS-dT, DSSP-T, DSSP-1, or DS-dTSC3.")
        return(NULL)
    }

    Sys.sleep(0.1)
    popUpInfo("The SESCA algorithm started. This may take a few minutes.")

    if (basisSet %in% c('DS-dTSC3','DS5-4SC1','DS6-1SC1','DS5-6SC','DSSP-TSC1','DSSP-1SC3','HBSS-3SC1')) {

        if (is.null(input$pdbFileSESCA_est)) {
            popUpWarning("Please provide a PDB file if you want to take into account the side chain contributions.")
            return(NULL)
        }
        sescaPyClass$bayes_estimate("sesca_reference.dat",basisSet,input$sescaIterations,input$pdbFileSESCA_est$datapath)
    } else {
        sescaPyClass$bayes_estimate("sesca_reference.dat",basisSet,input$sescaIterations)
    }

    popUpSuccess("The SESCA algorithm finished.")
    Sys.sleep(0.1)
    bayes_results <- format_bayes_results(extract_scaling_factor_and_ss_composition('Bayes_est_1.out'))

    bayes_results[,'Parameter'] <- fix_ss_labels(bayes_results[,'Parameter'],basisSet)

    append_record_to_logbook("Estimating the secondary structure composition of the reference spectrum using the SESCA bayesian method.")
    append_record_to_logbook(paste0("Selected basis set: ",basisSet))
    append_record_to_logbook(paste0("Number of iterations: ",input$sescaIterations))

    output$sesca_sec_str_est <- renderTable({bayes_results})

    bin_data <- extract_bin_data('Bayes_est_1.out')

    output$sesca_plot_bayes <- renderPlotly({plot_heatmap_bayes_posterior(bin_data,input$sesca_bayes_pp_x,input$sesca_bayes_pp_y,
    input$sesca_plot_sesca_plot_width,input$sesca_plot_height,input$sesca_plot_type,input$sesca_axis_size)})

    reactives$sesca_est_was_run <- TRUE

})
