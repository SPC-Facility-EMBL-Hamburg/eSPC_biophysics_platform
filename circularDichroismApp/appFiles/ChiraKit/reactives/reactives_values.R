# Set of initial reactiveValues that are useful to manipulate the UI elements

reactives <- reactiveValues(
  logbook                                = list(),# record data manipulation steps
  data_loaded                            = NULL,  # show plots/tables
  showVoltageThreshold                   = FALSE, # show the HT filter
  thermalWorkingUnits                    = NULL,  # To plot the y-axis with the correct units
  thermalDatasetCreated                  = FALSE, # allow thermal  analysis (thermal_ramp_reactives.R)
  melting_data_was_fitted                = FALSE, # show thermal   results  (thermal_ramp_reactives.R)
  melting_data_was_fitted_svd_or_pca     = FALSE, # show thermal   results  (thermal_ramp_reactives.R)
  spectra_was_decomposed                 = FALSE, # show thermal   results  (thermal_ramp_reactives.R),
  chemicalWorkingUnits                   = NULL,  # To plot the y-axis with the correct units
  chemicalDatasetCreated                 = FALSE, # allow chemical analysis (chemical_denaturation_reactives.R)
  chemical_data_was_fitted               = FALSE, # show chemical  results  (chemical_denaturation_reactives.R)
  chemical_data_was_fitted_svd_or_pca    = FALSE, # show chemical  results  (chemical_denaturation_reactives.R)
  spectra_was_decomposed_chemical        = FALSE, # show chemical  results  (chemical_denaturation_reactives.R)
  show_basis_change_option               = FALSE, # shared assuming the user will analyse one dataset at a time 
                                                  #  used in thermal_ramp_reactives.R, chemical_denaturation_reactives.R and custom_analysis_reactives.R
  secStrFittingWasDone                   = FALSE, # show export options - secondary str. from spectra
  secStrCalcWasDone                      = FALSE, # show export options - secondary str. from PDB file
  secStrCalculationTabsNames             = c(),   # panel 'menu_sec_structure' (secondary_structure_reactives.R)
  
  secStrRefMatrixC                       = NULL, # to use user reference sets, panel 'menu_sec_structure'
  secStrRefMatrixF                       = NULL, # to use user reference sets, panel 'menu_sec_structure'
  
  customWorkingUnits                     = NULL,  # To plot the y-axis with the correct units
  customDatasetCreated                   = FALSE, # allow custom analysis
  custom_data_was_fitted                 = FALSE, # show  custom analysis results
  custom_data_was_fitted_svd_or_pca      = FALSE, # show  custom analysis results
  spectra_was_decomposed_custom          = FALSE, # show  custom  results  (custom_analysis_reactives.R)
  
  compareDatasetCreated                   = FALSE, # allow compare spectra analysis
  
  spectra_decomposition_method_thermal   = 'None',# can be 'svd' or 'pca' (thermal_ramp_reactives.R)
  spectra_decomposition_method_chemical  = 'None',# can be 'svd' or 'pca' (chemical_denaturation_reactives.R)
  spectra_decomposition_method_custom    = 'None',# can be 'svd' or 'pca' (custom_analysis_reactives.R)
  fitted_coefficients_method_thermal     = 'None',# can be 'svd' or 'pca' (thermal_ramp_reactives.R)
  fitted_coefficients_method_chemical    = 'None',# can be 'svd' or 'pca' (chemical_denaturation_reactives.R)
  fitted_coefficients_method_custom      = 'None',# can be 'svd' or 'pca' (custom_analysis_reactives.R)

  sesca_pred_was_run                     = FALSE,# to show sesca results (sesca_reactives.R)
  sesca_est_was_run                      = FALSE,# to show sesca results (sesca_reactives.R)

  GQ_ref_load                            = FALSE, #(gQuadruplex_reactives.R)
  GQ_sample_load                         = FALSE, #(gQuadruplex_reactives.R)
  secondary_parameters                   = NULL,  #(gQuadruplex_reactives.R)
  tertiary_parameers                     = NULL   #(gQuadruplex_reactives.R)
  
)

# Allow the UI to know the values of these reactives

output$show_basis_change_option   <- reactive( { return( reactives$show_basis_change_option  ) } )

output$thermalDatasetCreated   <- reactive( { return( reactives$thermalDatasetCreated  ) } )
output$chemicalDatasetCreated  <- reactive( { return( reactives$chemicalDatasetCreated ) } )
output$customDatasetCreated    <- reactive( { return( reactives$customDatasetCreated ) } )
output$compareDatasetCreated   <- reactive( { return( reactives$compareDatasetCreated ) } )

output$melting_data_was_fitted  <- reactive( { return( reactives$melting_data_was_fitted  ) } )
output$chemical_data_was_fitted <- reactive( { return( reactives$chemical_data_was_fitted ) } )
output$custom_data_was_fitted   <- reactive( { return( reactives$custom_data_was_fitted ) } )

output$melting_data_was_fitted_svd_or_pca  <- reactive( { return( reactives$melting_data_was_fitted_svd_or_pca  ) } )
output$chemical_data_was_fitted_svd_or_pca <- reactive( { return( reactives$chemical_data_was_fitted_svd_or_pca ) } )
output$custom_data_was_fitted_svd_or_pca   <- reactive( { return( reactives$custom_data_was_fitted_svd_or_pca ) } )

output$secStrFittingWasDone   <- reactive( { return( reactives$secStrFittingWasDone  ) } )
output$secStrCalcWasDone      <- reactive( { return( reactives$secStrCalcWasDone     ) } )

output$sesca_pred_was_run     <- reactive( { return( reactives$sesca_pred_was_run    ) } )
output$sesca_est_was_run      <- reactive( { return( reactives$sesca_est_was_run     ) } )


outputOptions(output, "show_basis_change_option" , suspendWhenHidden = FALSE)

outputOptions(output, "melting_data_was_fitted" , suspendWhenHidden = FALSE)
outputOptions(output, "chemical_data_was_fitted", suspendWhenHidden = FALSE)
outputOptions(output, "custom_data_was_fitted",   suspendWhenHidden = FALSE)

outputOptions(output, "melting_data_was_fitted_svd_or_pca",  suspendWhenHidden = FALSE)
outputOptions(output, "chemical_data_was_fitted_svd_or_pca", suspendWhenHidden = FALSE)
outputOptions(output, "custom_data_was_fitted_svd_or_pca",   suspendWhenHidden = FALSE)

outputOptions(output, "thermalDatasetCreated",  suspendWhenHidden = FALSE)
outputOptions(output, "chemicalDatasetCreated", suspendWhenHidden = FALSE)
outputOptions(output, "customDatasetCreated",   suspendWhenHidden = FALSE)
outputOptions(output, "compareDatasetCreated",   suspendWhenHidden = FALSE)

outputOptions(output, "secStrFittingWasDone" , suspendWhenHidden = FALSE)
outputOptions(output, "secStrCalcWasDone"    , suspendWhenHidden = FALSE)

outputOptions(output, "sesca_pred_was_run" , suspendWhenHidden = FALSE)
outputOptions(output, "sesca_est_was_run"  , suspendWhenHidden = FALSE)