box(
  title = "Thermal data", width = 4, solidHeader = T, status = "primary", 
  fluidRow(
    
    column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_melting_data_all_spectra' ,'Melting spectra (all wavelengths)'))),
    column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_melting_data'         ,'Melting curves'))),
    
    conditionalPanel(
      "input.analysis_model_thermal != 'fixedWL'",
      
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_basis_spectra' ,
                            'Basis spectra'))),
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_explained_variance' ,
                            'Explained variance'))),
      
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_fitted_spectra_melting' ,
                            'Reconstructed spectra')))
      
    ),
    
    conditionalPanel(
      "output.melting_data_was_fitted || output.melting_data_was_fitted_svd_or_pca",
                     
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_melting_data_fit'     ,'Fitted melting curves'))),
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_melting_params'       ,'Fitted parameters'))),
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_melting_params_error' ,'Fitted parameters error (%)'))),
    
    )
    
))