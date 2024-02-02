box(
  title = "Custom data", width = 4, solidHeader = T, status = "primary", 
  fluidRow(
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_custom_data_all_spectra' ,'Custom spectra (all wavelengths)'))),
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_custom_data'             ,'Custom curves'))),
    
    conditionalPanel(
      "input.analysis_model_custom != 'fixedWL'",
      
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_basis_spectra_custom' ,
                            'Basis spectra'))),
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_explained_variance_custom' ,
                            'Explained variance'))),
      
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_fitted_spectra_custom' ,
                            'Reconstructed spectra')))
      
    ),
    
    conditionalPanel(
      "output.custom_data_was_fitted || output.custom_data_was_fitted_svd_or_pca",
      
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_custom_data_fit'     ,'Fitted unfolding curves'))),
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_custom_params'       ,'Fitted parameters'))),
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_custom_params_error' ,'Fitted parameters error (%)'))),
      
    )
    
  ))