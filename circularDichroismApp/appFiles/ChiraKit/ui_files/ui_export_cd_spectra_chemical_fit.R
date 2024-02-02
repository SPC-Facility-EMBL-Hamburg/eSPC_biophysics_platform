box(
  title = "Chemical data", width = 4, solidHeader = T, status = "primary", 
  fluidRow(
    
    column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_chemical_data_all_spectra' ,'Unfolding spectra (all wavelengths)'))),
    column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_chemical_data'             ,'Unfolding curves'))),
    
    conditionalPanel(
      "input.analysis_model_chemical != 'fixedWL'",
      
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_basis_spectra_chemical' ,
                            'Basis spectra'))),
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_explained_variance_chemical' ,
                            'Explained variance'))),
      
      column(12, 
             p(style = "font-size: 120%",HTML(""),
               downloadLink('download_fitted_spectra_chemical' ,
                            'Reconstructed spectra')))
      
    ),
    
    conditionalPanel(
      "output.chemical_data_was_fitted || output.chemical_data_was_fitted_svd_or_pca",
      
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_chemical_data_fit'     ,'Fitted unfolding curves'))),
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_chemical_params'       ,'Fitted parameters'))),
      column(12, p(style = "font-size: 120%",HTML(""),downloadLink('download_chemical_params_error' ,'Fitted parameters error (%)'))),
      
    )
    
  ))