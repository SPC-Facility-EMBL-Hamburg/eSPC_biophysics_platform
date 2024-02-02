box(title = "Fitting Information", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_table', 'Fitted Parameters Table'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_errors_table', 'Relative Errors Table'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_fitted_cond_table', 'Fitted Conditions Table')))
    ))