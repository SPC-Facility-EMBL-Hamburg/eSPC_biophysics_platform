box(title = "Fitting Information", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_table', 'Fitted Parameters (Fluorescence)'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_errors_table', 'Relative Errors (Fluorescence)'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_table_fu', 'Fitted Parameters (Fraction Unfolded)'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_errors_table_fu', 'Relative Errors (Fraction Unfolded)')))
    ))