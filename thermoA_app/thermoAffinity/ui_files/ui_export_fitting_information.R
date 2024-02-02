box(title = "Fitting information", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_params_table', 'Fitted parameters')))
    ))