box(title = "Fitting Information", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitting_params_table', 'Fitting parameters'))),
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_params_table', 'Fitted parameters'))),
      
      conditionalPanel(condition = "output.calibrationMethod == 'calibrationFile'",
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitting_params_table_calibration', 'Fitting parameters - calibration'))),
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_params_table_calibration', 'Fitted parameters - calibration')))
      
      )
      
    ))