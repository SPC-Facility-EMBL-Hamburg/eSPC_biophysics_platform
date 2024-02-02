box(title = "Session data", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
     column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_session' ,'JSON file (to re load the whole analysis later on)')))
    ))