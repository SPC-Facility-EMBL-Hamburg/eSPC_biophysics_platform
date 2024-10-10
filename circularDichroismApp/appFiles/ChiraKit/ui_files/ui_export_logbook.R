box(title = "Logbook", width = 2, solidHeader = T, status = "primary", 
    fluidRow(

      column(12, p(style = "font-size: 120%",HTML(""),
                   downloadLink('download_log_book'  ,'Logbook')))
    ))