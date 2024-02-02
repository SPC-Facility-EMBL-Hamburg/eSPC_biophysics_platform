box(title = "2.1 Model selection", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      column(6, p(HTML("<b></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu2-1"),
                  selectInput("model_selected", NULL,
                              c( "Local"     = "Local",
                                 "Global"    = "Global",
                                 "Global_CP"    = "Global_CP"
                                 )),
                  tippy::tippy_this(elementId = "info_uu2-1",
                                    tooltip = "Select Global_CP to fit CP and to fit the same value of the slope parameters for all curves.
                           Select Global to fix CP and to fit the same value of the slope parameters for all curves.
                           Select Local to fix CP and to fit each curve individually.",placement = "right"))),
      
      column(8,   withBusyIndicatorUI(actionButton("btn_cal","Run Fitting",class = "btn-primary")))
    ))