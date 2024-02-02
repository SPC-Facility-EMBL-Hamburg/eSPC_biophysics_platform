box(title = "Sort Fitted Parameters", width = 2, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(8, p(HTML("<b></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu2-4"),
                  
                  selectInput("sort_table_parameter", NULL,c("Tm"   = "Tm")),
                  
                  tippy::tippy_this(elementId = "info_uu2-4",
                                    tooltip = "Select which fitted parameter to sort tables",placement = "right")))
      
    ))