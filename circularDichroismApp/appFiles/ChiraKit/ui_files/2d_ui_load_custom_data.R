box(title = "1. Experimental parameters", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      # The hidden elements (shinyjs::hidden) are there
      # in case that we want to extend the code to fit functions
      # that depend on two experimental parameters, e.g., 
      # fit the CD(signal) as a function of both temperature and chemical agent concentration
      
      shinyjs::hidden(column(6, p(HTML("<b>Number of experimental parameters</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_numberOfExpParams"),
                  numericInput('numberOfExpParams',NULL, 1,min = 1, max = 2),
                  tippy::tippy_this(elementId = "info_uu_numberOfExpParams",
                                    tooltip = "Select how to analyse the CD signal.
                                    For example, if you want to analyse the CD signal as a function of time, 
                                    use one experimental parameter. If you want to analyse the CD signal 
                                    as a function of time and temperature, 
                                    use two experimental parameters.",placement = "right")))),
      
      column(6, p(HTML("<b>Parameter name</b>"),
                  textInput('firstParamName',NULL,'A'))),
      
      shinyjs::hidden(column(3, p(HTML("<b>2nd parameter name</b>"),
                  textInput('secondParamName',NULL,'B')))),
      
      column(width = 12,rHandsontableOutput('custom_parameters_data'))
      
    ))

