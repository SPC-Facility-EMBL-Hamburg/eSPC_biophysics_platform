box(title = "Ku Parameters", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(4, p(HTML("<b>ΔH</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu4-5"),
                  numericInput('dHu',NULL, 90,min = 0, max = 140),
                  tippy::tippy_this(elementId = "info_uu4-5",
                                    tooltip = "Enthalpy of unfolding. Units are kcal/mol.",placement = "right"))),
      
      column(4, p(HTML("<b>ΔCP</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu4-6"),
                  numericInput('dCPu',NULL,1 ,min = 0, max = 5),
                  tippy::tippy_this(elementId = "info_uu4-6",
                                    tooltip = "Heat capacity of unfolding. Units are kcal/°C.",placement = "right"))),
      
      #column(3, p(HTML("<b>DS</b>"),
      #           span(shiny::icon("info-circle"), id = "info_uu4-7"),
      #           numericInput('dSu',NULL, 0.19,min = 0, max = 2),
      #           tippy::tippy_this(elementId = "info_uu4-7",
      #                             tooltip = "Entropy of unfolding",placement = "right"))),
      
      column(4, p(HTML("<b>T<sub>m</sub></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu4-8"),
                  numericInput('Tu',NULL, 60,min = 5, max = 100),
                  tippy::tippy_this(elementId = "info_uu4-8",
                                    tooltip = "Reference temperature (in Celsius) were [F] = [U]",placement = "right")))
      
    ))