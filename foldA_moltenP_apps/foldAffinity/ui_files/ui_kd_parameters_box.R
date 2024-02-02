box(title = "Kb (= 1 / Kd) Parameters", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3, p(HTML("<b>Constant Kb</b>"),
                  checkboxInput("const_kb", "", TRUE))),
      
      conditionalPanel(condition = "(input.const_kb)", 
       column(3, p(HTML("<b>Kd in (uM)</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu4-0"),
                   numericInput('kd_simulation_const',NULL, 20,min = 0, max = 10000),
                   tippy::tippy_this(elementId = "info_uu4-0",
                                     tooltip = "Dissociation constant Kd",placement = "right")))
      ),
      
      conditionalPanel(condition = "!(input.const_kb)", 
      
        column(2, p(HTML("<b>ΔH</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu4-1"),
                    numericInput('dHb',NULL, -5,min = -10, max = 10),
                    tippy::tippy_this(elementId = "info_uu4-1",
                                      tooltip = "Enthalpy of binding. Units are kcal/mol.",placement = "right"))),
        
        column(2, p(HTML("<b>ΔCP</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu4-2"),
                    numericInput('dCPb',NULL, -0.2,min = -2, max = 2),
                    tippy::tippy_this(elementId = "info_uu4-2",
                                      tooltip = "Heat capacity of binding. Units are kcal/°C.",placement = "right"))),
        
        column(2, p(HTML("<b>ΔS</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu4-3"),
                    numericInput('dSb',NULL, 0.01,min = -1, max = 1),
                    tippy::tippy_this(elementId = "info_uu4-3",
                                      tooltip = "Entropy of binding. Units are kcal/mol.",placement = "right"))),
        
        column(2, p(HTML("<b>T<sub>b</sub></b>"),
                    span(shiny::icon("info-circle"), id = "info_uu4-4"),
                    numericInput('Tb',NULL, 37,min = 5, max = 100),
                    tippy::tippy_this(elementId = "info_uu4-4",
                                      tooltip = "Reference temperature (in Celsius)",placement = "right")))
      
      )
      
      
    ))