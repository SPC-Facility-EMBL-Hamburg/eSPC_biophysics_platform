box(title = "Simulation", width = 6, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(2, p(HTML("<b>Activate</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuSim1"),
                  
                  checkboxInput("runSimulation", "", FALSE)),
             tippy::tippy_this(elementId = "info_uuSim1",
                               tooltip = "Simulate a truncated gaussian and display it in the plot.",placement = "right")),
      
      conditionalPanel('input.runSimulation',
                       
      
        column(2, p(HTML("<b>Position</b>"),
               span(shiny::icon("info-circle"), id = "info_uuSim2"),
               
               numericInput("positionSimulate", label = NULL, 80, min = 0, max = 1e6)),
               tippy::tippy_this(elementId = "info_uuSim2",
                                 tooltip = "In kDa. Mean of the truncated gaussian.",placement = "right")),
        
        column(2, p(HTML("<b>Std.</b>"),
               span(shiny::icon("info-circle"), id = "info_uuSim3"),
               
               numericInput("stdSimulate", label = NULL, 20, min = 0, max = 1e6)),
               tippy::tippy_this(elementId = "info_uuSim3",
                                 tooltip = "In kDa. Standard deviation of the truncated gaussian.",placement = "right")),
        
      column(2, p(HTML("<b>Amplitude</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuSim4"),
                  
                  numericInput("amplitudeSimulate", label = NULL, 100, min = 0, max = 1e6)),
             tippy::tippy_this(elementId = "info_uuSim4",
                               tooltip = "In counts. Maximal height of the truncated gaussian.",placement = "right")),
      
      column(2, p(HTML("<b>Left limit</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuSim5"),
                  
                  numericInput("leftLimitSimulate", label = NULL, 30, min = 0, max = 1e6)),
             tippy::tippy_this(elementId = "info_uuSim5",
                               tooltip = "In kDa. Lower bound of the truncated gaussian. 
                               Same parameter as min. observed mass. (Input panel)",placement = "right"))
    
      )
      
    ))


