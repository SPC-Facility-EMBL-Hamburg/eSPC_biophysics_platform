box(title = "Simulation launcher", width = 2, solidHeader = T, status = "primary", 
    fluidRow(
      column(12, p(HTML("<b></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_sim_launcher_2-1"),
                  checkboxInput("runSimulation","Launch simulation",FALSE),
                  
                  tippy::tippy_this(elementId = "info_uu_sim_launcher_2-1",
                    tooltip = "Check this button to compute the results of the simulation.",
                    placement = "right")))
      
    ))


