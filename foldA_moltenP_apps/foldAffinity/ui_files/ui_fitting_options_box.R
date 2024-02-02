box(title = "2.2 Fitting options", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(6, p(HTML("<b>ΔCP (kcal/K/mol)</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu2-2"),
                  numericInput('delta_cp',NULL, 0,min = -30, max = 30),
                  tippy::tippy_this(elementId = "info_uu2-2",
                                    tooltip = "Specify the Heat capacity change upon unfolding. 
                                    If unknown, assume it is zero, or approximate it by the number of residues 'nRes' using 13.88 * 'nRes' / 1000.",placement = "right")))
    ))