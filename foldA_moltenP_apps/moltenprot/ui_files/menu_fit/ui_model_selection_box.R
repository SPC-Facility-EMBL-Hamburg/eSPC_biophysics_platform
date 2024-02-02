box(title = "2.1 Model selection", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      column(10, p(HTML("<b></b>"),
                   span(shiny::icon("info-circle"), id = "info_uu2-1"),
                   selectInput("model_selected", NULL,
                               c("Equilibrium Two State"      = "EquilibriumTwoState",
                                 "Empirical Two State"        = "EmpiricalTwoState",
                                 "Equilibrium Three State"    = "EquilibriumThreeState",
                                 "Empirical Three State"      = "EmpiricalThreeState",
                                 "Irreversible Two State"     = "IrreversibleTwoState")),
                   tippy::tippy_this(elementId = "info_uu2-1",
                                     tooltip = "More information about the models can be found in the User Guide section.
                                     The result of the fitting procedure is presented in the 'Fitted conditions' Table. 
                                     If the model fails to fit the data, we suggest first changing the temperature window,
                                     and second, changing the 'Temperature range for baseline estimation'.",placement = "right"))),
      
      column(6,   withBusyIndicatorUI(
        actionButton("btn_cal","Run Fitting",class = "btn-primary")))
    ))