tabBox(title = "", width = 10,id = "tabset1",
       tabPanel("Signal",withSpinner(plotlyOutput("signal"))),
       tabPanel("First Derivative",withSpinner(plotlyOutput("first_der"))),
       tabPanel("Tm versus [Ligand]",withSpinner(plotlyOutput("tm_vs_lig"))),
       tabPanel("Initial fluorescence versus [Ligand]",withSpinner(plotlyOutput("initialFluo_vs_lig"))))

