tabBox(title = "", width = 10,id = "tabset3",
       tabPanel("Score Table",tableOutput("score_table")),
       tabPanel("Results Plot", withSpinner(plotlyOutput("results_plot")))
)