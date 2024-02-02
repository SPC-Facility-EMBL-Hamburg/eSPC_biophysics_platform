tabBox(title = "", width = 10,id = "tabset2",
  tabPanel("Fitting Plot", plotlyOutput("fluo_fit_plot")),
  tabPanel("Fitted Parameters",tableOutput("params_table"))
)       


