tabBox(title = "", width = 10,id = "tabset2",
  tabPanel("Fitting plot", plotlyOutput("fluo_fit_plot")),
  tabPanel("Fitted parameters",tableOutput("params_table"))
)       


