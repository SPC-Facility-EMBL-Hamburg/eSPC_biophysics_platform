tabBox(title = "", width = 10,id = "tabset2_exploreFhot",
  tabPanel("Kd vs Fhot plot", plotlyOutput("fluo_fit_advanced_plot")),
  tabPanel("Parameters Fhot", tableOutput("params_table_advanced")),
  tabPanel("Fitting plot", plotlyOutput("fluo_fit_plot2")),
  tabPanel("Fitted parameters",tableOutput("params_table2"))

)       


