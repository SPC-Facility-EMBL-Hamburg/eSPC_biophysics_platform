tabBox(title = "", width = 11,id = "tabset2",
       tabPanel("Fitted Parameters",tableOutput("params_table")),
       tabPanel("Fitted Parameters Relative Errors",tableOutput("params_table_errors")),
       tabPanel("Fitting Plots", withSpinner(plotOutput("fluo_fit_plot")))
)