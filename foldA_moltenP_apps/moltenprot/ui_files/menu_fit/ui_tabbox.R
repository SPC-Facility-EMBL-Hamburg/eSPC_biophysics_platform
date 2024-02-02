tabBox(title = "", width = 11,id = "tabset2",
       tabPanel("Fitted Conditions",tableOutput("fitted_conditions_table")),
       tabPanel("Fitted Parameters (kcal °C mol)",tableOutput("params_table")),
       tabPanel("Parameters Relative Errors",tableOutput("params_table_errors")),
       tabPanel("Fitting", withSpinner(plotOutput("fluo_fit_plot"))),
       tabPanel("Fitting residuals", withSpinner(plotOutput("fluo_residuals_plot")))
)