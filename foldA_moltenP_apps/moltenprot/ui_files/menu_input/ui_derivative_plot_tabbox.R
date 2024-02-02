tabBox(title = "", width = 11,id = "tabset1",
       tabPanel("First derivative",withSpinner(plotlyOutput("signal_der1"))),
       tabPanel("Second derivative",withSpinner(plotlyOutput("signal_der2"))),
       tabPanel("Tm (using the 1st derivative)",withSpinner(plotlyOutput("tm_derivative"))))