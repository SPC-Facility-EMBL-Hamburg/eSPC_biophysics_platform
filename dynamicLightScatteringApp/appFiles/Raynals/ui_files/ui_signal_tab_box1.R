tabBox(title = "", width = 6,id = "tabset1",
       tabPanel("Autocorrelation",withSpinner(plotlyOutput("autocorrelation"))),
       tabPanel("Fitted autocorrelation",(plotlyOutput("fittedAutocorrelation1"))))