box(title = "Analysis", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3, p(HTML("<b><br></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-analysisBox1"),
                  actionButton("updateInfo", "Run fitting"),
                  tippy::tippy_this(elementId = "info_uu-analysisBox1",
                                    tooltip = "Apply changes done in the Samples Tables, fit the data again, 
                                    and reload the plots.",
                                    placement = "right"))),
      
      column(3, p(HTML("<b>Group by</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-analysisBox2"),
                  selectInput("splitFactor", NULL,
                              c("None"          =  "None",
                                "Condition"     =  "Condition",
                                "Read"          =  "Read",
                                "Scan"          =  "Scan",
                                "Experiment"    =  "Experiment")),
                  tippy::tippy_this(elementId = "info_uu-analysisBox2",
                                    tooltip = "Colour the plots according to different groups. 
                                    If no change is observed (but you expect it), it's because the Samples Tables was modified: 
                                    Press 'Update'",
                                    placement = "right"))),
      
      column(3, p(HTML("<b>Residuals filter</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-analysisBox3"),
                  numericInput("filterByResidualsValue",label = NULL, 1, min = 1, max = 10,step=0.02),
                  tippy::tippy_this(elementId = "info_uu-analysisBox3",
                                    tooltip = "Remove from the plots the autocorrelation curves 
                                    where the (absolut) maximum value of the residuals is higher than
                                    this value.",
                                    placement = "right"))),
      
      column(3, p(HTML("<b>Show Advanced</b>"),
                  checkboxInput("show_advanced_fitting_opts", "", FALSE)))
      
    ),
    
    conditionalPanel(
      condition = "input.show_advanced_fitting_opts",
      fluidRow(

        column(3, p(HTML("<b>Max time (μs)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-advanced_fitting_opts4"),
                    numericInput('maxTime',NULL, 1e7,min = 0, max = 1e12),
                    tippy::tippy_this(elementId = "info_uu1-advanced_fitting_opts4",
                                      tooltip = "Perform the fitting using only the autocorrelation 
                                       before certain time.",placement = "right"))),
        
      column(4, p(HTML("<b>Reg method</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-advanced_fitting_opts6"),
                  selectInput("regMethod", NULL,
                              c("L-curve"     =  "lCurve",
                                "Fixed alpha" =  "fixedAlpha"
                                )),
                  
                  tippy::tippy_this(elementId = "info_uu1-advanced_fitting_opts6",
                                    tooltip = "Select L-curve to determine the optimal
                                    value of alpha using an empiric criteria.",placement = "right"))),
      
      conditionalPanel("input.regMethod == 'fixedAlpha'",
      
        column(4, p(HTML("<b>Reg term</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-advanced_fitting_opts5"),
                    numericInput('alphaRegTerm',NULL, 0.05,min = 0.00, max = 1,step=0.001),
                    tippy::tippy_this(
                      elementId = "info_uu1-advanced_fitting_opts5",
                      tooltip = "Regularization term for smoothness.",
                      placement = "right"))))),
    
    fluidRow(
      conditionalPanel("input.regMethod != 'fixedAlpha'",
      
       column(4, p(HTML("<b>Reg terms Start</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-advanced_fitting_opts7"),
                   numericInput('alphaRegTermLcurve_start',NULL, -6,min = 0.00, max = 1,step=0.2),
                   tippy::tippy_this(
                     elementId = "info_uu1-advanced_fitting_opts7",
                     tooltip = "Initial value for the list of regularisation terms that we 
                     want to evaluate. This sequence is built using the following formula:
                     5**seq(start,stop,step))**2",
                     placement = "right"))),
       
       column(4, p(HTML("<b>Reg terms Stop</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-advanced_fitting_opts8"),
                   numericInput('alphaRegTermLcurve_stop',NULL, 2.5,min = 0.00, max = 10,step=0.2),
                   tippy::tippy_this(
                     elementId = "info_uu1-advanced_fitting_opts8",
                     tooltip = "Stop value for the list of regularisation terms that we 
                     want to evaluate. This sequence is built using the following formula:
                     5**seq(start,stop,step))**2",
                     placement = "right"))),
       
       column(4, p(HTML("<b>Reg terms Step</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-advanced_fitting_opts9"),
                   numericInput('alphaRegTermLcurve_step',NULL, 0.25,min = 0.00, max = 1,step=0.05),
                   tippy::tippy_this(
                     elementId = "info_uu1-advanced_fitting_opts9",
                     tooltip = "Step value for the list of regularisation terms that we 
                     want to evaluate. This sequence is built using the following formula:
                     5**seq(start,stop,step))**2",
                     placement = "right")))
       
       )                       
                       

        )
    
    )
)
    
    