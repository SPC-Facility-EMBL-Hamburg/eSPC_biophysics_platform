box(title = "Settings", width = 5, solidHeader = T, status = "primary", 
    fluidRow(

      column(4, p(HTML("<b>Model</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu2-1"),
                  selectInput("model_selected", NULL,
                              c( "1 Site  "                                         = "one_site",
                                 "2 Sites, 1 Kd, Shared Signal"                     = "2_Sites_1_Kd_Shared_Signal",
                                 "2 Sites, 1 Kd, Different Signal"                  = "2_Sites_1_Kd_Different_Signal",
                                 "2 Sites, 1 Kd, Shared Signal, Cooperativity"      = "2_Sites_1_Kd_Shared_Signal_coop",
                                 "2 Sites, 1 Kd, Different Signal, Cooperativity"   = "2_Sites_1_Kd_Different_Signal_coop",
                                 "2 Sites 2 Kd, Shared Signal"                      = "2_Sites_2_Kd_Shared_Signal",
                                 "2 Sites 2 Kd, Different Signal"                   = "2_Sites_2_Kd_Different_Signal")),
                  
                  
                  tippy::tippy_this(elementId = "info_uu2-1",
                                    tooltip = "If you use a model with 2 Kds, it is better to fix one based on previous knowledge.
                                    Different Signal models should be only used when the signal of P
                                    (free protein) could be the same as PL (Kd2) but different to LP (Kd1) and LPL (double bound).",placement = "right"))),
      
      conditionalPanel("!output.is_csv",
      
      column(4, p(HTML("<b>Signal</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu2-2"),
                  selectInput("signal_selected", NULL,
                              c( "F_hot / F_cold  "       = "f_norm",
                                 "Initial Fluorescence"   = "init_fluo")),

                  tippy::tippy_this(elementId = "info_uu2-2",
                                    tooltip = "If you use the initial fluorescence, 
                                    verify that the ligand alone is not contributing to the signal.",
                                    placement = "right")))
      
      ),
      
      column(4, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "btn_cal",label = "Run Fitting!",
                    icon("meteor"),
                    style="color: #fff; background-color: #337ab7; 
               border-color: #2e6da4")))
      
      ),
    
    fluidRow(
      
      column(width = 4,rHandsontableOutput('experiments2fit1',height="200px")),
      column(width = 4,rHandsontableOutput('experiments2fit2',height="200px")),
      column(width = 4,rHandsontableOutput('experiments2fit3',height="200px"))
      
      
    ))