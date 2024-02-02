box(title = "Report", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(10,p(HTML("<b></b>")),
             textInput("filename_report",
                       label="Report filename",value = paste0("MP_report_", Sys.Date()))),
      column(6,
              span(shiny::icon("info-circle"), id = "info_uu4-1"),
             withBusyIndicatorUI(actionButton("downloadReport","Generate",class = "btn-primary")),
             tippy::tippy_this(elementId = "info_uu4-1",
             tooltip = "Set the report name before creating it. The results plots are only available for 
             the Equilibrium Two State and Empirical Two State models.",placement = "right")),
      
      conditionalPanel(condition = "output.report_was_created",
                       column(10,p(HTML("<b><br></b>"),
                                  downloadButton("downloadReportHidden",label = "Download")))
      ),
    
      conditionalPanel(condition = "output.data_loaded",
      
        column(10, p(HTML("<b><br>Include</b>"), checkboxInput("report_include_signal"         , "Signal Plot"               , TRUE  ))),
        column(10, p(HTML("<b></b>")       ,     checkboxInput("report_include_der"            , "Derivative Plot"           , TRUE  ))),
        column(10, p(HTML("<b></b>")       ,     checkboxInput("report_include_max_der"        , "Tm from Derivative Plot" , FALSE )))
      
      ),
      
      conditionalPanel(condition = "output.data_was_fitted",
      
        column(10, p(HTML("<b></b>")       , checkboxInput("report_include_fitted_params_table", "Fitted Parameters Table"   , FALSE  ))),
        
        column(10, p(HTML("<b></b>")       , checkboxInput("report_include_fitted_params_errors_table", 
                                                           "Fitted Parameters Errors Table"           ,        FALSE  ))),
        
        column(10, p(HTML("<b></b>")       , checkboxInput("report_include_fitting_plots", "Fitting Plots"   , FALSE  ))),
        column(10, p(HTML("<b></b>")       , checkboxInput("report_include_results_plots", "Results Plots"   , FALSE  )))
      
      )

    ))

