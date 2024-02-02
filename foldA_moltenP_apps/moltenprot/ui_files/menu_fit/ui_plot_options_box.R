box(title = "Plot options", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(9,p(HTML("<b></b>"),
                 span(shiny::icon("info-circle"), id = "info_uu222-4"),
                 
                 selectInput("select_fitting_plot", NULL,c("Fitted conditions" = "Fit_1")),
                 
                 tippy::tippy_this(elementId = "info_uu222-4",
                                   tooltip = "Select which fittings to display. 
                                   This option affects only the 'Fitting' and 'Fitting residuals' windows.",placement = "right"))),
      
      conditionalPanel(condition = "output.data_was_fitted",
      column(9,p(HTML("<b></b>"),
                 span(shiny::icon("info-circle"), id = "info_uu222-5"),
                       withBusyIndicatorUI(downloadButton("download_fit_plots", "Download fitting plots")),
                 tippy::tippy_this(elementId = "info_uu222-5",
                  tooltip = "Create a zip file containing the fitting plots, 
                  including one PNG file every 16 fitted conditions (this function can take up to 30-60 seconds per PNG file).",
                  placement = "right"))))

    ))
