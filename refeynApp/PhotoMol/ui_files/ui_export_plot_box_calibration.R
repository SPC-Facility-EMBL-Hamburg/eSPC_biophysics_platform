box(title = "Plot download options", width = 2, solidHeader = T, status = "primary", 
    fluidRow(
      
      conditionalPanel(condition = "input.show_plot_download_options_inputCalibration",
                       
                       column(6, p(HTML("<b>Width</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalib1-13"),
                                   numericInput('plot_widthCalibration',NULL, 12,min = 1, max = 100),
                                   tippy::tippy_this(elementId = "info_uuCalib1-13",
                                                     tooltip = "Units are pixels * 50",placement = "right"))),
                       
                       column(6, p(HTML("<b>Height</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalib1-14"),
                                   numericInput('plot_heightCalibration',NULL, 10,min = 1, max = 100),
                                   tippy::tippy_this(elementId = "info_uuCalib1-14",
                                                     tooltip = "Units are pixels * 50",placement = "right"))),
                       
                       column(6, p(HTML("<b>File type</b>"),
                                   selectInput("plot_typeCalibration", NULL,
                                               c("PNG"                 = "png",
                                                 "SVG"    = "svg",
                                                 "JPEG"    = "jpeg")))),
                       
                       column(6, p(HTML("<b>Axis text size</b>"),
                                   numericInput('plot_axis_sizeCalibration',NULL, 18,min = 4, max = 40)))
                       
      ),
      
      column(6, p(HTML("<b>Show</b>"),
                  checkboxInput("show_plot_download_options_inputCalibration", "", FALSE)))
      
      
    ))