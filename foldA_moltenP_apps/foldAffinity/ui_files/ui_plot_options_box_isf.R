box(title = "Plot download options", width = 2, solidHeader = T, status = "primary", 
    fluidRow(
      
      conditionalPanel(condition = "input.show_plot_download_options_input_isf",
      
        column(6, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-13_isf"),
                    numericInput('plot_width_isf',NULL, 18,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uu1-13_isf",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(6, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-14_isf"),
                    numericInput('plot_height_isf',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uu1-14_isf",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(6, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_isf", NULL,
                                c("PNG"                 = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
      
        column(6, p(HTML("<b>Axis text size</b>"),
                    numericInput('plot_axis_size_isf',NULL, 14,min = 4, max = 40))),
        
        column(6, p(HTML("<b>Log Axis type</b>"),
                    selectInput("logScaleType_isf", NULL,
                                c("Micromolar"     =  "micromolar",
                                  "Milimolar"      =  "milimolar",
                                  "Molar"          =  "molar",
                                  "Default Plotly" =  "defaultPlotly"))))
      
      ),
      
      column(6, p(HTML("<b>Show</b>"),
                  checkboxInput("show_plot_download_options_input_isf", "", FALSE)))

      
    ))

