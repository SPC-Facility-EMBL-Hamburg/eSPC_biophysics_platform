box(title = "Plot options", width = 6, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(4,p(HTML("<b>Type</b>"),
                 span(shiny::icon("info-circle"), id = "info_uu3-1"),
                 selectInput("select_plot_type", NULL,c("Fraction Unfolded")),
                 tippy::tippy_this(elementId = "info_uu3-1",
                                   tooltip = "Select type of plot to display in the \"Results Plot\" panel",placement = "right"))),
      
      conditionalPanel(condition = "input.show_plot_download_options_analyze",
      
        column(2, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu3-2"),
                    numericInput('plot_width_results',NULL, 18,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uu3-2",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(2, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu3-3"),
                    numericInput('plot_height_results',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uu3-3",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(3, p(HTML("</b>File type</b>"),
                    selectInput("plot_type_results", NULL,
                                c("PNG"                 = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
        
        column(3, p(HTML("<b>Legend text size</b>"),
                    numericInput('plot_font_size_results',NULL, 13,min = 6, max = 40))),
                
        column(2, p(HTML("<b>Axis text size</b>"),
                    numericInput('plot_axis_size_results',NULL, 14,min = 4, max = 40)))
        ),
      
      column(6, p(HTML("<b>Show Download Options</b>"),
                  checkboxInput("show_plot_download_options_analyze", "", FALSE)))
      
    ))