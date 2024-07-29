box(title = "4. Plotting", width = 12, solidHeader = T, status = "primary",
    fluidRow(column(width = 12,rHandsontableOutput('legendInfo'))),
    fluidRow(                       
      column(width = 3,
             p(HTML("<b><br>Choose experiment</b>")),
             selectInput("mol2changeColor", label = NULL,c("X"))),
      
      column(width = 3,
             p(HTML("<b><br>Set colour</b>")),
             colourpicker::colourInput("colorForLegend",label=NULL, value = "#E41A1C")),
      
      column(4, p(HTML("<b><br>Show plot export options</b>")),
                  checkboxInput('showAdvancedPlottingOptions',NULL,FALSE))
      ),
      
    fluidRow( 
    
      conditionalPanel(
        "input.showAdvancedPlottingOptions",
        
        column(3, p(HTML("<b>File type</b>"),
                    selectInput("plot_type", NULL,
                                c("PNG"    = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
        
        column(3, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotWidth"),
                    numericInput('plot_width',NULL, 12,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotWidth",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(3, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotHeight"),
                    numericInput('plot_height',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotHeight",
                                      tooltip = "Units are pixels * 50",
                                      placement = "right"))),                     
        
        column(3, p(HTML("<b>Text size</b>"),
                    numericInput('plot_axis_size',NULL, 16,min = 4, max = 40)))
        
        )
      
    ))