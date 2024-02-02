box(title = "Plot download options", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(1, p(HTML("<b>Show</b>"),
                  checkboxInput("show_plot_download_options_input_fit", "", TRUE))),
      
      conditionalPanel(condition = "input.show_plot_download_options_input_fit",
        
         column(2, p(HTML("<b>Plot title</b>"),
                     span(shiny::icon("info-circle"), id = "info_uu1-15_fit"),
                     textInput("fittingPlotTitle", NULL, value = "", width = NULL),
                     tippy::tippy_this(elementId = "info_uu1-15_fit",
                                       tooltip = "Use a descriptive title. For example, 
                                       Ligand T1.",placement = "right"))),
         
         column(1, p(HTML("<b>Kd in title</b>"),
                     checkboxInput("show_kd_in_title_fitPlot", "", TRUE))),

         column(2, p(HTML("<b>X-axis label</b>"),
                     textInput("xAxisLabel", NULL,
                               value = "Binding Partner Concentration", width = NULL))),

         column(2, p(HTML("<b>Y-axis label</b>"),
                     textInput("yAxisLabel", NULL, value = "Signal (AU)", width = NULL))),
                                
        column(1, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-13_fit"),
                    numericInput('plot_width_fit',NULL, 14,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uu1-13_fit",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(1, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-14_fit"),
                    numericInput('plot_height_fit',NULL, 10,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uu1-14_fit",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(1, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_fit", NULL,
                                c("PNG"                 = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
      
        column(1, p(HTML("<b>Axis size</b>"),
                    numericInput('plot_axis_size_fit',NULL, 16,min = 4, max = 40)))
      
        )),
        
        conditionalPanel(condition = "input.show_plot_download_options_input_fit",
                     
        fluidRow(
        
        column(1, p(HTML("<b>Legend size</b>"),
                    numericInput('plot_legend_text_size_fit',NULL, 18,min = 4, max = 40))),
        
        column(1, p(HTML("<b>Log Axis type</b>"),
                    selectInput("logScaleType_fit", NULL,
                                c("Micromolar"     =  "micromolar",
                                  "Milimolar"      =  "milimolar",
                                  "Molar"          =  "molar")))),
          
        column(2, p(HTML("<b>X-axis upper limit</b>"),
                    numericInput('xAxis_upperLimit_fitPlot',NULL, 1e6))),
        
        column(2, p(HTML("<b>X-axis lower limit</b>"),
                    numericInput('xAxis_lowerLimit_fitPlot',NULL, 0))),
        
        column(1, p(HTML("<b># X-axis ticks</b>"),
                    numericInput('xAxis_ticksNumber_fitPlot',NULL, 2,min = 0, max = 100))),
        
        column(2, p(HTML("<b>Y-axis upper limit</b>"),
                    numericInput('yAxis_upperLimit_fitPlot',NULL, 1))),
        
        column(2, p(HTML("<b>Y-axis lower limit</b>"),
                    numericInput('yAxis_lowerLimit_fitPlot',NULL, 0))),
        
        column(1, p(HTML("<b># Y-axis ticks</b>"),
                    numericInput('yAxis_ticksNumber_fitPlot',NULL, 3,min = 0, max = 100)))
        
        ),
        
        fluidRow(
        
        column(1, p(HTML("<b>Autoscale</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuAutoScale"),
                    checkboxInput("autoscale", "", FALSE),
                    tippy::tippy_this(
                      elementId = "info_uuAutoScale",
                      tooltip = "Force all curves to start around 0. 
                      This is done by substracting the signal corresponding
                      to the uncomplexed protein (parameter RF1 * protein concentration). If the experiment 
                      was defined as a control experiment (wasn't fitted)
                      we subtract the mean value of the first 3 points.
                      Only useful for visualization purposes.",placement = "right"))),
        
        column(1, p(HTML("<b>Display sd</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuDisplaySD"),
                    checkboxInput("displaySD", "", TRUE),
                    tippy::tippy_this(
                      elementId = "info_uuDisplaySD",
                      tooltip = "Show data points as mean ± sd (standard deviation) 
                      in the 'Fitting Plot'.
                      Only useful for visualization purposes.",placement = "right")))
        
        
        )
      
      ))
      
    

