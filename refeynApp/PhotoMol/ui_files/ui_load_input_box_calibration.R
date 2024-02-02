box(title = "Input", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(6, p(HTML("<b>Activate calibration</b>"),
                  checkboxInput("activateCalibration", "", FALSE))),
      
      conditionalPanel(condition = "input.activateCalibration",
      
      column(6, p(HTML("<b>Calibration method</b>"),
                  selectInput("calibrationMethod",NULL,
                              c("Calibration file" = "calibrationFile",
                                "Custom parameters" = "custom"
                                )))),
      
      conditionalPanel(condition = "output.calibrationMethod == 'custom'",
                       
                       column(6, p(HTML("<b>Slope * 1e6</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalibCustom1"),
                                   
                                   numericInput("slopeCustom", label = NULL, -57, min = -1e6, max = 0)),
                              tippy::tippy_this(elementId = "info_uuCalibCustom1",
                                                tooltip = "Slope obtained by fitting a line to the 
                                                observed ratiometric contransts against known masses (kDa).
                                                ",placement = "right")),
                       
                       column(6, p(HTML("<b>Intercept * 1e6</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalibCustom2"),
                                   
                                   numericInput("interceptCustom", label = NULL, 47, min = 0, max = 1e6)),
                              tippy::tippy_this(elementId = "info_uuCalibCustom2",
                                                tooltip = "Intercept obtained by fitting a line to the 
                                                observed ratiometric contrasts against known masses (kDa).
                                                ",placement = "right"))
                       
                       ),
      
      conditionalPanel(condition = "output.calibrationMethod == 'calibrationFile'",
      
      column(10, p(HTML("<b>MassPhotometry file </b>"),
                   span(shiny::icon("info-circle"), id = "info_uuCalib1-1"),
                   fileInput("massPhotometryFileCalibration", NULL,accept = c(".h5")),
                   tippy::tippy_this(elementId = "info_uuCalib1-1",
                                     tooltip = ".h5 (Hierarchical Data Format) file with a 1D dataset
                                     called 'contrasts'"))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(HTML("<b><br></b>")),
             withBusyIndicatorUI(
               shinyjs::hidden(actionButton("Go2","2. Load data!",class = "btn-primary")))),
      
      conditionalPanel(condition = "output.data_loadedCalibration",
                       
                       column(6, p(HTML("<b>Bin width * 1e3</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalib1-2"),
                                   
                                   numericInput("bin_widthContrast", label = NULL, 0.0004*1e3, min = 0, max = 20)),
                              tippy::tippy_this(elementId = "info_uuCalib1-2",
                                                tooltip = "Used to group the data and build the histogram. 
                                                ",placement = "right")),
                       
      column(8, p(HTML("<b>Initial guesses * 1e3</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuCalib1-12"),
                  
                  textInput("starting_valuesContrast", label=NULL,value="50"),
                  tippy::tippy_this(elementId = "info_uuCalib1-12",
                                    tooltip = "Input one starting value for each gaussian that you want to fit. 
                                 Values should be separated by spaces. Units are contrasts. 
                                    ",placement = "right"))),
      
                       column(8, p(HTML("<b>Known masses</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalib1-5"),
                                   
                                   textInput("knownMasses", label=NULL,value="50"),
                                   tippy::tippy_this(elementId = "info_uuCalib1-5",
                                                     tooltip = "Input one starting value for each truncated gaussian that you want to fit. 
                                 Values should be separated by spaces. Units are kDa.",placement = "right"))),
                       
                       column(12, p(HTML("<b>Window range * 1e3</b>"),
                                    span(shiny::icon("info-circle"), id = "info_uuCalib1-7"),
                                    
                                    sliderInput("window_rangeContrast", NULL,min = -1e3, max = 0,value = c(-600,0))),
                              tippy::tippy_this(elementId = "info_uuCalib1-7",
                                                tooltip = "Set the limits to build the histogram. 
                                                ",placement = "right")),
                       
                       column(5, p(HTML("<b>Slider left limit</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalib1-10"),
                                   
                                   numericInput("leftLimitWindowRangeContrast", label = NULL, -600, min = -1e6, max = 0)),
                              tippy::tippy_this(elementId = "info_uuCalib1-10",
                                                tooltip = "Set the left limit for the window range slider. 
                                 Changing this value will automatically update the selected window range.",placement = "right")),
                       
                       column(5, p(HTML("<b>Slider right limit</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuCalib1-11"),
                                   
                                   numericInput("rightLimitWindowRangeContrast", label = NULL, 0, min = 0, max = 1e6)),
                              tippy::tippy_this(elementId = "info_uuCalib1-11",
                                                tooltip = "Set the right limit for the window range slider. 
                                 Changing this value will automatically update the selected window range.",placement = "right"))
                       
    )))))


