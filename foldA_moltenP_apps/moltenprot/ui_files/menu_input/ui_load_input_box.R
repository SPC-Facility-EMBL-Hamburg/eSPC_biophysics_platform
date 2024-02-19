box(title = "Input", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(9, p(HTML("<b>1. nDSF or DSF file </b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-1"),
                    fileInput("FLf", NULL,accept = c(".xlsx",".zip",".xls",".txt",'.csv')),
                   tippy::tippy_this(elementId = "info_uu1-1",
                                     tooltip = "Check the User Guide to understand the format of the input files. 
                   Hint: In the case of Prometheus, Panta, or Tycho spreadsheet files, you can combine them in a zip (.zip extension) 
                   and load them together. The zip should have only one type of file.",placement = "right"))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      # This button is not visible in the UI
      column(1, p(HTML("<b><br></b>")),
                  withBusyIndicatorUI(
        shinyjs::hidden(actionButton("Go","2. Load data!",class = "btn-primary")))),
      # End of Little hack to use the withBusyIndicatorUI function (loading spinner)
      
      column(8, p(HTML("<b>Load layout</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-12"),
                  fileInput("layout_file", NULL,accept = c(".xlsx")),
                  tippy::tippy_this(elementId = "info_uu1-12",
                                    tooltip = "If desired, you can upload a xlsx file to complete the condition versus position Table.",placement = "right"))),
      
      column(4, p(HTML("<b>Sort</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-7"),
                  checkboxInput("sort_conditions", "", FALSE),
                  tippy::tippy_this(elementId = "info_uu1-7",
                                    tooltip = "Rearrange the conditions order by the original condition names 
                                            (Useful when analyzing many nDSF files with labels present in the Overview sheet)",placement = "right"))),
      
      column(5, p(HTML("<b>Signal</b>"),
                  selectInput("which", NULL,c("Ratio"="Ratio","350nm" = "350nm",
                                              "330nm" = "330nm","Scattering"="Scattering")))),
      
      column(6, p(HTML("<b>Normalization</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-2"),
                  selectInput("normalization_type", NULL,
                              c("Raw Signal"                 = "Raw_signal",
                                "Divide by initial value"    = "Divide_by_init_value",
                                "Min-max normalization"            = "MC_Normalization",
                                "Divide by area value"       = "area_Normalization")),
                  tippy::tippy_this(elementId = "info_uu1-2",
                                    tooltip = "Rescale the fluorescence signal to facilitate comparison",placement = "right"))),
      
      column(12, p(HTML("<b>\nSignal Window Range (ºC)</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-4"),
                   
                   #Change colour of slider (this code should be re-written in a cleaner way. For now, it works)
                   tags$style(HTML('.js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                            
                                                            background: #00829c;
                                                            border-top: 1px solid #00829c ;
                                                            border-bottom: 1px solid #00829c ;}
          
                                      /* changes the colour of the number tags */
                                     .irs-from, .irs-to, .irs-single { background: #00829c }'
                   )),
                   sliderInput("sg_range", NULL,min = 5, max = 95,value = c(25,90)),
                   tippy::tippy_this(elementId = "info_uu1-4",
                                     tooltip = "To obtain a better fitting, select a temperature range close to the melting(s) transition(s).",placement = "right"))),
      
      column(5, p(HTML("<b>Median filter</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-5"),
                  numericInput('median_filter',NULL, 0,min = 0, max = 6),
                  tippy::tippy_this(elementId = "info_uu1-5",
                                    tooltip = "Apply a rolling median filter to smooth curves and remove spikes, the minimum and maximum values are respectively 0 and 6 degrees",placement = "right"))),
      
      column(5, p(HTML("<b>SG window size</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-6"),
                  numericInput('SG_window2',NULL, 10,min = 0, max = 15),
                  tippy::tippy_this(elementId = "info_uu1-6",
                    tooltip = "The Savitzky-Golay window size (in °C) is used to estimate the first and second derivative",placement = "right"))),
      
      column(4, p(HTML("<b>Select Series</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-8"),
                  selectInput("selected_cond_series", NULL,c("All")),
                  tippy::tippy_this(elementId = "info_uu1-8",
                    tooltip = "Change the Series column labels (from \'A\' to \'B\' for example) and display only the desired conditions",
                    placement = "right"))),
      
    ))