box(title = "Input", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(9, p(HTML("<b>1. nDSF or DSF file </b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-1"),
                   fileInput("FLf", NULL,accept = c(".xlsx",".zip",".xls",
                                                    ".csv",".txt",'.json')),
                   tippy::tippy_this(elementId = "info_uu1-1",
                   tooltip = "Check the User Guide to understand the format of the input files. 
                   Hint: Multiple Prometheus and Panta spreadsheet files can be combined in zip an loaded
                   together.",placement = "right"))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(HTML("<b><br></b>"),
             withBusyIndicatorUI(
               shinyjs::hidden(actionButton("Go","2. Load data!",class = "btn-primary"))))),
      
      column(12, p(HTML(""),
               actionButton("GoLoadExample","Load example data",class = "btn-primary"))),
      
      column(6, p(HTML("<b><br>Units</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-2"),
                  selectInput("conc_units", NULL,c(
                    "Molar"="Molar","Milimolar" = "Milimolar","Micromolar" = "Micromolar","Nanomolar"="Nanomolar")),
                  tippy::tippy_this(elementId = "info_uu1-2",
                                    tooltip = "Ligand concentration units (displayed in the 
                                            Position versus Concentration Table)",placement = "right"))),
      
      column(4, p(HTML("<b><br>Signal</b>"),
                  selectInput("which", NULL,
                              c("Ratio"="Ratio","350nm" = "350nm","330nm" = "330nm","Scattering"="Scattering")))),
      
      column(12, p(HTML("<b>\nSignal Window Range (ºC)</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-3"),
                   
                   #Change colour of slider (this code should be re-written in a cleaner way. For now, it works)
                   tags$style(HTML('.js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                            
                                                            background: #00829c;
                                                            border-top: 1px solid #00829c ;
                                                            border-bottom: 1px solid #00829c ;}
          
                                      /* changes the colour of the number tags */
                                     .irs-from, .irs-to, .irs-single { background: #00829c }'
                   )),
                   sliderInput("sg_range", NULL,min = 5, max = 95,value = c(25,90)),
                   tippy::tippy_this(elementId = "info_uu1-3",
                                     tooltip = "Remove temperature data far from the melting transition. 
                We recommend selecting a temperature range that covers not more than 45 degrees",placement = "right"))),
      
      column(6, p(HTML("<b>Median filter</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-4"),
                  numericInput('median_filter',NULL, 0,min = 0, max = 6),
                  tippy::tippy_this(elementId = "info_uu1-4",
                                    tooltip = "Apply a rolling median filter to smooth curves and remove spikes, the minimum and maximum values are respectively 0 and 6 degrees",placement = "right"))),
      
      conditionalPanel(condition = "input.fill_table",
                       column(4, p(HTML("<b># Replicates</b>"),numericInput('n_replicates',NULL, 2,min = 1, max = 6)))),
      
      conditionalPanel(condition = "input.fill_table", 
                       column(4, p(HTML("<b>Dilution Factor</b>"),numericInput('dil_factor',NULL, 2,min = 1, max = 10)))),
      
      conditionalPanel(condition = "input.fill_table",
                       column(5, p(HTML("<b>Initial Ligand Concentration</b>"),
                                   numericInput('initial_ligand',NULL, 500.0,min = 1, max = 5000)))),
      
      conditionalPanel(condition = "input.fill_table", 
                       column(4, p(HTML("<b>Reverse Order</b>"),checkboxInput("rev_order", "", FALSE)))),
      
      column(6, p(HTML("<b>Autofill Table</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-5"),
                  checkboxInput("fill_table", "", FALSE),
                  tippy::tippy_this(elementId = "info_uu1-5",
                                    tooltip = "Select this option to complete the \'Position versus Concentration\'
                                            Table using a constant dilution factor
                                            and an initial ligand concentration",placement = "right")))
      
      
    ))