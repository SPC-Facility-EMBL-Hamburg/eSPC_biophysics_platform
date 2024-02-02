box(title = "Input", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(9, p(HTML("<b>1. MST file </b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-1"),
                   fileInput("FLf", NULL,accept = c(".xlsx",".xls",".csv",".zip")),
                   tippy::tippy_this(elementId = "info_uu1-1",
                                     tooltip = "MST file from Nanothemper or csv with no header and two columns: Ligand Concentration | Signal.
                                     Hint: You can load a zip file with many MST files (.xlsx)",placement = "right"))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(HTML("<b><br></b>")),
             withBusyIndicatorUI(
               shinyjs::hidden(actionButton("Go","2. Load data!",class = "btn-primary")))),
      
      column(12, p(HTML(""),
                   actionButton("GoLoadExample","Load example data",class = "btn-primary"))),
      
      column(6, p(HTML("<b><br>Units</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-2"),
                  selectInput("conc_units", NULL,c(
                    "Molar"="molar","Milimolar" = "milimolar","Micromolar" = "micromolar","Nanomolar"="nanomolar")),
                  tippy::tippy_this(elementId = "info_uu1-2",
                                    tooltip = "Ligand & protein concentration units (displayed in the 
                                            Position versus Concentration Table)",placement = "right"))),
      
      conditionalPanel("!output.is_csv",
      
      column(6, p(HTML("<b><br>Normalization</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-33"),
                  selectInput("normalization_type", NULL,
                              c("Raw Signal"                 = "Raw_signal",
                                "Divide by the mean of the cold region (negative time)" = "Divide_by_init_value")),
                  tippy::tippy_this(elementId = "info_uu1-33",
                                    tooltip = "Normalize the raw fluorescence. Only useful for trace versus time datasets.",placement = "right"))),
      
      column(6, p(HTML("<b>Cold Region</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-3"),
                   
                   sliderInput("cold_range", NULL,min = -10, max = 0,value = c(-1,0)),
                   tippy::tippy_this(elementId = "info_uu1-3",
                                     tooltip = "Select the cold region",placement = "right"))),
      
      column(6, p(HTML("<b>Hot Region</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-34"),
                  
                  sliderInput("hot_range", NULL,min = 5, max = 30,value = c(20,21)),
                  tippy::tippy_this(elementId = "info_uu1-34",
                                    tooltip = "Select the hot region",placement = "right"))),
      
      column(5, p(HTML("<b>Median filter</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-4"),
                  numericInput('median_filter',NULL, 0,min = 0, max = 4),
                  tippy::tippy_this(elementId = "info_uu1-4",
                                    tooltip = "Apply a rolling median filter to smooth curves and remove spikes. 
                                    The minimum and maximum values are respectively 0 and 4 seconds",placement = "right"))),
      
      conditionalPanel(condition = "input.fill_table",
                       column(4, p(HTML("<b># Replicates</b>"),numericInput('n_replicates',NULL, 2,min = 1, max = 6)))),
      
      conditionalPanel(condition = "input.fill_table", 
                       column(4, p(HTML("<b>Dilution Factor</b>"),numericInput('dil_factor',NULL, 2,min = 1, max = 10)))),
      
      conditionalPanel(condition = "input.fill_table",
                       column(5, p(HTML("<b>Initial Ligand Concentration</b>"),
                                   numericInput('initial_ligand',NULL, 500.0,min = 1, max = 5000)))),
      
      conditionalPanel(condition = "input.fill_table", 
                       column(4, p(HTML("<b>Reverse Order</b>"),checkboxInput("rev_order", "", FALSE)))),
      
      column(5, p(HTML("<b>Autofill Table</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-5"),
                  checkboxInput("fill_table", "", FALSE),
                  tippy::tippy_this(elementId = "info_uu1-5",
                                    tooltip = "Select this option to complete the \'Position versus Concentration\'
                                            Table using a constant dilution factor
                                            and an initial ligand concentration",placement = "right")))
      
    )
      
    ))