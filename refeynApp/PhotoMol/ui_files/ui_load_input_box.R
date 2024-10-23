box(title = "Input", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(10, p(HTML("<b>Mass Photometry file(s) (up to five) </b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-1"),
                   fileInput("massPhotometryFile", NULL,accept = c(".h5",".csv"),multiple=TRUE),
                   tippy::tippy_this(elementId = "info_uu1-1",
                                     tooltip = ".h5 (Hierarchical Data Format) file with a 1D dataset called 'masses_kDa'. You can load
                                     up to 5 files simultaneously."))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(HTML("<b><br></b>")),
             withBusyIndicatorUI(
               shinyjs::hidden(actionButton("Go","2. Load data!",class = "btn-primary")))),
      
      column(12, p(HTML(""),
                   actionButton("GoLoadExample","Load example data",class = "btn-primary"))),
      
      conditionalPanel(condition = "output.dataLoaded", 
      
        column(4, p(HTML("<b>Bin width</b>"),
               span(shiny::icon("info-circle"), id = "info_uu1-2"),
               
               numericInput("bin_width", label = NULL, 6, min = 1, max = 20)),
               tippy::tippy_this(elementId = "info_uu1-2",
                                 tooltip = "In kDa. Used to group the data and build the histogram.",placement = "right")),
        
        column(6, p(HTML("<b>Min. observed mass</b>"),
               span(shiny::icon("info-circle"), id = "info_uu1-8"),
               
               numericInput("min_observed_mass", label = NULL, 30, min = -500, max = 500)),
               tippy::tippy_this(elementId = "info_uu1-8",
                                 tooltip = "Minimum mass (kDa) that could be theoretically observed.
                                 Lower (absolute) values are not used for the fitting.",placement = "right")),
        
        column(6, p(HTML("<b>Upper limit for the standard deviation</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-3"),
                    
                    sliderInput("upper_limit_std", NULL,min = 5, max = 200,value = 30)),
                    tippy::tippy_this(elementId = "info_uu1-3",
                                      tooltip = "The fitted standard deviations (kDa) will not be higher than this value.",placement = "right")),
        
        column(6, p(HTML("<b>Tolerance to the initial guesses</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-4"),
                    
                    sliderInput("position_tolerance", NULL,min = 1, max = 200,value = 30)),
               tippy::tippy_this(elementId = "info_uu1-4",
                                 tooltip = "Used to set the maximum accepted deviance for the fitted means (kDa). For example,
                                 if we set this value to 10, and as starting values we use '50', the fitted mean won't be lower
                                 than 40 or higher than 60.",placement = "right")),
        
        column(8, p(HTML("<b>Starting values</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-5"),
                    
                    textInput("starting_values1", label=NULL,value=""),
                    tippy::tippy_this(elementId = "info_uu1-5",
                                 tooltip = "Input one starting value for each truncated gaussian that you want to fit. 
                                 The starting values should be separated by spaces, 
                                 and, in absolute units, higher than the minimum observed mass. Units are kDa.",placement = "right"))),
      
        column(4, p(HTML("<b>Baseline</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-9"),
                    
                    numericInput("baseline", label = NULL, 0, min = 0, max = 20)),
               tippy::tippy_this(elementId = "info_uu1-9",
                                 tooltip = "Useful when there is a constant background noise. More info in the User guide.",placement = "right")),

        conditionalPanel(condition = "output.nFiles > 1",

        column(8, p(HTML("<b>Starting values (file 2)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-St2"),

                    textInput("starting_values2", label=NULL,value=""),
                    tippy::tippy_this(elementId = "info_uu1-St2",
                                 tooltip = "Input one starting value for each truncated gaussian that you want to fit.
                                 The starting values should be separated by spaces,
                                 and, in absolute units, higher than the minimum observed mass. Units are kDa.",placement = "right")))),

        conditionalPanel(condition = "output.nFiles > 2",

        column(8, p(HTML("<b>Starting values (file 3)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-St3"),

                    textInput("starting_values3", label=NULL,value=""),
                    tippy::tippy_this(elementId = "info_uu1-St3",
                                 tooltip = "Input one starting value for each truncated gaussian that you want to fit.
                                 The starting values should be separated by spaces,
                                 and, in absolute units, higher than the minimum observed mass. Units are kDa.",placement = "right")))),

        conditionalPanel(condition = "output.nFiles > 3",

        column(8, p(HTML("<b>Starting values (file 4)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-St4"),

                    textInput("starting_values4", label=NULL,value=""),
                    tippy::tippy_this(elementId = "info_uu1-St4",
                                 tooltip = "Input one starting value for each truncated gaussian that you want to fit.
                                 The starting values should be separated by spaces,
                                 and, in absolute units, higher than the minimum observed mass. Units are kDa.",placement = "right")))),

        conditionalPanel(condition = "output.nFiles > 4",

        column(8, p(HTML("<b>Starting values (file 5)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-St5"),

                    textInput("starting_values5", label=NULL,value=""),
                    tippy::tippy_this(elementId = "info_uu1-St5",
                                 tooltip = "Input one starting value for each truncated gaussian that you want to fit.
                                 The starting values should be separated by spaces,
                                 and, in absolute units, higher than the minimum observed mass. Units are kDa.",placement = "right")))),

        column(12, p(HTML("<b>Window range</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-7"),
                    
                    sliderInput("window_range", NULL,min = 0, max = 0,value = c(0,0))),
               tippy::tippy_this(elementId = "info_uu1-7",
                                 tooltip = "Set the limits to build the histogram. Units are kDa.",placement = "right")),
        
        column(4, p(HTML("<b>Slider left limit</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-10"),
                    
                    numericInput("leftLimitWindowRange", label = NULL, value=0, min = -1e6, max = 0)),
               tippy::tippy_this(elementId = "info_uu1-10",
                                 tooltip = "Set the left limit for the window range slider. 
                                 Changing this value will automatically update the selected window range.",placement = "right")),
        
        column(4, p(HTML("<b>Slider right limit</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-11"),
                    
                    numericInput("rightLimitWindowRange", label = NULL, value=0, min = 0, max = 1e6)),
               tippy::tippy_this(elementId = "info_uu1-11",
                                 tooltip = "Set the right limit for the window range slider. 
                                 Changing this value will automatically update the selected window range.",placement = "right")),

        column(4, p(HTML("<b>Automatic fitting</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu1-12"),

                    checkboxInput("automaticFit", label = NULL, value = TRUE)),
               tippy::tippy_this(
               elementId = "info_uu1-12",
               tooltip = "Re-fit the data automatically after changing any parameter (e.g., bin width, starting values, etc.).",
               placement = "right"))
           
      )   
    ))


