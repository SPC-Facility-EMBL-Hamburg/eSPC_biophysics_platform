box(title = "1. Input", width = 3, solidHeader = T, status = "primary",
    
    fluidRow(
      
      column(8, p(HTML("<b>Input mode</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-inputMode"),
                  selectInput("inputMode", NULL,
                              c('Custom'  = 'custom',
                                'Automatic baseline subtraction'  = 'automaticBaselineSub')),
                  tippy::tippy_this(
                    elementId = "info_uu-inputMode",
                    tooltip = "Select custom if you want to process the spectra manually.
                    Select 'Automatic baseline subtraction' if you want to load the CD scans of the sample,
                    the CD scans of the baseline  and process them automatically. 
                    This includes 1)
                    averaging the scans, 2) subtracting the baseline, 3) zeroing the spectrum, and 4)
                    normalising the final spectrum (e.g., by converting it to molar ellipticity), if desired.
                    ",placement = "right"))),
      
    ),
    
    conditionalPanel("input.inputMode == 'custom'",
                     fluidRow(
                       
                       column(9, p(HTML("<b>CD data files </b>"),
                                   span(shiny::icon("info-circle"), id = "info_uu1-1"),
                                   fileInput("cdFiles", NULL,multiple = TRUE),
                                   tippy::tippy_this(
                                     elementId = "info_uu1-1",
                                     tooltip = "Check the User Guide for an explanation of the input data format.
                     Hint: You can load a zip with any combination of accepted input files.
                     ",placement = "right"))),
                       
                       # Little hack to use the withBusyIndicatorUI function (loading spinner)
                       column(3, p(HTML("<b></b>"),
                                   withBusyIndicatorUI(
                                     shinyjs::hidden(actionButton("Go","",class = "btn-primary")))))
                       
                     )
                     ),

    conditionalPanel("input.inputMode == 'automaticBaselineSub'",
                     fluidRow(
                       
                       column(6, p(HTML("<b>Sample files </b>"),
                                   span(shiny::icon("info-circle"), id = "info_uu-cdFilesSample"),
                                   fileInput("cdFilesSample", NULL,multiple = TRUE),
                                   tippy::tippy_this(
                                     elementId = "info_uu-cdFilesSample",
                                     tooltip = "Check the User Guide for an explanation of the input data format.
                     Hint: You can load a zip with any combination of accepted input files.
                     ",placement = "right"))),

                       column(6, p(HTML("<b>Baseline files </b>"),
                                   span(shiny::icon("info-circle"), id = "info_uu-cdFilesBaseline"),
                                   fileInput("cdFilesBaseline", NULL,multiple = TRUE),
                                   tippy::tippy_this(
                                     elementId = "info_uu-cdFilesBaseline",
                                     tooltip = "Check the User Guide for an explanation of the input data format.
                     Hint: You can load a zip with any combination of accepted input files.
                     ",placement = "right"))),
                       
                       column(3, p(HTML("<b></b>"))),
                       
                       column(4, p(HTML("<b></b>"),
                              actionButton(inputId = "automaticProcess",
                                           label = "Process data!",
                                           icon("caret-right"), 
                                           style="color: #FFFFFF; background-color: #00829c; 
               border-color: #00829c")))  
                       
                     )
    ),
          
    fluidRow(
    
      column(7, p(HTML("<b>Select experiment to remove</b>"),
                  selectInput("experiment2delete", NULL,
                              c(	'None'        = 'None')))),
      
      column(4, p(HTML("<b><br></b>")),
             actionButton(inputId = "triggerDeletion",label = "",
                          icon("trash-can"), 
                          style="color: #0E090D; background-color: #DABFDF; 
               border-color: #6A4D71"))
      
      ),
    
    fluidRow(
      
      column(12, p(HTML("<b>Wavelength range</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu1-5"),
                   
                   #Change colour of slider (this code should be re-written in a cleaner way. For now, it works)
                   tags$style(HTML('.js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {
                                                            
                                                            background: #00829c;
                                                            border-top: 1px solid #00829c ;
                                                            border-bottom: 1px solid #00829c ;}
          
                                      /* changes the colour of the number tags */
                                     .irs-from, .irs-to, .irs-single { background: #00829c }'
                   )),
                   sliderInput("wavelengthRange", NULL,min = 190, max = 250,value = c(190,250)),
                   tippy::tippy_this(elementId = "info_uu1-5",
                                     tooltip = "Select the wavelength range taking into account the signal noise.
                                     For example, check the CD curves of the buffer alone and the 
                                     'High Tension (HT) Voltage' reads.",placement = "right"))),
      
      
      
    ))