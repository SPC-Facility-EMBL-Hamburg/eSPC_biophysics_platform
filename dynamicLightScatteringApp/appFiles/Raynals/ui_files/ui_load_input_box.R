box(title = "Input", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(9, p(HTML("<b>1. DLS files </b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-1"),
                  fileInput("dlsFiles", NULL,accept = c(".csv",".zip",".7z"),multiple = TRUE),
                  tippy::tippy_this(elementId = "info_uu1-1",
                                    tooltip = "Upload a CSV file (output format from Wyatt) 
                                    or a compressed file (output format from Nanotemper Panta) 
                                    containing multiple Excel files. For more details, refer to the User Guide.",
                                    placement = "right"))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(HTML("<b><br></b>")),
             withBusyIndicatorUI(
               shinyjs::hidden(actionButton("Go","2. Load data!",class = "btn-primary")))),
      
      conditionalPanel(condition = '!output.inputFileGiven',
      
      column(12, p(HTML(""),
                   actionButton("GoLoadExample","Load example data",class = "btn-primary")))
      
      )
      
    ),
    
    fluidRow(
      
      column(4, p(HTML("<b><br>Set experiment parameters</b>"),
                  checkboxInput("setExperimentParameters", NULL, value = TRUE))),
      
      conditionalPanel('input.setExperimentParameters',
                       
                       column(4, p(HTML("<b><br>Select experiment</b>"),
                                   selectInput("experiment2delete", NULL,
                                               c(	'None'        = 'None')))),
                       
                       column(3, p(HTML("<b><br><br></b>")),
                              actionButton(inputId = "triggerDeletion",label = "Remove"))
                       
      )
      
    )
    
    )