box(title = "1. Initialize and input the temperature data", width = 12, solidHeader = T, status = "primary", 
    
    fluidRow(
    
      column(6,p(
        HTML("<b><br></b>"),
        span(shiny::icon("info-circle"), id = "info_uu_create_peptide_dataset"),
        actionButton("btn_create_peptide_dataset","1a. Initialize the peptide dataset",class = "btn-primary"),
        tippy::tippy_this(
          elementId = "info_uu_create_peptide_dataset",
          tooltip = "Initiate the dataset for computing peptide helicity. You'll need to input the temperature data.
          ",placement = "right")
      ))
      
    ),
    
    fluidRow(
      
      column(width = 12,rHandsontableOutput('thermal_peptide_data'))
      
    ))