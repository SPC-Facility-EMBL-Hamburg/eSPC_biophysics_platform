box(title = "1. Reference dataset", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3, p(HTML("<b>1a. Use default references</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_useDefaultRefSetGQuad"),
                  checkboxInput("useDefaultReferenceSetGQuad", NULL,TRUE),
                  tippy::tippy_this(
                    elementId = "info_uu_useDefaultRefSetGQuad",
                    tooltip = "Change to False if you would like to use your 
                    own reference set of CD spectra, definition of secondary and tertiary
                    structure elements. You'll need three csv files containing 
                    three matrices: C (reference spectra), F 
                    (reference secondary structure elements). 
                    ",placement = "right"))),
      
      conditionalPanel(
        "!input.useDefaultReferenceSetGQuad",
        
        column(3, p(HTML("<b>Spectra</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_GQ_load_ref"),
                    fileInput("cd_ref_spectra_GQ", NULL,accept = '.csv'),
                    tippy::tippy_this(
                      elementId = "info_uu_GQ_load_ref",
                      tooltip = "Load a CSV file with the wavelength data in the first column and the 
                      spectral data in the subsequent columns (in molar extinction units). The file should have a header.
                      Proceed then to load secondary parameters and/or tertiary parameters.
                    ",placement = "right"))),
        
        column(3, p(HTML("<b>Secondary params</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_GQ_load_sec"),
                    fileInput("cd_ref_sec_params_GQ", NULL,accept = '.csv'),
                    tippy::tippy_this(
                      elementId = "info_uu_GQ_load_sec",
                      tooltip = "Load a CSV file with the spectra names in the first column and the 
                      associated secondary fractions in the subsequent columns. The file should have a header. 
                    ",placement = "right"))),
        
        column(3, p(HTML("<b>Tertiary params</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_GQ_load_ter"),
                    fileInput("cd_ref_ter_params_GQ", NULL,accept = '.csv'),
                    tippy::tippy_this(
                      elementId = "info_uu_GQ_load_ter",
                      tooltip = "Load a CSV file with the spectra names in the first column and the 
                      associated tertiary fractions in the subsequent columns. The file should have a header. 
                    ",placement = "right")))
        
        ), 
      
      column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "launchReferenceGquad",label = "1b. Load references",
                    icon("book"),
                    style="color: #fff; background-color: #337ab7; 
               border-color: #2e6da4")))
      
    )
)