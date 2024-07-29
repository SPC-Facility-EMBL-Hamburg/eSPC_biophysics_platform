box(title = "1. Reference dataset", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(4, p(HTML("<b>Use default reference set</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_useDefaultRefSetGQuad"),
                  checkboxInput("useDefaultReferenceSetGQuad", NULL,TRUE),
                  tippy::tippy_this(
                    elementId = "info_uu_useDefaultRefSetGQuad",
                    tooltip = "Change to False if you would like to use your 
                    own reference set of CD spectra, definition of secondary and tertiary
                    structure elements. You'll need three csv files containing 
                    three matrices: C (reference spectra), F 
                    (reference secondary structure elements). 
                    ",placement = "right")))
      
    )
)