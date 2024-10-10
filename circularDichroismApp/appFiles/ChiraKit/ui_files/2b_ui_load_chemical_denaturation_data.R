box(title = "1. [Denaturant agent]", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      conditionalPanel(
        "input.inputMode != 'chemicalUnfolding'",
      
        column(8,p(
           HTML("<b>Oligomeric state</b>"),
           span(shiny::icon("info-circle"), id = "info_uu_oligomeric_state_chem"),
           selectInput('oligomeric_state_chem',NULL,
                       choices = c('Monomer'      = 'monomer',
                                   'Homodimer'    = 'dimer',
                                   'Homotrimer'   = 'trimer',
                                   'Homotetramer' = 'tetramer')),
           tippy::tippy_this(
             elementId = "ui_oligomeric_state_chem",
             tooltip = "Select if you want to apply unfolding models for monomers,
            dimers, trimers, or tetramers. In the last three cases, you need to input the 
            concentration of the n-mer equivalent, in micromolar units.",
             placement = "right")
         
       )),

      column(width = 12,rHandsontableOutput('chemical_denaturation_conc')))
      
    ))