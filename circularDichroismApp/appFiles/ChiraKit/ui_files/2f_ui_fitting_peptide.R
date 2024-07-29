box(title = "2. Fitting", width = 12, solidHeader = T, status = "primary", 

    fluidRow(
      
      column(6,p(
        HTML("<b><br></b>"),
        span(shiny::icon("info-circle"), id = "info_uu_btn_fit_helicity"),
        actionButton("btn_fit_helicity","2a. Fit peptide helicity",class = "btn-primary"),
        tippy::tippy_this(
          elementId = "info_uu_btn_fit_helicity",
          tooltip = "Compute the peptide helicity for peptides that undergo helix to coil transitions.
          Otherwise, using this feature for non-helix-coil type peptides can result in misleading or 
          inaccurate estimations.
          ",placement = "right")
      ))
      
      )
    
)