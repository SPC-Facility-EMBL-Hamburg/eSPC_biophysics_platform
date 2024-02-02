box(title = "Explore Fhot selection", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
       column(3, p(HTML("<b></b>"),
        span(shiny::icon("info-circle"), id = "info_uu_adv_sett_Fhot"),
        checkboxInput("explore_Fhot", "", FALSE),
        tippy::tippy_this(elementId = "info_uu_adv_sett_Fhot",
        tooltip = "Estimate the Kd using different Fhot intervals (of 1 sec). 
                               The results displayed in the plot will be automatically filtered 
                               to leave only Kds that satisfy
                               (Standard_Error / Kd) < 0.5 and Lower_CI95 > 0.",
        placement = "right"))),
       
       conditionalPanel("input.explore_Fhot == 1",
       
       column(6,  
              withBusyIndicatorUI(actionButton("btn_cal_advance","Run Fitting",class = "btn-primary")))),

    )
)




