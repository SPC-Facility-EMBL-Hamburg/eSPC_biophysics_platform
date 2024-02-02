box(title = "Fit the isothermal shift", width = 3, solidHeader = T, status = "primary", 
    fluidRow(
      column(10, p(HTML("<b>Model</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuTm-0"),
                  selectInput("model_selected_tm_shift", NULL,
                              c( "One Site"    = "one_site_tm_shift",
                                 "Two Sites One Kd"    = "two_site_one_kd_tm_shift",
                                 #'Three Sites One Kd'  = 'three_site_one_kd_tm_shift',
                                 "Two Sites Two Kd"     = "two_site_two_kd_tm_shift")),
                  tippy::tippy_this(elementId = "info_uuTm-0",
                                    tooltip = "If the protein has two binding sites, try to fit first
                                    the two sites one kd model.",placement = "right"))),
      column(4,   
             withBusyIndicatorUI(actionButton("btn_cal_tm_fit","Run Fitting",class = "btn-primary")),
             ),
      column(4, p(HTML("<b></b>"),
                  span(shiny::icon("info-circle"), id = "info_uuTm-1"),
                  tippy::tippy_this(elementId = "info_uuTm-1",
                                    tooltip = "Use this  thermodynamically based 
                                    model if you can suppose L > P (where the signal changes)",placement = "right")))
    ))

