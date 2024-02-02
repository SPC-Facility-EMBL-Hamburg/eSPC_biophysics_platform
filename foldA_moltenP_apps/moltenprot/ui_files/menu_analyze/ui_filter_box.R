box(title = "Filters", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      column(12, p(HTML(""),
                   checkboxInput("sd_factor_bool",   "For all fitted parameters: (Standard_deviation / value) * 100 < 50", FALSE),
                   checkboxInput("bs_factor_bool",   "Baseline separation factor > 0.5 (Only for Two State models)", FALSE),
                   checkboxInput("far_from_bounds",  "Estimated parameters far from fitting bounds", FALSE)))
      
    ))