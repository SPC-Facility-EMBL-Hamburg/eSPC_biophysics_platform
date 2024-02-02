box(title = "Fluorescence dependence on temperature", width = 6, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(2, p(HTML("<b>b1</b>"), numericInput('b1', NULL, 80000,  min = 0,    max = 800000))),
      column(2, p(HTML("<b>b2</b>"), numericInput('b2', NULL, 80000,  min = 0,    max = 8000000))),
      column(2, p(HTML("<b>b3</b>"), numericInput('b3', NULL, 200000, min = 0,    max = 2000000))),
      column(2, p(HTML("<b>m1</b>"), numericInput('m1', NULL, -120,   min = -500, max = 500))),
      column(2, p(HTML("<b>m2</b>"), numericInput('m2', NULL, -120,   min = -500, max = 500))),
      column(2, p(HTML("<b>m3</b>"), numericInput('m3', NULL, -200,   min = -500, max = 500)))
    ))