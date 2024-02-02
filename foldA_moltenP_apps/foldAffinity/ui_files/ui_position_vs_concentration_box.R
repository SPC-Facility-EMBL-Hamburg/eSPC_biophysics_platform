box(title = "Position versus Concentration", width = 9, solidHeader = T, status = "primary",
    fluidRow(
      column(width = 3,rHandsontableOutput('table1')),
      column(width = 3,rHandsontableOutput('table2')),
      column(width = 3,rHandsontableOutput('table3')),
      column(width = 3,rHandsontableOutput('table4'))
      
    ))