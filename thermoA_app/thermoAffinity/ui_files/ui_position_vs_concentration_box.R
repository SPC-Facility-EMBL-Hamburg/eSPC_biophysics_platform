box(title = "Position versus Concentration", width = 9, solidHeader = T, status = "primary",
    fluidRow(
      column(width = 3,rHandsontableOutput('table1',height="500px")),
      column(width = 3,rHandsontableOutput('table2',height="500px")),
      column(width = 3,rHandsontableOutput('table3',height="500px")),
      column(width = 3,rHandsontableOutput('table4',height="500px"))
      
    ))