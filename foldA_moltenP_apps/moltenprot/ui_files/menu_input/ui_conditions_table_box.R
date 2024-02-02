box(title = "Conditions table", width = 9, solidHeader = T, status = "primary",
    fluidRow(
      #4. 
      column(width = 3,rHandsontableOutput('table1',height="600px")),
      column(width = 3,rHandsontableOutput('table2',height="600px")),
      column(width = 3,rHandsontableOutput('table3',height="600px")),
      column(width = 3,rHandsontableOutput('table4',height="600px"))
    ))
