box(title = "Preprocessing", width = 8, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(2, p(HTML("<b><br></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-filter0"),
                  actionButton("selectAll", "(Un)select all"),
                  tippy::tippy_this(elementId = "info_uu-filter0",
                                    tooltip = "Select all the samples.",
                                    placement = "right"))),
      
      column(2, p(HTML("<b><br></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-filter3"),
                  actionButton("applyChanges", "Update"),
                  tippy::tippy_this(elementId = "info_uu-filter3",
                                    tooltip = "Apply the selected filters and 
                                    the changes done in the conditions Table 
                                    (i.e, you changed a sample name).",
                                    placement = "right"))),
      
      column(2, p(HTML("<b>Filter by <br> initial value</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-filter1"),
                  numericInput("filterByInitialValue",label = NULL, 1, min = 1, max = 10,step=0.02),
                  tippy::tippy_this(elementId = "info_uu-filter1",
                                    tooltip = "Unselect the autocorrelation curves where the initial value is
                                    lower than this value.",
                                    placement = "right"))),

      column(2, p(HTML("<b>BQ filter</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-filter2"),
                  numericInput("bumpRemovalTolerance",label = "Tolerance", 1.2, min = 1, max = 10,step=0.05),
                  tippy::tippy_this(
                    elementId = "info_uu-filter2",
                    tooltip = "Baseline quality filter: Unselect the autocorrelation curves that contain (after certain time)
                    at least one point whose value is higher than this threshold",
                    placement = "right"))),
      
      column(2, p(HTML("<b><br></b>"),
                  numericInput("bumpRemovalTimeLimit",label = "Time limit (μs)", 100, min = 0, max = 100,step=1),
                  )),
      
      column(2, p(HTML("<b>Group by</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-filter4"),
                  selectInput("splitFactorPreproccess", NULL,
                              c("None"          =  "None",
                                "Condition"     =  "Condition",
                                "Read"          =  "Read",
                                "Scan"          =  "Scan",
                                "Experiment"    =  "Experiment")),
                  tippy::tippy_this(elementId = "info_uu-filter4",
                                    tooltip = "Colour the plots according to different groups. 
                                    If no change is observed (but you expect it), it's because the Samples Tables was modified: 
                                    Press 'Update'",
                                    placement = "right")))
      
      
    ))