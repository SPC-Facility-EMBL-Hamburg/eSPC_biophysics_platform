box(title = "2. Compare", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3,p(
        HTML("<b><br></b>"),
        span(shiny::icon("info-circle"), id = "info_uu_createCompareDataset"),
        actionButton("btn_create_compare_dataset","2a. Create dataset",class = "btn-primary"),
        tippy::tippy_this(
          elementId = "info_uu_createCompareDataset",
          tooltip = "Use the information from the 'Spectra labels' Table to create a suitable dataset
          for the comparison analysis.",placement = "right")
      )),
    
      column(3,p(
        HTML("<b>2b. Normalisation</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_normalise"),
        selectInput("normalise", NULL,
                    c("None","L2_norm")),
        tippy::tippy_this(
          elementId = "info_uu_normalise",
          tooltip = "Use the L2 normalisation to compare shapes.",placement = "right")
      )),
        
      column(3,p(
        HTML("<b>2c. Reference group</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_reference"),
        selectInput("comparison_reference", NULL,
                    c("None" = 'unknown')),
        tippy::tippy_this(
          elementId = "info_uu_reference",
          tooltip = "Select a reference to plot only the comparisons 
          against this group.",placement = "right")
      )),
      
      column(3,p(
        HTML("<b>Show plot export options</b>"),
        checkboxInput("showPlotExportOptionsCompare",NULL,FALSE)
        
      ))
      ),
  
    conditionalPanel(
      'input.showPlotExportOptionsCompare',
      fluidRow(
        column(3, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_compare", NULL,
                                c("PNG"    = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
        
        column(3, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotWidthCompare"),
                    numericInput('plot_width_compare',NULL, 12,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotWidthCompare",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(3, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotHeightCompare"),
                    numericInput('plot_height_compare',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotHeightCompare",
                                      tooltip = "Units are pixels * 50",
                                      placement = "right"))),                     
        
        column(3, p(HTML("<b>Text size</b>"),
                    numericInput('plot_axis_size_compare',NULL, 16,min = 4, max = 40))))
      
    )

)

