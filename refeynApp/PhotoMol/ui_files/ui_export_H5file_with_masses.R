box(title = "Export H5 dataset", width = 2, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(6, p(HTML(""),
                  span(shiny::icon("info-circle"), id = "info_uuExportH5"),
             downloadButton('downloadH5file', label = "Download")),
             tippy::tippy_this(elementId = "info_uuExportH5",
                               tooltip = "Download a h5 file with a dataset called masses_kDa.
                               Export this file to avoid doing the calibration again.",
                               placement = "right"))
      
      
    ))


