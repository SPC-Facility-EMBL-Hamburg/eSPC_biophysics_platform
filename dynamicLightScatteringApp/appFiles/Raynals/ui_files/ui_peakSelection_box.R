box(title = "Peak selection", width = 10, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(2, p(HTML("<b>Show</b>"),
                  checkboxInput("show_peaks_in_plot", "", TRUE))),
      
      column(2, p(HTML("<b>Estimation method</b>"),
                  span(shiny::icon("info-circle"), 
                       id = "info_uu_EstimationMethod"),
                  selectInput("hrEstimationMethod", NULL,
                              c("Harmonic mean"    =  "harmonicMean",
                                "Peak maximum"     =  "peakMax")),
                                #"Gaussian fit"     = "gaussFit")),
                  tippy::tippy_this(
                    elementId = "info_uu_EstimationMethod",
                    tooltip = "Select 'Harmonic mean' if you want to report
                    the intensity-weighted harmonic mean inside the selected Rh ranges.
                    Select 'Peak maximum' to report the maximum values.",
                    placement = "right"))),
      
      column(4, p(HTML("<b>Left bounds</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_peakSelectionLeftBound"),
                  textInput("peakLeftBounds", NULL, value = "0.1",),
                  tippy::tippy_this(elementId = "info_uu_peakSelectionLeftBound",
                                    tooltip = "Set the left limit for the peak selection area.
                                    Values should be separated by spaces (e.g., '0.1 10').
                                    The number of left limits should be the same as right limits.",
                                    placement = "right"))),
      
      column(4, p(HTML("<b>Right bounds</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_peakSelectionRightBound"),
                  textInput("peakRightBounds", NULL, value = "50"),
                  tippy::tippy_this(elementId = "info_uu_peakSelectionRightBound",
                                    tooltip = "Set the right limit for the peak selection area.
                                    Values should be separated by spaces (e.g., '50 500').
                                    The number of right limits should be the same as left limits.",
                                    placement = "right")))
      
    )
)
