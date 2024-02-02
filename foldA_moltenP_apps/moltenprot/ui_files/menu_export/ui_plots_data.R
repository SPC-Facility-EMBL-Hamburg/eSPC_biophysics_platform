box(title = "Plots data", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(6, p(HTML("</b>Data structure export format</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu4-2"),
                  selectInput("data_export_format", NULL,
                              c("Row-wise"      = "row_wise",
                                "Column-wise"    = "col_wise")),
                  tippy::tippy_this(
                    elementId = "info_uu4-2",
                    tooltip = "Select row-wise to export a csv file where the first column stores the temperature data and the
                    following columns (conditions) store the signal data (one column per condition). 
                    Select column-wise to export a csv file with 4 columns:
                    temperature, condition name used in moltenprot, signal value and original condition name",
                    placement = "right"))),
      
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_signal_plot'     , 'Signal'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_derivative_plot' , 'First Derivative'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_2derivative_plot', 'Second Derivative'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_max_derivative_plot', 'Maximum of Derivative')))#,
      #column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_fluo_fits_plot', 'Fittings')))
    ))