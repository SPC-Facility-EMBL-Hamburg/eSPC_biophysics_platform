box(title = "Plots data", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_signal_plot'     , 'Fluorescence signal'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_signalDerivative_plot'     , 'Fluorescence signal - derivative'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_initialSignal_plot' , 'Initial fluorescence signal'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_fluo_fit_plot'    , 'Fluorescence signal fitting'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_fraction_unfolded_plot' ,'Unfolded fraction'))),
      column(8, p(style = "font-size: 120%",HTML(""),downloadLink('download_fraction_unfolded_fit_plot' ,'Unfolded fraction fitting')))
    ))