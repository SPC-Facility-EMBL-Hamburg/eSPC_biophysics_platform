box(title = "Plots data", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_raw_signal_plot', 'Raw Signal'))),
      column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_lig_signal_plot', 'Ligand Signal'))),
      column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_init_fluo_plot', 'Initial Fluorescence versus [Ligand]'))),
      column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_fnorm_plot', 'Fnorm versus [Ligand]'))),
      column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_fitting_plot', 'Fitted Curve')))
    ))