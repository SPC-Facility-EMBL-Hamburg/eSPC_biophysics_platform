box(title = "Plots data", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      column(10, p(style = "font-size: 120%",HTML(""),
                   downloadLink('download_mass_histogram', 'Histogram of masses'))),
      column(10, p(style = "font-size: 120%",HTML(""),
                   downloadLink('download_fitted_gaussians', 'Fitted gaussians'))),
      column(10, p(style = "font-size: 120%",HTML(""),
                   downloadLink('download_mass_histogramNormalised', 'Normalised histogram of masses '))),
      column(10, p(style = "font-size: 120%",HTML(""),
                   downloadLink('download_fitted_gaussiansNormalised', 'Normalised fitted gaussians')))
      #column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_lig_signal_plot', 'Ligand Signal'))),
      #column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_init_fluo_plot', 'Initial Fluorescence versus [Ligand]'))),
      #column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_fnorm_plot', 'Fnorm versus [Ligand]'))),
      #column(10, p(style = "font-size: 120%",HTML(""),downloadLink('download_fitting_plot', 'Fitted Curve')))
    ))