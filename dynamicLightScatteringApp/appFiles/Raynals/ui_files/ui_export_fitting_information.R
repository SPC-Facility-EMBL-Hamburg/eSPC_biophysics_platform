box(title = "Fitted data", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      
      conditionalPanel('output.data_loaded',
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitted_hr_colWise', 'Hydrodynamic radius - column wise'))),
  
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitted_hr_rowWise', 'Hydrodynamic radius - row wise'))),

      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitted_dfs_colWise', 'Diffusion coefficients - column wise'))),
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitted_dfs_rowWise', 'Diffusion coefficients - row wise'))),
  
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitted_autocorrelation_colWise', 'Fitted autocorrelation - column wise'))),
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fitted_autocorrelation_rowWise', 'Fitted autocorrelation - row wise'))),
      
      conditionalPanel("input.regMethod != 'fixedAlpha'",
      
      column(8, p(style = "font-size: 120%",HTML(""),
                  downloadLink('download_fidelityPenaltyNorms', 'Fidelity and penalty norms (L-curve)'))))
      
      )
    ))