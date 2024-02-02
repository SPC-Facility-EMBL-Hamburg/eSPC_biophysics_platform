box(
  title = "Comparison data", width = 4, solidHeader = T, status = "primary", 
  fluidRow(
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_average_spectra' ,'Average spectra by label'))),
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_standard_deviation' ,'Standard deviation by label'))),
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_difference_spectra' ,'Difference spectra'))),
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_difference_std' ,"Standard deviation of the 'difference' spectra"))),
    
    column(12, p(style = "font-size: 120%",HTML(""),
                 downloadLink('download_distance_data','Normalised euclidean distance matrix')))
    
  ))