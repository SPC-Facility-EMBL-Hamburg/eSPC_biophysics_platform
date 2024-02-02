box(
  title = "Secondary structure data", width = 4, solidHeader = T, status = "primary", 
  fluidRow(
    
    conditionalPanel(
      'output.secStrFittingWasDone',
      
      column(12,p(
        style = "font-size: 120%",HTML(""),
        downloadLink('download_estimated_sec_str' ,'Estimated secondary structure (from spectra)'))),
      
      column(12,p(
        style = "font-size: 120%",HTML(""),
        downloadLink('download_sec_str_spectra_fittings' ,'Fitted spectra')))
    
    ),
    
    conditionalPanel(
      'output.secStrCalcWasDone', 
      
      column(12, p(
        style = "font-size: 120%",HTML(""),
        downloadLink('download_calculated_sec_str','Calculated secondary structure (from PDB)')))
      
    )
  ))