box(
  title = "SESCA predicted spectra", width = 4, solidHeader = T, status = "primary",
  fluidRow(

    column(12, p(
      style = "font-size: 120%",HTML(""),
      downloadLink('download_sesca_pred_data','Predicted spectra')))

  ))