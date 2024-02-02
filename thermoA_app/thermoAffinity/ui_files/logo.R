#Ref: https://github.com/nik01010/dashboardthemes/blob/master/R/dashboardthemes.R
shinyDashboardLogoDIY <- function(boldText, mainText, textSize = 20, badgeText, badgeTextColor,
                                  badgeTextSize = 2, badgeBackColor, badgeBorderRadius = 3) {
  
  htmltools::HTML(
    
    paste0(
      
      "<p style=\"font-size:", textSize, "px\">
      <b> ", boldText, " </b>",
      
      mainText ,"<span> &nbsp; </span>
      <span style=\"background-color: ", badgeBackColor, ";
      border-radius: ", badgeBorderRadius ,"px; \"> &nbsp;
      <font color=\"", badgeTextColor, "\" size=\"", badgeTextSize, "\">",
      
      badgeText ,"  </font> &nbsp; </span> </p>"
    )
  )
}

logo_grey_light <- shinyDashboardLogoDIY(
  boldText = "ThermoAffinity"
  ,mainText = ""
  ,textSize = 16
  ,badgeText = "" #"1.0"
  ,badgeTextColor = "white"
  ,badgeTextSize = 2
  ,badgeBackColor = "rgb(0,124,130)"#"rgb(0,176,155)" # default gray"rgb(150,150,150)"
  ,badgeBorderRadius = 3
)
