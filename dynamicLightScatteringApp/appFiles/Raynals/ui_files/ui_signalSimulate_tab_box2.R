tabBox(
  title = "", width = 6,id = "tabset4",
  tabPanel("Number    dist",withSpinner(plotOutput("hrDistributionNumberSim"))),
  tabPanel("Volume    dist",withSpinner(plotOutput("hrDistributionVolumeSim"))),
  tabPanel("Intensity dist",withSpinner(plotOutput("hrDistributionIntensitySim"))))
