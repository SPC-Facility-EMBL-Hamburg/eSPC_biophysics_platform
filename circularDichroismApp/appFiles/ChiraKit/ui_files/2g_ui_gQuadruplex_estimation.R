box(title = "2. Estimation", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "launchSamplesPCAGquad",label = "2a. Run PCA on samples",
                    icon("chart-line"),
                    style="color: #fff; background-color: #337ab7; 
               border-color: #2e6da4"))),
      
      column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "launchSVD_GQ",label = "2b. Run SVD structure estimation",
                    icon("tornado"),
                    style="color: #fff; background-color: #337ab7; 
               border-color: #2e6da4")))
    )
)