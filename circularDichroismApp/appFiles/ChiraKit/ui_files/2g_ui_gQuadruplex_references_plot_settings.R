box(title = "Reference dataset - Plot settings", width = 12, solidHeader = T, status = "primary", 
    
    fluidRow(
      
      column(3, p(HTML("<b>PC - Horizontal axis</b>"),
                  numericInput("PC_ha", NULL,1,min = 1,max = 8))),
      
      column(3, p(HTML("<b>PC - Vertical axis</b>"),
                  numericInput("PC_va", NULL,2,min = 2,max = 8))),            
      
      column(3, p(HTML("<b>Number of clusters</b>"),
                  selectInput("GQnClust", NULL,
                              paste0(3:10)))),
                              #c('Auto',paste0(3:10))))),
    
    column(3, p(HTML("<b>Show graph options</b>"),
                selectInput("showAdvanced_GQ", NULL,
                            choices = c('None','PCA','Cluster'))))),
    
    conditionalPanel(
      "input.showAdvanced_GQ == 'PCA'",
    
    fluidRow(
      
        column(3, p(HTML("<b>Graph title</b>"),
                    textInput("pca_graph_title", NULL,""))),
        
        column(3, p(HTML("<b>Graph label size</b>"),
                    numericInput("pca_graph_label", NULL,6,min = 0,max = 15))),
        
        column(3, p(HTML("<b>Axis label size</b>"),
                    numericInput("pca_axis_label", NULL,16,min = 0,max = 30))), 
        
        column(3, p(HTML("<b>Point size</b>"),
                    numericInput("pca_point_size", NULL,3,min = 0,max = 10))),
        
        column(3, p(HTML("<b>Arrow width</b>"),
                    numericInput("pca_arrow_width", NULL,1,min = 0,max = 10))),
        
        column(3, p(HTML("<b>Ellipse confidence level</b>"),
                    numericInput("pca_ellipse_confidence", NULL,85,min = 0,max = 100))),
        
        column(3, p(HTML("<b>Show ellipse</b>"),
                    checkboxInput("show_ellipse_pcaGQ", NULL,TRUE))),
        
        column(3, p(HTML("<b>Show legend</b>"),
                    checkboxInput("show_legend_pcaGQ", NULL,FALSE)))),
    
    fluidRow(
      
      column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  downloadButton("download_reference_png", "Export to PNG" ))),
      
      column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  downloadButton("download_reference_pdf", "Export to PDF" )))
    )
      
    ),
    
    conditionalPanel(
      "input.showAdvanced_GQ == 'Cluster'",
      
      fluidRow(
        
        column(3, p(HTML("<b>Graph title</b>"),
                    textInput("cluster_graph_title", NULL,""))),
        
        column(3, p(HTML("<b>Text size</b>"),
                    numericInput("cluster_text_size", NULL,9,min = 0,max = 20))),
      
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_cluster_png", "Export to PNG" ))),
        
        column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                    downloadButton("download_cluster_pdf", "Export to PDF" )))
      )
    )
)
