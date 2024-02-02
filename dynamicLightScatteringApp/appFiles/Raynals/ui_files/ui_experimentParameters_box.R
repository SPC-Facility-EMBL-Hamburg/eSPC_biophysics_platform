box(title = "Experiment parameters", width = 8, solidHeader = T, status = "primary", 
    fluidRow(
      
                       column(12,DT::dataTableOutput('dlsFilesInfo')),
                       tags$script(HTML("Shiny.addCustomMessageHandler('unbind-DT', function(id) {
          Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
        })"))
                       
      
    )
    
)