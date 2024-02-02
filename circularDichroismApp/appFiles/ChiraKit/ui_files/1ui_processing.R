box(title = "3. Processing", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(12,DT::dataTableOutput('proccesingInfo')),
      tags$script(HTML("Shiny.addCustomMessageHandler('unbind-DT', function(id) {
          Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
        })")),
      
      column(2, p(HTML("<b><br></b>"),
                  actionButton(
                    inputId = "triggerProcessing",label = "Apply",
                    icon("forward"), 
                    style="color: #fff; background-color: #337ab7; 
               border-color: #2e6da4"))),
      
      column(4, p(HTML("<b><br>Allow linear interpolation</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-processing1"),
             checkboxInput(inputId = 'allowLinearInterpolation',label = NULL,value = FALSE),
             tippy::tippy_this(
               elementId = "info_uu-processing1",
               tooltip = "Check to allow manipulation of spectra where the signal was measured 
               at different wavelength points and/or ranges. 
               This option is not advised for reference/baseline subtraction of a sample.
               ",placement = "right"))),
      
      column(4, p(HTML("<b><br>Working units</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu1-4"),
                  selectInput("workingUnits", NULL,global_cd_units_choice),
                  tippy::tippy_this(
                    elementId = "info_uu1-4",
                    tooltip = "Select the final units of analysis. 
                    The working units will be used to plot the CD signal, 
                    to create the thermal / chemical unfolding datasets, 
                    and to export the CD data. 
                    To calculate the molar ellipticity/extinction, 
                    please input the molecular weight, concentration and path length
                    ('2. Parameters for molar ellipticity / extinction' Box).
                    To calculate the mean residue molar ellipticity/extinction, 
                    please add the number of amino acids.",placement = "right"))), 
      conditionalPanel(
        'output.showVoltageThreshold',
                       
        column(2, p(HTML("<b><br>HT max value</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu-maxVoltage"),
                    numericInput('maxHTvalue',NULL,value = 0),
                    tippy::tippy_this(
                      elementId = "info_uu-maxVoltage",
                      tooltip = "Adapt the wavelength range based 
                      on the maximum accepted value for the high tension (HT)
                      voltage curves. As a result, all the CD curves will have 
                      associated HT values below the selected threshold.",
                      placement = "right")))
      )
))

