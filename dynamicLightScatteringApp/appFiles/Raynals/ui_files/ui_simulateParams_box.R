box(title = "Parameters", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(12,p(HTML(""),
                  rHandsontableOutput("simulation_params_table")))),
      
    fluidRow(
      
      column(2, p(HTML("<b> Number of simulations </b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-SimBox0"),
                  numericInput('SimNumber', NULL, value = 5,step = 5),
                  tippy::tippy_this(
                    elementId   =  "info_uu-SimBox0",
                    tooltip     =  "Select how many distributions of 
                    spheric particles you want to simulate. The
                    number & volume distributions depend on the population mean(s),
                    population sd(s) & sample size. The intensity 
                    distribution is based on the Mie Theory and depends also on the
                    wavelength, angle of detection &
                    refractive index. The autocorrelation curve depends on the Mie intensity
                    contributions and the diffusion coefficient (determined by the
                    temperature & viscosity).",
                                    placement   =  "right"))),
      
      column(1, p(HTML("<b><br></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-SimBox1"),
                  actionButton("launchSimulation", "Run"),
                  tippy::tippy_this(elementId = "info_uu-SimBox1",
                                    tooltip = "Launch the simulation. Press this button 
                                    everytime you change the parameters.",
                                    placement = "right"))),
      
      column(2, p(HTML("<b> Intercept </b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-SimBox3"),
                  numericInput('autocorrelationIntercept', NULL, value = 0.2,step = 0.1,
                               min=0, max=1),
                  tippy::tippy_this(
                    elementId   =  "info_uu-SimBox3",
                    tooltip     =  "Intercept of the autocorrelation curves.",
                    placement   =  "right"))),
      
      column(2, p(HTML("<b> Gaussian error </b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-SimBox4"),
                  numericInput('autocorrelationError', NULL, value = 0,step = 0.001),
                  tippy::tippy_this(
                    elementId   =  "info_uu-SimBox4",
                    tooltip     =  "Standard deviation of the gaussian error 
                    added to the autocorrelation curves.",
                    placement   =  "right"))),
      
      column(1, p(HTML("<b><br></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu-SimBox2"),
                  downloadButton("exportSimData", "Export"),
                  tippy::tippy_this(elementId = "info_uu-SimBox2",
                                    tooltip = "Export the simulated (meta)data.",
                                    placement = "right"))),
      
      )
      
)