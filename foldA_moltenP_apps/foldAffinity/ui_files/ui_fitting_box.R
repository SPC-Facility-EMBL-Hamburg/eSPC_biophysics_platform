box(title = "Fitting", width = 4, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(6, p(HTML("<b>Model</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu3-1"),
                  selectInput("model_sites", NULL,
                              c("1 Site" = "One_Site","2 Sites 1 Kd" = "Two_Sites_One_Kd")),#,
                  #"2 Sites 2 Kd" = "Two_Sites_Two_Kd")),
                  tippy::tippy_this(elementId = "info_uu3-1",
                                    tooltip = "Select the model according the number of binding sites",placement = "right"))),
      
      column(6, p(HTML("<b>Protein Concentration (µM)</b>"),numericInput('protein_conc',NULL, 6,min = 0, max = 100))),
      
      column(10, p(HTML("<b>Isothermal Fit range (ºC)</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu3-3"),
                   #Change colour of slider (this code should be re-written in a cleaner way. For now. it works)                                             
                   tags$style(make_css(list('.irs-bar', c('border-top', 'border-bottom', 'background'), 
                                            rep('#00829c', 3)),
                                       list('.irs-bar-edge',c('background', 'border'),c('#00829c', '0px !important')),
                                       list('.irs-single','background','#00829c'))),
                   
                   sliderInput("if_range", NULL,min = 20, max = 90,value = c(40,80)),
                   tippy::tippy_this(elementId = "info_uu3-3",
                                     tooltip = "Select the range to perform the isothermal fitting. 
                                     This model works best when using temperatures slightly above the 
                                     melting temperature of the protein alone.",placement = "right"))),
      column(6, actionButton("btn_cal_fit_its", "Run Fitting"))
    ))