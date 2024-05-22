box(title = "2. Fitting", width = 12, solidHeader = T, status = "primary", 
    fluidRow(

      column(3,p(
        HTML("<b><br></b>"),
        span(shiny::icon("info-circle"), id = "info_uu_createThermalDataset"),
        actionButton("btn_create_thermal_dataset","2a. Create dataset",class = "btn-primary"),
        tippy::tippy_this(
          elementId = "info_uu_createThermalDataset",
          tooltip = "Use the information from the Temperature Table to create a suitable dataset
          for the thermal unfolding analysis.",placement = "right")
        )),
      
      column(4,p(
        HTML("<b>2b. Analysis type</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_analysisModeThermal"),
        selectInput('analysis_model_thermal',NULL,
                    choices = c('Fixed wavelength'     = 'fixedWL',
                                'Spectra decomposition (SVD)' = 'spectraDecompositionSVD',
                                'Spectra decomposition (PCA)' = 'spectraDecompositionPCA')),
        tippy::tippy_this(
          elementId = "info_uu_analysisModeThermal",
          tooltip = "Select 'Fixed wavelength' to fit the CD signal at a certain 
          wavelength versus temperature. Select 'Spectra decomposition' to apply
          singular value decomposition (SVD), and then fit 
          the change in the coefficients from the basis spectra.",
          placement = "right")
        
        )),
      
      conditionalPanel(
        "input.analysis_model_thermal != 'fixedWL'",
        
        column(4,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuThermalFitting-Decompose"),
          actionButton("btn_decompose_spectra","2c. Decompose spectra",
                       class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuThermalFitting-Decompose",
            tooltip = "Apply the SVD/PCA algorithm.",placement = "right")
        ))
        
        )
      
      ),
    
    fluidRow(
      
      conditionalPanel(
        "input.analysis_model_thermal != 'fixedWL'",
        
        column(3,p(
          HTML("<b>2d. Variance threshold</b>"),
          span(shiny::icon("info-circle"), id = "info_uuVarianceThreshold"),
          numericInput("explained_variance_threshold",NULL,99,50,100),
          tippy::tippy_this(
            elementId = "info_uuVarianceThreshold",
            tooltip = "Modify the explained variance threshold to 
            determine how many basic spectra to use. 
            For SVD, we use the raw data variance. 
            For PCA, we use the variance of the centered data.",placement = "right")
          
        )),
       
        conditionalPanel(
          "output.show_basis_change_option",
                         
          column(2,p(
            
            HTML("<b><br></b>"),
            span(shiny::icon("info-circle"), id = "info_uuChangeBasis"),
            actionButton("btn_change_basis","Change basis",class = "btn-primary"),
            tippy::tippy_this(
              elementId = "info_uuChangeBasis",
              tooltip   = "Optional step. Combine the first and second basis spectra to create
              a new set of basis spectra. Afterwards, the new first basis spectrum will
              be similar to the first acquired spectra.",placement = "right")
            
          ))
        
        ),
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuThermalFitting-flipSpectrum"),
          actionButton("btn_flip_spectrum","Invert spectrum",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuThermalFitting-flipSpectrum",
            tooltip = "Optional step. Invert a selected basis spectrum.",placement = "right")
          
        )),
        
       column(3,p(
         HTML("<b>2e. Coefficients of interest</b>"),
         span(shiny::icon("info-circle"), id = "info_uuThermalFitting-selectK"),
         numericInput("selectedK",NULL,0,0,3),
         tippy::tippy_this(
           elementId = "info_uuThermalFitting-selectK",
           tooltip = "Select the number of coefficients to be
           fitted.",placement = "right")
         
         ))
       
      ),
      
      conditionalPanel(
        "input.analysis_model_thermal == 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu_filteringAlgorithm"),
          actionButton("btn_find_wl","2c. Find wavelength(s)",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu_filteringAlgorithm",
            tooltip = "Optional step. Apply a two-step filtering algorithm
            to select the wavelengths with the highest amplitude from the ones with 
            an adequate signal-to-noise ratio. ",
            placement = "right")
          
          )),
    
    column(3,p(
      HTML("<b>Selected wavelength(s)</b>"),
      span(shiny::icon("info-circle"), id = "info_uuSelectWLTherm"),
      textInput("selected_wavelength_thermal_unfolding", label=NULL,value="220"),
      tippy::tippy_this(
        elementId = "info_uuSelectWLTherm",
        tooltip = "Select the wavelength(s) for analyzing the CD signal against 
        temperature. The chosen curves will be fitted with shared thermodynamic 
        parameters, specifically the temperature of melting (Tm) and the 
        enthalpy of unfolding (dh). Please use spaces and/or dashes as 
        separators. For example, valid input formats include: '220 224', 
        '220-224', and '220-224 228'.",
        placement = "right")
      )),
    
    )),
    
    fluidRow(
    
      column(3,p(
        HTML("<b>Unfolding model</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_unfolding_model"),
        selectInput('thermal_unfolding_model',NULL,
                    choices = c('Reversible two-state'     = 'twoState',
                                'Reversible three-state'   = 'threeState')),
        tippy::tippy_this(
          elementId = "info_uu_unfolding_model",
          tooltip = "Apply an unfolding model where the protein unfolds in one step (reversible two-state), or 
          an unfolding model where where the proteins unfolds via an intermediate (reversible three-state).",
          placement = "right")
        
      )),
      
      conditionalPanel(
        "input.analysis_model_thermal != 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuThermalFitting-fit_SVD_coeff"),
          actionButton("btn_fit_melting_data_svd","2f. Fit data",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuThermalFitting-fit_SVD_coeff",
            tooltip = "Fit the CD melting curves.",placement = "right")
          
        )),
        
        ),
      
      conditionalPanel(
        "input.analysis_model_thermal == 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu-fitMeltingCurves"),
          actionButton("btn_fit_melting_data","2d. Fit data",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu-fitMeltingCurves",
            tooltip = "Fit the CD melting curves.",placement = "right")
        ))
        
        ),
      
    column(2,p(
      HTML("<b>Fit native slope</b>"),
      span(shiny::icon("info-circle"), id = "info_uu-fitSlopeNativeTherm"),
      checkboxInput("fitSlopeNative",NULL,FALSE),
      tippy::tippy_this(
        elementId = "info_uu-fitSlopeNativeTherm",
        tooltip = "Set to TRUE if there's a linear dependence
        between the CD signal (or PCA/SVD coefficients) and the temperature for the
        native state.",placement = "right")
      )),
    
    column(2,p(
      HTML("<b>Fit unfolded slope</b>"),
      span(shiny::icon("info-circle"), id = "info_uu-fitSlopeUnfoldedTherm"),
      checkboxInput("fitSlopeUnfolded",NULL,FALSE),
      tippy::tippy_this(
        elementId = "info_uu-fitSlopeUnfoldedTherm",
        tooltip = "Set to TRUE if there's a linear dependence
        between the CD signal (or PCA/SVD coefficients) and the temperature for the
        unfolded state.",placement = "right")
    
    )),
    
    # Little hack to use the withBusyIndicatorUI function (loading spinner)
    column(1,p(
      HTML("<b><br></b>"),
      withBusyIndicatorUI(
        shinyjs::hidden(actionButton("fitThermalHidden","",class = "btn-primary")))
    ))
  ),
  
  conditionalPanel(
    "input.thermal_unfolding_model == 'threeState'",
    
    fluidRow(
      
      column(3, p(HTML("<b>T1 Initial value (F <-> I)</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuT1_init"),
                  numericInput('T1_init',NULL, 0,min = 0, max = 10),
                  tippy::tippy_this(
                    elementId = "info_uuT1_init",
                    tooltip = "Initial value for the parameter T1. 
                      If zero, the parameter is constrained
                      within [minT + 4, maxT-7], where minT and maxT are 
                      respectively the minimum and maximum temperatures.
                      If non-zero, the parameter is constrained between [initial_guess - 15, initial_guess + 15],
                      Note: Ensure the fitted parameter is not too close to the boundaries 
                      after fitting.",
                    placement = "right"))),
      
      column(3, p(HTML("<b>T2 Initial value (I <-> U)</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuT2_init"),
                  numericInput('T2_init',NULL, 0,min = 0, max = 10),
                  tippy::tippy_this(
                    elementId = "info_uuT2_init",
                    tooltip = "Initial value for the parameter T2. 
                      If zero, the parameter is constrained
                      within [minT + 4, maxT], where minT and maxT are 
                      respectively the minimum and maximum temperatures.
                      If non-zero, the parameter is constrained between [initial_guess - 15, initial_guess + 15],
                      Note: Ensure the fitted parameter is not too close to the boundaries 
                      after fitting.",
                    placement = "right")))
      
    )),
  
  fluidRow(
    
    column(3,p(
      HTML("<b>Show plot export options</b>"),
      checkboxInput("showPlotExportOptionsMelting",NULL,FALSE)
      
    ))
    
    ),
  
  conditionalPanel(
    'input.showPlotExportOptionsMelting',
    fluidRow(
      column(3, p(HTML("<b>File type</b>"),
                  selectInput("plot_type_melt", NULL,
                              c("PNG"    = "png",
                                "SVG"    = "svg",
                                "JPEG"    = "jpeg")))),
      
      column(3, p(HTML("<b>Width</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuPlotWidthMelt"),
                  numericInput('plot_width_melt',NULL, 12,min = 1, max = 100),
                  tippy::tippy_this(elementId = "info_uuPlotWidthMelt",
                                    tooltip = "Units are pixels * 50",placement = "right"))),
      
      column(3, p(HTML("<b>Height</b>"),
                  span(shiny::icon("info-circle"), id = "info_uuPlotHeightMelt"),
                  numericInput('plot_height_melt',NULL, 11,min = 1, max = 100),
                  tippy::tippy_this(elementId = "info_uuPlotHeightMelt",
                                    tooltip = "Units are pixels * 50",
                                    placement = "right"))),                     
      
      column(3, p(HTML("<b>Text size</b>"),
                  numericInput('plot_axis_size_melt',NULL, 16,min = 4, max = 40)))),
    
    conditionalPanel(
      
      "input.analysis_model_thermal != 'fixedWL'",
      
    fluidRow(
      
      column(4, p(HTML("<b>Plot style (melting spectra)</b>"),
                  selectInput("plot_style_melt", NULL,
                              c("markers",
                                "lines"
                                ))))
    ))
  )
)

