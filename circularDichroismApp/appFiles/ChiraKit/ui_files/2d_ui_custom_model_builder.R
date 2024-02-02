box(title = "2. Fitting", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3,p(
        HTML("<b><br></b>"),
        span(shiny::icon("info-circle"), id = "info_uu_createCustomDataset"),
        actionButton("btn_create_custom_dataset","2a. Create dataset",class = "btn-primary"),
        tippy::tippy_this(
          elementId = "info_uu_createCustomDataset",
          tooltip = "Use the information from the Experimental Parameters Table to create a suitable dataset
          for the custom analysis.",placement = "right")
      )),
      
      column(4,p(
        HTML("<b>2b. Analysis type</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_analysisModeCustom"),
        selectInput('analysis_model_custom',NULL,
                    choices = c('Fixed wavelength'     = 'fixedWL',
                                'Spectra decomposition (SVD)' = 'spectraDecompositionSVD',
                                'Spectra decomposition (PCA)' = 'spectraDecompositionPCA')),
        tippy::tippy_this(
          elementId = "info_uu_analysisModeCustom",
          tooltip = "Select 'Fixed wavelength' to fit the CD signal at a certain 
          wavelength versus temperature. Select 'Spectra decomposition' to apply
          singular value decomposition (SVD), and then fit 
          the change in the coefficients from the basis spectra.",
          placement = "right")
        
      )),
      
      conditionalPanel(
        "input.analysis_model_custom != 'fixedWL'",
        
        column(4,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuCustomFitting-Decompose"),
          actionButton("btn_decompose_spectra_custom","2c. Decompose spectra",
                       class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuCustomFitting-Decompose",
            tooltip = "Apply the SVD/PCA algorithm.",placement = "right")
        ))
        
      )
      
    ),
    
    fluidRow(
      
      conditionalPanel(
        "input.analysis_model_custom != 'fixedWL'",
        
        column(3,p(
          HTML("<b>2d. Variance threshold</b>"),
          span(shiny::icon("info-circle"), id = "info_uuVarianceThresholdCustom"),
          numericInput("explained_variance_threshold_custom",NULL,99,50,100),
          tippy::tippy_this(
            elementId = "info_uuVarianceThresholdCustom",
            tooltip = "Modify the explained variance threshold to 
            determine how many basic spectra to use. 
            For SVD, we use the raw data variance. 
            For PCA, we use the variance of the centered data.",placement = "right")
          
        )),
        
        conditionalPanel(
          "output.show_basis_change_option",
          
          column(2,p(
            
            HTML("<b><br></b>"),
            span(shiny::icon("info-circle"), id = "info_uuChangeBasisCustom"),
            actionButton("btn_change_basis_custom","Change basis",class = "btn-primary"),
            tippy::tippy_this(
              elementId = "info_uuChangeBasisCustom",
              tooltip   = "Optional step. Combine the first and second basis spectra to create
              a new set of basis spectra. Afterwards, the new first basis spectrum will
              be similar to the first acquired spectra.",placement = "right")
            
          ))
          
        ),
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuCustomFitting-flipSpectrum"),
          actionButton("btn_flip_spectrum_custom","Invert spectrum",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuCustomFitting-flipSpectrum",
            tooltip = "Optional step. Invert a selected basis spectrum.",placement = "right")
          
        )),
        
        column(3,p(
          HTML("<b>2e. Coefficients of interest</b>"),
          span(shiny::icon("info-circle"), id = "info_uuThermalFitting-selectK"),
          numericInput("selectedK_custom",NULL,0,0,3),
          tippy::tippy_this(
            elementId = "info_uuCustomFitting-selectK",
            tooltip = "Select the number of coefficients to be
           fitted.",placement = "right")
          
        ))
        
      ),
      
      conditionalPanel(
        "input.analysis_model_custom == 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu_filteringAlgorithmCustom"),
          actionButton("btn_find_wl_custom","2c. Find wavelength(s)",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu_filteringAlgorithmCustom",
            tooltip = "Optional step. Apply a two-step filtering algorithm
            to select the wavelengths with the highest amplitude from the ones with 
            an adequate signal-to-noise ratio. 
            Only the first experimental parameter will be taken into account
            for this procedure.",
            placement = "right")
          
        )),
        
        column(3,p(
          HTML("<b>Selected wavelength(s)</b>"),
          span(shiny::icon("info-circle"), id = "info_uuSelectWLCustom"),
          textInput("selected_wavelength_custom", label=NULL,value="220"),
          tippy::tippy_this(
            elementId = "info_uuSelectWLCustom",
            tooltip = "Select the wavelength(s) for analyzing the CD signal against 
        the experimental parameters. The chosen curves will be fitted with shared fitting 
        parameters. Please use spaces and/or dashes as 
        separators. For example, valid input formats include: '220 224', 
        '220-224', and '220-224 228'.",
            placement = "right")
        ))
        
      )),
    
    fluidRow(
      
      column(12,p(
        HTML("<b>Fitting function</b>"),
        span(shiny::icon("info-circle"), id = "info_uuCustomFittingFunctionText"),
        textInput("custom_model",NULL,
        '( e^(-GlobalDhPos * (1 - (T+273) / (GlobalTmPos+273))  / (0.001987 * (T+273))) / (1 +  e^(-GlobalDhPos * (1 - (T+273) / (GlobalTmPos+273))  / (0.001987 * (T+273))))) * (kUNeg * (T+273) + bU) + (1 / (1 +  e^(-GlobalDhPos * (1 - (T+273) / (GlobalTmPos+273))  / (0.001987 * (T+273))))) * (0 * (T+273) + bN)'),
        tippy::tippy_this(
          elementId = "info_uuCustomFittingFunctionText",
          tooltip = "Provide a custom function to be fitted. 
          All curves will be fitted simultaneously.
          Parameters with the pattern 'global' in it, will have the same value in all curves.
          Parameters with the pattern 'pos' in it, will be constrained to be positive.
          Parameters with the pattern 'neg' in it, will be constrained to be negative
          ",placement = "right")
        
      )),
      
      column(3,p(
        HTML("<b><br></b>"),
        span(shiny::icon("info-circle"), id = "info_uuCustomGetInitialParams"),
        actionButton("btn_get_initial_params","2f. Get initial estimates",class = "btn-primary"),
        tippy::tippy_this(
          elementId = "info_uuCustomGetInitialParams",
          tooltip = "Use a log spaced grid search to obtain initial estimates of the function parameters.
          The initial values will also determine the fitting boundaries using a factor of 500. 
          If desired, use the 'Initial parameter estimates' Table to change either the initial values or the 
          fitting boundaries.",
          placement = "right")
        
      )),
      
      conditionalPanel(
        "input.analysis_model_custom != 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuCustomFitting-fit_SVD_coeff"),
          actionButton("btn_fit_custom_data_svd","2e. Fit data",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuCustomFitting-fit_SVD_coeff",
            tooltip = "Fit the CD custom curves. The fitted parameters will be constrained to lie
            inside the intervales delimited by the 'lower_limit' and 'high_limit' columns
            ('Initial parameter estimates Tab').",placement = "right")
          
        ))
        
      ),
      
      conditionalPanel(
        "input.analysis_model_custom == 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu-fitCustomCurves"),
          actionButton("btn_fit_custom_data","2e. Fit data",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu-fitCustomCurves",
            tooltip = "Fit the CD custom curves. The fitted parameters will be constrained to lie
            inside the intervales delimited by the 'lower_limit' and 'high_limit' columns
            ('Initial parameter estimates Tab').",placement = "right")
        ))
        
      ),
      
      column(3,p(
        HTML("<b>Show plot export options</b>"),
        checkboxInput("showPlotExportOptionsCustom",NULL,FALSE)
        
      )),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1,p(
        HTML("<b><br></b>"),
        withBusyIndicatorUI(
          shinyjs::hidden(actionButton("fitCustomHidden","",class = "btn-primary")))
      ))
    ),
    
    conditionalPanel(
      'input.showPlotExportOptionsCustom',
      fluidRow(
        column(3, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_custom", NULL,
                                c("PNG"    = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
        
        column(3, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotWidthCustom"),
                    numericInput('plot_width_custom',NULL, 12,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotWidthCustom",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(3, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotHeightCustom"),
                    numericInput('plot_height_custom',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotHeightCustom",
                                      tooltip = "Units are pixels * 50",
                                      placement = "right"))),                     
        
        column(3, p(HTML("<b>Text size</b>"),
                    numericInput('plot_axis_size_custom',NULL, 16,min = 4, max = 40)))),
      
      fluidRow(
        
        column(3, p(HTML("<b>Log scale for the x-axis</b>"),
                    checkboxInput('use_log_axis_custom',NULL,F)))
      )
    )
    
)

