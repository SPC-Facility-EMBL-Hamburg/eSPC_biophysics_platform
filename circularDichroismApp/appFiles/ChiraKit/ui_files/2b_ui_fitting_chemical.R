box(title = "2. Fitting", width = 12, solidHeader = T, status = "primary", 
    fluidRow(
      
      conditionalPanel(
        "input.inputMode != 'chemicalUnfolding'",
      
      column(3,p(
        HTML(
          "<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu-createChemDataset"),
          actionButton("btn_create_chemical_dataset","2a. Create dataset",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu-createChemDataset",
            tooltip = "Use the information from the denaturant Table to create 
            a suitable dataset for the chemical unfolding analysis.",
            placement = "right")
          
          )),
      
      column(4,p(
        HTML("<b>2b. Analysis type</b>"),
        span(shiny::icon("info-circle"), id = "info_uu-ChemAnalysisType"),
        selectInput('analysis_model_chemical',NULL,
                    choices = c('Fixed wavelength'     = 'fixedWL',
                                'Spectra decomposition (SVD)' = 'spectraDecompositionSVD',
                                'Spectra decomposition (PCA)' = 'spectraDecompositionPCA')),
        tippy::tippy_this(
          elementId = "info_uu-ChemAnalysisType",
          tooltip = "Select 'Fixed wavelength' to fit the CD signal at a certain 
          wavelength versus the denaturant agent concentration. 
          Select 'Spectra decomposition' to apply singular value decomposition (SVD) 
          and then fit the change in the coefficients from the basis spectra.",
          placement = "right")
        )),
      
      conditionalPanel(
        "input.analysis_model_chemical != 'fixedWL'",
        
        column(4,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu-ChemDecompose"),
          actionButton("btn_decompose_spectra_chemical","2c. Decompose spectra",
                       class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu-ChemDecompose",
            tooltip = "Apply the SVD/PCA algorithm.",placement = "right")
        ))
        
      ))),
    
    fluidRow(
      
      conditionalPanel(
        "input.inputMode != 'chemicalUnfolding'",
      
      conditionalPanel(
        "input.analysis_model_chemical != 'fixedWL'",
        
        column(3,p(
          HTML("<b>2d. Variance threshold</b>"),
          span(shiny::icon("info-circle"), id = "info_uuVarianceThreshold"),
          numericInput("explained_variance_threshold_chemical",NULL,99,50,100),
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
            span(shiny::icon("info-circle"), id = "info_uuChangeBasisChemical"),
            actionButton("btn_change_basis_chemical","Change basis",class = "btn-primary"),
            tippy::tippy_this(
              elementId = "info_uuChangeBasisChemical",
              tooltip   = "Optional step. Combine the first and second basis spectra to create
              a new set of basis spectra. Afterwards, the new first basis spectrum will
              be similar to the first acquired spectra.",placement = "right")
            
          ))
          
        ),
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuChemicalFitting-flipSpectrum"),
          actionButton("btn_flip_spectrum_chemical","Invert spectrum",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuChemicalFitting-flipSpectrum",
            tooltip = "Optional step. Invert a selected basis spectrum.",placement = "right")
          
        )),
        
        column(3, p(
          HTML("<b>2e. Coefficients of interest</b>"),
          span(shiny::icon("info-circle"), id = "info_uu-selectKchem"),
          numericInput("selectedK_chemical",NULL,0,0,3),
          tippy::tippy_this(
            elementId = "info_uu-selectKchem",
            tooltip = "Select the number of coefficients to be fitted.",
            placement = "right")
          ))
        
      ),
      
      conditionalPanel(
        "input.analysis_model_chemical == 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu-findWLchem"),
          actionButton("btn_find_wl_chemical","2c. Find wavelength(s)",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu-findWLchem",
            tooltip = "Optional step. Apply a two-step filtering algorithm
            to select the wavelengths with the highest amplitude from the ones with 
            an adequate signal-to-noise ratio.",placement = "right")
          )),
        
        column(3,p(
          HTML("<b>Selected wavelength(s)</b>"),
          span(shiny::icon("info-circle"), id = "info_uu-selectedWLchem"),
          textInput("selected_wavelength_chemical_unfolding", label=NULL,value="220"),
          tippy::tippy_this(
            elementId = "info_uu-selectedWLchem",
            tooltip = "Select the wavelength(s) for analyzing the CD signal against 
            the denaturant agent concentration. The chosen curves will be fitted 
            with shared parameters, specifically the concentration of denaturant
            where 50 % of the protein is denaturated (D50) and the slope
            of the transition (m). Please use spaces and/or dashes as 
            separators. For example, valid input formats include: '220 224', 
            '220-224', and '220-224 228'.",placement = "right")
          ))
        
        )
      )),
    
    fluidRow(
      
      column(3,p(
        HTML("<b>Unfolding model</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_unfolding_model_chem"),
        selectInput('chemical_unfolding_model',NULL,
                    choices = c('Create the dataset first')),
        tippy::tippy_this(
          elementId = "info_uu_unfolding_model_chem",
          tooltip = "Apply an unfolding model where the protein unfolds in one step or two steps.",
          placement = "right")
        
      )),
      
      conditionalPanel(
        "input.analysis_model_chemical != 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uuChemicalFitting-fit_SVD_coeff"),
          actionButton("btn_fit_chemical_data_svd","2f. Fit data",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uuChemicalFitting-fit_SVD_coeff",
            tooltip = "Fit the CD melting curves.",placement = "right")
          
        ))
        
      ),
      
      conditionalPanel(
        "input.analysis_model_chemical == 'fixedWL'",
        
        column(3,p(
          HTML("<b><br></b>"),
          span(shiny::icon("info-circle"), id = "info_uu-fitChemicalCurves"),
          actionButton("btn_fit_chemical_data","2d. Fit data",class = "btn-primary"),
          tippy::tippy_this(
            elementId = "info_uu-fitChemicalCurves",
            tooltip = "Fit the CD melting curves.",placement = "right")
        ))
        
      ),
      
      column(2,p(
        HTML("<b>Fit native slope</b>"),
        span(shiny::icon("info-circle"), id = "info_uu-fitSlopeNativeChem"),
        checkboxInput("fitSlopeNative_chemical",NULL,FALSE),
        tippy::tippy_this(
          elementId = "info_uu-fitSlopeNativeChem",
          tooltip = "Set to TRUE if there's a linear dependence
          between the CD signal (or PCA/SVD coefficients) and the denaturant concentration for the
          native state.",placement = "right")
        )),
      
      column(2,p(
        HTML("<b>Fit unfolded slope</b>"),
        span(shiny::icon("info-circle"), id = "info_uu-fitSlopeUnfoldedChem"),
        checkboxInput("fitSlopeUnfolded_chemical",NULL,FALSE),
        tippy::tippy_this(
          elementId = "info_uu-fitSlopeUnfoldedChem",
          tooltip = "Set to TRUE if there's a linear dependence
          between the CD signal (or PCA/SVD coefficients) and the denaturant concentration for the
          unfolded state.",placement = "right")
        )),
    
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(
        HTML("<b><br></b>"),
        withBusyIndicatorUI(
          shinyjs::hidden(actionButton("fitChemicalHidden","",class = "btn-primary")))
      ))
      
    ),
    
    conditionalPanel(
      "input.chemical_unfolding_model == 'threeState'           ||
       input.chemical_unfolding_model == 'threeStateDimerMI'    ||
       input.chemical_unfolding_model == 'threeStateDimerDI'    ||
       input.chemical_unfolding_model == 'threeStateTrimerMI'   ||
       input.chemical_unfolding_model == 'threeStateTrimerDI'   ||
       input.chemical_unfolding_model == 'threeStateTetramerMI'
      ",
  
      fluidRow(
        
        column(3, p(HTML("<b>D50 Initial guess (DG1 = 0)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuD50v1_init"),
                    numericInput('D50v1_init',NULL, 0,min = 0, max = 10),
                    tippy::tippy_this(
                      elementId = "info_uuD50v1_init",
                      tooltip = "Initial value for the parameter D50 (first step).
                      If zero, ChiraKit will try to find a good initial value.
                      If non-zero, the parameter is constrained
                      within [initial_guess - 1.75, initial_guess + 1.75].
                      Note: Ensure the fitted parameter is not too close to the boundaries 
                      after fitting.",
                      placement = "right"))),
        
        column(3, p(HTML("<b>D50 Initial guess (DG2 = 0)</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuD50v2_init"),
                    numericInput('D50v2_init',NULL, 0,min = 0, max = 10),
                    tippy::tippy_this(
                      elementId = "info_uuD50v2_init",
                      tooltip = "Initial value for the parameter D50 (second step).
                      If zero, ChiraKit will try to find a good initial value.
                      If non-zero, the parameter is constrained
                      within [initial_guess - 1.75, initial_guess + 1.75].
                      Note: Ensure the fitted parameter is not too close to the boundaries 
                      after fitting.",
                      placement = "right")))
      
    )),
    
    fluidRow(
      
      column(3,p(
        HTML("<b>Show plot export options</b>"),
        checkboxInput("showPlotExportOptionsChemical",NULL,FALSE)
        
      ))
      
    ),
    
    conditionalPanel(
      'input.showPlotExportOptionsChemical',
      fluidRow(
        column(3, p(HTML("<b>File type</b>"),
                    selectInput("plot_type_chem", NULL,
                                c("PNG"    = "png",
                                  "SVG"    = "svg",
                                  "JPEG"    = "jpeg")))),
        
        column(3, p(HTML("<b>Width</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotWidthChem"),
                    numericInput('plot_width_chem',NULL, 12,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotWidthChem",
                                      tooltip = "Units are pixels * 50",placement = "right"))),
        
        column(3, p(HTML("<b>Height</b>"),
                    span(shiny::icon("info-circle"), id = "info_uuPlotHeightChem"),
                    numericInput('plot_height_chem',NULL, 11,min = 1, max = 100),
                    tippy::tippy_this(elementId = "info_uuPlotHeightChem",
                                      tooltip = "Units are pixels * 50",
                                      placement = "right"))),                     
        
        column(3, p(HTML("<b>Text size</b>"),
                    numericInput('plot_axis_size_chem',NULL, 16,min = 4, max = 40)))),
      
      conditionalPanel(
        
        "input.analysis_model_chemical != 'fixedWL'",
        
        fluidRow(
          
          column(4, p(HTML("<b>Plot style (melting spectra)</b>"),
                      selectInput("plot_style_chem", NULL,
                                  c("markers",
                                    "lines"
                                  ))))
        ))
    )
)

