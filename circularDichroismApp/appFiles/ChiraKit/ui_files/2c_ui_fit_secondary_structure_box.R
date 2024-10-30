box(title = "1. Settings", width = 12, solidHeader = T, status = "primary",

    fluidRow(

      column(6, p(HTML("<b>Analysis mode</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_secStr_analysis_mode"),
                  selectInput('secStr_analysis_mode',NULL,
                              choices=c('SELCON (estimation)' = 'SELCON',
                                        'SESCA (prediction)' = 'SESCA_pred',
                                        'SESCA (estimation)' = 'SESCA_est')),
                  tippy::tippy_this(
                    elementId = "info_uu_secStr_analysis_mode",
                    tooltip =
                      "Use SELCON to estimate the protein secondary structure composition.
                      Use SESCA (prediction) to predict the spectra of provided PDB files.
                      Use SESCA (estimation) to determine the most likely secondary structure composition.",
                    placement = "right")))

    ),

    conditionalPanel('input.secStr_analysis_mode == "SESCA_est"',

    fluidRow(

            column(4, p(HTML("<b>a) Spectrum</b>"),
                  span(shiny::icon("info-circle"), id = "info_select_sesca_ref_est"),
                  selectInput("sescaReferenceEstimation", NULL,
                              choices = c('None')),
                  tippy::tippy_this(
                    elementId = "info_select_sesca_ref_est",
                    tooltip = "Select the reference spectrum to estimate the secondary structure composition."
                    ,placement = "right"))),

            column(4, p(HTML("<b>b) Basis Set </b>"),
                        span(shiny::icon("info-circle"), id = "info_select_sesca_basis_est"),
                        selectInput("sescaBasisSetEstimation", NULL,
                        # We leave only choices that have sense (6 or less labelled secondary structure elements & good performance)
                                    choices = c(
                                      "DS-dT","DSSP-T", "DSSP-1", "DS3-1", "DS6-1", "DS5-4",
                                      "HBSS-3",
                                      "DS-dTSC3", "DS5-4SC1",
                                      "DS6-1SC1", "DSSP-TSC1","DSSP-1SC3", "HBSS-3SC1"
                                    )),
                        tippy::tippy_this(
                          elementId = "info_select_sesca_basis_est",
                          tooltip = "Select the basis set for the SESCA Bayesian estimation."
                          ,placement = "right"))),

            conditionalPanel('input.sescaBasisSetEstimation == "DS-dTSC3" ||
                              input.sescaBasisSetEstimation == "DS5-4SC1" ||
                              input.sescaBasisSetEstimation == "DS6-1SC1" ||
                              input.sescaBasisSetEstimation == "DS5-6SC" ||
                              input.sescaBasisSetEstimation == "DSSP-TSC1" ||
                              input.sescaBasisSetEstimation == "DSSP-1SC3" ||
                              input.sescaBasisSetEstimation == "HBSS-3SC1"',

            column(4, p(HTML("<b>PDB file </b>"),
                        span(shiny::icon("info-circle"), id = "info_uu_pdbSESCA_est"),
                        fileInput("pdbFileSESCA_est", NULL,multiple = FALSE,accept = ".pdb"),
                        tippy::tippy_this(
                          elementId = "info_uu_pdbSESCA_est",
                          tooltip = "Import the PDB file (extension .pdb) so the side chain correction can be applied.
                          Only the sequence information is used!"
                          ,placement = "right")))

                           )),

        fluidRow(

              column(4, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "runSESCA_est",label = "Run estimation!",
                    icon("jedi"),
                    style="color: #fff; background-color: #337ab7;
               border-color: #2e6da4"))),

              column(4, p(HTML("<b># Iterations</b>"),
        span(shiny::icon("info-circle"), id = "info_select_sesca_iterations"),
          numericInput("sescaIterations", NULL,  250, 50, 1000, step = 5),
          tippy::tippy_this(
            elementId = "info_select_sesca_iterations",
            tooltip = "Select the number of iterations for the bayesian estimation of the secondary structure composition.
            We recommend at least 250. Increasing the number of iterations improves the accuracy,
            but also extends the computation time.",
            placement = "right")))

        )


    ),

    conditionalPanel('input.secStr_analysis_mode == "SESCA_pred"',

    fluidRow(

      column(4, p(HTML("<b>PDB file(s) </b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_pdbSESCA"),
                  fileInput("pdbFilesSESCA", NULL,multiple = TRUE,accept = ".pdb"),
                  tippy::tippy_this(
                    elementId = "info_uu_pdbSESCA",
                    tooltip = "Import the PDB file(s) (extension .pdb) to be used for predicting a CD spectrum."
                    ,placement = "right"))),

      column(4, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "runSESCA",label = "Run prediction!",
                    icon("jedi"),
                    style="color: #fff; background-color: #337ab7;
               border-color: #2e6da4"))),

      column(4, p(HTML("<b>Reference spectrum</b>"),
                  span(shiny::icon("info-circle"), id = "info_select_sesca_ref"),
                  selectInput("sescaReference", NULL,
                              choices = c('None')),
                  tippy::tippy_this(
                    elementId = "info_select_sesca_ref",
                    tooltip = "If you select a reference spectrum, the SESCA algorithm will
                    compare the structure-based predicted CD spectra with the reference spectrum."
                    ,placement = "right")))

      ),

     fluidRow(

      column(4, p(HTML("<b>Basis Set </b>"),
                  span(shiny::icon("info-circle"), id = "info_select_sesca_basis"),
                  selectInput("sescaBasisSet", NULL,
                              choices = c(
                                "DS-dT","DSSP-T", "DSSP-1", "DS3-1", "DS6-1", "DS5-4",
                                "HBSS-3",
                                "DS-dTSC3", "DS5-4SC1",
                                "DS6-1SC1", "DSSP-TSC1","DSSP-1SC3", "HBSS-3SC1"
                              )),
                  tippy::tippy_this(
                    elementId = "info_select_sesca_basis",
                    tooltip = "Select the basis set for the SESCA algorithm."
                    ,placement = "right"))),

      column(4, p(HTML("<b>Scaling factor</b>"),
        span(shiny::icon("info-circle"), id = "info_select_sesca_scaling_factor"),
          numericInput("sescaScalingFactor", NULL,  1, 0.5, 1.5, step = 0.01),
          tippy::tippy_this(
            elementId = "info_select_sesca_scaling_factor",
            tooltip = "Select the scaling factor for the reference spectrum.
            Can be applied only before running the prediction.",
            placement = "right"))),

      column(4, p(HTML("<b>Average ensemble</b>"),
        span(shiny::icon("info-circle"), id = "info_select_sesca_average_ensemble"),
          checkboxInput("sescaAverageEnsemble", NULL,FALSE),
          tippy::tippy_this(
            elementId = "info_select_sesca_average_ensemble",
            tooltip = "Average the ensemble of the predicted CD spectra to produce only one
            predicted spectrum. Can be applied before or after running the prediction.",
            placement = "right")))

    )

    ),

    conditionalPanel('input.secStr_analysis_mode == "SELCON"',

    fluidRow(

      column(3, p(HTML("<b>Lower WL</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_LowerWL"),
                  numericInput('lower_wl_secStr',NULL, 170,min = 170, max = 190),
                  tippy::tippy_this(elementId = "info_uu_LowerWL",
                                    tooltip = "Lower wavelength for the secondary 
                                    structure fitting. When using the default reference sets: 1) it 
                                    can't be lower than 175 and higher than 190. 
                                    2) For values between 180 (inclusive) and 190, the AU-SMP180 reference 
                                    set will be used. 
                                    3) For values between 175 and 180, the AU-SP175 reference set will
                                    be used.",placement = "right"))),
      
      column(3, p(HTML('<p style="margin-bottom:0px;"><br></p>'),
                  actionButton(
                    inputId = "runSecStrEstimation",label = "Run estimation!",
                    icon("jedi"),
                    style="color: #fff; background-color: #337ab7; 
               border-color: #2e6da4"))),
      
      column(1, 
             p(HTML("<b><br></b>"),
               withBusyIndicatorUI(
                 shinyjs::hidden(actionButton("hiddenBtnSecStr","",class = "btn-primary"))))),
      
      column(4, p(HTML("<b>PDB files </b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_pdbToSecStr"),
                  fileInput("pdbFiles", NULL,multiple = TRUE),
                  tippy::tippy_this(
                    elementId = "info_uu_pdbToSecStr",
                    tooltip = "Compute the secondary structure fractions from PDB/mmCIF files.
                    The following elements will be detected:
                    Regular alpha-helix, distorded alpha-helix, regular beta-sheet, distorted beta sheet, 
                    turns and 'other'. The total alpha-helix content corresponds to the DSSP fraction 'H' and 'G'.
                    The total beta-sheet content corresponds to the DSSP fraction 'E'.
                    Turns are derived from the DSSP fraction 'T'. 
                    The first and last two residues of an alpha-helix are considered 'distorted'. 
                    The first and last residue of a beta-sheet are considered 'distorted'.
                    ",placement = "right"))),
      
      # Little hack to use the withBusyIndicatorUI function (loading spinner)
      column(1, p(HTML("<b><br></b>"),
                  withBusyIndicatorUI(
                    shinyjs::hidden(actionButton("hiddenBtnPDBfile","",class = "btn-primary")))))
      
    ),
    
    fluidRow(
      
      column(4, p(HTML("<b>Use default reference set</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_useDefaultRefSet"),
                  checkboxInput("useDefaultReferenceSet", NULL,TRUE),
                  tippy::tippy_this(
                    elementId = "info_uu_useDefaultRefSet",
                    tooltip = "Change to False if you would like to use your 
                    own reference set of CD spectra and definition of secondary 
                    structure elements. You'll need two text files containing 
                    two matrices: C (reference spectra) and F 
                    (reference secondary structure elements). 
                    ",placement = "right"))),
      
      conditionalPanel(
        '!input.useDefaultReferenceSet',
                       
        column(4, p(HTML("<b>a) Matrix C</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_matrixC"),
                    fileInput("matrixC", NULL,accept = c('.txt','.csv','.dat')),
                    tippy::tippy_this(
                      elementId = "info_uu_matrixC",
                      tooltip = "m × n matrix containing the CD spectra in 
                      delta epsilon units (mean residue molar extinction).
                      m and n are respectively the number of wavelengths and proteins.
                    ",placement = "right"))),
        
        column(4, p(HTML("<b>b) Matrix F</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_matrixF"),
                    fileInput("matrixF", NULL,accept = c('.txt','.csv','.dat')),
                    tippy::tippy_this(
                      elementId = "info_uu_matrixF",
                      tooltip = "l × n matrix containing the fractions of the secondary 
                      structure elements. l and n are respectively the number of 
                      secondary structure elements and proteins.",
                      placement = "right"))),
        
        column(4, p(HTML("<b>c) Max wavelength</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_maxWL_refSet"),
                    numericInput("maxWL_refSet", NULL,240,0,1e3),
                    tippy::tippy_this(
                      elementId = "info_uu_maxWL_refSet",
                      tooltip = "Maximum wavelength of the reference set.
                      For example, 240 nm.",
                      placement = "right"))),
        
        column(4, p(HTML("<b>d) Wavelength step</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_stepWL_refSet"),
                    numericInput("stepWL_refSet",  NULL,1,0,5),
                    tippy::tippy_this(
                      elementId = "info_uu_stepWL_refSet",
                      tooltip = "Wavelength step of the reference set. 
                      For example, 1 nm (the reference spectra were measured
                      every 1 nm).",
                      placement = "right"))),
        
        column(4, p(HTML("<b>e) Secondary structure elements names</b>"),
                    
                    rHandsontableOutput('secondary_structure_elements_names')))
        
        )
      )
    )
)



