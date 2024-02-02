box(title = "1. Estimation", width = 12, solidHeader = T, status = "primary", 
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



