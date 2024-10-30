source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

shinyUI(dashboardPage(

  title = paste0(appName),
  
  dashboardHeader(  title = logo_grey_light, titleWidth = 200), #logo_grey_light is described in logo.R
  dashboardSidebar( collapsed = F,width = 200,
    
    sidebarMenu(style = "white-space: normal;",
      
      menuItem(div(tags$img(src = "cd_spectrum.svg", width="20px"),
                   "1. Import data"), tabName = "menu_input"),
      
      menuItem(
        '2. Analysis',icon = icon("magnifying-glass-chart"),
        
        menuSubItem("2a. Thermal  unfolding",   icon = icon("temperature-high"),  tabName = "menu_thermal"),
        menuSubItem("2b. Chemical  unfolding",  icon = icon("flask"),             tabName = "menu_chemical"),
        
        menuSubItem(div(tags$img(src = "protein_helix_icon.svg", width="20px"),
      "2c. Protein secondary structure"), tabName = "menu_sec_structure"),
      
        menuSubItem("2d. Custom  analysis",              icon = icon("chart-line"),     tabName = "menu_custom"),
        menuSubItem("2e. Spectra comparison",            icon = icon("scale-balanced"), tabName = "menu_spectra_comparison"),
        menuSubItem("2f. Peptide helicity",              icon = icon("percent"),        tabName = "menu_peptide")

        #menuSubItem("2g. G-Quadruplex structure",        icon = icon("dna"),            tabName = "menu_gQuadruplex")
      
      ),
      
      menuItem("3. Export data",  icon = icon("file-export"),      tabName = "menu_export"),
      menuItem("User guide",      icon = icon("user-astronaut"),   tabName = "menu_user_guide"),
      menuItem("Tutorial",        icon = icon("book-open"),        tabName = "menu_tutorial"),
      menuItem("About",           icon = icon("circle-info"),      tabName = "menu_about"))
    
      ),
  
  dashboardBody(theme_grey_light,
                
                tags$style(HTML("
      iframe {
        width: 80%;
        height: 100vh;
        border: none;
      }
    ")), 
                
    tabItems(
      tabItem(
        tabName = "menu_input",
        
        fluidRow(
          
          source("ui_files/1ui_load_input_box.R",local = TRUE)$value,
          
          #Custom CSS to increase plot height
          tags$head(tags$style("
          #cdSpectra{height:500px !important;}
          #htSpectra{height:500px !important;}
          #cd_ht_spectra{height:500px !important;}
          #cdSpectraMiliDeg{height:500px !important;}"
                               )),
          
          # TabBox to plot the CD spectra and the associated voltage
          tabBox(title = "", width = 9,id = "tabBoxcdSpectra",
                 tabPanel("CD signal / working units",plotlyOutput("cdSpectra")),
                 tabPanel("HT signal",plotlyOutput("htSpectra")),
                 tabPanel("CD (working units) and HT",plotlyOutput("cd_ht_spectra")),
                 tabPanel("CD signal / milidegrees",plotlyOutput("cdSpectraMiliDeg"))
                 
                 )
          ),
          
        fluidRow(
        
          column(5,
                 source("ui_files/1ui_load_experiment_parameters.R",local = TRUE)$value,
                 tabBox(title = "", width = 12,id = "metadata")),
          
          column(7,
                 source("ui_files/1ui_processing.R",local = TRUE)$value,
                 source("ui_files/1ui_plotting_box.R",local = TRUE)$value))
        ),
      
      tabItem(
        tabName = "menu_thermal",
        fluidRow(
          
          column(4,
                 source("ui_files/2a_ui_load_thermal_denaturation_data.R", local = TRUE)$value),
          
          column(8,
                 source("ui_files/2a_ui_fitting_thermal.R",local = TRUE)$value,
                 
                 conditionalPanel(
                   "input.analysis_model_thermal == 'fixedWL'",
                   
                   #Custom CSS to increase plots height
                   tags$head(tags$style("
                   #meltingCurves{height:600px !important;}
                   #fittedMeltingCurves{height:600px !important;}
                   #residualsMeltingCurves{height:600px !important;}
                   #fractions_melting{height:600px !important;}"
                                        
                                        )),
                   
                   # TabBox to plot the thermal unfolding curves and the fitted parameters
                   tabBox(title = "", width = 12,id = "tabBoxcdSpectra1",
                          tabPanel("Melting curves",                      plotlyOutput("meltingCurves")),
                          tabPanel("Fitted melting curves",               plotlyOutput("fittedMeltingCurves")),
                          tabPanel("Residuals",                           plotOutput("residualsMeltingCurves")),
                          tabPanel("Fitted parameters (kcal Celsius mol)",tableOutput("fittedParams_melting")),
                          tabPanel("Relative errors (%)",                 tableOutput("fittedErrors_melting")),
                          tabPanel("Fitting bounds",                      tableOutput("fittingBounds_melting")),
                          tabPanel("Fractions",                           plotlyOutput("fractions_melting"))
                          
                          )
                   ),
                 
                 conditionalPanel(
                   "input.analysis_model_thermal != 'fixedWL'",
                   
                   #Custom CSS to increase plots height
                   tags$head(tags$style("
                   #meltingSpectra{height:600px !important;}
                   #basisSpectra{height:600px !important;}
                   #fittedSpectra{height:600px !important;}
                   #explainedVariance{height:600px !important;}
                   #svdCoefficients{height:600px !important;}
                   #fittedSVDCoefficients{height:600px !important;}
                   #residualsSVDCoefficients{height:600px !important;}
                   #fractions_melting_svd{height:600px !important;}"
                                        
                                        )),
                   
                   tabBox(title = "", width = 12,id = "tabBoxcdSpectra2",
                          
                          tabPanel("Melting spectra",        plotlyOutput("meltingSpectra")),
                          tabPanel("Basis spectra",          plotlyOutput("basisSpectra")),
                          tabPanel("Reconstructed spectra",  plotlyOutput("fittedSpectra")),
                          tabPanel("Explained variance",     plotlyOutput("explainedVariance")),
                          tabPanel("Coefficients",           plotlyOutput("svdCoefficients")),
                          tabPanel("Fitted coefficients",    plotlyOutput("fittedSVDCoefficients")),
                          tabPanel("Residuals",              plotOutput("residualsSVDCoefficients")),
                          tabPanel("Fractions",              plotlyOutput("fractions_melting_svd")),
                          
                          tabPanel("Fitted parameters (kcal Celsius mol)",tableOutput("fittedParams_meltingSVD")),
                          tabPanel("Relative errors (%)",                 tableOutput("fittedErrors_meltingSVD")),
                          tabPanel("Fitting bounds",                      tableOutput("fittingBoundsSVD_melting"))
                          
                      )))

            )),
      
      tabItem(tabName = "menu_chemical",
              fluidRow(
                
                column(4,
                       source("ui_files/2b_ui_load_chemical_denaturation_data.R", local = TRUE)$value),
                
                column(8,
                       source("ui_files/2b_ui_fitting_chemical.R"               , local = TRUE)$value,
                       
                       conditionalPanel(
                         "input.analysis_model_chemical == 'fixedWL'",
                         
                         #Custom CSS
                         tags$head(tags$style("
                         #chemicalCurves{height:600px !important;}
                         #fittedChemicalCurves{height:600px !important;}
                         #residualsChemicalCurves{height:600px !important;}
                         #fractions_chemical{height:600px !important;}"
                                              
                                              )),
                         
                         tabBox(title = "", width = 12,id = "tabBoxcdSpectraChem",
                                tabPanel("Chemical unfolding curves",plotlyOutput("chemicalCurves")),
                                tabPanel("Fitted parameters (kcal Celsius mol)",tableOutput("fittedParams_chemical")),
                                tabPanel("Relative errors (%)",tableOutput("fittedErrors_chemical")),
                                tabPanel("Fitting bounds",     tableOutput("fittingBounds_chemical")),
                                tabPanel("Fitted curves",plotlyOutput("fittedChemicalCurves")),
                                tabPanel("Residuals",    plotOutput("residualsChemicalCurves")),
                                tabPanel("Fractions",    plotlyOutput("fractions_chemical"))
                                
                                )),
                       
                       conditionalPanel(
                         "input.analysis_model_chemical != 'fixedWL'",
                         
                         #Custom CSS
                         tags$head(tags$style("
                         #chemUnfoldingSpectra{height:600px !important;}
                         #chemBasisSpectra{height:600px !important;}
                         #chemFittedSpectra{height:600px !important;}
                         #chemExplainedVariance{height:600px !important;}
                         #chemSVDCoefficients{height:600px !important;}
                         #chemFittedSVDCoefficients{height:600px !important;}
                         #chemResidualsSVDCoefficients{height:600px !important;}
                         #fractions_SVDchemical{height:600px !important;}"
                                              
                                              )),
                         
                         tabBox(title = "", width = 12,id = "tabBoxcdSpectraChemSVD",
                                tabPanel("Unfolding spectra",       plotlyOutput("chemUnfoldingSpectra")),
                                tabPanel("Basis spectra",         plotlyOutput("chemBasisSpectra")),
                                tabPanel("Reconstructed spectra", plotlyOutput("chemFittedSpectra")),
                                tabPanel("Explained variance",    plotlyOutput("chemExplainedVariance")),
                                tabPanel("Coefficients",          plotlyOutput("chemSVDCoefficients")),
                                tabPanel("Fitted coefficients",   plotlyOutput("chemFittedSVDCoefficients")),
                                tabPanel("Residuals",             plotOutput("chemResidualsSVDCoefficients")),
                                tabPanel("Fractions",             plotlyOutput("fractions_SVDchemical")),
                                
                                tabPanel("Fitted parameters (kcal Celsius mol)",tableOutput("fittedParams_chemicalSVD")),
                                tabPanel("Relative errors (%)",                 tableOutput("fittedErrors_chemicalSVD")),
                                tabPanel("Fitting bounds",                      tableOutput("fittingBounds_chemicalSVD"))
                                
                                ))
                )
                
              )),
      
      tabItem(
        tabName = "menu_sec_structure",
        fluidRow(
          
          column(5,
                 source("ui_files/2c_ui_fit_secondary_structure_box.R",local=T)$value,
                 conditionalPanel('input.secStr_analysis_mode == "SELCON"',
                 tabBox(title = "", width = 12,id = "secondary_structure_tabBox",
                        tabPanel("Secondary structure estimation",tableOutput("secondary_structure"))))
          ),
          
          #Custom CSS
          tags$head(tags$style("
                         #fitted_CD_spectra_Sec_Str{height:600px !important;}
                         "
          )),

          conditionalPanel('input.secStr_analysis_mode == "SELCON"',

          column(7,
                 tabBox(title = "", width = 12,id = "tabBoxcdSpectraSecStr",
                        tabPanel("Fitted spectra / mean unit molar extinction", plotlyOutput("fitted_CD_spectra_Sec_Str"))),
                        
                 tabBox(title = "", width = 12,id = "secondary_structure_calc_tabBox")
          )),

          conditionalPanel('input.secStr_analysis_mode == "SESCA_pred" || input.secStr_analysis_mode == "SESCA_est"',

                    column(7,

          conditionalPanel('input.secStr_analysis_mode == "SESCA_pred"',

          #Custom CSS to increase plot height
          tags$head(tags$style("
          #sesca_plot{height:600px !important;}"
          )),

          tabBox(title = "", width = 12,id = "sesca_results",
          tabPanel("Spectra / mean unit molar extinction" ,    plotlyOutput("sesca_plot" )),
          tabPanel("PDB Sec. Str.",                                tableOutput( "sesca_sec_str")),
          tabPanel("Model stats",                              tableOutput( "sesca_comparison_stats"))
          )),

          conditionalPanel('input.secStr_analysis_mode == "SESCA_est"',

          #Custom CSS to increase plot height
          tags$head(tags$style("
          #sesca_plot_bayes{height:600px !important;}"
          )),

          tabBox(title = "", width = 12,id = "sesca_results_est",
          tabPanel("Estimated Sec. Str.",       tableOutput( "sesca_sec_str_est")),
          tabPanel("Posterior probability" ,    plotlyOutput("sesca_plot_bayes" ))
          )),

          source("ui_files/2c_ui_sesca_plot_box.R",local=T)$value

          ))
        
        )),
      
      tabItem(
        tabName = "menu_custom",
        fluidRow(
          
          column(4,
                 source("ui_files/2d_ui_load_custom_data.R", local = TRUE)$value),
          
          column(8,
                 source("ui_files/2d_ui_custom_model_builder.R", local = TRUE)$value,
                 
                 conditionalPanel(
                   "input.analysis_model_custom == 'fixedWL'",
                   
                   #Custom CSS
                   tags$head(tags$style("
                         #customCurves{height:600px !important;}
                         #fittedCustomCurves{height:600px !important;}
                         #residualsCustomCurves{height:600px !important;}"
                   )),
                   
                   tabBox(title = "", width = 12,id = "tabBoxcdSpectraCustom",
                          tabPanel("Custom analysis curves",plotlyOutput("customCurves")),
                          tabPanel("Initial parameter estimates",rHandsontableOutput("initialParamsValues")),
                          tabPanel("Fitted parameters",tableOutput("fittedParams_custom")),
                          tabPanel("Relative errors (%)",tableOutput("fittedErrors_custom")),
                          tabPanel("Fitted curves",plotlyOutput("fittedCustomCurves")),
                          tabPanel("Residuals",plotOutput("residualsCustomCurves"))
                          
                   )),
                 
                 conditionalPanel(
                   "input.analysis_model_custom != 'fixedWL'",
                   
                   #Custom CSS
                   tags$head(tags$style("
                         #customSpectra{height:600px !important;}
                         #customBasisSpectra{height:600px !important;}
                         #customFittedSpectra{height:600px !important;}
                         #customExplainedVariance{height:600px !important;}
                         #customSVDCoefficients{height:600px !important;}
                         #customFittedSVDCoefficients{height:600px !important;}
                         #customResidualsSVDCoefficients{height:600px !important;}"
                   )),
                   
                   tabBox(title = "", width = 12,id = "tabBoxcdSpectraCustomSVD",
                          tabPanel("Custom analysis spectra", plotlyOutput("customSpectra")),
                          tabPanel("Basis spectra",         plotlyOutput("customBasisSpectra")),
                          tabPanel("Reconstructed spectra", plotlyOutput("customFittedSpectra")),
                          tabPanel("Explained variance",    plotlyOutput("customExplainedVariance")),
                          tabPanel("Coefficients",          plotlyOutput("customSVDCoefficients")),
                          tabPanel("Initial parameter estimates",rHandsontableOutput("initialParamsValuesSVD")),
                          tabPanel("Fitted coefficients",   plotlyOutput("customFittedSVDCoefficients")),
                          tabPanel("Residuals",   plotOutput("customResidualsSVDCoefficients")),
                          tabPanel("Fitted parameters",tableOutput("fittedParams_customSVD")),
                          tabPanel("Relative errors (%)",                 tableOutput("fittedErrors_customSVD"))
                   ))             
                 
                 )
          
        )),
      
      tabItem(
        tabName = "menu_spectra_comparison",
        fluidRow(
          
          column(4,
                 source("ui_files/2e_ui_load_spectra_comparison_data.R", local = TRUE)$value),
          
          column(8,
                 source("ui_files/2e_ui_spectra_comparison_workflow.R", local = TRUE)$value,
          
                #Custom CSS to increase plot height
                tags$head(tags$style("
                #cdSpectraAvg{height:600px !important;}
                #cdSpectraDiff{height:600px !important;}
                #cdSpectraDist{height:700px !important;}
                #cdSpectraTree{height:700px !important;}
                #cdSpectraSim{height:800px !important;}"
                )),
                
                # TabBox to plot the CD spectra and the associated voltage
                tabBox(title = "", width = 12,id = "tabBoxCompareSpectra",
                       tabPanel("CD avg ± sd",         plotlyOutput("cdSpectraAvg")),
                       tabPanel("Difference spectra",  plotlyOutput("cdSpectraDiff")),
                       tabPanel("Distances",           plotlyOutput("cdSpectraDist")),
                       tabPanel("Dendogram",           plotlyOutput("cdSpectraTree")),
                       tabPanel("Spectral similarity plots",plotlyOutput("cdSpectraSim"))
                )
          )
          
        )),
      
      tabItem(
        tabName = "menu_peptide",
        fluidRow(
          
          column(6,
                 source("ui_files/2f_ui_load_peptide_data.R", local = TRUE)$value),
          
          column(6,
                 source("ui_files/2f_ui_fitting_peptide.R" , local = TRUE)$value,
          
          tabBox(title = "", width = 12,id = "helicityResults",
                 tabPanel("Helicity Table", tableOutput( "helicityTable")),
                 tabPanel("Helicity Plot" , plotlyOutput("helicityPlot" ))
                 ))

          
        )),

      tabItem(
        tabName = "menu_gQuadruplex",
        
        column(6,
               fluidRow(
                 
                 column(12,
                        source("ui_files/2g_ui_gQuadruplex_references.R",local=T)$value,
                        
                        #Custom CSS to increase plot height
                        tags$head(tags$style("
                        #cdSpectraGQ{height:470px !important;}
                        #pca_results_GQ{height:470px !important;}
                        #pca_clustering_GQ{height:470px !important;}
                        "
                        )),
                        
                        # TabBox to plot the CD spectra and the associated voltage
                        tabBox(title = "", width = 12,id = "tabBoxRefSpectraGQuadruplex1",
                               tabPanel("Spectra - Ref",     plotlyOutput("cdSpectraGQ")),
                               tabPanel("PCA - Ref",         plotOutput("pca_results_GQ")),
                               tabPanel("Clusters - Ref",    plotOutput("pca_clustering_GQ")),
                               tabPanel("Secondary params",  tableOutput("secondary_params_GQ")),
                               tabPanel("Tertiary params",   tableOutput("tertiary_params_GQ"))
                               
                        ),
                        
                        source("ui_files/2g_ui_gQuadruplex_references_plot_settings.R",local=T)$value
                        
                        )
               )),
        
        column(6,
               fluidRow(
                 
                 column(12,
                        source("ui_files/2g_ui_gQuadruplex_estimation.R",local=T)$value,
                        
                        #Custom CSS to increase plot height
                        tags$head(tags$style("
                        #cdSpectraGQ_samples{height:470px !important;}
                        #pca_results_GQ_samples{height:470px !important;}
                        #pca_results_GQ_combined{height:470px !important;}
                        #pca_clustering_GQ_samples{height:470px !important;}
                        #pca_clustering_GQ_combined{height:470px !important;}
                        "
                        )),
                        
                        # TabBox to plot the CD spectra and the associated voltage
                        tabBox(title = "", width = 12,id = "tabBoxRefSpectraGQuadruplex2",
                               tabPanel("Spectra - Samples",      plotlyOutput("cdSpectraGQ_samples")),
                               tabPanel("PCA - Samples",          plotOutput("pca_results_GQ_samples")),
                               tabPanel("Clusters - Samples",     plotOutput("pca_clustering_GQ_samples")),
                               tabPanel("PCA - Ref+Samples",      plotOutput("pca_results_GQ_combined")),
                               tabPanel("Clusters - Ref+Samples", plotOutput("pca_clustering_GQ_combined")),
                               tabPanel("Secondary str.",         tableOutput("fitted_secondary_str_GQ")),
                               tabPanel("Tertiary str.",          tableOutput("fitted_tertiary_str_GQ"))
                               
                               
                        ),
                        source("ui_files/2g_ui_gQuadruplex_samples_plot_settings.R",local=T)$value
                        )
               ))
                
        ),      
            
      tabItem(tabName = "menu_export",
              fluidRow(
                
                source("ui_files/ui_export_cd_spectra.R"      , local=T)$value,
                source("ui_files/ui_export_logbook.R"         , local=T)$value,
                
                conditionalPanel(
                  "output.thermalDatasetCreated",
                  source("ui_files/ui_export_cd_spectra_thermal_fit.R"  , local=T)$value
                  
                  ),
                
                conditionalPanel(
                  "output.chemicalDatasetCreated",
                  source("ui_files/ui_export_cd_spectra_chemical_fit.R"  , local=T)$value
                  
                ),
                
                conditionalPanel(
                  "output.secStrFittingWasDone || output.secStrCalcWasDone",
                  source("ui_files/ui_export_cd_spectra_sec_str_fit.R"  , local=T)$value
                  
                ),

                conditionalPanel(
                  "output.sesca_pred_was_run",
                  source("ui_files/ui_export_sesca_pred.R"  , local=T)$value

                ),


                conditionalPanel(
                  "output.compareDatasetCreated",
                  source("ui_files/ui_export_cd_spectra_comparison.R"  , local=T)$value
                  
                ),
                
                conditionalPanel(
                  "output.customDatasetCreated",
                  source("ui_files/ui_export_cd_spectra_custom_fit.R"  , local=T)$value
                  
                )
              )),
      
      tabItem(tabName = "menu_user_guide", includeHTML("www/docs/user_guide.html")),
      tabItem(tabName = "menu_tutorial"  , includeHTML("www/docs/tutorial.html"  )),
      tabItem(tabName = "menu_about"     , includeHTML("www/docs/about.html"     ))
            
      ))))
