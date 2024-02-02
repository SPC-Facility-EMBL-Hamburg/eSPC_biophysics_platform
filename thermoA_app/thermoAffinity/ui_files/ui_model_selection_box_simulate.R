box(title = "Settings", width = 5, solidHeader = T, status = "primary", 
    fluidRow(
      column(4, p(HTML("<b>Model</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_sim_2-1"),
                  selectInput("model_selected_sim", NULL,
                              c( "1 Site  "                                         = "one_site",
                                 "2 Sites, 1 Kd, Shared Signal"                     = "2_Sites_1_Kd_Shared_Signal",
                                 "2 Sites, 1 Kd, Different Signal"                  = "2_Sites_1_Kd_Different_Signal",
                                 "2 Sites, 1 Kd, Shared Signal, Cooperativity"      = "2_Sites_1_Kd_Shared_Signal_coop",
                                 "2 Sites, 1 Kd, Different Signal, Cooperativity"   = "2_Sites_1_Kd_Different_Signal_coop",
                                 "2 Sites 2 Kd, Shared Signal"                      = "2_Sites_2_Kd_Shared_Signal",
                                 "2 Sites 2 Kd, Different Signal"                   = "2_Sites_2_Kd_Different_Signal")),
                  
                  tippy::tippy_this(elementId = "info_uu_sim_2-1",
                                    tooltip = "Select the model for the simulation. 
                                    More information in the User Guide.",placement = "right"))),
      
      column(4, p(HTML("<b>[Protein] (uM)</b>"),
                  numericInput('protein_conc_sim',NULL, 5,min = 0, max = 1e6))),
      
      column(4, p(HTML("<b>Initial [Ligand]  (uM)</b>"),
                  numericInput('init_lig_sim',NULL, 2000,min = 0, max = 1e6)))
      ),
      
    fluidRow(
      
      column(4, p(HTML("<b># Dilutions </b>"),
                  numericInput('numb_dil_sim',NULL, 12,min = 0, max = 1e6))),
      
      column(4, p(HTML("<b>Ligand Dilution Factor </b>"),
                  numericInput('lig_dil_factor_sim',NULL, 2,min = 0, max = 1e6))),
      
      column(4, p(HTML("<b></b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_sim_2-2"),
                  checkboxInput("runSimulation","Launch simulation",FALSE),
                  tippy::tippy_this(elementId = "info_uu_sim_2-2",
                                    tooltip = "Check this button to compute the results of the simulation.",placement = "right")))
      
    ))