tabBox(title = "", width = 10,id = "tabset1legend",
       tabPanel("Signal plot -  Figure caption", 
                p(HTML("Suggestion 1 - (Normalized) thermophoresis traces of protein X
                       incubated with different concentrations of
                       ligand Y. Plot created by ThermoAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - MST assay of protein X and ligand Y ( 
                       © eSPC, spc.embl-hamburg.de)."))),
       
       tabPanel("Ligand fluorescence -  Figure caption", 
                p(HTML("Suggestion 1 - Fluorescence at the maximum and minimum ligand concentrations.
                       Plot created by ThermoAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Fluorescence at the maximum and minimum ligand concentrations (
                © eSPC, spc.embl-hamburg.de).
                       "))),
       
       tabPanel("Initial fluorescence versus [ligand] -  Figure caption", 
                p(HTML("Suggestion 1 - Initial fluorescence versus ligand concentration. 
                Plot created by ThermoAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Initial fluorescence versus ligand concentration
                       (© eSPC, spc.embl-hamburg.de)."))),
       
       tabPanel("F_hot / F_cold versus [ligand] -  Figure caption", 
                p(HTML("Suggestion 1 - Normalized fluorescence versus ligand concentration.
                       Plot created by ThermoAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Normalized fluorescence versus ligand concentration (
                       © eSPC, spc.embl-hamburg.de).")))
       
)