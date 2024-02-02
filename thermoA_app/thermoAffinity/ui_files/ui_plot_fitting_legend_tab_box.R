tabBox(title = "", width = 9,id = "tabset1legend",
       tabPanel("Fitting plot -  Figure caption", 
                p(HTML("Suggestion 1 - Protein X binds ligand Y with high/low 
                micromolar/nanomolar binding affinity. 
                Fitting of the (normalized) fluorescence traces / initial fluorescence.
                The confidence interval (marginal asymmetric confidence interval 
                at a 95 % confidence level) was estimated as suggested by 
                Vaida Paketuryte, et al., 2021. 
                The analysis and plot was done with ThermoAffinity (spc.embl-hamburg.de).")),
               
                p(HTML("Suggestion 2 - Protein X binds ligand Y with high/low 
                micromolar/nanomolar binding affinity. 
                Fitting of the MST traces / initial fluorescence.
                The confidence interval (marginal asymmetric confidence interval 
                at a 95 % confidence level) was estimated as suggested by 
                Vaida Paketuryte, et al., 2021 (
                © eSPC, spc.embl-hamburg.de). 
                ")))
)