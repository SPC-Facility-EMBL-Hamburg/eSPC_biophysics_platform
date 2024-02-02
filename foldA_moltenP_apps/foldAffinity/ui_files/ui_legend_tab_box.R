tabBox(title = "", width = 10,id = "tabset1legend",
       tabPanel("Signal plot -  Figure caption", 
                p(HTML("Suggestion 1 - Protein Y stability induced by ligand X binding. 
                       Fluorescence-based melting curves of protein Y at different concentrations of
                       ligand X. Plot created by FoldAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Thermal shift assay of protein X and ligand Y ( 
                       © eSPC, spc.embl-hamburg.de)."))),
       

       tabPanel("First derivative plot -  Figure caption", 
                p(HTML("Suggestion 1 - Protein Y stability induced by ligand X binding. 
                       First derivative of the fluorescence-based melting curves of protein Y at 
                       different concentrations of
                       ligand X. Plot created by FoldAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Thermal shift assay of protein X and ligand Y. 
                The first derivative of the fluorescence signal is shown (
                © eSPC, spc.embl-hamburg.de).
                       "))),
       
       tabPanel("Tm versus [ligand] -  Figure caption", 
                p(HTML("Suggestion 1 - Protein Y stability induced by ligand X binding. 
                       Melting temperature estimation (using the first derivative) at 
                       different concentrations of
                       ligand X. Plot created by FoldAffinity (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Ligand Y increases the melting temperature of protein X. 
                Tm estimation was done using the first derivative (
                       © eSPC, spc.embl-hamburg.de).")))
       
       )