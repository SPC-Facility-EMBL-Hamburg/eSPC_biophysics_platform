tabBox(title = "", width = 10,id = "tabset1legend",
       tabPanel("Signal plot -  Figure caption", 
                p(HTML("Suggestion 1 - Conditions screen of protein Y.
                       Fluorescence based melting curves were plotted 
                       using the MoltenProt tool (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Thermal shift assay of protein X  ( 
                       © eSPC, spc.embl-hamburg.de)."))),
       
       tabPanel("First derivative plot -  Figure caption", 
                p(HTML("Suggestion 1 - Conditions screen of protein Y.  
                       The first derivative curves 
                       were plotted using the MoltenProt tool (spc.embl-hamburg.de)")),
                p(HTML("Suggestion 2 - Thermal shift assay of protein X. 
                The first derivative of the fluorescence signal is shown (
                © eSPC, spc.embl-hamburg.de).
                       "))),
       
       tabPanel("Tm plot -  Figure caption", 
                p(HTML("Suggestion 1 - Conditions screen of protein Y. 
                       Melting temperature estimation (using the first derivative) 
                       was done with the MoltenProt tool (spc.embl-hamburg.de).")),
                p(HTML("Suggestion 2 - Condition X/Y increases the stability of protein X. 
                Tm estimation was done using the first derivative (
                       © eSPC, spc.embl-hamburg.de).")))
       
       )