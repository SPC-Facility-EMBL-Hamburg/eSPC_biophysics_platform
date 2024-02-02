tabBox(title = "", width = 9,id = "tabset1legend",
       
       tabPanel("Counts vs Mass (kDa) plot -  Figure caption",
                
                p(HTML("Option 1 - MP mass distribution histogram and Gaussian fit of the
                       sample XX. The plot and analysis was done with the
                       PhotoMol tool (spc.embl-hamburg.de).")),
                p(HTML("Option 2 - Mass distribution fitting of sample Y ( 
                       © eSPC, spc.embl-hamburg.de)."))),
       
       tabPanel("All (unbinding events) -  Figure caption", 
                p(HTML("Option 1 - MP binding
                and unbinding events histogram.
                       Plot was generated using the PhotoMol tool (spc.embl-hamburg.de)")),
                p(HTML("Option 2 - MP binding
                and unbinding events histogram (
                © eSPC, spc.embl-hamburg.de).")))
       
       
)