box(title = "Parameters", width = 6, solidHeader = T, status = "primary", 
    fluidRow(
  
        column(3, p(HTML("<b>RF1</b>"),numericInput('RF1_sim',NULL, 10,min = 0, max = 5000))),
        column(3, p(HTML("<b>RF2</b>"),numericInput('RF2_sim',NULL, 20,min = 0, max = 5000))),
      
       conditionalPanel(
         condition = "input.model_selected_sim != '2_Sites_2_Kd_Shared_Signal' && 
                      input.model_selected_sim != '2_Sites_2_Kd_Different_Signal'",
         
          column(3, p(HTML("<b>Kd (uM)</b>"),numericInput('kd_sim',NULL, 1,min = 0, max = 5000)))),
         
       conditionalPanel(
         condition = "input.model_selected_sim == '2_Sites_1_Kd_Shared_Signal_coop' || 
                      input.model_selected_sim == '2_Sites_1_Kd_Different_Signal_coop'",
         
         column(3, p(HTML("<b>Cooperativity factor</b>"),numericInput('c_factor_sim',NULL, 1,min = 0, max = 5000)))),
       
       conditionalPanel(
         condition = "input.model_selected_sim == '2_Sites_2_Kd_Shared_Signal' || 
                      input.model_selected_sim == '2_Sites_2_Kd_Different_Signal'",
         
          column(3, p(HTML("<b>Kd1 (uM)</b>"),numericInput('kd1_sim',NULL, 1,min = 0, max = 5000))),
          column(3, p(HTML("<b>Kd2 (uM)</b>"),numericInput('kd2_sim',NULL, 100,min = 0, max = 5000)))
         
         )
    )
)




