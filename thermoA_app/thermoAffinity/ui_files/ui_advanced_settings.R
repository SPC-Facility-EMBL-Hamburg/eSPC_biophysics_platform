box(title = "Advanced Settings", width = 6, solidHeader = T, status = "primary", 
    fluidRow(
  
      conditionalPanel(condition = "input.show_advanced",
               
       column(3, p(HTML("<b>RF1</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_adv_sett_1"),
        checkboxInput("fix_RF1", "", FALSE),
        tippy::tippy_this(elementId = "info_uu_adv_sett_1",
        tooltip = "Fix the value of RF1",
        placement = "right"))),                       
                
       column(3, p(HTML("<b>RF2</b>"),
                   span(shiny::icon("info-circle"), id = "info_uu_adv_sett_2"),
                   checkboxInput("fix_RF2", "", FALSE),
                   tippy::tippy_this(elementId = "info_uu_adv_sett_2",
                                     tooltip = "Fix the value of RF2",
                                     placement = "right"))),
                      
        conditionalPanel(
          condition = "input.model_selected != '2_Sites_2_Kd_Shared_Signal' && 
                      input.model_selected  != '2_Sites_2_Kd_Different_Signal'",
                       
            column(3, p(HTML("<b>Kd</b>"),
              span(shiny::icon("info-circle"), id = "info_uu_adv_sett_3"),
              checkboxInput("fix_kd", "", FALSE),
              tippy::tippy_this(elementId = "info_uu_adv_sett_3",
              tooltip = "Fix the value of Kd",
              placement = "right"))),

        ),
       
       conditionalPanel(
         condition = "input.model_selected == '2_Sites_2_Kd_Shared_Signal' || 
                      input.model_selected == '2_Sites_2_Kd_Different_Signal'",
         
         column(3, p(HTML("<b>Kd1</b>"),
                     span(shiny::icon("info-circle"), id = "info_uu_adv_sett_4"),
                     checkboxInput("fix_kd1", "", FALSE),
                     tippy::tippy_this(elementId = "info_uu_adv_sett_4",
                                       tooltip = "Fix the value of Kd1",
                                       placement = "right"))),
         
         column(3, p(HTML("<b>Kd2</b>"),
                     span(shiny::icon("info-circle"), id = "info_uu_adv_sett_5"),
                     checkboxInput("fix_kd2", "", FALSE),
                     tippy::tippy_this(elementId = "info_uu_adv_sett_5",
                                       tooltip = "Fix the value of Kd2",
                                       placement = "right"))),
       ),
       
       
       conditionalPanel(condition = "input.fix_RF1",
                        column(3, p(HTML("<b>RF1</b>"),numericInput('RF1',NULL, 1,min = 0, max = 5000)))),
       
       conditionalPanel(condition = "!input.fix_RF1",
                        column(3, p(HTML("<b>Min RF1</b>"),numericInput('min_RF1',NULL, 1,min = 0, max = 5000))),
                        column(3, p(HTML("<b>Max RF1</b>"),numericInput('max_RF1',NULL, 1,min = 0, max = 5000)))),
       
       conditionalPanel(condition = "input.fix_RF2",
                        column(3, p(HTML("<b>RF2</b>"),numericInput('RF2',NULL, 1,min = 0, max = 5000)))),
       
       conditionalPanel(condition = "!input.fix_RF2",
                        column(3, p(HTML("<b>Min RF2</b>"),numericInput('min_RF2',NULL, 1,min = 0, max = 5000))),
                        column(3, p(HTML("<b>Max RF2</b>"),numericInput('max_RF2',NULL, 1,min = 0, max = 5000)))),
       
       conditionalPanel(
         condition = "input.model_selected != '2_Sites_2_Kd_Shared_Signal' && 
                      input.model_selected  != '2_Sites_2_Kd_Different_Signal'",
         
         conditionalPanel(condition = "input.fix_kd",
                          column(3, p(HTML("<b>Kd (uM)</b>"),numericInput('kd',NULL, 1,min = 0, max = 5000)))),
         
         conditionalPanel(condition = "!input.fix_kd",
                          column(3, p(HTML("<b>Min Kd (uM)</b>"),numericInput('min_kd',NULL, 1,min = 0, max = 5000))),
                          column(3, p(HTML("<b>Max Kd (uM)</b>"),numericInput('max_kd',NULL, 1,min = 0, max = 5000)))),
         
       ),
       
       conditionalPanel(
         condition = "input.model_selected == '2_Sites_2_Kd_Shared_Signal' || 
                      input.model_selected == '2_Sites_2_Kd_Different_Signal'",
         
         conditionalPanel(condition = "input.fix_kd1",
                          column(3, p(HTML("<b>Kd1 (uM)</b>"),numericInput('kd1',NULL, 1,min = 0, max = 5000)))),
         
         conditionalPanel(condition = "!input.fix_kd1",
                          column(3, p(HTML("<b>Min Kd1 (uM)</b>"),numericInput('min_kd1',NULL, 1,min = 0, max = 5000))),
                          column(3, p(HTML("<b>Max Kd1 (uM)</b>"),numericInput('max_kd1',NULL, 1,min = 0, max = 5000)))),
         
         conditionalPanel(condition = "input.fix_kd2",
                          column(3, p(HTML("<b>Kd2 (uM)</b>"),numericInput('kd2',NULL, 1,min = 0, max = 5000)))),
         
         conditionalPanel(condition = "!input.fix_kd2",
                          column(3, p(HTML("<b>Min Kd2 (uM)</b>"),numericInput('min_kd2',NULL, 1,min = 0, max = 5000))),
                          column(3, p(HTML("<b>Max Kd2 (uM)</b>"),numericInput('max_kd2',NULL, 1,min = 0, max = 5000))))
         
         )

      ),
      
      column(5, p(HTML("<b>Show</b>"),
        span(shiny::icon("info-circle"), id = "info_uu_adv_sett_6"),
        checkboxInput("show_advanced", "", FALSE),
        tippy::tippy_this(elementId = "info_uu_adv_sett_6",
        tooltip = "Fix parameters and change fitting bounds.",
        placement = "right")))
      
    )
)




