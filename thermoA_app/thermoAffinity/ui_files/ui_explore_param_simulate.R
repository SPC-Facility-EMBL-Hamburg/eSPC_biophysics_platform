box(title = "Explore Parameters", width = 6, solidHeader = T, status = "primary", 
    fluidRow(
      
      column(3, p(HTML("<b>RF1</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_sim_explore_1"),
                  checkboxInput("explore_rf1", "", FALSE),
                  tippy::tippy_this(elementId = "info_uu_sim_explore_1",
                                    tooltip = "RF1 ± 10 %",
                                    placement = "right"))), 
      
      column(3, p(HTML("<b>RF2</b>"),
                  span(shiny::icon("info-circle"), id = "info_uu_sim_explore_2"),
                  checkboxInput("explore_rf2", "", FALSE),
                  tippy::tippy_this(elementId = "info_uu_sim_explore_2",
                                    tooltip = "RF2 ± 10 %",
                                    placement = "right"))), 
      conditionalPanel(
        condition = "input.model_selected_sim != '2_Sites_2_Kd_Shared_Signal' && 
                      input.model_selected_sim != '2_Sites_2_Kd_Different_Signal'",
        
        column(3, p(HTML("<b>Kd</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_sim_explore_3"),
                    checkboxInput("explore_kd", "", FALSE),
                    tippy::tippy_this(elementId = "info_uu_sim_explore_3",
                                      tooltip = "[ 0.1*kd - 10*kd ]",
                                      placement = "right")))),       
      conditionalPanel(
        condition = "input.model_selected_sim == '2_Sites_1_Kd_Shared_Signal_coop' || 
                      input.model_selected_sim == '2_Sites_1_Kd_Different_Signal_coop'",
        
        column(3, p(HTML("<b>C_factor</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_sim_explore_4"),
                    checkboxInput("explore_c_factor", "", FALSE),
                    tippy::tippy_this(elementId = "info_uu_sim_explore_4",
                                      tooltip = "[ 0.1*c_factor - 10*kd*c_factor ]",
                                      placement = "right")))),       
      conditionalPanel(
        condition = "input.model_selected_sim == '2_Sites_2_Kd_Shared_Signal' || 
                      input.model_selected_sim == '2_Sites_2_Kd_Different_Signal'",
        
        column(3, p(HTML("<b>Kd1</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_sim_explore_5"),
                    checkboxInput("explore_kd1", "", FALSE),
                    tippy::tippy_this(elementId = "info_uu_sim_explore_5",
                                      tooltip = "[ 0.1*kd - 10*kd ]",
                                      placement = "right"))),
        
        column(3, p(HTML("<b>Kd2</b>"),
                    span(shiny::icon("info-circle"), id = "info_uu_sim_explore_6"),
                    checkboxInput("explore_kd2", "", FALSE),
                    tippy::tippy_this(elementId = "info_uu_sim_explore_6",
                                      tooltip = "[ 0.1*kd - 10*kd ]",
                                      placement = "right"))) 
        
      )
    )
)




