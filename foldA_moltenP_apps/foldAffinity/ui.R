source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

packages <- c("shinydashboard","shinycssloaders","rhandsontable","tableHTML","plotly","shinyalert")
invisible(lapply(packages, library, character.only = TRUE))

shinyUI(dashboardPage(title = "FoldAffinity",
  
  dashboardHeader(  title = logo_grey_light, titleWidth = 200), #logo_grey_light is described in logo.R
  dashboardSidebar( collapsed = F,width = 200,
                   
    sidebarMenu(
      
      menuItem("1. Load input",            icon = icon("file-circle-plus"),       tabName = "menu_input"),
      menuItem("2. Fit fluorescence",      icon = icon("chart-line"),             tabName = "menu_fit"),
      menuItem("3. Fit unfolded fraction", icon = icon("boxes-stacked"),          tabName = "menu_fit2"),
      menuItem("4. Export results",        icon = icon("file-export"),            tabName = "menu_export"),
      menuItem("Tm fitting",               icon = icon("temperature-half"),       tabName = "menu_tm_fit"),
      menuItem("Simulate data",            icon = icon("magnifying-glass-chart"), tabName = "menu_simulate"),
      menuItem("User guide",               icon = icon("user-astronaut"),         tabName = "menu_user_guide"),
      menuItem("Tutorial",                 icon = icon("book-open"),              tabName = "menu_tutorial"),
      menuItem("About",                    icon = icon("circle-info"),            tabName = "menu_about"))
    
      ),
  
  dashboardBody(theme_grey_light,
    tabItems(
      tabItem(tabName = "menu_input",
        fluidRow(

          # input$FLf, conc_units, which, sg_range, median_filter, n_replicates, dil_factor, initial_ligand, rev_order, fill_table
          source("ui_files/ui_load_input_box.R",local = TRUE)$value, 
          source("ui_files/ui_position_vs_concentration_box.R",local = TRUE)$value, # input$table1, table2, table3, table4
          source("ui_files/ui_signal_tab_box.R",local = TRUE)$value,
          source("ui_files/ui_plot_options_box.R",local = TRUE)$value,
          source("ui_files/ui_legend_tab_box.R",local = TRUE)$value
          
            )),
          
      tabItem(tabName = "menu_fit",
        fluidRow(
           
          source("ui_files/ui_model_selection_box.R",local = TRUE)$value, # input$model_selected, btn_cal
              
          conditionalPanel(
            condition = "input.model_selected != 'Global_CP'", # input$cp
            source("ui_files/ui_fitting_options_box.R",local = TRUE)$value),
              
          source("ui_files/ui_plot_options_fitting_box.R",    local = TRUE)$value, # input$select_fitting_plot
          source("ui_files/ui_fluorescence_fitting_tabbox.R", local = TRUE)$value,
          source("ui_files/ui_legend_fluorescence_fitting_tabbox.R",local = TRUE)$value
          
          )),
      
      tabItem(tabName = "menu_fit2",
              fluidRow(
                
                # input$model_sites, protein_conc, if_range, btn_cal_fit_its
                source("ui_files/ui_fitting_box.R",local = TRUE)$value, 
                
                tabBox(
                  title = "", width = 10,id = "tabset3",tabPanel("Isothermal Fitting", withSpinner(plotlyOutput("isf")))),
                source("ui_files/ui_plot_options_box_isf.R",local = TRUE)$value,
                source("ui_files/ui_legend_isothermal_fitting_tabbox.R",local = TRUE)$value
                
              )),
    
      tabItem(tabName = "menu_tm_fit",
              fluidRow(
                
                source("ui_files/ui_tm_model.R",local = TRUE)$value, # input$model_selected, btn_cal
                tabBox(title = "", width = 10,id = "tabset3",
                       tabPanel("Thermal Shift Fitting", withSpinner(plotlyOutput("tm_shift_fitting"))),
                       tabPanel("Fitted Parameters",tableOutput("params_table_tm_shift"))),
                source("ui_files/ui_plot_options_box_tm_shift.R",local = TRUE)$value,
                source("ui_files/ui_tm_fit_equations.R",local = TRUE)$value # input$model_selected, btn_ca
 
              )),
      
      tabItem(tabName = "menu_simulate",
              fluidRow(
                 
                source("ui_files/ui_kd_parameters_box.R",local=T)$value,                  # input$dHb dCPb dSb Tb dHu
                source("ui_files/ui_ku_parameters_box.R",local=T)$value,                  # input$dHu dCPu Tm  

                box(title = "Protein and Ligand Params.", width = 2, solidHeader = T, status = "primary", 
                    fluidRow(
                      column(6, p(HTML("<b>[Protein] in (uM)</b>"),numericInput('pconc_sim',NULL, 10,min = 1, max = 100))),
                      column(6, p(HTML("<b>Max. [Ligand] in (uM)</b>"),numericInput('maxLigandConcSimulation',NULL, 500,min = 1, max = 1e6)))
                      
                    )),
                                
              
                tabBox(title = "", width = 10,id = "tabset1",
                       tabPanel("Kd vs temperature",withSpinner(plotlyOutput("kd_vs_temp"))),
                       tabPanel("Ku vs temperature",withSpinner(plotlyOutput("ku_vs_temp"))),
                       tabPanel("L0 (Total ligand concentration) versus F, FL and U",withSpinner(plotOutput("fractions_vs_total_ligand"))),
                       tabPanel("Fluorescence vs Temperature",withSpinner(plotOutput("fluo_vs_ligand")))
                ),
                
                source("ui_files/ui_launch_simulation.R",local = TRUE)$value,
                
                source("ui_files/ui_fluo_dependence_on_temperature_box.R",local=T)$value, # input$m1 m2 m3 b1 b2 b3

                source("ui_files/ui_equations.R",local=T)$value
                
              )),
      
      tabItem(tabName = "menu_export",
              fluidRow(
                
                source("ui_files/ui_export_fitting_information.R",local=T)$value,
                source("ui_files/ui_export_plots_data.R",local=T)$value,
                source("ui_files/ui_export_session.R",local=T)$value
                
              )),
      
      tabItem(tabName = "menu_user_guide",includeHTML("docs/user_guide.html")),
      tabItem(tabName = "menu_tutorial",includeHTML("docs/tutorial.html")),
      #tabItem(tabName = "menu_tutorial",includeMarkdown("docs/tutorial.md")),
      tabItem(tabName = "menu_about",includeHTML("docs/about.html"))
      
      ))))
