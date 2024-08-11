source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

packages <- c("shinydashboard","shinycssloaders","rhandsontable","shinyalert","plotly")
invisible(lapply(packages, library, character.only = TRUE))

shinyUI(dashboardPage(title = "ThermoAffinity",
  
  dashboardHeader(  title = logo_grey_light, titleWidth = 200), #logo_grey_light is described in logo.R
  dashboardSidebar( collapsed = F,width = 200,
                   
    sidebarMenu(
      
      menuItem("1. Load input",            icon = icon("file-circle-plus"),       tabName = "menu_input"),
      menuItem("2. Fitting",               icon = icon("chart-line"),             tabName = "menu_fit"),
      menuItem("3. Export results",        icon = icon("file-export"),            tabName = "menu_export"),
      menuItem("Simulate data",            icon = icon("magnifying-glass-chart"), tabName = "menu_simulate"),
      menuItem("User guide",               icon = icon("user-astronaut"),         tabName = "menu_user_guide"),
      menuItem("Tutorial",                 icon = icon("book-open"),              tabName = "menu_tutorial"),
      menuItem("About",                    icon = icon("circle-info"),            tabName = "menu_about"))
    
      ),
  
  dashboardBody(theme_grey_light,
    tabItems(
      tabItem(tabName = "menu_input",
        fluidRow(
          
          # input$FLf, conc_units, which, cold_range,hot_range, median_filter, n_replicates, dil_factor, initial_ligand, rev_order, fill_table
          source("ui_files/ui_load_input_box.R",local = TRUE)$value, 
          source("ui_files/ui_position_vs_concentration_box.R",local = TRUE)$value, # input$table1, table2, table3
          source("ui_files/ui_signal_tab_box.R",local = TRUE)$value,
          
          conditionalPanel(condition = "output.ui_signal_tab_box",
                           source("ui_files/ui_plot_options_box.R",local = TRUE)$value,
                           source("ui_files/ui_plot_legend_tab_box.R",local = TRUE)$value)
                                            
            )),
      
      tabItem(tabName = "menu_fit",
        fluidRow(
          
          source("ui_files/ui_model_selection_box.R",         local = TRUE)$value, 
          source("ui_files/ui_advanced_settings.R",           local = TRUE)$value,
          
          conditionalPanel("input.signal_selected == 'f_norm'",
          
          source("ui_files/ui_advanced_settings_explore_fhot.R",local = TRUE)$value),
          
          conditionalPanel("input.explore_Fhot != 1 || output.is_csv",
            source("ui_files/ui_fluorescence_fitting_tabbox_simple.R", local = TRUE)$value),
          
          conditionalPanel("input.explore_Fhot == 1 && !output.is_csv",
            source("ui_files/ui_fluorescence_fitting_tabbox_explore.R", local = TRUE)$value)
          
          ),
          
          conditionalPanel(
            condition = "output.data_loaded",
            
            fluidRow(source("ui_files/ui_plot_options_box_fit.R",        local = TRUE)$value)
            
          ),
          
        fluidRow(
          
          box(title = "Legends", width = 3, solidHeader = T, status = "primary",
              fluidRow(
                column(width = 8,rHandsontableOutput('legendInfo')),
                column(width = 4,
                       
                  colourpicker::colourInput("colorForLegend", NULL, value = "#E41A1C"),
                  selectInput("mol2changeColor", "Set colour",c("X")))
                
              )),
          
          source("ui_files/ui_plot_fitting_legend_tab_box.R",local = TRUE)$value),
        
        fluidRow(
          column(width = 12,style='padding-left:0px; padding-right:1px; padding-top:5px; padding-bottom:15px')
        )
        
        ),
  
      tabItem(tabName = "menu_export",
              fluidRow(
                
                source("ui_files/ui_export_fitting_information.R",local=T)$value,
                source("ui_files/ui_export_plots_data.R",local=T)$value
                
              )),
      
      tabItem(tabName = "menu_simulate",
              fluidRow(
                
                source("ui_files/ui_model_selection_box_simulate.R",local=T)$value,                  
                source("ui_files/ui_advanced_settings_simulate.R"  ,local=T)$value,  
                source("ui_files/ui_explore_param_simulate.R"      ,local=T)$value,  

                tabBox(title = "", width = 10,id = "tabset_sim",
                       tabPanel("Fraction of occupied binding sites",plotlyOutput("fractionOccupied_sim_plot")),
                       tabPanel("Signal",plotlyOutput("signal_sim_plot"))
                       
                ),
                source("ui_files/ui_plot_options_box_sim.R",        local = TRUE)$value
                
              )),
      
      tabItem(tabName = "menu_user_guide",includeHTML("docs/user_guide.html")),
      tabItem(tabName = "menu_tutorial",includeHTML("docs/tutorial.html")),
      tabItem(tabName = "menu_about",includeHTML("docs/about.html"))
            
      ))))
