source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

packages <- c("shinydashboard","shinycssloaders","rhandsontable","tableHTML","plotly","shinyalert")
invisible(lapply(packages, library, character.only = TRUE))

shinyUI(dashboardPage(title = "MoltenProt",
  #logo_grey_light is described in logo.R
  dashboardHeader(title = logo_grey_light, titleWidth = 200),
  dashboardSidebar(collapsed = F,width = 200,
    sidebarMenu(
      
      menuItem("1. Load input",       icon = icon("chart-simple"),                tabName = "menu_input"),
      menuItem("2. Fitting",          icon = icon("chart-simple"),                tabName = "menu_fit"),
      menuItem("3. Analyze",          icon = icon("chart-simple"),                tabName = "menu_analyze"),
      menuItem("4. Export",           icon = icon("file-export"),                tabName = "menu_export"),
      menuItem("User guide",          icon = icon("user-circle"),       tabName = "menu_user_guide"),
      menuItem("Tutorial",            icon = icon("user-circle"),       tabName = "menu_tutorial"),
      menuItem("About",               icon = icon("circle-info"), tabName = "menu_about"))
    
      ),
  #theme_grey_light is described in theme.R
  dashboardBody(theme_grey_light,
    tabItems(
      
      tabItem(tabName = "menu_input",
        fluidRow(
          useShinyalert(),
          
          source("ui_files/menu_input/ui_load_input_box.R",       local = TRUE)$value,
          
          source("ui_files/menu_input/ui_conditions_table_box.R", local = TRUE)$value,
    

          tabBox(title = "", width = 10,id = "tabset1",
                  tabPanel("Signal",withSpinner(plotlyOutput("signal")))),
          
          
          source("ui_files/menu_input/ui_plot_options_box.R",       local = TRUE)$value,
          source("ui_files/menu_input/ui_derivative_plot_tabbox.R", local = TRUE)$value,
          source("ui_files/menu_input/ui_legend_tab_box.R",local = TRUE)$value
          
            )),
          
      tabItem(tabName = "menu_fit",
        fluidRow(
          
          source("ui_files/menu_fit/ui_model_selection_box.R",    local = TRUE)$value,
          source("ui_files/menu_fit/ui_fitting_options_box.R",    local = TRUE)$value,
          source("ui_files/menu_fit/ui_sort_fitted_params_box.R", local = TRUE)$value,
          source("ui_files/menu_fit/ui_plot_options_box.R",       local = TRUE)$value,
          source("ui_files/menu_fit/ui_tabbox.R",                 local = TRUE)$value,
          source("ui_files/menu_fit/ui_legend_fitting_tab_box.R",local = TRUE)$value

          )),
      
      tabItem(tabName = "menu_analyze",
        fluidRow(
          
          source("ui_files/menu_analyze/ui_filter_box.R",                    local = TRUE)$value,
          source("ui_files/menu_analyze/ui_plot_options_box.R",              local = TRUE)$value,
          source("ui_files/menu_analyze/ui_score_table_resultPlot_tabbox.R", local = TRUE)$value

          )),
      
      tabItem(tabName = "menu_export",
              fluidRow(
                
                source("ui_files/menu_export/ui_report.R",              local = TRUE)$value,
                source("ui_files/menu_export/ui_fitting_information.R", local = TRUE)$value,
                source("ui_files/menu_export/ui_plots_data.R",          local = TRUE)$value
                
                )),
                                  
      tabItem(tabName = "menu_user_guide", includeHTML("docs/user_guide.html")),
      tabItem(tabName = "menu_tutorial",   includeHTML("docs/tutorial.html")),
      tabItem(tabName = "menu_about",      includeHTML("docs/about.html"   ))

      ))))

