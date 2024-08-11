source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

shinyUI(dashboardPage(title = paste0(appName),
  
  dashboardHeader(  title = logo_grey_light, titleWidth = 200), #logo_grey_light is described in logo.R
  dashboardSidebar( collapsed = F,width = 200,
                   
    sidebarMenu(
      
      menuItem("1. Load input",       icon = icon("file-circle-plus"),       tabName = "menu_input"     ),
      menuItem("2. Analysis",         icon = icon("chart-line"),             tabName = "menu_analysis"  ),
      menuItem("3. Export",           icon = icon("file-export"),            tabName = "menu_export"    ),
      menuItem("Simulation",          icon = icon("magnifying-glass-chart"), tabName = "menu_simulate"  ),
      menuItem("User guide",          icon = icon("user-astronaut"),         tabName = "menu_user_guide"),
      menuItem("About",               icon = icon("circle-info"),            tabName = "menu_about"     ))
    
      ),
  
  dashboardBody(theme_grey_light,
    tabItems(
      tabItem(tabName = "menu_input",
        fluidRow(
          
          source("ui_files/ui_load_input_box.R",local = TRUE)$value,
          
          tabBox(title = "", width = 9,id = "tabBoxDlsWellsInfo")),
          
        fluidRow(
        
          conditionalPanel('output.data_loaded',
            conditionalPanel('input.setExperimentParameters',
            source("ui_files/ui_experimentParameters_box.R",local = TRUE)$value),                           
           source("ui_files/ui_preprocessing_box.R",local = TRUE)$value,
           source("ui_files/ui_signal_tab_box0.R",local = TRUE)$value))
        
         ),

      tabItem(tabName = "menu_analysis",
                
              fluidRow(
                
                source("ui_files/ui_analysis_box.R",local = TRUE)$value,
                source("ui_files/ui_plotDownloadOptions_box.R",local = TRUE)$value),
              
              fluidRow(
                
                #Custom CSS
                tags$head(tags$style("
                     
                     #autocorrelation{height:600px !important;}
                     #fittedAutocorrelation1{height:600px !important;}
                     #intensityDistribution{height:600px !important;}
                     #intensityDistributionPlotly{height:600px !important;}
                     #residuals{height:600px !important;}
                     #lCurve{height:600px !important;}
                     #intensityDistributionAvg{height:600px !important;}

                     ")),
                
                source("ui_files/ui_signal_tab_box1.R",   local = TRUE)$value,
                source("ui_files/ui_signal_tab_box2.R",   local = TRUE)$value,
                source("ui_files/ui_peakSelection_box.R", local = TRUE)$value,
                source("ui_files/ui_parameters_tabBox.R", local = TRUE)$value
                
              )),
      
      tabItem(tabName = "menu_export",
              fluidRow(
                
                source("ui_files/ui_export_fitting_information.R",local=T)$value
                
              )),
      
      tabItem(tabName = "menu_simulate",
              
              fluidRow(
                
                source("ui_files/ui_simulateParams_box.R",local = TRUE)$value
                
              ),
                #source("ui_files/ui_Simulate_plotDownloadOptions_box.R",local = TRUE)$value),
              
              fluidRow(
                
                #Custom CSS
                tags$head(tags$style("
                     
                     #autocorrelationSimulate{height:600px !important;}
                     #hrDistributionNumberSim{height:600px !important;}
                     #hrDistributionVolumeSim{height:600px !important;}
                     #hrDistributionIntensitySim{height:600px !important;}
                    
                     ")),
                
                source("ui_files/ui_signalSimulate_tab_box2.R",local = TRUE)$value,
                source("ui_files/ui_signalSimulate_tab_box1.R",local = TRUE)$value)
                
              ),
      
      tabItem(tabName = "menu_user_guide",includeHTML("docs/user_guide.html")),
      tabItem(tabName = "menu_tutorial",includeHTML("docs/tutorial.html")),
      tabItem(tabName = "menu_about",includeHTML("docs/about.html"))
            
      ))))



