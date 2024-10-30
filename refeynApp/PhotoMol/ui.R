source("ui_files/theme.R")
source("ui_files/logo.R")
source("ui_files/busy_indicator.R")

shinyUI(dashboardPage(title = paste0(appName),
  
  dashboardHeader(  title = logo_grey_light, titleWidth = 200), #logo_grey_light is described in logo.R
  dashboardSidebar( collapsed = F,width = 200,
                   
    sidebarMenu(
      
      menuItem("Analyze",                  icon = icon("chart-simple"),                tabName = "menu_input"),
      menuItem("Calibration",              icon = icon("chart-simple"),                tabName = "menu_calibration"),
      menuItem("Export",                   icon = icon("file-export"),                tabName = "menu_export"),
      menuItem("User guide",               icon = icon("user-circle"),       tabName = "menu_user_guide"),
      menuItem("About",                    icon = icon("circle-info"), tabName = "menu_about"))
    
      ),
  
  dashboardBody(theme_grey_light,
#    tags$head(includeHTML(("ui_files/google-analytics.html"))),                
    tabItems(
      tabItem(tabName = "menu_input",
        fluidRow(

          source("ui_files/ui_load_input_box.R",local = TRUE)$value, 

                    #Custom CSS to increase plot height
          tags$head(tags$style("
          #counts_plot{height:600px !important;}
          #counts_plotNormalized{height:600px !important;}
          #binding_plot{height:600px !important;}
          #counts_plot_stacked{height:700px !important;}
          #counts_plotNormalized_stacked{height:700px !important;}
          "
                               )),

          tabBox(title = "", width = 9,id = "tabset_sim",
                   tabPanel("Mass histogram",plotlyOutput("counts_plot")),
                   tabPanel("Normalised mass histogram",plotlyOutput("counts_plotNormalized")),
                   tabPanel("All (un)binding events",plotlyOutput("binding_plot")),
                   tabPanel("Mass histogram (subplots)",plotlyOutput("counts_plot_stacked")),
                   tabPanel("Normalised mass histogram (subplots)",plotlyOutput("counts_plotNormalized_stacked"))),

          conditionalPanel(condition = "output.dataLoaded", 
                           source('ui_files/ui_load_input_legend.R',local=T)$value)),

          fluidRow(

          box(title = "Fitted parameters and counts", width = 5, solidHeader = T, status = "primary",
              fluidRow(
                column(width = 12,tableOutput('fittedParams'))
                
              )),
          
          box(title = "Legends", width = 4, solidHeader = T, status = "primary",
              fluidRow(
                conditionalPanel(condition = "output.dataLoaded",

                                 column(width = 8,
                                        rHandsontableOutput('legendInfo',    width = "100%"),
                                        rHandsontableOutput('legendInfoHist',width = "100%")),
                                 
                                 column(width = 4,
                                        selectInput("mol2changeColor", "Set colour",c("X")),
                                        colourInput("colorForLegend", NULL, value = "#E41A1C"),
                                 p(HTML("<b>Show estimated masses</b>"),
                                             span(shiny::icon("info-circle"), id = "info_uuL-1"),
                                             checkboxInput("show_massesPlot", "In the plot", TRUE),
                                             checkboxInput("show_massesLegend", "In the legend", FALSE),
                                             tippy::tippy_this(elementId = "info_uuL-1",
                                                               tooltip = "Display the means of the fitted gaussians
                                                               .",placement = "right")),
                                 p(HTML("<b>Show counts percentage</b>"),
                                   span(shiny::icon("info-circle"), id = "info_uuL-2"),
                                   checkboxInput("show_percentagePlot", "In the plot", TRUE),
                                   checkboxInput("show_percentageLegend", "In the legend", FALSE),
                                   tippy::tippy_this(elementId = "info_uuL-2",
                                                     tooltip = 
                                                     "Display the percentage of counts for each peak.",
                                                     placement = "right")))
                                 
                )
              )),
          
          source("ui_files/ui_export_plot_box.R",local = TRUE)$value),
        fluidRow(
          source("ui_files/ui_simulate_box.R",local = TRUE)$value,
          conditionalPanel('input.activateCalibration',
            source("ui_files/ui_export_H5file_with_masses.R",local = TRUE)$value)
        
        )),
      tabItem(tabName = "menu_calibration",
              fluidRow(
                source("ui_files/ui_load_input_box_calibration.R",local = TRUE)$value, 
                
                conditionalPanel(condition = "output.calibrationMethod == 'calibrationFile' 
                                 && input.activateCalibration",
                
                tabBox(title = "", width = 9,id = "tabset_calib",
                       tabPanel("Counts vs Contrast",     plotlyOutput("contrast_plot_calib")),
                       tabPanel("Masses vs Contrast",     plotlyOutput("mass_vs_contrast")),
                       tabPanel("Calibration parameters", tableOutput("calibParams"))),
                
                box(title = "Fitted parameters and counts", width = 4, solidHeader = T, status = "primary",
                    fluidRow(
                      column(width = 12,tableOutput('fittedParamsCalibration'))
                      
                    )),
                
                box(title = "Legends", width = 3, solidHeader = T, status = "primary",
                    fluidRow(
                      conditionalPanel(condition = "output.data_loadedCalibration", 
                         column(width = 8,rHandsontableOutput('legendInfoCalibration')),
                         column(width = 4,
                                selectInput("mol2changeColorCalibration", "Set colour",c("X")),
                                colourInput("colorForLegendCalibration", NULL, value = "#E41A1C")),
                         column(8, p(HTML("<b>Show estimated contrasts</b>"),
                                     span(shiny::icon("info-circle"), id = "info_uuCalib-1"),
                                     checkboxInput("show_contrastPlot", "In the plot", TRUE),
                                     checkboxInput("show_contrastLegend", "In the legend", TRUE),
                                     tippy::tippy_this(elementId = "info_uuCalib-1",
                                                       tooltip = "Display the means of the fitted gaussians
                                                 .",placement = "right")))
                                       
                      )
                    )),
                
                source("ui_files/ui_export_plot_box_calibration.R",local = TRUE)$value
                
                )
                
                
              )),
      
      tabItem(tabName = "menu_export",
              fluidRow(
                
                source("ui_files/ui_export_fitting_information.R",local=T)$value,
                source("ui_files/ui_export_plots_data.R",local=T)$value
                
              )),
      
      tabItem(tabName = "menu_user_guide",includeHTML("docs/user_guide.html")),
      tabItem(tabName = "menu_about",includeHTML("docs/about.html"))
            
      ))))