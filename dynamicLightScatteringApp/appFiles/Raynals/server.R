options(shiny.maxRequestSize=30*1024^2)
options(stringsAsFactors = F)

source("server_files/helpers.R")
source("server_files/plot_functions.R")
source_python("./dlsAnalyzer.py")
source_python("./helpers.py")
source_python("./simulation_helpers.py")
source("server_files/simulation_helpers.R")

function(input, output, session) {
  
  welcomeMessage() 

  source(paste0(base_dir,"reactives/reactives.R"          ), local = T)
  source(paste0(base_dir,"reactives/plot_reactives.R"     ), local = T)
  source(paste0(base_dir,"reactives/simulate_reactives.R" ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R" ), local = T)
  
  dlsAnalyzer  <- dlsAnalyzer()
  
}
