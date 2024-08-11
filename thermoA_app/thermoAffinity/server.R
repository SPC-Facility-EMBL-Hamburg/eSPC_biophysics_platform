options(shiny.maxRequestSize=30*1024^2)
options(stringsAsFactors = F)

source("server_files/global_variables.R")

source_python("mst.py")
source_python("helpers.py")

scripts <- list.files("server_files/fitting_helpers/")
lapply(paste0("server_files/fitting_helpers/",scripts), source)

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/plot_functions.R")
source("server_files/simulation_helpers.R")
source("server_files/explore_parameters_helpers.R")
### End of variables to change

function(input, output, session) {
  
  # ... Initialize python class ...
  mst <- MST_fit()
  
  welcome_message() 
  
  source(paste0(base_dir,"reactives/load_input_reactives.R"), local = T)
  source(paste0(base_dir,"reactives/reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"  ), local = T)
  
}
