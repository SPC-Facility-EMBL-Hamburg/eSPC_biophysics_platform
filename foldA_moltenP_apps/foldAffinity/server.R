options(shiny.maxRequestSize=100*1024^2)
options(stringsAsFactors = F)

source("server_files/global_variables.R")

source_python("fitting_helpers_thermal_unfolding.py")
source_python("fitting_helpers_unfolded_fraction.py")
source_python("foldAffinity.py")

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/simulation_helpers.R")
source("server_files/plot_functions.R")
source("server_files/fitting_helpers.R")
### End of variables to change

function(input, output, session) {
  
  # ... Initialize python_nanoDSF class ...
  dsf <- DSF_binding()
  
  welcome_message()
  
  source(paste0(base_dir,"reactives/load_input_reactives.R"), local = T)
  source(paste0(base_dir,"reactives/reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"  ), local = T)
  
}


