options(shiny.maxRequestSize=100*1024^2)
options(stringsAsFactors = F)

source("server_files/global_variables.R")

source_python("helpers.py")
source_python("moltenprot_shiny.py")

source("server_files/load_input_helpers.R")
source("server_files/helpers.R")
source("server_files/plot_functions.R")

function(input, output, session) {
  
  # ... Initialize python_nanoDSF class ...
  dsf <- DSF_molten_prot_fit()
  
  welcome_message()
  
  source(paste0(base_dir,"reactives/reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"  ), local = T)
  
}
