options(shiny.maxRequestSize=100*1024^2)
options(stringsAsFactors = F)
options(shiny.useragg = FALSE)

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
  
  cwd <- getwd()
  
  setwd(users_dir);   Sys.sleep(1)
  total_folders <- paste(count_folders(".")+1)
  
  # Create folder for this specific user # Useful when using Rstudio server.
  dir.create((total_folders)); Sys.sleep(1)
  setwd(total_folders)
  
  file.copy(paste0(cwd,"/www/testFluo.npy"),"testFluo.npy")
  file.copy(paste0(cwd,"/www/testTemp.npy"),"testTemp.npy")
  file.copy(paste0(cwd,"/www/testConcs.txt"),"testConcs.txt")
  
  # Delete all the files inside the folder we have created during this Session. 
  session$onSessionEnded(function() {
  #  setwd("..")
    system(paste0("rm -f ",users_dir,total_folders,"/*"))
  #  stopApp()
  })
  
}


