options(shiny.maxRequestSize=30*1024^2)
options(stringsAsFactors = F)

source("server_files/helpers.R")
source("server_files/plot_functions.R")
source_python("./dlsAnalyzer.py")
source_python("./helpers.py")
source_python("./simulation_helpers.py")
source("server_files/simulation_helpers.R")
### End of variables to change  

function(input, output, session) {
  
  welcomeMessage() # helpers.R

  source(paste0(base_dir,"reactives/reactives.R"          ), local = T)
  source(paste0(base_dir,"reactives/plot_reactives.R"     ), local = T)
  source(paste0(base_dir,"reactives/simulate_reactives.R" ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R" ), local = T)
  
  cwd <- getwd()
  
  # Create folder for this specific user
  setwd(users_dir);   Sys.sleep(1)
  total_folders <- paste(count_folders(".")+1)
  dir.create((total_folders)); Sys.sleep(1)
  setwd(total_folders)
  
  file.copy(paste0(cwd,"/www/test.csv"),"test.csv")
  
  dlsAnalyzer  <- dlsAnalyzer()
  
  # Delete all the files inside the folder we have created during this Session. 
  session$onSessionEnded(function() {
    #  setwd("..")
    system(paste0("rm -f ",users_dir,total_folders,"/*"))
    #  stopApp()
  })
}
