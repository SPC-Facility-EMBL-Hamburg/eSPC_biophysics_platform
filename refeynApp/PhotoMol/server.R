options(shiny.maxRequestSize=30*1024^2)
options(stringsAsFactors = F)

source("server_files/helpers.R")
source_python("refeyn.py")
source_python("refeynCalibration.py")
source_python("helpers.py")
source("server_files/plot_functions.R")
### End of variables to change

function(input, output, session) {
  
  welcomeMessage() # helpers.R
  
  refeyn      <- Refeyn()
  refeynCalib <- RefeynCalib()
  
  source(paste0(base_dir,"reactives/reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/reactivesCalibration.R"), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"  ), local = T)

  cwd <- getwd()
  
  # Create folder for this specific user
  setwd(users_dir);   Sys.sleep(1)
  total_folders <- paste(count_folders(".")+1)
  dir.create((total_folders)); Sys.sleep(1)
  setwd(total_folders)
  
  file.copy(paste0(cwd,"/www/demo.h5"),"demoTestFile")
  
  # Delete all the files inside the folder we have created during this Session. 
  session$onSessionEnded(function() {
    #  setwd("..")
    system(paste0("rm -f ",users_dir,total_folders,"/*"))
    #  stopApp()
  })
}
