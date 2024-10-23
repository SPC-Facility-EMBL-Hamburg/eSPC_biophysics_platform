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

  refeynCalib <- RefeynCalib()

  photoMolModels <- PhotoMolModels()

  source(paste0(base_dir,"reactives/reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/reactivesCalibration.R"), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"  ), local = T)

}
