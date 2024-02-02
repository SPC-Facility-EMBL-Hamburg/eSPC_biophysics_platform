options(shiny.maxRequestSize=100*1024^2)
options(stringsAsFactors = F)

source("server_files/helpers.R")
source("server_files/helpers_unfolding.R")
source("server_files/helpers_plotting.R")
source("server_files/plotFunctions.R")
source("server_files/plotFunctionsSpectraComparison.R")

source_python("./cdAnalyzer.py")
source_python("./cdUnitsConverter.py")
source_python("./helpers.py")
source_python('./loadCDfilesHelpers.py')
source_python("./get_dssp_summary.py")
source_python("./read_references_to_matrices.py")
### End of variables to change

function(input, output, session) {
    
  welcomeMessage() # helpers.R
  
  # To handle the general processing, the unfolding models, 
  # the secondary structure calculation, and the custom models
  cdAnalyzer                <- cdAnalyzer()     
  
  # To handle the spectra comparison module
  compareSpectraPyClass     <- cd_experiment_comparison()
  
  source(paste0(base_dir,"reactives/reactives_values.R"                 ), local = T)
  source(paste0(base_dir,"reactives/reactives.R"                        ), local = T)
  source(paste0(base_dir,"reactives/plot_reactives_load_input.R"        ), local = T)
  source(paste0(base_dir,"reactives/thermal_ramp_reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/chemical_denaturation_reactives.R"  ), local = T)
  source(paste0(base_dir,"reactives/secondary_structure_reactives.R"    ), local = T)
  source(paste0(base_dir,"reactives/custom_analysis_reactives.R"        ), local = T)
  source(paste0(base_dir,"reactives/spectra_comparison_reactives.R"     ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"               ), local = T)
  
  # Create folder for this specific user
  setwd(users_dir);   Sys.sleep(1)
  total_folders <- paste(count_folders(".")+1)
  dir.create((total_folders)); Sys.sleep(1)
  setwd(total_folders)
  
  # Delete all the files inside the folder we have created during this Session. 
  session$onSessionEnded(function() {
    #  setwd("..")
    system(paste0("rm -f ",users_dir,total_folders,"/*"))
    #  stopApp()
  })
}

