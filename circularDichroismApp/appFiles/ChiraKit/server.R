options(shiny.maxRequestSize=100*1024^2)
options(stringsAsFactors = F)

# List of Python script files to source
py_scripts <- c("cdAnalyzer.py","helpers.py","decomposition_helpers.py",
                "fitting_helpers.py","get_dssp_summary.py")

# Source the Python helper functions 
for (script in py_scripts) {source_python(paste0('python_src/', script))}

# List of R script files to source
r_scripts <- c("helpers.R","helpers_unfolding.R","helpers_plotting.R",
  "helpers_G-Quadruplex.R","plotFunctions.R","plotFunctionsSpectraComparison.R")

# Source the R helper functions 
for (script in r_scripts) {source(paste0('server_files/', script))}

### End of variables to change

function(input, output, session) {
    
  welcomeMessage() # helpers.R
  
  # To handle the general processing, the unfolding models, 
  # the secondary structure calculation, and the custom models
  cdAnalyzer                <- CdAnalyzer()
  
  # To handle the spectra comparison module
  compareSpectraPyClass     <- CdExperimentComparison()
  
  # To handle the G-Quadruplex module
  gQuadRefPyClass     <- CdExperimentGeneral()
  gQuadSamplePyClass  <- CdExperimentGeneral()
  
  source(paste0(base_dir,"reactives/reactives_values.R"                 ), local = T)
  source(paste0(base_dir,"reactives/reactives.R"                        ), local = T)
  source(paste0(base_dir,"reactives/plot_reactives_load_input.R"        ), local = T)
  source(paste0(base_dir,"reactives/thermal_ramp_reactives.R"           ), local = T)
  source(paste0(base_dir,"reactives/chemical_denaturation_reactives.R"  ), local = T)
  source(paste0(base_dir,"reactives/secondary_structure_reactives.R"    ), local = T)
  source(paste0(base_dir,"reactives/custom_analysis_reactives.R"        ), local = T)
  source(paste0(base_dir,"reactives/spectra_comparison_reactives.R"     ), local = T)
  source(paste0(base_dir,"reactives/peptide_helix_content_reactives.R"  ), local = T)
  source(paste0(base_dir,"reactives/gQuadruplex_reactives.R"            ), local = T)
  source(paste0(base_dir,"reactives/download_reactives.R"               ), local = T)
  
}

