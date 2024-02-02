# The ChiraKit app

Last time updated: January 2024

## Introduction

ChiraKit is an online multi-purpose tool developed to investigate circular dichroism (CD) data. 

The workflow of ChiraKit is divided into four steps:

1) Data Importing: Upload raw CD data from various instruments (e.g.,'Applied Photophysics' and 'JASCO') and file formats (e.g., .dat, .pcd, .gen). 

2) Preprocessing: Add, subtract, zero, smooth, and/or average the spectra. Choose between different working units, such as millidegrees, differential absorbance, molar extinction, and mean residue molar ellipticity.

3) Analysis: Calculate the protein secondary structure, fit user-defined models, or explore the CD data as a function of a certain experimental parameter. Noteworthy functionalities include two-state models for thermal and chemical unfolding.

4) Export: Export the finalised spectra, fitted parameters, and fitted curves.

Key functionalities:

Secondary Structure Prediction: Apply the Selcon3 algorithm to predict the protein secondary structure by comparing the CD spectrum of interest to a reference set.

Spectra Decomposition: Employ Singular Value Decomposition (SVD) or Principal Component Analysis (PCA) to decompose multiple spectra into a set of basis spectra.

Thermal and Chemical Unfolding: Fit unfolding curves with a reversible two-state model.

Customized Data Analysis: Study the CD signal as a function of any desired experimental parameter.

## Getting started

To run the app locally you need to 

1) Install R (tested with version 4.2.2)
2) Install python (tested with version 3.10.9) 
3) Install the required R packages (it may take a long time):

``` R 
Rscript install_r_packages.R
```

4) Install conda/miniconda:

``` bash 
user=$(whoami) 
wget --no-verbose https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh  -P /home/${user}/               
bash  /home/${user}/Miniconda3-latest-Linux-x86_64.sh -b                                                                
rm -f /home/${user}/Miniconda3-latest-Linux-x86_64.sh
```

5) Install the required python packages in the conda environment 'r-reticulate':

```bash
/home/${user}/miniconda3/bin/conda  create    --name r-reticulate                                                      
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda nomkl               
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda numpy     
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda pandas    
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda scipy     
```

6) Set the correct path for the app

``` bash 
if [ "$(basename "$(pwd)")" = "circularDichroismApp" ]; then
    sed -i "s|base_dir <- paste0.*|base_dir <- paste0('$PWD', '/appFiles/ChiraKit/')|" appFiles/ChiraKit/global.R
else
    echo "Change the working directory to circularDichroismApp"
fi
```

7) Create a folder to store temporary files:

``` R 
mkdir -p "/home/${user}/data_users/" 
```

8) Run the shiny app with R:

``` R 
cd appFiles/ChiraKit
shiny::runApp()
```

## General code structure

The following tree represents how the files are organised and their purpose:
```
circularDichroismApp
   |-- appFiles                     # Folder with the shiny app - used to create the docker image
   |   |-- Rprofile.site            # R configuration file (used inside the docker image)
   |   |-- install_r_packages.R     # R packages required by the shiny app - to be installed inside the docker image
   |   |-- mkdssp-4.4.0-linux-x64   # linux executable file of the DSSP algorithm 
   |   |-- ChiraKit 
   |   |   |-- global.R             # Global parameters to be loaded at the start of the app
   |   |   |-- server.R             # Server file 
   |   |   |-- ui.R                 # Main UI file 
   |   |   |-- cdAnalyzer.py        # Core script to analyse CD data. It works as a standalone programm too
   |   |   |-- cdUnitsConverter.py  # Helpers functions to convert between CD units. E.g. milidegrees to delta epsilon
   |   |   |-- fitting_helpers.py   # Helpers functions to apply PCA/SVD, fit the unfolding datasets, fit the custom datasets
   |   |   |-- get_dssp_summary.py  # Run the DSSP algorithm given a certain protein input file. 
   |   |   |-- helpers.py           # General helper functions
   |   |   |-- loadCDfilesHelpers.py # Helper functions to parse the CD data (in different formats)
   |   |   |-- read_references_to_matrices.py # Helper functions to load custom reference datasets for the sec. str. calculation
   |   |   |-- runApp.r                       # Used in the dev process to run the shiny app
   |   |   |-- test_reticulate.R              # Used in the dev process to test new functionalities
   |   |   |-- reactives                                 # Contains the required scripts to make the app interactive
   |   |   |   |-- chemical_denaturation_reactives.R     # Scripts to allow analysis of chemical unfolding curves
   |   |   |   |-- custom_analysis_reactives.R           # Scripts to allow analysis of custom curves
   |   |   |   |-- download_reactives.R                  # Scripts to allow exporting the data and results
   |   |   |   |-- plot_reactives_load_input.R           # Scripts with plots from the 1. Load Input Tab
   |   |   |   |-- reactives.R                           # Scripts to import and preprocess the CD spectra 
   |   |   |   |-- reactives_values.R                    # General reactive values that are used throughout the app
   |   |   |   |-- secondary_structure_reactives.R       # Scripts to calculate the secondary structure from the spectra 
   |   |   |   |-- thermal_ramp_reactives.R              # Scripts to allow analysis of thermal unfolding curves
   |   |   |-- server_files                              
   |   |   |   |-- helpers.R                          # helper functions used inside the server side
   |   |   |   |-- helpers_plotting.R                 # helper functions for plotting the data
   |   |   |   |-- helpers_unfolding.R                # helper functions to process the unfoldin and 'custom' data
   |   |   |   |-- plotFunctions.R                    # collection of plotting functions
   |   |   |-- docs                                   
   |   |   |   |-- about.html                         # HTML text for the About section
   |   |   |   |-- tutorial.html                      # HTML text for the Tutorial section (not used for now!!!!)
   |   |   |   |-- user_guide.html                    # HTML text for the User Guide section
   |   |   |-- secondary_structure_estimation_files   # Files provided by Søren Vrønning Hoffmann and Nykola Jones (Aarhus University, Denmark)
   |   |   |   |-- AU-A128_PCDDB-Nov22.txt            # A (also called 'C') matrix with the reference spectra (AU-SP175 and AUSMP180)
   |   |   |   |-- AU-F128_T-Nov22.txt                # F matrix with the reference secondary structure elements (AU-SP175 and AUSMP180)
   |   |   |   |-- Labels-SMP180_PCDDBOrder.txt       # Labels of the A matrix
   |   |   |   |-- SelconsFunction.py                 # Python implementation of the Selcon algorithm
   |   |   |   |-- TestProteins                       # Set of test proteins for SelconsFunction.py - you can use Lyzozyme.txt directly in the shiny app
   |   |   |-- ui_files                                     # UI elements
   |   |   |   |-- 1ui_load_experiment_parameters.R         # To handle the associated experimental parameters (e.g., molecular weight, concentration)
   |   |   |   |-- 1ui_load_input_box.R                     # To import the data and set the wavelength range
   |   |   |   |-- 1ui_plotting_box.R                       # To plot the CD spectra and the HT curves
   |   |   |   |-- 1ui_processing.R                         # To allow processing the spectra (addition, subtraction, average, ...)
   |   |   |   |-- 2a_ui_fitting_thermal.R                  # To fit  the thermal unfolding curves
   |   |   |   |-- 2a_ui_load_thermal_denaturation_data.R   # To load the thermal unfolding data
   |   |   |   |-- 2b_ui_fitting_chemical.R                 # To fit  the chemical unfolding curves
   |   |   |   |-- 2b_ui_load_chemical_denaturation_data.R  # To load the chemical unfolding data
   |   |   |   |-- 2c_ui_fit_secondary_structure_box.R      # To calculate the sec. str. with default or user provided reference sets
   |   |   |   |-- 2d_ui_custom_model_builder.R             # To fit the custom curves
   |   |   |   |-- 2d_ui_load_custom_data.R                 # To load the custom data
   |   |   |   |-- busy_indicator.R                         # Function to show the user a nice ...loading... image
   |   |   |   |-- logo.R                                   # Configuration for the app logo
   |   |   |   |-- theme.R                                  # Theme of the app
   |   |   |   |-- ui_export_cd_spectra.R                   # To export the finalised spectra
   |   |   |   |-- ui_export_cd_spectra_chemical_fit.R      # To export the chemical unfolding data
   |   |   |   |-- ui_export_cd_spectra_custom_fit.R        # To export the 'custom' data
   |   |   |   |-- ui_export_cd_spectra_sec_str_fit.R       # To export the sec. str. data
   |   |   |   |-- ui_export_cd_spectra_thermal_fit.R       # To export the thermal unfolding data
   |   |   |   |-- ui_export_logbook.R                      # To export a text file in order to reproduce the analysis
   |   |   |-- www                                          # files available to the shiny app 
```

## Acknowledgments

The ChiraKit app is possible thanks to:

R language: R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

R package shiny:   Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2020). shiny: Web Application Framework for R. R package version 1.4.0.2. https://CRAN.R-project.org/package=shiny

R package tidyverse: Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

R package shinydashboard:   Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard

R package ggplot2:   H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

R package reshape2:   Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

R package signal: signal developers (2013). signal: Signal processing. URL:
  http://r-forge.r-project.org/projects/signal/.

R package tippy:   John Coene (2018). tippy: Add Tooltips to 'R markdown' Documents or 'Shiny' Apps. R package version 0.0.1. https://CRAN.R-project.org/package=tippy

R package shinyalert:   Pretty Popup Messages (Modals) in 'Shiny'. R package version 1.1. https://CRAN.R-project.org/package=shinyalert

R package plotly:   C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.

R package rhandsontable:   Jonathan Owen (2018). rhandsontable: Interface to the 'Handsontable.js' Library. R package version 0.3.7. https://CRAN.R-project.org/package=rhandsontable

R package shinyjs:   Dean Attali (2020). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 1.1. https://CRAN.R-project.org/package=shinyjs

R package reticulate:   Kevin Ushey, JJ Allaire and Yuan Tang (2020). reticulate: Interface to 'Python'. R package version 1.16. https://CRAN.R-project.org/package=reticulate

R package shinycssloaders:   Andras Sali and Dean Attali (2020). shinycssloaders: Add CSS Loading Animations to 'shiny' Outputs. R package version 0.3. https://CRAN.R-project.org/package=shinycssloaders

R package DT: Xie Y, Cheng J, Tan X (2022). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.25, https://CRAN.R-project.org/package=DT.

R package colourpicker: Attali D (2021). _colourpicker: A Colour Picker Tool for Shiny and for Selecting Colours in Plots_. R package version 1.1.1, <https://CRAN.R-project.org/package=colourpicker>.

Python3.7 language: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.