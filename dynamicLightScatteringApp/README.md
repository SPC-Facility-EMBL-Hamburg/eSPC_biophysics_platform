# The Raynals app

Last time updated: March 2023

## Introduction

This folder contains a shiny app to analyse dynamic light scattering (DLS) data.
The fitting is done using the Tikhonov Philips regularisation and the L-curve criteria.

The input data for Raynals is the second-order autocorrelation and the output data 
are the estimated intensity weighted contributions of the hydrodynamic radii. 
Example data can be loaded directly when running the app.

## Getting started

To run the app locally you need to 

1) Install R (tested with version 4.2.2)
2) Install python (tested with version 3.10.9) 
3) Install the required R packages (it may take a long time):

``` R 
Rscript ./appFiles/install_r_packages.R
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
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda pip     
/home/${user}/miniconda3/envs/r-reticulate/bin/pip install miepython        
```

6) Set the correct path for the app

``` bash 
if [ "$(basename "$(pwd)")" = "dynamicLightScatteringApp" ]; then
    sed -i "s|base_dir <- paste0.*|base_dir <- paste0('$PWD', '/appFiles/Raynals/')|" appFiles/Raynals/global.R
else
    echo "Change the working directory to dynamicLightScatteringApp"
fi
```

7) Create a folder to store temporary files:

``` R 
mkdir -p "/home/${user}/data_users/" 
```

8) Run the shiny app with R:

``` R 
cd appFiles/Raynals
shiny::runApp()
```

## General code structure

The following tree represents how the files are organised and their purpose:
```
   dynamicLightScatteringApp
      |-- appFiles                           # Folder with the shiny app - used to create the docker image
      |   |-- Rprofile.site                  # R configuration file (used inside the docker image)
      |   |-- install_r_packages.R           # R packages required by the shiny app - to be installed inside the docker image
      |   |-- Raynals
      |   |   |-- global.R                   # Global parameters to be loaded at the start of the app
      |   |   |-- server.R                   # Server file 
      |   |   |-- ui.R                       # Main UI file
      |   |   |-- dlsAnalyzer.py             # Core script to analyse DLS data. It works as a standalone command line programm too
      |   |   |-- helpers.py                 # Helpers function to run dlsAnalyzer.py
      |   |   |-- loadDLSdataHelpers.py      # Functions to load DLS input files (second-order autocorrelation curves)
      |   |   |-- reactives                  # Contains the required scripts to make the app interactive
      |   |   |   |-- download_reactives.R   # Scripts to download the data
      |   |   |   |-- plot_reactives.R       # Scripts to plot the data
      |   |   |   |-- reactives.R            # Scripts to load the input data and create tables
      |   |   |   |-- simulate_reactives.R   # Scripts to run the simulations
      |   |   |-- server_files
      |   |   |   |-- helpers.R              # helper functions used inside the server side
      |   |   |   |-- plot_functions.R       # plot functions, called by plot_reactives.R
      |   |   |   |-- simulation_helpers.R   # simulation functions, called by simulate_reactives.R
      |   |   |-- simulation_helpers.py      # simulation helpers
      |   |   |-- test_reticulate.R          # Used in the dev process to test new functionalities
      |   |   |-- ui_files                   # UI elements
      |   |   |   |-- busy_indicator.R                # function to show the user a nice ...loading... image
      |   |   |   |-- logo.R                          # Configuration for the app logo
      |   |   |   |-- theme.R                         # Theme of the app
      |   |   |   |-- ui_analysis_box.R               # To run the fitting of the data
      |   |   |   |-- ui_experimentParameters_box.R   # To define the instrument setup and experimental conditions
      |   |   |   |-- ui_export_fitting_information.R # To export the fitted data
      |   |   |   |-- ui_load_input_box.R             # To load the input file (or example data)
      |   |   |   |-- ui_parameters_tabBox.R          # To show a Table with the information about the estimated hydrodynamic radii
      |   |   |   |-- ui_peakSelection_box.R          # To define the hydrodynamic radii regions of interest 
      |   |   |   |-- ui_plotDownloadOptions_box.R    # To change the aspect of the plot
      |   |   |   |-- ui_preprocessing_box.R          # To apply some preprocessing to the data
      |   |   |   |-- ui_signalSimulate_tab_box1.R    # To plot the simulated autocorrelation curve
      |   |   |   |-- ui_signalSimulate_tab_box2.R    # To plot the simulated distributions (intensity, volume, number)
      |   |   |   |-- ui_signal_tab_box0.R            # To plot the loaded autocorrelation data
      |   |   |   |-- ui_signal_tab_box1.R            # To plot the loaded and predicted autocorrelation data
      |   |   |   |-- ui_signal_tab_box2.R            # To plot the estimated contributions
      |   |   |   |-- ui_simulateParams_box.R         # To define the simulation parameters
      |   |   |-- www                                 # Contains images for the about section and example input files
      |   |   |-- docs
      |   |   |   |-- about.html             # HTML text for the about section 
      |   |   |   |-- user_guide.html        # HTML text for the user guide section 
```

## Future development 

Would be great to add support for multi-angle DLS data and temperature ramp analyses

## Acknowledgments

The Raynals app is possible thanks to:

R language: R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

R package shiny:   Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2020). shiny: Web Application Framework for R. R package version 1.4.0.2. https://CRAN.R-project.org/package=shiny

R package viridis: Simon Garnier (2018). viridis: Default Color Maps from 'matplotlib'. R package version 0.5.1. https://CRAN.R-project.org/package=viridis

R package tidyverse: Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

R package pracma: Hans W. Borchers (2019). pracma: Practical Numerical Math Functions. R package version 2.2.9. https://CRAN.R-project.org/package=pracma

R package shinydashboard:   Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard

R package ggplot2:   H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

R package reshape2:   Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

R package tippy:   John Coene (2018). tippy: Add Tooltips to 'R markdown' Documents or 'Shiny' Apps. R package version 0.0.1. https://CRAN.R-project.org/package=tippy

R package shinyalert:   Pretty Popup Messages (Modals) in 'Shiny'. R package version 1.1. https://CRAN.R-project.org/package=shinyalert

R package plotly:   C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.

R package rhandsontable:   Jonathan Owen (2018). rhandsontable: Interface to the 'Handsontable.js' Library. R package version 0.3.7. https://CRAN.R-project.org/package=rhandsontable

R package shinyjs:   Dean Attali (2020). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 1.1. https://CRAN.R-project.org/package=shinyjs

R package reticulate:   Kevin Ushey, JJ Allaire and Yuan Tang (2020). reticulate: Interface to 'Python'. R package version 1.16. https://CRAN.R-project.org/package=reticulate

R package shinycssloaders:   Andras Sali and Dean Attali (2020). shinycssloaders: Add CSS Loading Animations to 'shiny' Outputs. R package version 0.3. https://CRAN.R-project.org/package=shinycssloaders

R package stringr: Wickham H (2022). stringr: Simple, Consistent Wrappers for Common String Operations. R package version 1.4.1, https://CRAN.R-project.org/package=stringr.

R package Ciaro: Urbanek S, Horner J (2022). Cairo: R Graphics Device using Cairo Graphics Library for Creating High-Quality Bitmap (PNG, JPEG, TIFF), Vector (PDF, SVG, PostScript) and Display (X11 and Win32) Output. R package version 1.6-0, https://CRAN.R-project.org/package=Cairo.

R package svglite: Wickham H, Henry L, Pedersen T, Luciani T, Decorde M, Lise V (2022). svglite: An 'SVG' Graphics Device. R package version 2.1.0, https://CRAN.R-project.org/package=svglite.

R package RColorBrewer: Neuwirth E (2022). RColorBrewer: ColorBrewer Palettes. R package version 1.1-3, https://CRAN.R-project.org/package=RColorBrewer.

R package DT: Xie Y, Cheng J, Tan X (2022). DT: A Wrapper of the JavaScript Library 'DataTables'. R package version 0.25, https://CRAN.R-project.org/package=DT.

R package scales: Wickham H, Seidel D (2022). scales: Scale Functions for Visualization. R package version 1.2.1, https://CRAN.R-project.org/package=scales.

Python3.7 language: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

Python package miepython: Prahl, S. "miepython v1. 3.0." (2017).