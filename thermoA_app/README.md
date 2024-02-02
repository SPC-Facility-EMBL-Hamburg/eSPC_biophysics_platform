# The ThermoAffinity app

Last time updated: March 2023

## Introduction

This folder contains a shiny app developed for analyzing microscale thermophoresis (MST) experiments. 
The fitting is done using non-linear least squares and the estimated equilibrium dissociation constand Kd is reported
together with asymmetric confidence intervals.

The input data for ThermoAffinity is the thermophoresis curve and the output data 
are the estimated parameters (Kd and baselines). ThermoAffinity can also be used with any kind of data where the 
signal of the complex differs from the signal of the protein alone (e.g., fluorescence quenching). 
Example data is available when running the app.

## Getting started

To run the apps locally you need to 

1) Install R (tested with version 4.2.2)
2) Install Python (tested with version 3.10.9) 
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

5) Install the python packages numpy, pandas, xlrd and openpyxl in the conda environment 'r-reticulate':

```bash
/home/${user}/miniconda3/bin/conda  create    --name r-reticulate                                                      
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda nomkl           
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda numpy   
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c conda-forge pandas   
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda xlrd      
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda openpyxl        
```

6) Set the correct path for the app

``` bash 
if [ "$(basename "$(pwd)")" = "thermoA_app" ]; then
    sed -i "s|base_dir <- paste0.*|base_dir <- paste0('$PWD', '/thermoAffinity/')|" thermoAffinity/server_files/global_variables.R
else
    echo "Change the working directory to thermoA_app"
fi
```

7) Create a folder to store temporary files:

``` R 
mkdir -p "/home/${user}/data_users/" 
```

8) Run the shiny app with R:

``` R 
cd thermoAffinity
shiny::runApp()
```

## General code structure

```
thermoA_app
   |-- install_r_packages.R     # R packages required by the shiny app - to be installed inside the docker image
   |-- Rprofile.site            # R configuration file (used inside the docker image)
   |-- thermoAffinity
   |   |-- ui.R                                 # Main UI file
   |   |-- server.R                             # Server file
   |   |-- README.md
   |   |-- mst.py                               # Core script to analyse MST data. It works as a standalone command line programm too
   |   |-- helpers.py                           # Helpers function to run mst.py
   |   |-- reactives                            # Contains the required scripts to make the app interactive
   |   |   |-- download_reactives.R             # Scripts to download the data
   |   |   |-- load_input_reactives.R           # Scripts to load the data
   |   |   |-- reactives.R                      # Scripts to process and plot the data
   |   |-- server_files 
   |   |   |-- explore_parameters_helpers.R     # helpers function for the simulation panel 
   |   |   |-- fitting_helpers/                 # Functions to fit the binding curves
   |   |   |-- global_variables.R               # Global parameters to be loaded at the start of the app
   |   |   |-- helpers.R                        # helper functions used inside the server side
   |   |   |-- load_input_helpers.R             # helpers for loading the data
   |   |   |-- plot_functions.R                 # plotting functions, called by reactives.R
   |   |   |-- simulation_helpers.R             # helpers to simulate ligand binding curves
   |   |-- test_reticulate.R                    # Used in the dev process to test new functionalities
   |   |-- ui_files                             # UI elements 
   |   |   |-- busy_indicator.R                 # function to show the user a nice ...loading... image
   |   |   |-- logo.R                           # Configuration for the app logo 
   |   |   |-- theme.R                          # theme of the app
   |   |   |-- ui_advanced_settings.R                       # Bounds for the estimated parameters
   |   |   |-- ui_advanced_settings_explore_fhot.R          # To fit the Kd for different Fnorm intervals
   |   |   |-- ui_advanced_settings_simulate.R              # To set the parameters for the simulaton
   |   |   |-- ui_explore_param_simulate.R                  # To expand the simulated curves, e.g., simulate Kd*1.1 and Kd*0.9
   |   |   |-- ui_export_fitting_information.R              # To export the fitted data
   |   |   |-- ui_export_plots_data.R                       # To export the data for the plots
   |   |   |-- ui_fluorescence_fitting_tabbox_explore.R     # To show the result of fitting different Fnorm intervals
   |   |   |-- ui_fluorescence_fitting_tabbox_simple.R      # To show the result of the fitting (one Fnorm)
   |   |   |-- ui_load_input_box.R                          # To load the input data
   |   |   |-- ui_model_selection_box.R                     # To select the binding model
   |   |   |-- ui_model_selection_box_simulate.R            # To select the binding model for the simulation
   |   |   |-- ui_plot_fitting_legend_tab_box.R             # Suggested legend for the fitted binding curves
   |   |   |-- ui_plot_legend_tab_box.R                     # Suggested legend for the raw data
   |   |   |-- ui_plot_options_box.R                        # To customise the raw data plots
   |   |   |-- ui_plot_options_box_fit.R                    # To customise the fitted data plots
   |   |   |-- ui_plot_options_box_sim.R                    # To customise the simulation plots
   |   |   |-- ui_position_vs_concentration_box.R           # To set the ligand and protein concentration
   |   |   |-- ui_signal_tab_box.R
   |   |-- www                                    # Contains images for the about section and example input files
   |   |-- docs
   |   |   |-- about.html        # HTML text for the about section                 
   |   |   |-- tutorial.html     # HTML text for the tutorial section 
   |   |   |-- user_guide.html   # HTML text for the user guide section 
```

## Future development 

It would be great to add the possibility of exporting a PDF report.

## References

Burastero, Osvaldo, et al. "eSPC: an online data-analysis platform for molecular biophysics." Acta Crystallographica Section D: Structural Biology 77.10 (2021): 1241-1250.

## Acknowledgments

ThermoAffinity is possible thanks to: 

R language: R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

R package shiny:   Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2020). shiny: Web Application Framework for R. R package version 1.4.0.2. https://CRAN.R-project.org/package=shiny

R package viridis: Simon Garnier (2018). viridis: Default Color Maps from 'matplotlib'. R package version 0.5.1. https://CRAN.R-project.org/package=viridis

R package tidyverse: Wickham et al., (2019). Welcome to the tidyverse. Journal of Open Source Software, 4(43), 1686, https://doi.org/10.21105/joss.01686

R package pracma: Hans W. Borchers (2019). pracma: Practical Numerical Math Functions. R package version 2.2.9. https://CRAN.R-project.org/package=pracma

R package shinydashboard:   Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard

R package ggplot2:   H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

R package xlsx:   Adrian Dragulescu and Cole Arendt (2020). xlsx: Read, Write, Format Excel 2007 and Excel 97/2000/XP/2003 Files. R package version 0.6.3. https://CRAN.R-project.org/package=xlsx

R package reshape2:   Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

R package tippy:   John Coene (2018). tippy: Add Tooltips to 'R markdown' Documents or 'Shiny' Apps. R package version 0.0.1. https://CRAN.R-project.org/package=tippy

R package shinyalert:   Pretty Popup Messages (Modals) in 'Shiny'. R package version 1.1. https://CRAN.R-project.org/package=shinyalert

R package plotly:   C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.

R package tableHTML:   Theo Boutaris, Clemens Zauchner and Dana Jomar (2019). tableHTML: A Tool to Create HTML Tables. R package version 2.0.0. https://CRAN.R-project.org/package=tableHTML

R package rhandsontable:   Jonathan Owen (2018). rhandsontable: Interface to the 'Handsontable.js' Library. R package version 0.3.7. https://CRAN.R-project.org/package=rhandsontable

R package remotes:   Jim Hester, Gábor Csárdi, Hadley Wickham, Winston Chang, Martin Morgan and Dan Tenenbaum (2020). remotes: R Package Installation from Remote Repositories, Including 'GitHub'. R package version 2.1.1. https://CRAN.R-project.org/package=remotes

R package devtools:   Hadley Wickham, Jim Hester and Winston Chang (2020). devtools: Tools to Make Developing R Packages Easier. R package version 2.3.0. https://CRAN.R-project.org/package=devtools

R package shinyjs:   Dean Attali (2020). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 1.1. https://CRAN.R-project.org/package=shinyjs

R package data.table:   Matt Dowle and Arun Srinivasan (2019). data.table: Extension of data.frame. R package version 1.12.8. https://CRAN.R-project.org/package=data.table

R package reticulate:   Kevin Ushey, JJ Allaire and Yuan Tang (2020). reticulate: Interface to 'Python'. R package version 1.16. https://CRAN.R-project.org/package=reticulate

R package shinycssloaders:   Andras Sali and Dean Attali (2020). shinycssloaders: Add CSS Loading Animations to 'shiny' Outputs. R package version 0.3. https://CRAN.R-project.org/package=shinycssloaders

R package nlstools: Florent Baty, Christian Ritz, Sandrine Charles, Martin Brutsche, Jean-Pierre Flandrois, Marie-Laure Delignette-Muller (2015). A Toolbox for Nonlinear Regression in R: The Package nlstools. Journal of Statistical Software, 66(5), 1-21. URL  http://www.jstatsoft.org/v66/i05/

R package minpack.lm: Timur V. Elzhov, Katharine M. Mullen, Andrej-Nikolai Spiess and Ben Bolker (2016). minpack.lm: R Interface to the Levenberg-Marquardt Nonlinear Least-Squares Algorithm Found in MINPACK, Plus Support for Bounds. R package version 1.2-1.  https://CRAN.R-project.org/package=minpack.lm

R package broom: David Robinson, Alex Hayes and Simon Couch (2020). broom: Convert Statistical Objects into Tidy Tibbles. R package version 0.7.1. https://CRAN.R-project.org/package=broom

R pacakge data.table: Matt Dowle and Arun Srinivasan (2021). data.table: Extension of `data.frame`. R package version 1.14.2. https://CRAN.R-project.org/package=data.table

Python3.7 language: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

Python package xlrd: https://xlrd.readthedocs.io/en/latest/index.html

Python package natsort: https://natsort.readthedocs.io/en/master/
