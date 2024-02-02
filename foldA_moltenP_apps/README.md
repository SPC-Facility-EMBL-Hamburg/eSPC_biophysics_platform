# The FoldAffinity and MoltenProt apps

Last time updated: March 2023

## Introduction

This folder contains two shiny apps developed for analyzing differential scanning fluorimetry (DSF) data. The first app, FoldAffinity, estimates binding affinities using a two-state unfolding model coupled with ligand binding. It provides one model based on fitting the isotherms and one model based on fitting the observed melting temperatures. The second app, MoltenProt, provides a user-friendly and reliable method for fitting thermal unfolding curves to different models.

The input data for both apps consists of fluorescence-based melting curves, while the output data includes the estimated parameters. Example data can be found in the 'www' folders of each app. In the case of FoldAffinity, the dissociation constant Kd is of particular interest, whereas for MoltenProt, the melting temperature and enthalpy of unfolding are the key parameters. 

## Getting started

To run the apps locally you need to 

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
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c conda-forge pandas   
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda scipy     
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda xlrd      
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda openpyxl        
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda natsort   
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c conda-forge python-kaleido  
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c plotly plotly    
```

6) Set the correct path for both apps

``` bash 
if [ "$(basename "$(pwd)")" = "foldA_moltenP_apps" ]; then
   sed -i "s|base_dir <- paste0.*|base_dir <- paste0('$PWD', '/foldAffinity/')|" foldAffinity/server_files/global_variables.R
   sed -i "s|base_dir <- paste0.*|base_dir <- paste0('$PWD', '/moltenprot/')|" moltenprot/server_files/global_variables.R
else
   echo "Change the working directory to foldA_moltenP_apps"
fi
```
7) Create a folder to store temporary files:

``` R 
mkdir -p "/home/${user}/data_users/" 
```

8a) Run the shiny app with R (FoldAffinity):

``` R 
cd foldAffinity
shiny::runApp()
```

8b) Run the shiny app with R (MoltenProt):

``` R 
cd moltenprot
shiny::runApp()
```

## General code structure

The following tree represents how the files of FoldAffinity are organised and their purpose:
```
foldA_moltenP_apps
   |-- install_r_packages.R     # R packages required by the shiny app - to be installed inside the docker image
   |-- Rprofile.site            # R configuration file (used inside the docker image)
   |-- FoldAffinity                             
   |   |-- README.md
   |   |-- ui.R                                 # Main UI file
   |   |-- server.R                             # Server file
   |   |-- foldAffinity.py                      # Core script to analyse DSF data. It works as a standalone command line programm too
   |   |-- fitting_helpers_thermal_unfolding.py # Helpers to fit the fluorescence based melting curves
   |   |-- fitting_helpers_unfolded_fraction.py # Helpers to fit the unfolded fraction versus ligand concentration curves
   |   |-- helpers.py                           # Helpers to run foldAffinity.py
   |   |-- reactives                            # Contains the required scripts to make the app interactive
   |   |   |-- download_reactives.R             # Scripts to download the data
   |   |   |-- load_input_reactives.R           # Scripts to load the data
   |   |   |-- reactives.R                      # Scripts to process and plot the data
   |   |-- server_files
   |   |   |-- fitting_helpers.R                # fitting helpers for the unfolding model (of the observerd melting temperatures)
   |   |   |-- global_variables.R               # Global parameters to be loaded at the start of the app
   |   |   |-- helpers.R                        # helper functions used inside the server side
   |   |   |-- load_input_helpers.R             # helpers for loading the data
   |   |   |-- plot_functions.R                 # plotting functions, called by reactives.R
   |   |   |-- simulation_helpers.R             # helpers to simulate ligand induced thermal shift
   |   |-- test_reticulate.R                    # Used in the dev process to test new functionalities
   |   |-- ui_files                             # UI elements 
   |   |   |-- busy_indicator.R                 # function to show the user a nice ...loading... image
   |   |   |-- logo.R                           # Configuration for the app logo 
   |   |   |-- theme.R                          # theme of the app
   |   |   |-- ui_equations.R                   # equations shown in the simulation panel
   |   |   |-- ui_export_fitting_information.R  # to export the fitted parameters
   |   |   |-- ui_export_plots_data.R           # to export the the plots data
   |   |   |-- ui_fitting_box.R                 # to select the binding model and the temperature range of analysis
   |   |   |-- ui_fitting_options_box.R         # to fix the value of the heat capacity of change 
   |   |   |-- ui_fluo_dependence_on_temperature_box.R # parameters for the simulation panel
   |   |   |-- ui_fluorescence_fitting_tabbox.R        # Plot and table with the estimated parameters of the fitting of the unfolding curves
   |   |   |-- ui_kd_parameters_box.R                  # parameters for the simulation panel
   |   |   |-- ui_ku_parameters_box.R                  # parameters for the simulation panel
   |   |   |-- ui_launch_simulation.R                  # to run or not the simulation
   |   |   |-- ui_legend_fluorescence_fitting_tabbox.R # suggested legend for the fitted fluorescence curves
   |   |   |-- ui_legend_isothermal_fitting_tabbox.R   # suggested legend for the fitted isothermal curves
   |   |   |-- ui_legend_tab_box.R                     # suggested legend for the fitted isothermal curves
   |   |   |-- ui_load_input_box.R                     # suggested legend for the initial plots (raw curves and first derivative)
   |   |   |-- ui_model_selection_box.R           # to select the binding model for the unfolding curves (isothermal model)      
   |   |   |-- ui_plot_options_box.R              # to change the aspect of the plot -  raw data
   |   |   |-- ui_plot_options_box_isf.R          # to change the aspect of the plot -  isothermal fitting 
   |   |   |-- ui_plot_options_box_tm_shift.R     # to change the aspect of the plot -  observed Tms fitting 
   |   |   |-- ui_plot_options_fitting_box.R      # to download the fitted unfolding curves (isothermal model)
   |   |   |-- ui_position_vs_concentration_box.R # tables with the ligand concentration
   |   |   |-- ui_signal_tab_box.R                # plots with the raw data and 1st derivative
   |   |   |-- ui_tm_fit_equations.R              # equations for the Tm fit
   |   |   |-- ui_tm_model.R                      # to select the binding model for the fitting of the observed Tms 
   |   |-- www                                    # Contains images for the about section and example input files
   |   |-- docs
   |   |   |-- about.html        # HTML text for the about section                 
   |   |   |-- tutorial.html     # HTML text for the tutorial section 
   |   |   |-- user_guide.html   # HTML text for the user guide section 
```

The following tree represents how the files of MoltenProt are organised and their purpose:

```
foldA_moltenP_apps 
   |-- install_r_packages.R     # R packages required by the shiny app - to be installed inside the docker image
   |-- Rprofile.site            # R configuration file (used inside the docker image)                           
   |-- moltenprot
   |   |-- README.md
   |   |-- ui.R                         # Main UI file
   |   |-- server.R                     # Server file
   |   |-- moltenprot_shiny.py          # Core script to analyse DSF data. It works as a standalone command line programm too
   |   |-- helpers.py                   # Helpers to run moltenprot_shiny.py
   |   |-- reactives                    # Contains the required scripts to make the app interactive
   |   |   |-- download_reactives.R     # Scripts to download the data
   |   |   |-- reactives.R              # Scripts to load, process and plot the data
   |   |-- report_template              # folder with the report template 
   |   |   |-- header.tex               # latex configuration
   |   |   |-- report.Rmd               # report file - will be knited to a PDF
   |   |-- server_files
   |   |   |-- global_variables.R       # Global parameters to be loaded at the start of the app
   |   |   |-- helpers.R                # helper functions used inside the server side
   |   |   |-- load_input_helpers.R     # helpers for loading the data
   |   |   |-- plot_functions.R         # plotting functions, called by reactives.R
   |   |-- ui_files                     # UI elements 
   |   |   |-- busy_indicator.R         # function to show the user a nice ...loading... image
   |   |   |-- logo.R                   # Configuration for the app logo 
   |   |   |-- theme.R                  # theme of the app
   |   |   |-- menu_analyze
   |   |   |   |-- ui_filter_box.R                      # To filter the fitted curves
   |   |   |   |-- ui_plot_options_box.R                # To change the aspect of the plots
   |   |   |   |-- ui_score_table_resultPlot_tabbox.R   # Table with 'scores'
   |   |   |-- menu_export
   |   |   |   |-- ui_fitting_information.R             # To export the fitted data
   |   |   |   |-- ui_plots_data.R                      # To export the plots data
   |   |   |   |-- ui_report.R                          # To create a PDF report
   |   |   |-- menu_fit
   |   |   |   |-- ui_fitting_options_box.R             # To input parameters of the selected unfolding model
   |   |   |   |-- ui_legend_fitting_tab_box.R          # Example legend for the fitted curves
   |   |   |   |-- ui_model_selection_box.R             # To select the unfolding model
   |   |   |   |-- ui_plot_options_box.R                # To filter the fitted curves
   |   |   |   |-- ui_sort_fitted_params_box.R          # To sort the estimated parameters
   |   |   |   |-- ui_tabbox.R                          # To plot the estimated unfolding curves
   |   |   |-- menu_input
   |   |   |   |-- ui_conditions_table_box.R            # To load the sample names
   |   |   |   |-- ui_derivative_plot_tabbox.R          # To plot the first and second derivatives
   |   |   |   |-- ui_legend_tab_box.R                  # Example legend for the raw DSF curves (and derivative)
   |   |   |   |-- ui_load_input_box.R                  # To load the input data
   |   |   |   |-- ui_plot_options_box.R                # to change the aspect of the plot -  raw data
   |   |-- www
   |   |-- docs
   |   |   |-- about.html        # HTML text for the about section                 
   |   |   |-- tutorial.html     # HTML text for the tutorial section 
   |   |   |-- user_guide.html   # HTML text for the user guide section 
```

## Future development 

For MoltenProt, it would be great to add support for best-model selection, e.g., to automatically tell the user 
if they should use a two or three states model. For FoldAffinity, to add a new feature such that users can 
download a PDF report.

## References

Kotov, Vadim, et al. "In‐depth interrogation of protein thermal unfolding data with MoltenProt." Protein Science 30.1 (2021): 201-217.

Bai, Nan, et al. "Isothermal analysis of ThermoFluor data can readily provide quantitative binding affinities." Scientific reports 9.1 (2019): 1-15.

Niebling, Stephan, et al. "FoldAffinity: binding affinities from nDSF experiments." Scientific Reports 11.1 (2021): 9572.

Burastero, Osvaldo, et al. "eSPC: an online data-analysis platform for molecular biophysics." Acta Crystallographica Section D: Structural Biology 77.10 (2021): 1241-1250.

## Acknowledgments

Packages

FoldAffinity is possible thanks to: 

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

Python3.7 language: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

Python package xlrd: https://xlrd.readthedocs.io/en/latest/index.html

Python package natsort: https://natsort.readthedocs.io/en/master/

###########################################################################################
###########################################################################################

MoltenProt is possible thanks to: 

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

 Baptiste Auguie (2019). egg: Extensions for 'ggplot2': Custom Geom, Custom Themes, Plot Alignment, Labelled Panels, Symmetric Scales, and Fixed Panel Size. R package version 0.4.5. https://CRAN.R-project.org/package=egg

Python3.7 language: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

Python package xlrd: https://xlrd.readthedocs.io/en/latest/index.html

Python package natsort: https://natsort.readthedocs.io/en/master/

