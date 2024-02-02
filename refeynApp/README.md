# The PhotoMol app

Last time updated: March 2023

## Introduction

This folder contains a shiny app developed for analyzing mass photometry (MP) experiments. 
The fitting is done using a multi-gaussian model with a user defined number of gaussians.

The input data for PhotoMol is the distribution of events (counts versus masses) and the output data 
are the estimated peaks and widths of the gaussians. Contrasts can be converted to masses by loading a second MP experiment 
with known species. Example data is available when running the app.

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
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda numpy    
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda pandas    
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda scipy     
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda h5py     
/home/${user}/miniconda3/bin/conda  install   --freeze-installed -n r-reticulate -c anaconda xlrd        
```

6) Set the correct path for the app

``` bash 
if [ "$(basename "$(pwd)")" = "refeynApp" ]; then
    sed -i "s|base_dir <- paste0.*|base_dir <- paste0('$PWD', '/PhotoMol/')|" PhotoMol/global.R
else
    echo "Change the working directory to refeynApp"
fi
```

7) Create a folder to store temporary files:

``` R 
mkdir -p "/home/${user}/data_users/" 
```

8) Run the shiny app with R:

``` R 
cd PhotoMol
shiny::runApp()
```

## General code structure

The following tree represents how the files of PhotoMol are organised and their purpose:
```
refeynApp
   |-- PhotoMol
   |-- install_r_packages.R     # R packages required by the shiny app - to be installed inside the docker image
   |-- Rprofile.site            # R configuration file (used inside the docker image)
   |   |-- ui.R
   |   |-- server.R         # Server file
   |   |-- global.R         # Main UI file
   |   |-- refeyn.py        # Core script to analyse MP data. It works as a standalone command line programm too
   |   |-- helpers.py       # Helpers function to run refeyn.py
   |   |-- refeynCalibration.py # Core script to calibrate MP data. It works as a standalone command line programm too
   |   |-- reactives                        # Contains the required scripts to make the app interactive
   |   |   |-- download_reactives.R         # Scripts to download the data
   |   |   |-- reactives.R                  # Scripts to load, analyse and plot the data
   |   |   |-- reactivesCalibration.R       # Scripts to load, analyse and plot the data for calibration
   |   |-- server_files
   |   |   |-- helpers.R                    # helper functions used inside the server side
   |   |   |-- plot_functions.R             # plot functions, called by reactives.R and reactivesCalibration.R
   |   |-- test_reticulate.R                # Used in the dev process to test new functionalities
   |   |-- ui_files                                 # UI elements 
   |   |   |-- busy_indicator.R                     # function to show the user a nice ...loading... image
   |   |   |-- logo.R                               # Configuration for the app logo
   |   |   |-- theme.R                              # Theme of the app
   |   |   |-- ui_export_H5file_with_masses.R       # To export the calibrated h5 file
   |   |   |-- ui_export_fitting_information.R      # To export the fitted data
   |   |   |-- ui_export_plot_box.R                 # To customiset the distribution plot
   |   |   |-- ui_export_plot_box_calibration.R     # To customiset the distribution plot (calibration panel)
   |   |   |-- ui_export_plots_data.R               # To explort the data for the plots
   |   |   |-- ui_load_input_box.R                  # To load the input file
   |   |   |-- ui_load_input_box_calibration.R      # To load the input file (calibration panel)
   |   |   |-- ui_load_input_legend.R               # Suggested legend for the distribution plot
   |   |   |-- ui_simulate_box.R                    # To simulate a gaussian 
   |   |-- www                   # Contains images for the about section and example input files
   |   |-- docs
   |   |   |-- about.html        # HTML text for the about section 
   |   |   |-- user_guide.html   # HTML text for the user guide section 
```

## Future development 

It would be great to add support for loading multiple input files.

## References

Niebling, Stephan, et al. "Biophysical Screening Pipeline for Cryo-EM Grid Preparation of Membrane Proteins." Frontiers in Molecular Biosciences (2022): 535.

## Acknowledgments

Packages

R language: R Core Team (2020). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.

R package shiny:   Winston Chang, Joe Cheng, JJ Allaire, Yihui Xie and Jonathan McPherson (2020). shiny: Web Application Framework for R. R package version 1.4.0.2. https://CRAN.R-project.org/package=shiny

R package shinydashboard:   Winston Chang and Barbara Borges Ribeiro (2018). shinydashboard: Create Dashboards with 'Shiny'. R package version 0.7.1. https://CRAN.R-project.org/package=shinydashboard

R package ggplot2:   H. Wickham. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

R package reshape2:   Hadley Wickham (2007). Reshaping Data with the reshape Package. Journal of Statistical Software, 21(12), 1-20. URL http://www.jstatsoft.org/v21/i12/.

R package tippy:   John Coene (2018). tippy: Add Tooltips to 'R markdown' Documents or 'Shiny' Apps. R package version 0.0.1. https://CRAN.R-project.org/package=tippy

R package shinyalert:   Pretty Popup Messages (Modals) in 'Shiny'. R package version 1.1. https://CRAN.R-project.org/package=shinyalert

R package plotly:   C. Sievert. Interactive Web-Based Data Visualization with R, plotly, and shiny. Chapman and Hall/CRC Florida, 2020.

R package shinyjs:   Dean Attali (2020). shinyjs: Easily Improve the User Experience of Your Shiny Apps in Seconds. R package version 1.1. https://CRAN.R-project.org/package=shinyjs

R package reticulate:   Kevin Ushey, JJ Allaire and Yuan Tang (2020). reticulate: Interface to 'Python'. R package version 1.16. https://CRAN.R-project.org/package=reticulate

R package shinycssloaders:   Andras Sali and Dean Attali (2020). shinycssloaders: Add CSS Loading Animations to 'shiny' Outputs. R package version 0.3. https://CRAN.R-project.org/package=shinycssloaders

Python3.7 language: Van Rossum, G., & Drake, F. L. (2009). Python 3 Reference Manual. Scotts Valley, CA: CreateSpace.

Python package numpy: Travis E, Oliphant. A guide to NumPy, USA: Trelgol Publishing, (2006). Stéfan van der Walt, S. Chris Colbert, and Gaël Varoquaux. The NumPy Array: A Structure for Efficient Numerical Computation, Computing in Science & Engineering, 13, 22-30 (2011), DOI:10.1109/MCSE.2011.37

Python package pandas: Wes McKinney. Data Structures for Statistical Computing in Python, Proceedings of the 9th Python in Science Conference, 51-56 (2010)

Python package scipy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser, Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K. Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern, Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas, Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero, Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa, Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.

