SESCA: Semi-Empirical CD Spectrum Calculation Algorithm version 0.9.7

Installation guide:
Most of the scripts in this package are written as python scripts, and require an operating PYTHON environment.
SESCA and HBSS modules found in the 'scripts' and 'HBSS' subdirectories were tested and are compatible with Python 2.7
and Python 3.4 environments. DISICL scripts found in the DISICL subdirectory were compatible with Python 2.7 only,
but since version 0.9.5, they were also updated to function in Python 3.4 as well. Note that the deconvolution 
(SESCA_deconv.py) and scaling modules (SESCA_scale.py) of SESCA are dependent on the numeric python (numpy) and 
scientific python (scipy.optimize) packages, and may be disabled in the absence of those packages. Note that SESCA 
relies on previously published secondary structure analysis tools (DISICL and DSSP) to determine the secondary
structure of proteins from their three-dimensional structure, stored in a standard protein databank (PDB) file format.
SESCA can also process structural ensembles and simulation trajctories in PDB format, if the structural models are
presented between lines starting with the 'MODEL' and 'ENDMDL' expressions.

Automatic SESCA installation:
1) copy and unpack the zip file in your preferred directory (e.g.: '/home/usr/Programs/')
2) run setup.py in from the main SESCA directory 
3) to use SESCA, run modules from the scripts directory as stand-alone Python programs
4) if the modules do not work optimally see the instructions given for manual installation (below)

To leave the source files unmodified and install SESCA to a different directory, provide setup.py with a preffered
installation path using the @prefix flag (example: setup.py @prefix /home/usr/Programs/CDtools/SESCA).

If you are installing SESCA on windows-based operating systems, provide setup.py with the @win 1 arguments.
For UNIX-based operating systems, use @win 0 instead (default behaviour). Please note that system paths under Windows
-based operating systems start with a volume label (eg. "C:") and separated by "\" instead of "/". Because Python
denotes special characters also with a "\", system paths should be denoted using "\\" separators under Windows. 
(eg: "C:\\Programs\\SESCA")

use the arguments @help or --help to print usage for setup.py 


Manual SESCA installation:
1) copy and unpack the zip file in your preferred directory (e.g.: '/home/usr/Programs/')
2) enter the 'scripts' subdirectory, and open the SESCA_main.py file with your favorite text editor
3) modify the variable 'SESCA_dir' (in line 35) to point to the main SESCA directory (e.g.: SESCA_dir = "/home/usr/Programs/SESCA")
4) if you plan to use SESCA on a computer with a Windows operating system, change the value of the 'win' variable to 1 (line 34)
5) save the changes and close SESCA_main.py
6) repeat steps 3 to 5 for the files SESCA_bayes.py, SESCA_deconv.py, SESCA_scale.py, and SESCA_solver.py as well.

We included compatible versions of a pre-compiled DSSP program and ready-to-use DISICL and HBSS modules into this package. Please
follow the installation instructions of these algorithms if you wish to use them. Please cite the corresponding algorithms if you
used them to calculate CD spectra from their output.

To use SESCA under UNIX-based operating systems, you may need to provide permissions to be executed for the pyhton
scripts (files with the .py extension) in the 'scripts', 'DISICL', and 'HBSS' subdirectories (e.g.: chmod +rwx scripts/*.py),
as well as the the binary files (dssp-2.0.4-linux-amd64 and dssp-2.0.4-linux-i386) in the subdirectory 'DSSP' 

if you want to use the SESCA_dssp.py module as a stand-alone script, open it with you text editor, and change the 'win' ad 'dssp_main'
variables. The 'win' varible should be 0 on Unix based systems, and 1 on Windows systems. The 'dssp_main' variable should point to an
operational DSSP binary file. Note that both variables are overwritten by the ones specified in SESCA_main.py module if SESCA_dssp is 
called in automatic mode.

The output from DISICL is compatible with SESCA, and the output files can be used directly with SESCA_pred.py and
the appropriate basis sets (see libs directory) to generate the CD spectra.

The main script to calculate and compare CD spectra of proteins is SESCA_main.py in the "scripts" directory. It requires a basis set
(several basis sets can be found in the "libs" subdirectory), a protein structure file in PDB (protein data bank) format or a 
secondary structure summary file (from DSSP, DISICL or HBSS) compatible with the basis set, as well as a reference CD spectrum (optional). 
The standard output of SESCA and the basis sets are defined in mean residue ellipticity units (1000 degrees*cm^2/dmol), the wavelength 
information by default is in nanometers. 

Please take a look at the exercises provided in the "examples" directory on how to operate SESCA and compare the calculated
CD spectra with experimental CD measurements.

It is recommended to add the "scripts" directory of SESCA to your $PATH environmental variable or define aliases
for the scripts if you wish to use them. It is also recommended to add the absoulte paths to your 'scripts', 'DISICL' and 'HBSS'
subdirectories to the PYTHONPATH environmental varibale, if you wish to import them into your python scripts and programs.


License: 
Note that all files in the 'scripts', 'libs', 'examples' and 'HBSS' directories are considered as part of the SESCA package, 
and fall under the the GNU general public licese, version 3 (provded in the license.txt file in the main directory). 
For further details, see CD_calc.py in the scripts directory. 

This software is distributed in the hope that it will be useful but comes without any warranty. The author does not accept
responsibility to anyone for the consequences of using it, for whether it serves any particular purpose, or that it works at all.
No warranty is made about the software or its performance. 

If you found SESCA useful during your work please cite the following atricle(s):

Nagy, G.; Igaev, M.; Hoffmann, S. V.; Jones, N. C.; Grubmueller, H.;
SESCA: Predicting circular dichroism spectra from protein molecular structures;
Journal of Chemical Theory and Computation 15 (9); pp. 5087 - 5102 (2019);
doi: 10.1021/acs.jctc.9b00203 

Nagy G.; Grubmueller H.;
How accurate is circular dichroism-based  model  validation?;
European Biophysics Journal 49(6); pp. 497-510 (2020);
doi: 10.1007/s00249-020-01457-6

Nagy G.; Grubmueller H.;
Implementation of a Bayesian SecondaryStructure Estimation Method for the SESCA Circular Dichroism Analysis Package;
bioRxiv preprint doi: https://doi.org/10.1101/2020.12.02.408302
