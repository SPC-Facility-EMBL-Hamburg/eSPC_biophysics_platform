import sys
import os
import re

import numpy  as np
import pandas as pd

# Module to be imported from server.R !!!

#Step 1: import basic modules
Workdir = os.getcwd() # get current working directory, where server.R is located

#Step 2: import SESCA modules
#define SESCA_dir based on current directory:
SESCA_dir = os.path.join(Workdir,"./SESCA_v097/")
SESCA_dir = os.path.normpath(SESCA_dir)
Script_dir = os.path.join(SESCA_dir,"scripts")
#make sure that the SESCA scripts directory is in the system path:
sys.path.append(Script_dir)

#import SESCA_main.py and SESCA_pred.py from the scripts directory:
import SESCA_main  as Main
import SESCA_bayes as Bayes

#get all imports for SESCA main as well:
Imports = Main.Import_Custom()

mixed_basis_sets_names = [
    "DS-dTSC3",
    "DS5-4SC1",
    "DS6-1SC1",
    "DS5-6SC",
    "DSSP-TSC1",
    "DSSP-1SC3",
    "HBSS-3SC1"
]

"""
Description of the basis set
Name        Type		Number of SS    
DSSP-T		Dssp		3
DSSP-1		Dssp		4
DS3-1		DS_sim		5
DS6-1       DS_det      6    
DS5-4		DS_det		6
DS-dT		DS_det		3
HBSS-3		Hbss_ext	5
DSSP-F		Dssp		8
DS-simF		DS_sim		8 
Per-1		DS_det		X
Seq2-2		Seq		    X
DS-B4R1		DS_det		4   (no labels!)

#with side chain contributions:
DS-dTSC3	DS_det		3
DS5-4SC1	DS_det		6
DS6-1SC1	DS_det		6
DS5-6SC		DS_det		6
DSSP-TSC1	Dssp		6
DSSP-1SC3	Dssp		6
HBSS-3SC1	Hbss_ext	5
"""

def read_sesca_contributions_file(file_path):
    """
    Reads a SESCA output file and extracts wavelength and signal data.

    :param file_path: Path to the SESCA output file
    :return: Two lists containing wavelength and signal data
    """
    wavelength = []
    signal     = []

    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Find the start of the data section
        data_start = False
        for line in lines:
            if line.strip().startswith("#") or not line.strip():
                continue
            if not data_start:
                data_start = True

            if data_start:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        wl  = float(parts[0])
                        sig = float(parts[1])
                        wavelength.append(wl)
                        signal.append(sig)
                    except ValueError:
                        continue

    return np.array(wavelength), np.array(signal)


def extract_scaling_factor_and_ss_composition(file_path):
    """
    Extracts the average scaling factor, its standard deviation, and SS composition with their standard deviations from the given file.

    :param file_path: Path to the file
    :return: Dictionary containing the average scaling factor, its standard deviation, and SS composition with their standard deviations
    """
    avg_scaling_factor = None
    scaling_factor_sd = None
    ss_composition = {}
    ss_composition_sd = {}

    with open(file_path, 'r') as file:
        lines = file.readlines()
        ss_section = False

        for line in lines:
            # Extract average scaling factor and its standard deviation
            if "Avg. Scaling factor" in line:
                match = re.search(r'Avg\. Scaling factor\s*:\s*([\d.]+)\s*\+/-\s*([\d.]+)', line)
                if match:
                    avg_scaling_factor = float(match.group(1))
                    scaling_factor_sd = float(match.group(2))

            # Detect the start of the SS composition section
            if "Weighting factors for the most likely SS composition" in line:
                ss_section = True
                continue

            # Extract SS composition and their standard deviations
            if ss_section:
                match = re.search(r'(\w+)\s*:\s*([\d.]+)\s*\+/-\s*([\d.]+)', line)
                if match:
                    ss_composition[match.group(1)]    = float(match.group(2))
                    ss_composition_sd[match.group(1)] = float(match.group(3))
                else:
                    # End of SS composition section
                    ss_section = False

    return {
        "avg_scaling_factor": avg_scaling_factor,
        "scaling_factor_sd": scaling_factor_sd,
        "ss_composition": ss_composition,
        "ss_composition_sd": ss_composition_sd
    }

def format_bayes_results(results):
    """
    Formats the extracted Bayesian results into a pandas DataFrame for better readability.

    :param results: Dictionary containing the average scaling factor, its standard deviation, and SS composition with their standard deviations
    :return: Formatted pandas DataFrame
    """
    data = {
        "Parameter": ["Avg. Scaling Factor"] + list(results["ss_composition"].keys()),
        "Value (± SD)": [f"{results['avg_scaling_factor']} ± {results['scaling_factor_sd']}"] +
                 [f"{results['ss_composition'][key]} ± {results['ss_composition_sd'][key]}" for key in results["ss_composition"].keys()]
    }

    df = pd.DataFrame(data)

    return df


def extract_bin_data(file_path):

    """
    Extracts the SS bin data from the given file.
    :param file_path:
    :return: dataframe with the SS bin data. One column per SS element and one column for the probability of the bin
    """

    with open(file_path, 'r') as file:

        text_data = file.read()
        lines     = text_data.splitlines()

        for i,line in enumerate(lines):

            if '#Discretized SS dsitribution map:' in line:

                classes_start = i + 2

            if  '#bin:' in line and 'P(bin)' in line:

                classes_end = i

                break

        # Extract the SS classes
        ss_classes = [line.split()[-1] for line in lines[classes_start:classes_end]]

        bin_data = []

        read = True
        while read:

            i = i + 1

            info     = lines[i].strip().split()
            bin_data.append(np.array(info))

            bin_data.append(info)
            # stop reading if the next line is empty
            read = not len(lines[i+1].strip()) == 0

        columns = ss_classes + ['P(bin)']
        df      = pd.DataFrame(bin_data, columns=columns).astype(float)

        return df


class CdSpectraPredictor:

    def __init__(self):

        self.spectraNames = None
        self.CD_Spectra   = None
        self.SS_Comp      = None
        self.wavelength   = None

        self.comparison_stats = None

    def predict(self, files,names,basis,compare_ref=None):

        """
        Predict CD spectra from a list of PDB files and a given basis set
        :param files: list of PDB files
        :param names: list of PDB names
        :param basis: basis set to use
        :param compare_ref: reference file to compare the spectra with
        """

        self.comparison_stats = None

        if isinstance(files, str):
            files = [files]
            names = [names]

        try:

            spectra_names     = []
            CD_Spectra        = []
            SS_dfs            = []
            comparison_dfs    = []

            for cnt,file in enumerate(files):

                #define SESCA arguments:
                if compare_ref is None:
                    SESCA_args = " @pdb %1s @lib %1s @verb 0" % (file, basis)
                else:
                    SESCA_args = " @pdb %1s @ref %1s @lib %1s @verb 0" % (file, compare_ref, basis)

                if basis in mixed_basis_sets_names:

                    SESCA_args += " @write test.out"

                Processed_Args = Main.Read_Args(SESCA_args.split())
                #call SESCA_Main function to execute the program:
                Data = Main.SESCA_Main(Processed_Args)
                #Data array 0: SS composition, CD spectra

                SS_elements = [x[0]     for x in Data[0][0]]
                SS_values   = [x[1]*100 for x in Data[0][0]]

                # Create a DataFrame
                df        = pd.DataFrame(columns=SS_elements)
                df.loc[0] = SS_values
                df['pdb'] = names[cnt].replace(".pdb", "")

                SS_dfs.append(df)

                spectra_names.append(names[cnt])

                name = names[cnt].replace(".pdb", "")

                if compare_ref is not None:
                    df        = pd.DataFrame(columns=['RMSD (delta epsilon units)','MUSE (delta epsilon units)',
                                                      'Iref_norm_RMSD (no units)','dSS-est (percentage)','dSS-sd',
                                                      'dSS-conf.interval-up','dSS-conf.interval-low'])
                    stats     = np.array(Data[1][1:3]) / 3.298 # / 3.298 to convert to delta epsilon units
                    stats     = np.append(stats,Data[1][3])

                    stats     = np.concatenate((stats,np.array(Data[1][5][:4])))

                    df.loc[0] = stats
                    df['pdb'] = name

                    comparison_dfs.append(df)

                if basis in mixed_basis_sets_names:

                    _, bb_cont_signal  = read_sesca_contributions_file('BB_calc.out')
                    wl, sc_cont_signal = read_sesca_contributions_file('SC_calc.out')

                    spectra_names.append(name + ' BB contr.')
                    spectra_names.append(name + ' SC contr.')

                    CD_Spectra.append((bb_cont_signal + sc_cont_signal))
                    CD_Spectra.append(bb_cont_signal)
                    CD_Spectra.append(sc_cont_signal)
                    self.wavelength = wl

                else:

                    # Get the position with the predicted CD spectrum
                    pos = 2 if len(Data[0][1][0]) > 4 else 1

                    # Get the spectrum data from the SESCA output
                    # / 3.298 to convert to delta epsilon units, aka mean unit molar extinction
                    spectrum = np.array([x[pos] for x in Data[0][1]])
                    CD_Spectra.append(spectrum)
                    self.wavelength = np.array([x[0] for x in Data[0][1]])

            self.CD_Spectra   = np.column_stack(CD_Spectra) / 3.298 # / 3.298 to convert to delta epsilon units

            self.SS_Comp      = pd.concat(SS_dfs, axis=0).reset_index(drop=True)

            self.spectraNames = spectra_names

            if compare_ref is not None:
                self.comparison_stats = pd.concat(comparison_dfs, axis=0).reset_index(drop=True)

        except:

            pass

        # Remove files generated by SESCA
        os.system('rm -f *DISICL* HBSS*.out Hbss*.out *.dih *.log *.stat *.dssp BB_calc.out Dssp_0.out SC_calc.out Seq_*.out test.out' )

        return None

    def bayes_estimate(self,ref_file,basis,iterations,pdb_file_for_sc_corr=None):

        """
        Perform Bayesian estimation of the SS composition and scaling factor
        :param ref_file: space separated file with the reference CD spectra
        :param basis:   basis set to use
        :param iterations:  number of iterations to perform
        :param pdb_file_for_sc_corr: pdb file to extract side chain contributions (only the sequence information is used)
        """

        os.system('rm -f Bayes_est_1.out Bayes_map_1.out')

        if pdb_file_for_sc_corr is not None:

            SESCA_args = " @pdb %1s @lib %1s @verb 0 @write test.out" % (pdb_file_for_sc_corr, basis)

            # create the file 'SC_calc.out' with the SC contributions
            Processed_Args = Main.Read_Args(SESCA_args.split())
            # call SESCA_Main function to execute the program:
            Main.SESCA_Main(Processed_Args)

        SESCA_args = " @spect %1s @lib %1s @write Bayes_est_1.out @proj Bayes_map_1.out  @iter %1s  @verb 0 " % (ref_file, basis, iterations)

        if pdb_file_for_sc_corr is not None:

            SESCA_args += " @corr SC_calc.out" # include the correction for the SC contributions

        Processed_Args = Bayes.Read_Args(SESCA_args.split())
        Bayes.SSbayes_Main(Processed_Args)

        os.system('rm -f *DISICL* *.dih *.log *.stat *.dssp BB_calc.out Dssp_0.out SC_calc.out Seq_0.out test.out')

        return None
