import copy, sys, warnings

import pandas as pd

from loadCDfilesHelpers    import *
from helpers               import *
from cdUnitsConverter      import *
from decomposition_helpers import *

from fitting_helpers_thermal import *
from fitting_helpers_chemical import *
from fitting_helpers_thermochemical import *

sys.path.append('./secondary_structure_estimation_files')
from SelconsFunction import *

sys.path.append('./Dichroic-CD-model-main')
from single_fit_dichroic_estimator import *
from HC_thermal_fit                import *

from scipy.optimize  import minimize
"""
Classes for the analysis of circular dichroism data
Code written by Osvaldo Burastero 

If you use this script please use a reference similar to the following:

Osvaldo Burastero, Nykola C. Jones, Søren Vrønning Hoffmann, & Maria M. Garcia-Alai (2024). 
ChiraKit (Version 1.0). Manuscript in preparation. https://spc.embl-hamburg.de/app/chirakit

No warranty whatsoever
If you have questions please contact me:    oburastero@gmail.com
"""

# Matrix of protein secondary structure components (six components)
F1 = np.loadtxt("./secondary_structure_estimation_files/AU-F128_T-Nov22.txt", dtype='f', delimiter='\t')
# Matrix of known (reference) CD spectra - Proteins    
A1 = np.loadtxt("./secondary_structure_estimation_files/AU-A128_PCDDB-Nov22.txt", dtype='f', delimiter='\t')

# Load the req. data for the peptide helicity model
# Polynomials (partition function) for all helices (up to v**5) and double helices v**4
###################################
poly_total_all  = []
poly_double_all = []

# Iterate over different peptide lengths. The variable 'i' refers to the number of peptide bonds
for i_pep_bonds in range(4, 41):
    with open("./Dichroic-CD-model-main/Q_total/Q_total_" + str(i_pep_bonds) + ".txt", 'r') as l:
        poly_total = l.readlines()[0]
        poly_total_all.append(poly_total)

    with open("./Dichroic-CD-model-main/Q_double_H/Q_doubleH_" + str(i_pep_bonds) + ".txt", 'r') as l:
        poly_double = l.readlines()[0]
        poly_double_all.append(poly_double)

# F1 and A1 contain the AU-SP175 and AU-SMP180 reference datasets
# Labels of the known CD spectra
#Lbl = np.loadtxt("./secondary_structure_estimation_files/Labels-SMP180_PCDDBOrder.txt", dtype='str', delimiter='\t')

class CdExperimentGeneral:
    """
    Basic class for all CD spectra. It allows converting between different CD units
    """

    def __init__(self):

        self.wavelength = None  # 1D numpy array
        self.signalInput = None  # 2D dataframe with one column per measurement
        self.signalAbs = None  # 2D dataframe with one column per measurement, absorbance units
        self.signalDesiredUnit = None  # 2D dataframe with one column per measurement, signal in user selected units
        self.internalID = None  # Internal ID based on the experiment name and spectra names
        self.spectraNames = None  # Names of the measured spectra, e.g.,

        self.signalHT = None  # 2D dataframe with one column per measurement, unknown units

        # Metadata of the experiment, e.g, concentration, path length ... 
        self.metadata = None

        # Input units of the experiment 
        self.units = 'millidegrees'  # Default Input units

        # Working units, selected by the user
        # String, one of: 'milliabsorbance', 'absorbance', molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity', 
        #                 'meanUnitMolarExtinction' or 'meanUnitMolarEllipticity'
        self.desiredUnits = 'millidegrees'

        # To change back to previous working units
        self.currentUnits = None

        # Start with an impossible temperature (in degree Celsius), so the user has to change it.
        self.temperature = np.nan

        self.numberOfCroms = 0  # Integer. If protein sample, number of peptide bonds
        self.concentration = 0  # Float, in mg/ml
        self.pathLength = 0  # Float, in centimeters
        self.molecularWeight = 0  # Float, in Dalton

        # Boolean to decide if the user needs to modify the parameters numberOfCroms, concentration, pathLength and molecularWeight
        # if 'isFakeExperiment' is set to True, then the matrix signalDesiredUnit won't be changed once it is created for the first time 
        self.isFakeExperiment = False

        # String, to know how the 'fake experiment' matrix 'signalDesiredUnit' was created
        # For example, it should be either 'molarEllipticity',  'molarExtinction', 'meanUnitMolarExtinction' or 'meanUnitMolarEllipticity'
        self.fakeExperimentSignal = 'unknown'

        self.isGenerated = False  # Boolean to help exporting only the generated data

        # List of dataframes containing the secondary structure content, one df per CD spectrum
        self.secondary_structure_content = None

        # Name of the experiment
        self.name = None

        # For the helicity method. Mean residue ellipticity at 222 nm
        self.MRE_222nm = None

        # For the secondary structure method.
        self.secondary_structure_content = None
        self.fitted_spectra_sec_str = None
        self.query_spectra_sec_str = None
        self.lowerWL_sec_str = None
        self.F_sel = None
        self.A_sel = None
        self.SStruct_labels = None
        self.wl_even_sequence = None

    def experiment_to_absorbance_units(self, units_begin):

        """
        unitsBegin should be: 
            'absorbance' milliabsorbance', 'molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity',
            'meanUnitMolarExtinction' or 'meanUnitMolarEllipticity'
        """

        self.units = units_begin

        signal_in_abs_units = convert2absorbance(
            self.signalInput, units_begin, self.concentration,
            self.pathLength, self.molecularWeight, self.numberOfCroms)

        self.signalAbs = signal_in_abs_units

        return None

    def experiment_from_abs_to_other_units(self, units_end):

        """
        unitsEND should be: 
            'absorbance' milliabsorbance', 'molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity',
            'meanUnitMolarExtinction' or meanUnitMolarEllipticity

        """

        self.desiredUnits = units_end

        # Just copy the matrix in case the input and desired units are the same
        if self.desiredUnits == self.units:

            self.signalDesiredUnit = self.signalInput

        # Otherwise, compute the signal using the desired CD units
        else:

            signal_new = absorbance2desiredUnits(self.signalAbs,
                                                 units_end, self.concentration, self.pathLength, self.molecularWeight,
                                                 self.numberOfCroms)

            self.signalDesiredUnit = signal_new

        return None

    def load_data(self, file, name=''):

        file_type = detect_file_type(file)

        # Define the dictionaries to call the data loading function 

        read_data_function_dict = {

            "ChirascanFile": read_chirascan_file_data,
            "PCCDBFile": read_pccdb_file_data,
            "GenCDfile": read_gen_file_data,
            "DatFile": read_dat_file_data,
            "d0xFile": read_d0x_file_data,
            "plain_csv": read_custom_csv,
            "jasco_simple": read_jasco_single_sample_csv,
            "chirakit_txt_with_header": read_chirakit_txt_data,
            'ChirascanFileTemperatureRamp': read_chirascan_file_data_thermal_ramp
        }

        read_metadata_function_dict = {

            "ChirascanFile": read_chirascan_file_meta_data,
            "PCCDBFile": read_pccdb_file_meta_data,
            "GenCDfile": read_gen_file_meta_data,
            "DatFile": read_dat_file_meta_data,
            "d0xFile": read_d0x_file_meta_data,
            "plain_csv": read_chirakit_txt_meta_data,
            "jasco_simple": read_jasco_single_meta_data,
            "chirakit_txt_with_header": read_chirakit_txt_meta_data,
            "ChirascanFileTemperatureRamp": read_chirascan_file_meta_data
        }

        # Read the metadata
        self.metadata = read_metadata_function_dict[file_type](file)

        # If the file is of type Chirascan temperature ramp, assign spectra names based on temperature
        if file_type == 'ChirascanFileTemperatureRamp':

            self.wavelength, self.signalInput, self.temperature, self.signalHT = read_data_function_dict[file_type](
                file)

            self.internalID = [name + ' ' + str(round(t, 1)) for t in self.temperature]

        else:

            self.wavelength, self.signalInput, self.spectraNames, self.signalHT = read_data_function_dict[file_type](
                file)

            # If the file is of type d0x, load the temperature data directly from the file
            if file_type == 'd0xFile':

                self.temperature = read_temperature_d0x(file)

            # If the file is not of type d0x and not of type Chirascan temperature ramp, 
            # try to load the temperature data from the metadata 
            else:

                try:

                    for key, value in self.metadata.items():

                        if 'temp' in key.lower():

                            value = value.replace(",", " ")

                            digits = ''.join(ch for ch in value if ch.isdigit() or ch == '.')

                            # In case we have a numeric like string o less than 6 characters, we assume we have only one temperature value
                            if is_float(digits) and len(digits) < 6:
                                self.temperature = float(digits)
                                break

                            # Try to split based on spaces
                            if are_all_strings_numeric(value.split()):
                                self.temperature = np.array([float(x) for x in value.split()])
                                break

                except:

                    pass  # Keep default temperature

            # Assign internalID if we only have one spectrum

            if self.signalInput.shape[1] == 1:
                self.internalID = [name]

            else:

                # Expand the temperature float into a vector, if required. 
                if isinstance(self.temperature, float):
                    self.temperature = np.array([self.temperature for _ in self.spectraNames])

                # Assign internalID if we have many spectra
                self.internalID = [name + ' ' + x for x in self.spectraNames]

        self.spectraNames = self.internalID

        # Verify that the wavelength is numeric
        is_arr_numeric = np.issubdtype(self.wavelength.dtype, np.number)
        if not is_arr_numeric:
            raise TypeError("The wavelength array is not numeric. Please provide a numeric array.")

        self.units = guess_input_units_from_metadata_dictionary(self.metadata)

        self.numberOfCroms = guess_parameter_from_metadata_dictionary(self.metadata, ['chromophore', 'ncrom'])
        self.concentration = guess_parameter_from_metadata_dictionary(self.metadata, ['concentration', 'conc', 'mg'])
        self.pathLength = guess_parameter_from_metadata_dictionary(self.metadata, ['path length', 'path', ' length'])
        self.molecularWeight = guess_parameter_from_metadata_dictionary(self.metadata,
                                                                        ['molecular weight', 'mw', 'dalton'])

        self.name = name

        # Sort in increasing order
        # Get the indices that would sort the wavelength vector
        sorted_indices = np.argsort(self.wavelength)

        # Use the sorted_indices to rearrange the rows of the matrix
        self.signalInput = self.signalInput[sorted_indices]
        self.wavelength = self.wavelength[sorted_indices]
        self.signalHT = self.signalHT[sorted_indices]

        return None

    def init_and_check_helicity_method(self):

        # Verify that we have the protein concentration, path length, number of chromophore units, and weight (if needed)

        # To convert back the CD data to the current working units
        self.currentUnits = self.desiredUnits

        should_stop1 = any([
            self.numberOfCroms == 0, self.concentration == 0,
            self.pathLength == 0, self.molecularWeight == 0
        ])

        # Check that we actually need the parameters: protein concentration, path length, number of chromophores, and weight
        should_stop2 = self.units != 'meanUnitMolarEllipticity'

        should_stop3 = not np.any(np.logical_and(self.wavelength >= 221, self.wavelength <= 223))

        should_stop = (should_stop1 and should_stop2) or should_stop3 or self.numberOfCroms > 40 or self.numberOfCroms < 4

        return not should_stop

    def get_mre_222nm(self):

        self.experiment_from_abs_to_other_units('meanUnitMolarEllipticity')

        self.MRE_222nm = signal_at_222nm(self.signalDesiredUnit, self.wavelength)

        self.experiment_from_abs_to_other_units(self.currentUnits)

        return None

    def init_and_check_secondary_str_method(self, lower_wl):

        # Verify that we have the protein concentration, path length, number of chromophores, and weight (if needed)
        # Verify that the higher WL >= 240 nm and lower WL <= 190 nm

        self.secondary_structure_content = []
        self.fitted_spectra_sec_str = {}
        self.query_spectra_sec_str = {}
        self.lowerWL_sec_str = lower_wl

        self.currentUnits = self.desiredUnits

        should_stop1 = any([
            self.numberOfCroms == 0, self.concentration == 0,
            self.pathLength == 0, self.molecularWeight == 0
        ])

        # Check that we actually need the parameters: protein concentration, path length, number of chromophores, and weight
        should_stop2 = self.units != 'meanUnitMolarExtinction'

        should_stop3 = any([
            np.max(self.wavelength) < 240,
            np.min(self.wavelength) > 190
        ])

        should_stop = (should_stop1 and should_stop2) or should_stop3

        return not should_stop

    def set_secondary_structure_method_references_default(self):

        wavelength_temp = self.wavelength[self.wavelength >= self.lowerWL_sec_str]

        # Obtain the lower wavelength limit (can't be lower than 175) 
        lim_wl = int(np.max([np.min(wavelength_temp), 175]))

        # Select the most adequate spectra reference database
        if lim_wl >= 180:

            f = F1[:, :]
            self.A_sel = A1[:(240 - lim_wl + 1), :]  # At most, the first 61 WL (240-180), all 128 spectra (AU-SMP180)

            # Create a new F matrix with only four secondary structure components
            f_alpha = f[0, :] + f[1, :]
            f_beta = f[2, :] + f[3, :]
            self.F_sel = np.array([f_alpha, f_beta, f[4, :], f[5, :]])
            self.SStruct_labels = ['Alpha', 'Beta', 'Turns', 'Unord']

        else:

            self.F_sel = F1[:, :71]
            self.A_sel = A1[:(240 - lim_wl + 1), :71]  # At most, all 66 WL (240-175), first 71 spectra (AU-SP175)
            self.SStruct_labels = ['Alpha-r', 'Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Unord']

        # Wavelength even sequence with 1 nm step
        self.wl_even_sequence = np.arange(240, lim_wl - 1, -1)

        return None

    def set_secondary_structure_method_references_user(self, matrix_f, matrix_a, max_wl_ref, wl_ref_step,
                                                       ss_struct_labels):

        """
        Requires:

            - 'matrixF'         : numpy matrix, known secondary structure elements (each column should sum 1)
            - 'matrixA'         : numpy matrix, known CD spectra (in delta epsilon units)
            - 'max_wl_ref'      : float,        maximum wavelength of the known CD spectra
            - 'wl_ref_step'     : float,        wavelength step used to measure the known spectra
            - 'SStruct_labels'  : list,         list of strings with the names of the secondary structure elements

        Important detail:

            the matrix 'A' of the known CD spectra is assumed to have the wavelength in decreasing order.
            For example, the first row of 'A' has the CD signal at 240 nm
                     and the last  row of 'A' has the CD signal at 175 nm
        """

        wl_even_sequence = max_wl_ref - (np.arange(0, matrix_a.shape[0]) * wl_ref_step)

        cond1 = wl_even_sequence >= np.min(self.wavelength)
        cond2 = wl_even_sequence <= np.max(self.wavelength)
        cond3 = wl_even_sequence >= self.lowerWL_sec_str

        selected_idx = np.logical_and(cond1, np.logical_and(cond2, cond3))

        self.F_sel = matrix_f
        self.A_sel = matrix_a[selected_idx, :]
        self.wl_even_sequence = wl_even_sequence[selected_idx]
        self.SStruct_labels = ss_struct_labels

        return None

    def estimate_secondary_structure(self):

        """
        Run the Selcon3 algorithm to estimate the percentage of different secondary structure components
        """

        self.experiment_from_abs_to_other_units('meanUnitMolarExtinction')

        # Do secondary structure fitting
        signal_temp = self.signalDesiredUnit[self.wavelength >= self.lowerWL_sec_str]
        wavelength_temp = self.wavelength[self.wavelength >= self.lowerWL_sec_str]

        # Sort in increasing order, in case that the data  was not loaded with the self.load_data() function
        # Get the indices that would sort the wavelength vector
        sorted_indices = np.argsort(wavelength_temp)

        # Use the sorted_indices to rearrange the rows of the matrix and vector
        signal_temp = signal_temp[sorted_indices]
        wavelength_temp = wavelength_temp[sorted_indices]

        # Total CD spectra of this experiment
        n_spectra = self.signalDesiredUnit.shape[1]

        # Initalize list of empty lists
        # Each sublist will contain the results of the secondary structure fitting
        self.secondary_structure_content = [None for _ in range(n_spectra)]

        # Iterate over the columns 
        for j in range(n_spectra):

            # Interpolate 
            delta_epsilon = np.interp(self.wl_even_sequence, wavelength_temp, signal_temp[:, j], left=None, right=None,
                                      period=None)

            try:

                return_list, sec_str_df, mean_refit_prot, query_spectrum = SelconsPy(self.A_sel, self.F_sel,
                                                                                     delta_epsilon, self.SStruct_labels)

                self.secondary_structure_content[j] = sec_str_df
                self.fitted_spectra_sec_str[self.spectraNames[j]] = mean_refit_prot
                self.query_spectra_sec_str[self.spectraNames[j]] = query_spectrum

            except:

                pass

        # End of - Do secondary structure fitting

        self.experiment_from_abs_to_other_units(self.currentUnits)

        return None


class CdExperimentFittingModel(CdExperimentGeneral):
    """
    Advanced class for CD spectra with more measurements dimensions. For example, to load spectra at different temperatures or urea concentration
    This class is used as a starting point to build the class for thermal ramp analysis and for chemical denaturation analysis
    """

    def __init__(self):

        super().__init__()
        self.fit_params = None  # Values of the fitted parameters
        self.fit_rel_errors = None  # Relative error of the fitted parameters
        self.signal_predicted = None  # Predicted signal
        self.name = None  # String, experiment name
        self.decompositionDone = False
        self.last_col_name = 'WL (nm) / Basis spectrum'

        self.signal_useful, self.wavelength_useful = None, None

        self.oligo_conc_molar = None  # 1D numpy array, length 'n'
        self.oligo_conc_lst = None  # 1D list

        self.basis_spectra, self.coefficients = None, None
        self.minX, self.maxX = None, None

        self.basis_spectra_all, self.coefficients_all = None, None
        self.pca_based, self.explained_variance = None, None

        self.k = None  # Number of basis spectra
        self.fitted_spectra = None  # Reconstructed spectra from the basis spectra and coefficients
        self.wavelength_filtered = None  # Wavelengths that are useful for the analysis

        self.chem_concentration_ori, self.temperature_ori = None, None

        self.chem_concentration = None  # 1D numpy array, length 'n'. It could be for example, the concentration of Urea, Guanidim chloride or pH
        self.fractions = None
        self.bounds_df = None

        self.xAxis_lst_simple = None
        self.signal_lst2d = None  # List of signals, each element is a 2D numpy array
        self.xAxis_lst = None  # List of x values, each element is a 1D numpy array
        self.oligo_conc_lst = None  # List with the concentration of the n-mer equivalent
        self.signal_lst = None  # List of signals, each element is a 1D numpy array
        self.n = None  # Number of datasets to fit

        # Baselines and slopes for the folded and unfolded states
        self.bNs = None
        self.bUs = None
        self.kNs = None
        self.kUs = None
        self.bStart = None
        self.bEnd = None
        self.bIs = None

    def decompose_spectra_pca(self):

        try:

            x = self.signalDesiredUnit.T

            explained_variance, basis_spectra, coefficients = apply_pca(x)

            self.coefficients_all = coefficients
            self.basis_spectra_all = basis_spectra

            # Cumulated explained variance of the components
            self.explained_variance = explained_variance

            self.pca_based = True
            self.decompositionDone = True

        except:

            self.decompositionDone = False
            pass

        return None

    def decompose_spectra_svd(self):

        try:

            explained_variance, basis_spectra, coefficients = apply_svd(self.signalDesiredUnit)

            # Basis spectra and associated coefficients   
            self.basis_spectra_all = basis_spectra
            self.coefficients_all = coefficients

            # Cumulated explained variance of the components
            self.explained_variance = explained_variance

            self.pca_based = False
            self.decompositionDone = True

        except:

            self.decompositionDone = False
            pass

        return None

    def filter_basis_spectra(self, explained_variance_threshold=99):

        """
        Should be always run after decompose_spectra_svd or decompose_spectra_pca 
        so we create the matrices self.basis_spectra and self.coefficients
        based on the selected variance threshold
        """

        # Find the number of components (k) that capture at least threshold*100 percent of the variance or correlation
        k = np.sum(self.explained_variance < explained_variance_threshold) + 1

        self.k = k

        self.basis_spectra = self.basis_spectra_all[:, :k]
        self.coefficients = self.coefficients_all[:k, :]

        return None

    def align_basis_spectra_and_coefficients(self):

        #Try to align the n selected basis spectra and the associated coefficients to the observed peak

        self.basis_spectra, self.coefficients = align_basis_spectra_and_coefficients(self.signalDesiredUnit,
                                                                                     self.basis_spectra,
                                                                                     self.coefficients)

        return None

    def invert_selected_spectrum(self, n):

        n = int(n)

        if n <= self.basis_spectra.shape[1]:
            self.basis_spectra[:, n] = - self.basis_spectra[:, n]
            self.coefficients[n, :] = - self.coefficients[n, :]

        return None

    def rotate_basis_spectra(self):

        n_basis_spectra = self.basis_spectra.shape[1]

        if n_basis_spectra == 2:
            self.basis_spectra, self.coefficients = rotate_two_basis_spectra(self.signalDesiredUnit,
                                                                             self.basis_spectra_all, self.pca_based)

        if n_basis_spectra == 3:
            self.basis_spectra, self.coefficients = rotate_three_basis_spectra(self.signalDesiredUnit,
                                                                               self.basis_spectra_all, self.pca_based)

        # Compute again the explained variance after rotation!

        data = self.signalDesiredUnit

        if self.pca_based:
            data_mean = np.mean(data, axis=1, keepdims=True)
            data = data - data_mean

        total_variance = np.linalg.norm(data) ** 2  # Total variance in the data

        self.explained_variance = explained_variance_from_orthogonal_vectors(self.basis_spectra, self.coefficients,
                                                                             total_variance)

        return None

    def reconstruct_spectra(self):

        self.fitted_spectra = reconstruct_spectra(self.basis_spectra, self.coefficients, self.signalDesiredUnit,
                                                  self.pca_based)

        return None

    def assign_useful_signal_svd(self, relevant_k=1):
        """
        Trick to be able to fit the svd/pca coefficients
        This method will add the properties
            signal_useful and wavelength_useful

            signal_useful       contains the change in the coefficients
            wavelength_useful   contains the number of relevant coefficients
        """

        self.signal_useful = self.coefficients[0:relevant_k, :]
        self.wavelength_useful = np.arange(1, relevant_k + 1)

        return None

    def assign_useful_signal(self, user_selected_wavelength=None):

        # Default method, use the filtered wavelength
        if user_selected_wavelength is None:

            self.signal_useful     = self.signalDesiredUnit[np.isin(self.wavelength, self.wavelength_filtered), :]
            self.wavelength_useful = self.wavelength_filtered

        # Alternative, subset the CD data based on user given wavelengths    
        else:

            # Conver to list if required
            if not isinstance(user_selected_wavelength, list):
                user_selected_wavelength = [user_selected_wavelength]

            # It may happen that the user selected wavelengths were not measured
            # Therefore, we'll find the closest ones

            available_wavelength = []

            for wl in user_selected_wavelength:
                # Calculate absolute differences
                abs_diff = np.abs(self.wavelength - wl)

                # Find index of closest value
                closest_index = np.argmin(abs_diff)

                # Get the closest value
                closest_value = self.wavelength[closest_index]

                available_wavelength.append(closest_value)

            available_wavelength = np.array(available_wavelength)

            self.signal_useful = self.signalDesiredUnit[np.isin(self.wavelength, available_wavelength), :]
            # Remove wavelengths that are not in the data
            self.wavelength_useful = available_wavelength[np.isin(available_wavelength, self.wavelength)]

        return None

    def estimate_useful_signal_based_on_snr_and_amplitude(self, measurement_factor):

        if len(self.wavelength) < 5:
            self.wavelength_filtered = self.wavelength
            return None

        # 1) Apply the noise filtering algorithm
        wavelength_filtered = estimate_useful_wavelength_based_on_snr(self.wavelength, self.signalDesiredUnit)

        signal_useful = self.signalDesiredUnit[np.isin(self.wavelength, wavelength_filtered), :]

        # 2) Apply the amplitude filtering algorithm (on the already filtered data!)
        wavelength_filtered = estimate_useful_wavelength_based_on_amplitude(wavelength_filtered, signal_useful,
                                                                            measurement_factor)

        self.wavelength_filtered = wavelength_filtered

        return None

    def reshape_signal_oligomer(self, mode='Chemical'):

        """
        This is a method to treat differently the data collected at the same wavelength,
        but different protein concentration
        We need to have already self.signal_useful assigned
        """

        x_values = self.chem_concentration_ori if mode == 'Chemical' else self.temperature_ori

        self.minX, self.maxX = np.min(x_values), np.max(x_values)

        if self.oligo_conc_molar is None:

            na_masks = np.apply_along_axis(np.isnan, axis=1, arr=self.signal_useful)

            self.signal_lst = [self.signal_useful[i,:][~na_masks[i,:]] for i in range(self.signal_useful.shape[0])]
            self.xAxis_lst  = [x_values[~na_masks[i,:]] for i in range(self.signal_useful.shape[0])]

            self.n = len(self.xAxis_lst)

            self.oligo_conc_lst = [None for _ in self.xAxis_lst]
            return None

        unique_dimer_concs = np.unique(self.oligo_conc_molar)

        signal_lst, wavelength_lst, x_axis_lst, oligo_conc_lst = [], [], [], []

        signal_lst2d, x_axis_lst_simple, x_val_reorder = [], [], []

        for unc in unique_dimer_concs:
            selected = self.oligo_conc_molar == unc
            selected_signal = self.signal_useful[:, selected]

            wavelength_lst.append([str(wl) + ' [' + str(round(unc * 1e6, 2)) + ' μM]' for wl in self.wavelength_useful])

            signal_lst += selected_signal.tolist()

            x_vals = x_values[selected]

            x_val_reorder.append(x_vals)

            x_axis_lst += [x_vals for _ in self.wavelength_useful]
            oligo_conc_lst += [unc for _ in self.wavelength_useful]

            signal_lst2d.append(selected_signal)
            x_axis_lst_simple.append(x_vals)

        self.wavelength_useful = np.concatenate(wavelength_lst)

        are_all_same_length = len(np.unique([len(x) for x in x_axis_lst])) == 1

        if are_all_same_length:

            if all(np.allclose(x, x_axis_lst[0], atol=1e-3) for x in x_axis_lst):
                self.signal_useful = np.concatenate(signal_lst2d, axis=0)
                x_values = x_axis_lst[0]

        else:

            x_values = np.concatenate(x_val_reorder)
            self.signal_useful = signal_2d_matrices_to_full_matrix(signal_lst2d, len(self.wavelength_useful),
                                                                   len(self.oligo_conc_molar))

        if mode == 'Chemical':
            self.chem_concentration = x_values
        else:
            self.temperature = x_values

        self.xAxis_lst_simple = x_axis_lst_simple
        self.signal_lst2d = signal_lst2d
        self.xAxis_lst = x_axis_lst
        self.oligo_conc_lst = oligo_conc_lst
        self.signal_lst = signal_lst
        self.n = len(self.xAxis_lst)

        return None

    def baselines(self, model, fit_slope_native, fit_slope_unfolded, mode='Chemical'):

        """
        To treat same wavelength data with different protein concentration differently
        We need to have already self.signal_useful assigned
        """

        threshold = 2 if mode == 'Chemical' else 15

        stc = model_n(model)

        if model == 'Monomer':

            x_values = self.chem_concentration if mode == 'Chemical' else temperature_to_kelvin(self.temperature)

            self.bNs, self.kNs, self.bUs, self.kUs = fit_baselines(self.signal_useful, x_values, threshold)

            argmin_x  = [np.argmin(x) for x in self.xAxis_lst]
            argmax_x  = [np.argmax(x) for x in self.xAxis_lst]

            self.bStart = np.array([y[idx] for idx,y in zip(argmin_x,self.signal_lst)])
            self.bEnd   = np.array([y[idx] for idx,y in zip(argmax_x,self.signal_lst)])

        else:

            bns_all, kns_all, bus_all, kus_all = [], [], [], []
            b_start_all, b_end_all = [], []

            for selected_signal, x_values in zip(self.signal_lst2d, self.xAxis_lst_simple):

                if mode == 'Thermal': x_values = temperature_to_kelvin(x_values)

                bNs, kNs, bUs, kUs = fit_baselines(selected_signal, x_values, threshold)

                bStart = selected_signal[:, np.argmin(x_values)]
                bEnd   = selected_signal[:, np.argmax(x_values)] / stc

                b_start_all.append(bStart)
                b_end_all.append(bEnd)

                bns_all.append(bNs)
                bus_all.append(bUs)
                kns_all.append(kNs)
                kus_all.append(kUs)

            self.bNs = np.concatenate(bns_all)
            self.bUs = np.concatenate(bus_all)
            self.kNs = np.concatenate(kns_all)
            self.kUs = np.concatenate(kus_all)
            self.bStart = np.concatenate(b_start_all)
            self.bEnd = np.concatenate(b_end_all)

        # Initial signal value for the intermediate state
        self.bIs = (self.bStart + self.bEnd) / 2
        # Readjust intercept
        if not fit_slope_native:   self.bNs = self.bStart
        if not fit_slope_unfolded: self.bUs = self.bEnd

        return None


class CdExperimentComparison(CdExperimentFittingModel):
    """
    Advanced class to handle CD spectra with categorical labels. For example, 'control' and 'treatment'.

    We inherit the class 'cd_experiment_fitting_model' to take advantage of the PCA functions (for future developments)

    Some lines are commented with the pattern 'HW_Method'. I implemented the method 
    proposed by Hristova, Kalina, and William C. Wimley. 
    "Determining the statistical significance of the difference between arbitrary curves: A spreadsheet method." 
    Plos one 18.10 (2023): e0289619.

    The code works, but the results were not satisfactory so it is not used in production.
    To try it, uncomment the lines with that pattern
    """

    def __init__(self):
        super().__init__()
        self.signalDesiredUnitOri = None  # copy of self.signalDesiredUnit, to undo normalisation changes

        self.labels = None  # 1D (length m)   np array, labels of the CD signal matrix - one label per column
        self.labelsOri = None  # copy of self.labels, to undo normalisation changes

        self.workingUnits = None  # string, selected working units to create the dataset for comparisons

        self.labels_unique = None  # 1D (length z)   np array, unique labels
        self.labels_unique_N = None  # 1D (length z)   np array, number of curves per label

        self.means = None  # 2D (size n x z) np array, CD signal average            (per wavelength) of each category
        self.sds = None  # 2D (size n x z) np array, CD signal standard deviation (per wavelength) of each category
        #HW_Method self.standard_err          = None  # 2D (size n x z) np array, CD signal standard error     (per wavelength) of each category

        self.distance_matrix = None  # 2D (size m x m) np array, containing the all versus all comparisons
        self.distances = None  # list of lists containing the all versus all comparisons. One sublist per comparison
        self.comparison_labels = None  # list containing the labels of the comparisons.   One element per comparison

        # Let 'h' be the total number of comparisons to be done,
        self.difference_spectra = None  # 2D (size n x h) np array,
        self.difference_spectra_sd = None  # 2D (size n x h) np array,
        self.difference_spectra_lbl = None  # length h, np array, labels of the difference spectra

        self.lengths = None  # 1D (length m) np array, norms of the spectra

    '''
    def normalise_by_area(self):

        Divide each CD spectrum using the area under the curve (in absolute value, to keep the sign of the peaks)
        The area is calculated with the simpson rule

        self.integrals         = np.abs(simpson(self.signalDesiredUnitOri, x=self.wavelength,axis=0))
        self.signalDesiredUnit = self.signalDesiredUnitOri / self.integrals
        self.labels            = np.array(['Area norm. ' + str(x) for x in self.labelsOri])

        return None
    '''

    def normalise_by_L2norm(self):

        """
        Divide each CD spectrum using the length of the vector (in absolute value, to keep the sign of the peaks)
        The area is calculated with the simpson rule
        """

        self.lengths = np.linalg.norm(self.signalDesiredUnitOri, axis=0)
        self.signalDesiredUnit = self.signalDesiredUnitOri / self.lengths
        self.labels = np.array(['Norm. ' + str(x) for x in self.labelsOri])

        return None

    def undo_normalise(self):

        self.signalDesiredUnit = self.signalDesiredUnitOri
        self.labels = self.labelsOri

        return None

    def summarise_signal_per_label(self):

        """
        Obtain the average and standard deviation per unique label
        """

        u, ind = np.unique(self.labels, return_index=True)

        self.labels_unique = np.array(u[np.argsort(ind)])

        labels_unique_N = []
        means_by_label = []
        sds_by_label = []
        #HW_Method ses_by_label    = []

        for label in self.labels_unique:

            sel_idx = self.labels == label

            mean_by_wl = np.mean(self.signalDesiredUnit[:, sel_idx], axis=1)

            labels_unique_N.append(sum(sel_idx))

            if np.sum(sel_idx) > 1:

                sd_by_wl = np.std(self.signalDesiredUnit[:, sel_idx], axis=1,
                                  ddof=1)  # ddof=1 to use (N - 1) in the denominator
                #HW_Method se_by_wl   = sd_by_wl / np.sqrt(len(sel_idx))

            else:

                sd_by_wl = np.full(self.signalDesiredUnit.shape[0], np.nan)  # maybe np.nan would be better ...
                #HW_Method se_by_wl   = np.full(self.signalDesiredUnit.shape[0], np.nan) # maybe np.nan would be better ...

            means_by_label.append(mean_by_wl)
            sds_by_label.append(sd_by_wl)
            #ses_by_label.append(se_by_wl)

        self.means = np.array(means_by_label).T
        self.sds = np.array(sds_by_label).T
        #HW_Method self.standard_err    = np.array(ses_by_label).T
        self.labels_unique_N = np.array(labels_unique_N)

        return None

    def generate_comparison_labels(self):

        """
        Create the comparison labels, for example, if we have the unique labels 'Treatment' and 'Control'
        Then, the 'comparison_labels' will be 'Treatment versus Treatment', 'Treatment versus Control' and 'Control versus Control'
        """

        comparison_labels = []

        for j, label1 in enumerate(self.labels_unique):

            for label2 in self.labels_unique[j:]:
                comparison_labels.append(label1 + ' versus ' + label2)

        self.comparison_labels = np.array(comparison_labels)

        return None

    def generate_difference_spectra(self):

        """
        Generate all possible 'difference' spectra. For instance, if we have the unique labels
        'Control' and 'Treatment', we will create the spectrum 'Control - Treatment'.   
        """

        difference_spectra = []
        difference_spectra_sd = []
        difference_spectra_lbl = []

        # Uncomment the next lines to use test the method proposed by Hristova, Kalina, and William C. Wimley (2023). Plos one
        # Uncomment 
        #HW_Method difference_spectra_se     = []
        #HW_Method difference_spectra_nsmall = []
        #HW_Method difference_dfreedom       = []

        for i, label1 in enumerate(self.labels_unique[:-1]):

            for ii, label2 in enumerate(self.labels_unique[(i + 1):]):
                difference_spectra.append(self.means[:, i] - self.means[:, (ii + i + 1)])
                difference_spectra_sd.append(np.sqrt(self.sds[:, i] ** 2 + self.sds[:, (ii + i + 1)] ** 2))

                difference_spectra_lbl.append(label1 + ' - ' + label2)

                # Uncomment the next lines to use test the method proposed by Hristova, Kalina, and William C. Wimley (2023). Plos one
                #HW_Method difference_spectra_se.append(np.sqrt(self.standard_err[:,i]**2 + self.standard_err[:,(ii+i+1)]**2))
                #HW_Method nsmall = np.min([self.labels_unique_N[i],self.labels_unique_N[i+ii+1]])
                #HW_Method df     = nsmall*2 - 2

                #HW_Method difference_spectra_nsmall.append(nsmall)
                #HW_Method difference_dfreedom.append(df)

        self.difference_spectra = np.array(difference_spectra).T
        self.difference_spectra_sd = np.array(difference_spectra_sd).T
        self.difference_spectra_lbl = np.array(difference_spectra_lbl)

        return None

    def find_distances(self):

        """
        For each spectrum, compute the normalised euclidean distance against all other spectra
        
        As a result of running this method, we'll add to the class 

            self.distances  -   list of numpy arrays. Each numpy array contains the computed distances of one category. 
                                For instance, if the unique labels are 'Treatment' and 'Control', 
                                self.distances will have length 3, 
                                one array will have the distances of all the 'Treatment' spectra against the 'Treatment' spectra,   
                                one array will have the distances of all the 'Treatment' spectra against the 'Control'   spectra, and
                                one array will have the distances of all the 'Control'   spectra against the 'Control'   spectra.

            self.distance_matrix - matrix of the euclidean distance for all the spectra
        """

        n = len(self.labels)
        distance_matrix = np.full((n, n), np.nan)

        n_columns = self.signalDesiredUnit.shape[1]

        distance_matrix[n_columns - 1, n_columns - 1] = 0  # Fill the last element of the diagonal
        distances = [[] for _ in self.comparison_labels]

        comparison_labels_lst = self.comparison_labels.tolist()

        for i in range(n_columns - 1):

            distance_matrix[i, i] = 0  # Fill the diagonal

            label1 = self.labels[i]

            for ii in range(i + 1, n_columns):
                label2 = self.labels[ii]
                match_id = comparison_labels_lst.index(label1 + ' versus ' + label2)
                distance = np.linalg.norm(self.signalDesiredUnit[:, i] - self.signalDesiredUnit[:, ii])
                distances[match_id].append(distance)
                distance_matrix[i, ii] = distance

        distances = [np.array(x) for x in distances]
        self.distances = distances
        self.distance_matrix = distance_matrix

        return None


class CdExperimentChemicalUnfolding(CdExperimentFittingModel):
    """
    Advanced class to analyse chemical unfolding
    """

    def __init__(self):

        super().__init__()

    def load_unfolding_data_monomer(self, file_path):

        """
        Assign self.wavelength_useful, self.signal_useful directly from file
        """

        try:

            chem_concentration, signal_useful, wavelength_useful = read_unfolding_data_monomer(file_path)

            self.chem_concentration_ori, self.chem_concentration = chem_concentration, chem_concentration

            self.signal_useful     = signal_useful
            self.wavelength_useful = wavelength_useful

            self.last_col_name = 'Curve'
            self.temperature   = read_unfolding_file_temperature(file_path)

            return True

        except:

            return False

    def load_unfolding_data_oligomer(self, file_path):

        """
        Assign self.wavelength_useful, self.signal_useful directly from file
        """

        try:

            chem_concentration, signal, concentrations, self.wavelength_useful = read_unfolding_data_oligomer(file_path)

            self.chem_concentration_ori, self.chem_concentration = chem_concentration, chem_concentration

            if len(concentrations) == len(chem_concentration):
                self.oligo_conc_molar = concentrations

            self.signal_useful = signal
            self.last_col_name = 'N-mer concentration'

            self.temperature = read_unfolding_file_temperature(file_path)

            return True

        except:

            return False

    def fit_signal(self, fit_slope_native=True, fit_slope_unfolded=True, model='Monomer'):

        self.baselines(model, fit_slope_native, fit_slope_unfolded, 'Chemical')

        signal_fx    = map_two_state_chem_model_to_signal_fx(model)
        fractions_fx = map_two_state_chem_model_to_fractions_fx(model)

        # Initial parameters have to be in order: 
        # Global M, Global D50
        # Single intercepts folded, Single intercepts unfolded
        # Single slopes folded    , Single slopes unfolded

        # Set initial parameters
        p0 = np.concatenate(((1, np.median(self.chem_concentration)), self.bNs, self.bUs, self.kNs, self.kUs))

        d50_min = 0.05
        d50_max = self.maxX if model == 'Monomer' else self.maxX + 5

        if model == 'Trimer':
            d50_min, d50_max, p0[0] = d50_min + 1, d50_max + 3, p0[0] + 2

        if model == 'Tetramer':
            d50_min, d50_max, p0[0] = d50_min + 2, d50_max + 5, p0[0] + 3

        low_bound  = np.array([0.1, d50_min] + [x / 40 if x > 0 else 40 * x for x in p0[2:]])
        high_bound = np.array([20, d50_max]  + [40 * x if x > 0 else x / 40 for x in p0[2:]])

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 2, 2, self.n)

        global_fit_params, cov = fit_chemical_unfolding(self.xAxis_lst, self.signal_lst, self.temperature,
                                                        p0, low_bound, high_bound, fit_slope_native, fit_slope_unfolded,
                                                        signal_fx, self.oligo_conc_lst)

        ## Try to re fit the slopes and/or baselines if the fitted parameters are close to the boundaries
        re_fit, low_bound, high_bound = extend_bounds(global_fit_params, 2, low_bound, high_bound, 1, 250)

        if re_fit:
            global_fit_params, cov = fit_chemical_unfolding(self.xAxis_lst, self.signal_lst, self.temperature,
                                                            global_fit_params, low_bound, high_bound, fit_slope_native,
                                                            fit_slope_unfolded, signal_fx, self.oligo_conc_lst)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / global_fit_params * 100)

        params_splt = split_all_params_two_state(global_fit_params, self.n, fit_slope_native, fit_slope_unfolded)
        errors_splt = split_all_params_two_state(errors, self.n, fit_slope_native, fit_slope_unfolded)

        m, d50, blFold, slopeFold, blUnfold, slopeUnfold = params_splt
        MErr, D50Err, blFoldErr, slopeFoldErr, blUnfoldErr, slopeUnfoldErr = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        for i in range(self.n):
            bN, kN, bU, kU = blFold[i], slopeFold[i], blUnfold[i], slopeUnfold[i]
            bNe, kNe, bUe, kUe = blFoldErr[i], slopeFoldErr[i], blUnfoldErr[i], slopeUnfoldErr[i]

            Y = signal_fx(self.chem_concentration, self.temperature, d50, m, bN, kN, bU, kU, self.oligo_conc_lst[i])

            predicted.append(Y)

            fit_params.append([kN, bN, kU, bU, d50, m, str(self.wavelength_useful[i]), self.name])
            fit_errors.append([kNe, bNe, kUe, bUe, D50Err, MErr, str(self.wavelength_useful[i]), self.name])

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'D50', 'M', self.last_col_name, 'Dataset']

        # Convert the NumPy array to a Pandas DataFrame with column names
        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = generate_predicted_matrix(predicted, self.signal_useful)

        if model == 'Monomer':
            self.fractions = fractions_fx(self.temperature, self.chem_concentration, m, d50)

        else:

            fractions_lst = []
            unique_oligo_concs = np.unique(self.oligo_conc_lst)

            for unc in unique_oligo_concs:
                fractions_lst.append(fractions_fx(self.temperature, self.chem_concentration, m, d50, unc))

            self.fractions = merge_fractions_dictionaries(fractions_lst, unique_oligo_concs)

        # In the same order as given to the fitting function
        shared_params_names, local_params_names = ['M', 'D50'], ['bN', 'bU']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None

    def fit_signal_three_state(self, fit_slope_native=True, fit_slope_unfolded=True,
                               d50v1_init=0, d50v2_init=0, model='Monomer'):

        self.baselines(model, fit_slope_native, fit_slope_unfolded, 'Chemical')

        signal_fx = map_three_state_chem_model_to_signal_fx(model)
        fractions_fx = map_three_state_chem_model_to_fractions_fx(model)

        # Initial parameters have to be in order: 
        # Global M, Global D50
        # Single intercepts folded, Single intercepts unfolded
        # Single slopes folded    , Single slopes unfolded

        # Set initial parameters
        p0 = np.concatenate(
            ((0.25, self.minX + 0.5, 0.25, self.maxX - 0.5), self.bNs, self.bUs, self.bIs, self.kNs, self.kUs))

        D50min = 0.02
        D50max = self.maxX - 0.5 if model == 'Monomer' else self.maxX + 1

        increments = {"dimer": 2, "trimer": 3, "tetramer": 4}

        for key, value in increments.items():
            if key in model.lower():
                D50max += value

        low_bound = np.array([0.2, D50min, 0.2, D50min + 1.5] + [x / 30 if x > 0 else 30 * x for x in p0[4:]])
        high_bound = np.array([20, D50max, 20, D50max] + [30 * x if x > 0 else x / 30 for x in p0[4:]])

        # If required, change the initial guess and fitting limits
        if d50v1_init != 0:
            p0[1], low_bound[1], high_bound[1] = d50v1_init, d50v1_init - 1.75, d50v1_init + 1.75

        if d50v2_init != 0:
            p0[3], low_bound[3], high_bound[3] = d50v2_init, d50v2_init - 1.75, d50v2_init + 1.75

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 4, 3, self.n)

        # Expand boundaries for the intermediate signal
        start_id = 4 + 2 * self.n
        end_id = start_id + self.n

        low_bound[start_id:end_id] = low_bound[start_id:end_id] - np.abs(self.bStart) - np.abs(self.bEnd)
        high_bound[start_id:end_id] = high_bound[start_id:end_id] + np.abs(self.bStart) + np.abs(self.bEnd)

        # End of - Expand boundaries for the intermediate signal

        # Try to find good initia values for D50.1 and D50.2
        if d50v1_init == 0 and d50v2_init == 0:

            test_D1 = np.arange(D50min, D50max - 1.5, 1)
            test_D2 = np.arange(D50min + 2, D50max, 1)

            combinations = [(d1, d2) for d1 in test_D1 for d2 in test_D2]
            combinations = [(d1, d2) for d1, d2 in combinations if d1 < d2]

            df = pd.DataFrame(combinations, columns=['d1', 'd2'])

            rss_all = []

            for index, row in df.iterrows():
                fit_params, cov = fit_chemical_unfolding_three_species(self.xAxis_lst, self.signal_lst,
                                                                       self.temperature,
                                                                       p0, low_bound, high_bound, fit_slope_native,
                                                                       fit_slope_unfolded, signal_fx,
                                                                       self.oligo_conc_lst,
                                                                       True, row['d1'], row['d2'])

                # Insert fixed dh and ea 
                fit_params = np.insert(fit_params, 1, row['d1'])
                fit_params = np.insert(fit_params, 3, row['d2'])

                params_splt = split_all_params_three_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded)

                pred = predict_all_signal_nmer_with_intermediate_chemical(self.chem_concentration, self.temperature,
                                                                          *params_splt,
                                                                          self.oligo_conc_lst, signal_fx)

                rss = np.nansum((pred - self.signal_useful) ** 2)
                rss_all.append(rss)

            idx = np.argmin(rss_all)

            d1_init, d2_init = df['d1'][idx], df['d2'][idx]
            p0[1], p0[3] = d1_init, d2_init

            low_bound[1], low_bound[3] = d1_init - 3, d2_init - 3
            high_bound[1], high_bound[3] = d1_init + 3, d2_init + 3

        # End of - Try to find good initia values for D50.1 and D50.2

        global_fit_params, cov = fit_chemical_unfolding_three_species(self.xAxis_lst, self.signal_lst, self.temperature,
                                                                      p0, low_bound, high_bound, fit_slope_native,
                                                                      fit_slope_unfolded, signal_fx,
                                                                      self.oligo_conc_lst)

        ## Try to re fit the slopes and/or baselines if the fitted parameters are close to the boundaries
        re_fit, low_bound, high_bound = extend_bounds(global_fit_params, 4, low_bound, high_bound, 1, 40)

        if re_fit:
            global_fit_params, cov = fit_chemical_unfolding_three_species(self.xAxis_lst, self.signal_lst,
                                                                          self.temperature,
                                                                          global_fit_params, low_bound, high_bound,
                                                                          fit_slope_native, fit_slope_unfolded,
                                                                          signal_fx, self.oligo_conc_lst)

        ## Try a second re fit the slopes and/or baselines if the fitted parameters are still close to the boundaries
        re_fit2, low_bound, high_bound = extend_bounds(global_fit_params, 4, low_bound, high_bound, 5, 2e3)

        if re_fit and re_fit2:
            global_fit_params, cov = fit_chemical_unfolding_three_species(self.xAxis_lst, self.signal_lst,
                                                                          self.temperature,
                                                                          global_fit_params, low_bound, high_bound,
                                                                          fit_slope_native, fit_slope_unfolded,
                                                                          signal_fx, self.oligo_conc_lst)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / global_fit_params * 100)

        params_splt = split_all_params_three_state(global_fit_params, self.n, fit_slope_native, fit_slope_unfolded)
        errors_splt = split_all_params_three_state(errors, self.n, fit_slope_native, fit_slope_unfolded)

        m1, d50v1, m2, d50v2, bl_fold, slope_fold, bl_unfold, slope_unfold, bl_interm = params_splt
        m1_err, d50v1_err, m2_err, d50v2_err, bl_fold_err, slope_fold_err, bl_unfold_err, slope_unfold_err, bl_interm_err = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        # Fill the parameters dataframe
        for i in range(self.n):
            bN, kN, bU, kU, bI = bl_fold[i], slope_fold[i], bl_unfold[i], slope_unfold[i], bl_interm[i]
            bNe, kNe, bUe, kUe, bIe = bl_fold_err[i], slope_fold_err[i], bl_unfold_err[i], slope_unfold_err[i], \
            bl_interm_err[i]

            Y = signal_fx(self.chem_concentration, self.temperature, d50v1, m1, d50v2, m2, bN, kN, bU, kU, bI,
                          self.oligo_conc_lst[i])

            predicted.append(Y)
            fit_params.append([kN, bN, kU, bU, bI, d50v1, m1, d50v2, m2, str(self.wavelength_useful[i]), self.name])
            fit_errors.append(
                [kNe, bNe, kUe, bUe, bIe, d50v1_err, m1_err, d50v2_err, m2_err, str(self.wavelength_useful[i]),
                 self.name])

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'bI', 'D50_step1', 'M1', 'D50_step2', 'M2', self.last_col_name,
                        'Dataset']

        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = generate_predicted_matrix(predicted, self.signal_useful)

        if model == 'Monomer':
            self.fractions = fractions_fx(self.chem_concentration, self.temperature, d50v1, m1, d50v2, m2)

        else:

            fractions_lst = []
            unique_dimer_concs = np.unique(self.oligo_conc_lst)

            for unc in unique_dimer_concs:
                fr = fractions_fx(self.chem_concentration, self.temperature, d50v1, m1, d50v2, m2, unc)

                fractions_lst.append(fr)

            self.fractions = merge_fractions_dictionaries(fractions_lst, unique_dimer_concs)

        # In the same order as given to the fitting function 
        shared_params_names = ['M1', 'D50_step1', 'M2', 'D50_step2']
        local_params_names = ['bN', 'bU', 'bI']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None


class CdExperimentThermalRamp(CdExperimentFittingModel):

    def __init__(self):

        super().__init__()  # Load attributes from the parent class

        self.temperature_ori = None  # Temperature vector, one element per column in the original self.signalInput matrix

    def load_unfolding_data_monomer(self, file_path):

        """
        Assign self.wavelength_useful, self.signal_useful directly from file
        """

        try:

            temperature, signal_useful, wavelength_useful = read_unfolding_data_monomer(file_path)

            self.temperature_ori, self.temperature = temperature, temperature

            self.signal_useful     = signal_useful
            self.wavelength_useful = wavelength_useful

            self.last_col_name = 'Curve'

            return  True

        except:

                return False

    def load_unfolding_data_oligomer(self, file_path):

        """
        Assign self.wavelength_useful, self.signal_useful directly from file
         """

        try:

            temperature, signal, concentrations, self.wavelength_useful = read_unfolding_data_oligomer(file_path)

            self.temperature_ori, self.temperature = temperature, temperature

            if len(concentrations) == len(temperature):
                self.oligo_conc_molar = concentrations

            self.signal_useful = signal
            self.last_col_name = 'N-mer concentration'

            return True

        except:

            return False

    def fit_signal_three_state_irrev(self, scan_rate, fit_slope_native=False, fit_slope_unfolded=False, t1_init=0):

        self.baselines('Monomer', fit_slope_native, fit_slope_unfolded, 'Thermal')

        p0 = np.concatenate(((self.minX + 20, 120, 90, 20), self.bNs,  self.bUs, self.bIs, self.kNs, self.kUs))

        low_bound, high_bound = [self.minX + 5, 20, 65, 5], [self.maxX - 15, 400, 150, 60]

        for x in p0[4:]:

            if np.abs(x) < 0.1:
                low_bound.append(-0.25)
                high_bound.append(0.25)
            else:
                fct = 5

                low_bound.append(x / fct if x > 0 else x * fct)
                high_bound.append(x * fct if x > 0 else x / fct)

        if t1_init != 0:
            p0[0], low_bound[0], high_bound[0] = t1_init, t1_init - 15, t1_init + 15

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 4, 3, self.n)

        # Expand boundaries for the intermediate signal
        start_id = 4        + self.n*2
        end_id   = start_id + self.n

        low_bound[start_id:end_id]  = low_bound[start_id:end_id] - np.abs(self.bStart) - np.abs(self.bEnd)
        high_bound[start_id:end_id] = high_bound[start_id:end_id] + np.abs(self.bStart) + np.abs(self.bEnd)

        test_Tfs = np.arange(9, 90, 12) + t1_init
        test_Eas = 2.7 ** (np.arange(1, 6))

        combinations = [(tf, ea) for tf in test_Tfs for ea in test_Eas]
        df           = pd.DataFrame(combinations, columns=['Tfs', 'Eas'])

        rss_all = []

        for index, row in df.iterrows():

            fit_params, cov = fit_1rev_2irrev_unfolding(self.xAxis_lst, self.signal_lst, p0,
                                                        low_bound, high_bound, fit_slope_native, fit_slope_unfolded,
                                                        scan_rate, fixed_tf_ea=True, Tf_fix=row['Tfs'],
                                                        Ea_fix=row['Eas'])

            # Insert fixed dh and ea 
            fit_params = np.insert(fit_params, 2, row['Tfs'])
            fit_params = np.insert(fit_params, 3, row['Eas'])

            params_splt = split_all_params_three_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded)

            pred = predict_all_signal_irrev_three_state(self.temperature, *params_splt, scan_rate)
            rss = np.nansum((pred - self.signal_useful) ** 2)
            rss_all.append(rss)

        idx = np.argmin(rss_all)

        tf_init, ea_init = df['Tfs'][idx], df['Eas'][idx]

        p0[2], p0[3] = tf_init, ea_init

        low_bound[2], low_bound[3] = np.max((tf_init - 30, t1_init + 2)), ea_init / 2.6
        high_bound[2], high_bound[3] = tf_init + 50, ea_init * 2.6

        fit_params, cov = fit_1rev_2irrev_unfolding(self.xAxis_lst, self.signal_lst, p0,
                                                    low_bound, high_bound, fit_slope_native, fit_slope_unfolded,
                                                    scan_rate)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / fit_params * 100)

        params_splt = split_all_params_three_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded)
        errors_splt = split_all_params_three_state(errors, self.n, fit_slope_native, fit_slope_unfolded)

        t1, dh1, tf, ea, bl_fold, slope_fold, bl_unfold, slope_unfold, bl_interm = params_splt
        t1_err, dh1_err, tf_err, ea_err, bl_fold_err, slope_fold_err, bl_unfold_err, slope_unfold_err, bl_interm_err = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        # Fill the parameters dataframe
        for i in range(self.n):
            bN, kN, bU, kU, bI = bl_fold[i], slope_fold[i], bl_unfold[i], slope_unfold[i], bl_interm[i]
            bNe, kNe, bUe, kUe, bIe = bl_fold_err[i], slope_fold_err[i], bl_unfold_err[i], slope_unfold_err[i], \
            bl_interm_err[i]

            Y = rev_irrev_signal(self.temperature, dh1, t1, tf, ea, bU, bN, bI, kN, kU, scan_rate)

            predicted.append(Y)

            fit_params.append([kN, bN, kU, bU, bI, dh1, t1, tf, ea, str(self.wavelength_useful[i]), self.name])
            fit_errors.append(
                [kNe, bNe, kUe, bUe, bIe, dh1_err, t1_err, tf_err, ea_err, str(self.wavelength_useful[i]), self.name])

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'bI', 'DH1', 'T1', 'Tf', 'Ea', self.last_col_name, 'Dataset']

        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = np.array(predicted)
        self.fractions = rev_irrev_fractions(self.temperature, dh1, t1, tf, ea, scan_rate)

        # In the same order as given to the fitting function
        shared_params_names = ['T1', 'DH1', 'Tf', 'Ea']
        local_params_names = ['bN', 'bU', 'bI']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None

    def fit_signal_two_state_irrev(self, scan_rate, fit_slope_native=False, fit_slope_unfolded=False):

        self.baselines('Monomer', fit_slope_native, fit_slope_unfolded, 'Thermal')

        p0 = np.concatenate(((90, 20), self.bNs, self.bUs, self.kNs, self.kUs))

        low_bound = [65, 5]
        high_bound = [150, 60]

        for x in p0[2:]:

            if np.abs(x) < 0.1:
                low_bound.append(-0.25)
                high_bound.append(0.25)
            else:
                low_bound.append(x / 6 if x > 0 else x * 6)
                high_bound.append(x * 6 if x > 0 else x / 6)

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 2, 2, self.n)

        test_Tfs = np.arange(35, 125, 10)
        test_Eas = 2.7 ** (np.arange(1, 6))

        combinations = [(tf, ea) for tf in test_Tfs for ea in test_Eas]
        df = pd.DataFrame(combinations, columns=['Tfs', 'Eas'])

        rss_all = []

        for index, row in df.iterrows():
            fit_params, cov = fit_irrev_unfolding(self.xAxis_lst, self.signal_lst, p0,
                                                  low_bound, high_bound, fit_slope_native, fit_slope_unfolded,
                                                  scan_rate, fixed_tf_ea=True, Tf_fix=row['Tfs'], Ea_fix=row['Eas'])

            # Insert fixed dh and ea 
            fit_params = np.insert(fit_params, 0, row['Tfs'])
            fit_params = np.insert(fit_params, 1, row['Eas'])

            params_splt = split_all_params_two_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded)

            pred = predict_all_signal_irrev_two_state(self.temperature, *params_splt, scan_rate)
            rss = np.nansum((pred - self.signal_useful) ** 2)
            rss_all.append(rss)

        idx = np.argmin(rss_all)

        tf_init, ea_init = df['Tfs'][idx], df['Eas'][idx]
        p0[0], p0[1] = tf_init, ea_init

        low_bound[0], low_bound[1] = tf_init - 20, ea_init / 3
        high_bound[0], high_bound[1] = tf_init + 20, ea_init * 3

        fit_params, cov = fit_irrev_unfolding(self.xAxis_lst, self.signal_lst, p0,
                                              low_bound, high_bound, fit_slope_native, fit_slope_unfolded, scan_rate)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / fit_params * 100)

        params_splt = split_all_params_two_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded)
        errors_splt = split_all_params_two_state(errors, self.n, fit_slope_native, fit_slope_unfolded)

        Tf, Ea, blFold, slopeFold, blUnfold, slopeUnfold = params_splt
        TfErr, EaErr, blFoldErr, slopeFoldErr, blUnfoldErr, slopeUnfoldErr = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        # Fill the parameters dataframe
        for i in range(self.n):
            bN, kN, bU, kU = blFold[i], slopeFold[i], blUnfold[i], slopeUnfold[i],
            bNe, kNe, bUe, kUe = blFoldErr[i], slopeFoldErr[i], blUnfoldErr[i], slopeUnfoldErr[i]

            Y = irrev_signal(self.temperature, Tf, Ea, bU, bN, kN, kU, scan_rate)

            predicted.append(Y)

            fit_params.append([kN, bN, kU, bU, Tf, Ea, str(self.wavelength_useful[i]), self.name])
            fit_errors.append([kNe, bNe, kUe, bUe, TfErr, EaErr, str(self.wavelength_useful[i]), self.name])

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'Tf', 'Ea', self.last_col_name, 'Dataset']

        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = np.array(predicted)
        self.fractions = irrev_fractions(self.temperature, Tf, Ea, scan_rate)

        # In the same order as given to the fitting function
        shared_params_names, local_params_names = ['Tf', 'Ea'], ['bN', 'bU']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None

    def fit_signal_three_state(self, fit_slope_native=False, fit_slope_unfolded=False, t1_init=0, t2_init=0,
                               model='Monomer'):

        # Initial parameters have to be in order: 
        # Global melting temperature 1, Global enthalpy of unfolding 1,
        # Global melting temperature 2, Global enthalpy of unfolding 2,
        # Single intercepts folded, Single slopes folded,
        # Single intercepts unfolded, Single slopes unfolded

        # Set initial parameters

        self.baselines(model, fit_slope_native, fit_slope_unfolded, 'Thermal')

        signal_fx = map_three_state_model_to_signal_fx(model)
        fractions_fx = map_three_state_model_to_fractions_fx(model)

        p0 = np.concatenate(
            ((self.minX + 10, 120, self.maxX - 20, 120), self.bNs, self.bUs, self.bIs, self.kNs, self.kUs))

        t1_min = self.minX + 5
        t1_max = self.maxX - 5

        increments = {"dimer": 25, "trimer": 40, "tetramer": 55}

        for key, value in increments.items():
            if key in model.lower():
                t1_max += value

        low_bound = np.array([t1_min, 10, t1_min + 10, 10] + [x / 25 if x > 0 else x * 25 for x in p0[4:]])
        high_bound = np.array([t1_max, 1000, t1_max, 1000] + [x * 25 if x > 0 else x / 25 for x in p0[4:]])

        # Expand boundaries for the intermediate signal
        start_id = 4 + 2 * self.n
        end_id = start_id + self.n

        low_bound[start_id:end_id] = low_bound[start_id:end_id] - np.abs(self.bStart) - np.abs(self.bEnd)
        high_bound[start_id:end_id] = high_bound[start_id:end_id] + np.abs(self.bStart) + np.abs(self.bEnd)

        # End of - Expand boundaries for the intermediate signal

        if t1_init != 0:
            p0[0], low_bound[0], high_bound[0] = t1_init, t1_init - 15, t1_init + 15
            p0[2], low_bound[2] = t1_init + 5, t1_init - 15

        if t2_init != 0:
            p0[2], low_bound[2], high_bound[2] = t2_init, t2_init - 15, t2_init + 15

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 4, 3, self.n)

        step = 6
        num_rows = self.signal_useful.shape[0]

        if num_rows > 3:
            step += 2
        if num_rows > 4:
            step += 2

        if t1_init == 0 and t2_init == 0:

            test_T1s = np.arange(np.max([self.minX + 10, 20]), self.maxX - 25, step)
            test_T2s = np.arange(np.max([self.minX + 10, 20]) + step, self.maxX + 5, step)

            if 'trimer' in model:
                test_T1s = test_T1s + 10
                test_T2s = test_T1s + 10

            combinations = [(t1, t2) for t1 in test_T1s for t2 in test_T2s]
            combinations = [(t1, t2) for t1, t2 in combinations if t1 < t2]

            df = pd.DataFrame(combinations, columns=['t1', 't2'])

            rss_all = []

            for index, row in df.iterrows():
                fit_params, cov = fit_thermal_unfolding_three_species(self.xAxis_lst, self.signal_lst, p0,
                                                                      low_bound, high_bound, fit_slope_native,
                                                                      fit_slope_unfolded,
                                                                      signal_fx, listOfOligomerConc=self.oligo_conc_lst,
                                                                      fixed_t=True, t1=row['t1'], t2=row['t2'])

                # Insert fixed dh and ea 
                fit_params = np.insert(fit_params, 0, row['t1'])
                fit_params = np.insert(fit_params, 2, row['t2'])

                params_splt = split_all_params_three_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded)

                pred = predict_all_signal_nmer_with_intermediate(self.temperature, *params_splt, self.oligo_conc_lst,
                                                                 signal_fx)

                rss = np.nansum((pred - self.signal_useful) ** 2)
                rss_all.append(rss)

            idx = np.argmin(rss_all)

            t1_init, t2_init = df['t1'][idx], df['t2'][idx]
            p0[0], p0[2] = t1_init, t2_init

            low_bound[0], low_bound[2] = t1_init - 30, t2_init - 30
            high_bound[0], high_bound[2] = t1_init + 30, t2_init + 30

        global_fit_params, cov = fit_thermal_unfolding_three_species(self.xAxis_lst, self.signal_lst, p0,
                                                                     low_bound, high_bound, fit_slope_native,
                                                                     fit_slope_unfolded, signal_fx, self.oligo_conc_lst)

        ## Re fit the slopes and/or baselines if the fitted parameters are close to the boundaries
        re_fit, low_bound, high_bound = extend_bounds(global_fit_params, 4, low_bound, high_bound, 1, 500)

        if re_fit:
            global_fit_params, cov = fit_thermal_unfolding_three_species(self.xAxis_lst, self.signal_lst,
                                                                         global_fit_params,
                                                                         low_bound, high_bound, fit_slope_native,
                                                                         fit_slope_unfolded, signal_fx,
                                                                         self.oligo_conc_lst)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / global_fit_params * 100)

        params_splt = split_all_params_three_state(global_fit_params, self.n, fit_slope_native, fit_slope_unfolded)
        errors_splt = split_all_params_three_state(errors, self.n, fit_slope_native, fit_slope_unfolded)

        T1, DH1, T2, DH2, blFold, slopeFold, blUnfold, slopeUnfold, blInterm = params_splt
        T1Err, DH1Err, T2Err, DH2Err, blFoldErr, slopeFoldErr, blUnfoldErr, slopeUnfoldErr, blIntermErr = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        # Fill the parameters dataframe
        for i in range(self.n):
            bN, kN, bU, kU, bI = blFold[i], slopeFold[i], blUnfold[i], slopeUnfold[i], blInterm[i]
            bNe, kNe, bUe, kUe, bIe = blFoldErr[i], slopeFoldErr[i], blUnfoldErr[i], slopeUnfoldErr[i], blIntermErr[i]

            Y = signal_fx(self.temperature, T1, DH1, T2, DH2, bN, kN, bU, kU, bI, self.oligo_conc_lst[i], 0, 0)

            predicted.append(Y)

            fit_params.append([kN, bN, kU, bU, bI, DH1, T1, DH2, T2, str(self.wavelength_useful[i]), self.name])
            fit_errors.append(
                [kNe, bNe, kUe, bUe, bIe, DH1Err, T1Err, DH2Err, T2Err, str(self.wavelength_useful[i]), self.name])

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'bI', 'DH1', 'T1', 'DH2', 'T2', self.last_col_name, 'Dataset']

        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = generate_predicted_matrix(predicted, self.signal_useful)

        if model == 'Monomer':
            self.fractions = fractions_fx(self.temperature, DH1, DH2, T1, T2)

        else:

            fractions_lst = []
            unique_dimer_concs = np.unique(self.oligo_conc_lst)

            for unc in unique_dimer_concs:
                fr = fractions_fx(self.temperature, DH1, T1, DH2, T2, unc)
                fractions_lst.append(fr)

            self.fractions = merge_fractions_dictionaries(fractions_lst, unique_dimer_concs)

        # In the same order as given to the fitting function!
        shared_params_names = ['T1', 'DH1', 'T2', 'DH2']
        local_params_names = ['bN', 'bU', 'bI']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None

    def fit_signal_three_state_with_cp(self, fit_slope_native=False, fit_slope_unfolded=False, t1_init=0, t2_init=0,
                                       cp_th=0, model='Monomer'):

        # Initial parameters have to be in order: 
        # Global melting temperature 1, Global enthalpy of unfolding 1,
        # Global melting temperature 2, Global enthalpy of unfolding 2,
        # Single intercepts folded, Single slopes folded,
        # Single intercepts unfolded, Single slopes unfolded

        # Set initial parameters

        self.baselines(model, fit_slope_native, fit_slope_unfolded, 'Thermal')

        signal_fx = map_three_state_model_to_signal_fx(model)
        fractions_fx = map_three_state_model_to_fractions_fx(model)

        p0 = np.concatenate(
            ((self.minX + 10, 40, self.maxX - 20, 40, 0.1), self.bNs, self.bUs, self.bIs, self.kNs, self.kUs))

        t1_min = self.minX + 5
        t1_max = self.maxX - 5 if model == 'Monomer' else self.maxX + 25

        low_bound = np.array([t1_min, 10, t1_min + 10, 10, 0] + [x / 50 if x > 0 else x * 50 for x in p0[5:]])
        high_bound = np.array([t1_max, 1000, t1_max, 1000, cp_th] + [x * 50 if x > 0 else x / 50 for x in p0[5:]])

        # Expand boundaries for the intermediate signal
        start_id = 5 + 2 * self.n
        end_id = start_id + self.n

        low_bound[start_id:end_id] = low_bound[start_id:end_id] - np.abs(self.bStart) - np.abs(self.bEnd)
        high_bound[start_id:end_id] = high_bound[start_id:end_id] + np.abs(self.bStart) + np.abs(self.bEnd)

        # End of - Expand boundaries for the intermediate signal

        if t1_init != 0:
            p0[0], low_bound[0], high_bound[0] = t1_init, t1_init - 15, t1_init + 15
            p0[2], low_bound[2] = t1_init + 5, t1_init - 15

        if t2_init != 0:
            p0[2], low_bound[2], high_bound[2] = t2_init, t2_init - 15, t2_init + 15

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 5, 3, self.n)

        step = 7
        num_rows = self.signal_useful.shape[0]

        if num_rows > 3:    step += 2
        if num_rows > 4:    step += 2

        if t1_init == 0 and t2_init == 0:

            test_T1s = np.arange(np.max([self.minX + 10, 20]), self.maxX - 25, step)
            test_T2s = np.arange(np.max([self.minX + 10, 20]) + step, self.maxX + 8, step)

            combinations = [(t1, t2) for t1 in test_T1s for t2 in test_T2s]
            combinations = [(t1, t2) for t1, t2 in combinations if t1 < t2]

            df = pd.DataFrame(combinations, columns=['t1', 't2'])

            rss_all = []

            for index, row in df.iterrows():
                fit_params, cov = fit_thermal_unfolding_three_species_with_cp(self.xAxis_lst, self.signal_lst, p0,
                                                                              low_bound, high_bound, fit_slope_native,
                                                                              fit_slope_unfolded, cp_th,
                                                                              signal_fx,
                                                                              listOfOligomerConc=self.oligo_conc_lst,
                                                                              fixed_t=True, t1=row['t1'], t2=row['t2'])

                # Insert fixed dh and ea 
                fit_params = np.insert(fit_params, 0, row['t1'])
                fit_params = np.insert(fit_params, 2, row['t2'])

                params_splt = split_all_params_three_state(fit_params, self.n, fit_slope_native, fit_slope_unfolded,
                                                           withCp=True)

                T1, DH1, T2, DH2, Cp1, blFold, slopeFold, blUnfold, slopeUnfold, blInterm = params_splt

                pred = predict_all_signal_nmer_with_intermediate(
                    self.temperature, T1, DH1, T2, DH2, blFold, slopeFold, blUnfold, slopeUnfold, blInterm,
                    self.oligo_conc_lst, signal_fx, Cp1, cp_th)

                rss = np.nansum((pred - self.signal_useful) ** 2)
                rss_all.append(rss)

            idx = np.argmin(rss_all)

            t1_init, t2_init = df['t1'][idx], df['t2'][idx]
            p0[0], p0[2] = t1_init, t2_init

            threshold = 35 if 'trimer' in model else 20

            low_bound[0], low_bound[2] = t1_init - threshold, t2_init - threshold
            high_bound[0], high_bound[2] = t1_init + threshold, t2_init + threshold

        global_fit_params, cov = fit_thermal_unfolding_three_species_with_cp(self.xAxis_lst, self.signal_lst, p0,
                                                                             low_bound, high_bound, fit_slope_native,
                                                                             fit_slope_unfolded, cp_th, signal_fx,
                                                                             self.oligo_conc_lst)

        ## Re fit the slopes and/or baselines if the fitted parameters are close to the boundaries
        re_fit, low_bound, high_bound = extend_bounds(global_fit_params, 5, low_bound, high_bound, 1, 500)

        if re_fit:
            global_fit_params, cov = fit_thermal_unfolding_three_species_with_cp(self.xAxis_lst, self.signal_lst,
                                                                                 global_fit_params,
                                                                                 low_bound, high_bound,
                                                                                 fit_slope_native, fit_slope_unfolded,
                                                                                 cp_th, signal_fx, self.oligo_conc_lst)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / global_fit_params * 100)

        params_splt = split_all_params_three_state(global_fit_params, self.n, fit_slope_native, fit_slope_unfolded,
                                                   withCp=True)
        errors_splt = split_all_params_three_state(errors, self.n, fit_slope_native, fit_slope_unfolded, withCp=True)

        T1, DH1, T2, DH2, Cp1, blFold, slopeFold, blUnfold, slopeUnfold, blInterm = params_splt
        T1Err, DH1Err, T2Err, DH2Err, Cp1Err, blFoldErr, slopeFoldErr, blUnfoldErr, slopeUnfoldErr, blIntermErr = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        # Fill the parameters dataframe
        for i in range(self.n):
            bN, kN, bU, kU, bI = blFold[i], slopeFold[i], blUnfold[i], slopeUnfold[i], blInterm[i]
            bNe, kNe, bUe, kUe, bIe = blFoldErr[i], slopeFoldErr[i], blUnfoldErr[i], slopeUnfoldErr[i], blIntermErr[i]

            Y = signal_fx(self.temperature, T1, DH1, T2, DH2, bN, kN, bU, kU, bI, self.oligo_conc_lst[i], Cp1, cp_th)

            predicted.append(Y)

            fit_params.append([kN, bN, kU, bU, bI, DH1, T1, DH2, T2, Cp1, str(self.wavelength_useful[i]), self.name])
            fit_errors.append(
                [kNe, bNe, kUe, bUe, bIe, DH1Err, T1Err, DH2Err, T2Err, Cp1Err, str(self.wavelength_useful[i]),
                 self.name])

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'bI', 'DH1', 'T1', 'DH2', 'T2', 'Cp1', self.last_col_name, 'Dataset']

        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = generate_predicted_matrix(predicted, self.signal_useful)

        if model == 'Monomer':
            self.fractions = fractions_fx(self.temperature, DH1, DH2, T1, T2, None, Cp1, cp_th)

        else:

            fractions_lst = []
            unique_dimer_concs = np.unique(self.oligo_conc_lst)

            for unc in unique_dimer_concs:
                fr = fractions_fx(self.temperature, DH1, T1, DH2, T2, unc, Cp1, cp_th)

                fractions_lst.append(fr)

            self.fractions = merge_fractions_dictionaries(fractions_lst, unique_dimer_concs)

        # In the same order as given to the fitting function
        shared_params_names = ['T1', 'DH1', 'T2', 'DH2', 'Cp1']
        local_params_names = ['bN', 'bU', 'bI']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None

    def fit_signal(self, fit_slope_native=True, fit_slope_unfolded=True, cp=0, model='Monomer'):

        # Initial parameters have to be in order: 
        # Global melting temperature, Global enthalpy of unfolding,
        # Single intercepts folded, Single slopes folded,
        # Single intercepts unfolded, Single slopes unfolded

        self.baselines(model, fit_slope_native, fit_slope_unfolded, 'Thermal')

        signal_fx    = map_two_state_model_to_signal_fx(model)
        fractions_fx = map_two_state_model_to_fractions_fx(model)

        p0 = np.concatenate((((self.minX + self.maxX) / 2, 100), self.bNs, self.bUs, self.kNs, self.kUs))

        t1_min = self.minX + 5 if model == 'Monomer' else self.minX + 15
        t1_max = self.maxX - 5 if model == 'Monomer' else self.maxX + 25

        if model in ['Trimer', 'Tetramer']:

            t1_min, t1_max, p0[0] = t1_min + 20, t1_max + 30, p0[0] + 20

            if model == 'Tetramer':
                t1_min, t1_max, p0[0] = t1_min + 10, t1_max + 30, p0[0] + 20

        low_bound = [t1_min, 10] + [x / 50 if x > 0 else 50 * x for x in p0[2:]]
        high_bound = [t1_max, 600] + [50 * x if x > 0 else x / 50 for x in p0[2:]]

        p0, low_bound, high_bound = update_params_and_bounds(p0, low_bound, high_bound, fit_slope_native,
                                                             fit_slope_unfolded, 2, 2, self.n)

        global_fit_params, cov = fit_thermal_unfolding(self.xAxis_lst, self.signal_lst, p0,
                                                       low_bound, high_bound, fit_slope_native, fit_slope_unfolded,
                                                       signal_fx, cp, self.oligo_conc_lst)

        ## Re fit the slopes and/or baselines if the fitted parameters are close to the boundaries
        re_fit, low_bound, high_bound = extend_bounds(global_fit_params, 2, low_bound, high_bound, 1, 500)

        if re_fit:
            global_fit_params, cov = fit_thermal_unfolding(self.xAxis_lst, self.signal_lst, global_fit_params,
                                                           low_bound, high_bound, fit_slope_native, fit_slope_unfolded,
                                                           signal_fx, cp, self.oligo_conc_lst)

        errors = np.sqrt(np.diag(cov))
        errors = np.abs(errors / global_fit_params * 100)

        params_splt = split_all_params_two_state(global_fit_params, self.n, fit_slope_native, fit_slope_unfolded)
        errors_splt = split_all_params_two_state(errors, self.n, fit_slope_native, fit_slope_unfolded)

        Tm, dh, blFold, slopeFold, blUnfold, slopeUnfold = params_splt
        TmErr, dhErr, blFoldErr, slopeFoldErr, blUnfoldErr, slopeUnfoldErr = errors_splt

        predicted, fit_params, fit_errors = [], [], []

        # Fill the parameters dataframe
        for i in range(self.n):
            bN, kN, bU, kU = blFold[i], slopeFold[i], blUnfold[i], slopeUnfold[i],
            bNe, kNe, bUe, kUe = blFoldErr[i], slopeFoldErr[i], blUnfoldErr[i], slopeUnfoldErr[i]

            Y = signal_fx(self.temperature, Tm, dh, bN, kN, bU, kU, self.oligo_conc_lst[i], cp)

            predicted.append(Y)

            fit_params.append([kN, bN, kU, bU, Tm, dh, str(self.wavelength_useful[i]), self.name])
            fit_errors.append([kNe, bNe, kUe, bUe, TmErr, dhErr, str(self.wavelength_useful[i]), self.name])

        column_names = ['kN', 'bN', 'kU', 'bU', 'Tm', 'dH', self.last_col_name, 'Dataset']

        self.fit_params, self.fit_rel_errors = generate_params_dfs(column_names, fit_params, fit_errors,
                                                                   fit_slope_native, fit_slope_unfolded)

        self.signal_predicted = generate_predicted_matrix(predicted, self.signal_useful)

        if model == 'Monomer':
            self.fractions = fractions_fx(self.temperature, dh, Tm, None, cp)

        else:

            fractions_lst = []
            unique_oligo_concs = np.unique(self.oligo_conc_lst)

            for unc in unique_oligo_concs:
                fractions_lst.append(fractions_fx(self.temperature, dh, Tm, unc, cp))

            self.fractions = merge_fractions_dictionaries(fractions_lst, unique_oligo_concs)

        # In the same order as given to the fitting function!
        shared_params_names, local_params_names = ['Tm', 'dH'], ['bN', 'bU']

        self.bounds_df = generate_bounds_df(
            shared_params_names, local_params_names, self.wavelength_useful, fit_slope_native, fit_slope_unfolded,
            low_bound, high_bound, params_splt, errors_splt,
            self.last_col_name, self.name, self.fit_params, self.fit_rel_errors)

        return None

    def fit_signal_helix_peptide(self,numberOfCroms=32):

        # Caution! The method has sense only if self.signalDesiredUnit was created using mean unit molar ellipticity as working units

        p0         =  np.array([-0.22, -1.3, 0.002]) # Initial guess for [dg, dh, dcp]
        lowBounds  =  np.array([-5, -12, 0])
        highBounds =  np.array([5, 12, 0.008])

        if self.signalDesiredUnit is None:

            mre = self.signal_useful[0,:]

        else:

            mre     = signal_at_222nm(self.signalDesiredUnit,self.wavelength)

        fx_to_opt = get_objective_fx_peptide(numberOfCroms,poly_double_all,poly_total_all)

        params, cov = curve_fit(fx_to_opt,self.temperature, mre/(10**3), p0=p0, bounds=(lowBounds,highBounds)) # Scaling / 10**3 required for fitting

        self.fit_params = pd.DataFrame({
            'dG (kcal/mol/peptide_bond)':[params[0]],
            'dH (kcal/mol/peptide_bond)':[params[1]],
            'dCp (kcal/(mol K)/peptide_bond)':[params[2]],
            'Dataset':self.name})

        errors = np.abs(np.sqrt(np.diag(cov)) * 100 / params)

        self.fit_rel_errors = pd.DataFrame({
            'dG_rel_error '    : [errors[0]],
            'dH_rel_error '    : [errors[1]],
            'dCp_rel_error '   : [errors[2]],
            'Dataset':self.name})

        fh = calculate_cd_signal(params[0], params[1], params[2], self.temperature, numberOfCroms,poly_double_all,poly_total_all)[1]

        mre_pred = calculate_cd_signal(params[0], params[1], params[2], self.temperature, numberOfCroms,poly_double_all,poly_total_all)[0]
        mre_pred = mre_pred * 10**3 # Scaling back to original units

        self.signal_useful     = mre.reshape(1,-1)
        self.signal_predicted  = mre_pred.reshape(1, -1)
        self.wavelength_useful = ['222']

        self.fractions = {'Helix':fh,'Disordered':1-fh},

        self.bounds_df = pd.DataFrame({
            'Param': ['dG_per_peptide_bond','dH_per_peptide_bond','dCp_per_peptide_bond'],
            'Lower bound':lowBounds,
            'Fitted value':params,
            'Higher bound':highBounds,
            'Wavelength (nm)':222,
            'Dataset':self.name})

        return None


class CdExperimentTcUnfolding(CdExperimentFittingModel):
    """
    Advanced class to analyse chemical and thermal unfolding simultaneously
    """

    def __init__(self):

        super().__init__()

        self.chem_concentration = None  # 1D numpy array, length 'n'. It could be for example, the concentration of Urea, Guanidim chloride or pH
        self.chem_concentration_ori = None  # Same as self.chem_concentration, but will stay without modifications
        self.temperature_ori = None  # Same as self.temperature,        but will stay without modifications

    def tc_baselines(self, model='Monomer',
                     fitSlopeNativeT=True, fitSlopeUnfoldedT=True,
                     fitSlopeNativeD=True, fitSlopeUnfoldedD=True):

        '''
        To treat same wavelength data with different protein concentration differently
        We need to have already self.signal_useful assigned
        '''

        if model == 'Monomer':

            # Corner (Dmin,Tmin) to (Dmin,Tmax)
            selected = self.chem_concentration < np.min(self.chem_concentration) + 0.1

            signal_useful = self.signal_useful[:, selected]
            x_values = temperature_to_kelvin(self.temperature[selected])

            self.bNs, self.kNsT, self.bUs, self.kUsT = fit_baselines(signal_useful, x_values, 15)

            self.bStart = self.signal_useful[:, np.argmin(x_values)]
            self.bEnd   = self.signal_useful[:, np.argmax(x_values)]

            # Corner (Dmin,Tmin) to (Dmax,Tmin)
            selected = self.temperature < np.min(self.temperature) + 0.1

            signal_useful = self.signal_useful[:, selected]
            x_values = self.chem_concentration[selected]

            _, self.kNsD, _, self.kUsD = fit_baselines(signal_useful, x_values, 2)

        else:

            pass

    def fit_signal(self, fitSlopeNativeT=True, fitSlopeUnfoldedT=True,
                   fitSlopeNativeD=True, fitSlopeUnfoldedD=True, model='Monomer'):

        self.tc_baselines()

        signal_fx = map_two_state_tc_model_to_signal_fx(model)
        #fractions_fx = map_two_state_tc_model_to_fractions_fx(model)

        # Initial parameters have to be in order: 
        # Global M, Global D50
        # Single intercepts folded, Single intercepts unfolded
        # Single slopes folded    , Single slopes unfolded

        T50min = np.min(self.temperature) + 7
        T50max = np.max(self.temperature) - 7

        Tinit = 0.5 * (T50max + T50min)

        # Set initial parameters
        p0 = np.concatenate(((80, 50, 0, 1, 0.01), self.bNs, self.bUs, self.kNsT, self.kUsT, self.kNsD, self.kUsD))

        low_bound = np.array([10, T50min, 0, 0.01, 0] + [x / 50 if x > 0 else 50 * x for x in p0[5:]])
        high_bound = np.array([1000, T50max, 10, 20, 5] + [50 * x if x > 0 else x / 50 for x in p0[5:]])

        #p0, low_bound, high_bound = update_params_and_bounds(p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded,5,2,self.n)  

        global_fit_params, cov = fit_tc_unfolding([self.chem_concentration], self.signal_useful.tolist(),
                                                  [self.temperature], p0, low_bound, high_bound, True, True,
                                                  True, True, signal_fx, listOfOligomerConc=None)

        return None


class CdExperimentCustomAnalysis(CdExperimentFittingModel):
    '''
    Advanced class to analyse the CD data using custom models. 
    For example, binding or TFE titrations
    '''

    def __init__(self):

        super().__init__()  # Load attributes from the parent class
        self.model_function = None

        # Numpy array containing a certain measurement dimension
        self.first_measurement_dimension = None  # E.g., temperature or ligand concentration
        self.second_measurement_dimension = None

        # experimental parameters names
        self.first_exp_param_name = None
        self.second_exp_param_name = None

        self.cleaned_text_function = None
        
        self.global_params_names         = None
        self.local_params_names          = None
        self.local_params_names_extended = None

        self.n_global_params = None
        self.n_local_params = None

        self.p0 = None
        self.paramsSign = None
        self.all_fitted_params = None
        self.boundaries_updated = None
        self.low_bound = None
        self.high_bound = None
        

    def generate_model_function(self, text):

        cleaned_text = clean_function_text(text)

        self.cleaned_text_function = cleaned_text

        fit_params_names = extract_words(cleaned_text)

        fit_params_names = [x for x in fit_params_names if
                            x not in ['sqrt', 'log', 'np', 'exp', self.first_exp_param_name,
                                      self.second_exp_param_name]]

        # remove duplicates and sort it
        fit_params_names = sorted(list(set(fit_params_names)))

        # find global and local parameters
        global_params_names = [x for x in fit_params_names if 'global' in x.lower()]
        local_params_names = [x for x in fit_params_names if 'global' not in x.lower()]

        local_params_names_extended = []

        for wl in self.wavelength_useful:
            local_params_names_extended += [str(int(wl)) + ' ' + x for x in local_params_names]

        paramsSign = []

        for param in global_params_names + local_params_names:

            if 'Pos' in param:

                paramsSign.append('positive')

            elif 'Neg' in param:

                paramsSign.append('negative')

            else:

                paramsSign.append('unknown')

        self.paramsSign = paramsSign

        wavelengths = self.wavelength_useful

        self.global_params_names = global_params_names
        self.local_params_names = local_params_names
        self.local_params_names_extended = local_params_names_extended

        self.n_global_params = len(global_params_names)
        self.n_local_params = len(local_params_names)

        # Create a Python function using eval
        def custom_function(dummy_variable, *args):

            """
            args should be in this order: 
                global param 1, ... , global param N, 
                local param 1 (first wavelength), local param 1 (second wavelength), 
                ...  local param Z (last wavelength).

            The total number of arguments should equal 

                #global_params + (#local_params * #wavelengths)

            """

            signal = []

            dInit = {'np': np,
                     self.first_exp_param_name: self.first_measurement_dimension,
                     self.second_exp_param_name: self.second_measurement_dimension
                     }

            # initialize the i counter, just in case there are no global params
            i = -1
            for i, fit_param in enumerate(global_params_names):
                dInit[fit_param] = args[i]

            counter = 0

            for _ in wavelengths:

                d = dInit.copy()

                for fit_param in local_params_names:
                    d[fit_param] = args[i + counter + 1]

                    counter += 1

                signal.append(eval(cleaned_text, d))

            signal = np.concatenate(signal, axis=0)

            return signal

        self.model_function = custom_function

        return None

    def get_initial_values(self, left_limit_log_search=-3, right_limit_log_search=3):

        signal = self.signal_useful
        wls = self.wavelength_useful

        nParamsPerCurve = len(self.local_params_names) + len(self.global_params_names)

        # Generate a list of list containing a
        # combination of log-spaced values
        combis = generate_initial_params_combinations(nParamsPerCurve, self.paramsSign,
                                                      left_limit_log_search, right_limit_log_search)

        # Try to find acceptable initial values by running the user defined function with different parameter values

        dCombInit = {'np': np,
                     self.first_exp_param_name: self.first_measurement_dimension,
                     self.second_exp_param_name: self.second_measurement_dimension
                     }

        predicted = []

        # Suppress overflow warnings only for this block
        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)

            for comb in combis:

                dComb = dCombInit.copy()

                i = -1

                for i, fit_param in enumerate(self.global_params_names):
                    dComb[fit_param] = comb[i]

                for ii, fit_param in enumerate(self.local_params_names):
                    dComb[fit_param] = comb[i + ii + 1]

                predicted.append(eval(self.cleaned_text_function, dComb))

        all_wl_rss = []
        sel_combi = []

        for i, wl in enumerate(wls):

            wl_rss = []

            signalSlice = signal[i, :]

            for prediction in predicted:

                if np.isnan(prediction).any():
                    wl_rss.append(np.inf)

                else:

                    residuals = signalSlice - prediction
                    wl_rss.append(np.sum(residuals ** 2))

            all_wl_rss.append(np.min(wl_rss))
            sel_combi.append(combis[np.argmin(wl_rss)])

        p0 = []

        if self.n_global_params > 0:
            # find the overall best combination of initial parameters
            best_combi = sel_combi[np.argmin(all_wl_rss)]

            p0 = best_combi[:self.n_global_params]

            # Remove the initial guess of the global parameter
            sel_combi = [x[self.n_global_params:] for x in sel_combi]

        sel_combi = [item for sublist in sel_combi for item in sublist]

        p0 += sel_combi

        # End of trying to find initial acceptable values

        tol_factor = 500

        low_bound = [x / tol_factor if x > 0 else tol_factor * x for x in p0]
        high_bound = [tol_factor * x if x > 0 else x / tol_factor for x in p0]

        self.p0 = round_to_significant_digits(p0, 2)
        self.low_bound = round_to_significant_digits(low_bound, 2)
        self.high_bound = round_to_significant_digits(high_bound, 2)

        return None

    def fit_signal(self):

        """
        Find the fitting parameters by using the function generated through self.generate_model_function
        """

        # Reset the fitting params and the predicted signal
        self.fit_params = None
        self.fit_rel_errors = None
        self.signal_predicted = None
        self.boundaries_updated = False

        if self.p0 is None:
            self.get_initial_values()

        signal = self.signal_useful
        wls = self.wavelength_useful

        c1 = np.all(self.p0 <= self.high_bound)
        c2 = np.all(self.p0 >= self.low_bound)

        if not c1 or not c2:
            return None

        try:

            all_fitted_params, cov = curve_fit(self.model_function, 1, signal.flatten(),
                                               p0=self.p0,
                                               bounds=(tuple(self.low_bound), tuple(self.high_bound)))

            # Try to re-fit if the values are close to the boundaries +/- 500
            re_fit, self.low_bound, self.high_bound = extend_bounds(all_fitted_params, 0, self.low_bound,
                                                                    self.high_bound, 1, 500)

            if re_fit:
                self.boundaries_updated = True

                self.p0 = all_fitted_params

                all_fitted_params, cov = curve_fit(self.model_function, 1, signal.flatten(),
                                                   p0=self.p0,
                                                   bounds=(tuple(self.low_bound), tuple(self.high_bound)))

            self.all_fitted_params = all_fitted_params
            errors = np.sqrt(np.diag(cov))

            predicted = self.model_function(None, *all_fitted_params)
            predicted = predicted.reshape((self.signal_useful.shape[0], -1))

            self.signal_predicted = predicted

            # initialize the dataframe that will have the fitted parameters values
            local_params_df = pd.DataFrame({'Condition': wls})
            global_params_df = pd.DataFrame({'Condition': wls})

            local_params_df_error = pd.DataFrame({'Condition': wls})
            global_params_df_error = pd.DataFrame({'Condition': wls})

            if self.n_local_params > 0:
                column_names = self.local_params_names
                local_params1 = all_fitted_params[self.n_global_params:]
                local_params2 = local_params1.reshape((-1, self.n_local_params))
                local_params3 = pd.DataFrame(local_params2, columns=column_names)
                local_params_df = pd.concat([local_params3, local_params_df], axis=1)

                local_params_err1 = errors[self.n_global_params:]
                local_params_err2 = local_params_err1.reshape((-1, self.n_local_params))
                local_params_err3 = np.abs(local_params_err2 / local_params2 * 100)
                local_params_err4 = pd.DataFrame(local_params_err3, columns=column_names)
                local_params_df_error = pd.concat([local_params_err4, local_params_df_error], axis=1)

            if self.n_global_params > 0:

                for i, param in enumerate(all_fitted_params[:self.n_global_params]):
                    global_params_df[self.global_params_names[i]] = param

                    rel_err = np.abs(errors[i] / param * 100)

                    global_params_df_error[self.global_params_names[i]] = rel_err

            all_params_df = pd.merge(global_params_df, local_params_df)
            all_params_df_err = pd.merge(global_params_df_error, local_params_df_error)

            self.fit_params = all_params_df
            self.fit_rel_errors = all_params_df_err

        except:

            pass

        return None


class CdAnalyzer:
    """
    Useful to work with many different cd experiments
    """

    def __init__(self):

        """
        Create 
        """

        self.experimentsOri = {}  # Dictionary where each key value pair corresponds to one CD experiment
        self.experimentsModif = {}  # To subset the wavelength range of the CD experiments inside experimentsOri
        self.experimentNames = []
        self.sharedParameters = False  # To modify all the experiment at the same time with shared parameters

        self.experimentsThermal = {}  # Dictionary where each key value pair corresponds to one CD based thermal unfolding experiment
        self.experimentNamesThermal = []

        self.experimentsChemical = {}  # Dictionary where each key value pair corresponds to one CD based chemical unfolding experiment
        self.experimentNamesChemical = []

        self.experimentsCustom = {}  # Dictionary where each key value pair corresponds to one CD experiment with a custom formula for analysis
        self.experimentNamesCustom = []

    def load_experiment(self, file, name):

        """
        Append one experiment to experimentsOri 
        """

        if name in self.experimentNames:
            return "Experiment name already selected!"

        try:

            self.experimentsOri[name] = CdExperimentGeneral()
            self.experimentsOri[name].load_data(file, name)

            self.experimentNames.append(name)

            # Adapt the number of significant digits

            return "Data loaded successfully!!!"

        except:

            pass

        return "Data could not be loaded"

    def delete_experiment(self, names):

        # Convert to list if required
        if not isinstance(names, list):
            names = [names]

        for name in names:

            self.experimentNames.remove(name)
            del self.experimentsOri[name]

            try:
                del self.experimentsModif[name]
            except:
                pass

        return None

    def clean_experiments(self, experiment_type):

        if experiment_type == 'custom':
            self.experimentsCustom = {}
            self.experimentNamesCustom = []

        if experiment_type == 'thermal':
            self.experimentsThermal = {}
            self.experimentNamesThermal = []

        if experiment_type == 'chemical':
            self.experimentsChemical = {}
            self.experimentNamesChemical = []

        return None

    def experiments_to_absorbance_units(self, unitsDic):

        """
        Sets the proper units for each experiment in experimentsOri
            e.g. of unitsList - {'exp' : 'milliabsorbance','blank' : 'degrees', ... }
        """

        lst_of_experiments_with_new_units = []

        for experimentName, experimentUnits in unitsDic.items():

            if not self.experimentsOri[experimentName].isFakeExperiment:

                if experimentUnits != self.experimentsOri[experimentName].units:
                    lst_of_experiments_with_new_units.append(experimentName + ' : ' + experimentUnits)

                self.experimentsOri[experimentName].experiment_to_absorbance_units(experimentUnits)

        # Return which experiments were modified
        return lst_of_experiments_with_new_units

    def all_experiments_to_absorbance_units(self, experiment_units):

        for experimentName in self.experimentNames:

            if not self.experimentsOri[experimentName].isFakeExperiment:
                self.experimentsOri[experimentName].experiment_to_absorbance_units(experiment_units)

        return None

    def experiments_absorbance_units_to_other_units(self, units_dic):

        """
        Transform absorbance units to desired units
            e.g. of unitsList - {'exp' : 'milliabsorbance','blank' : 'degrees', ... }

        This function must be always called after running setExperimentUnits !!!
        """

        for experimentName, experimentUnits in units_dic.items():

            if not self.experimentsOri[experimentName].isFakeExperiment:

                self.experimentsOri[experimentName].experiment_from_abs_to_other_units(experimentUnits)

            else:

                self.experimentsOri[experimentName].signalDesiredUnit = self.experimentsOri[experimentName].signalInput

        return None

    def all_experiments_absorbance_units_to_other_units(self, experiment_units):

        for experimentName in self.experimentNames:

            if not self.experimentsOri[experimentName].isFakeExperiment:

                self.experimentsOri[experimentName].experiment_from_abs_to_other_units(experiment_units)

            else:

                self.experimentsOri[experimentName].signalDesiredUnit = self.experimentsOri[experimentName].signalInput

        return None

    def initialize_experiment_modif(self):

        """
        Required to store the original data and allow post-processing
        """

        self.experimentsModif = copy.deepcopy(self.experimentsOri)

        return None

    def set_experiment_properties(self, experiment_name, variable, value):

        """
        experimentName must be in self.experimentNames
        """

        setattr(self.experimentsOri[experiment_name], variable, value)
        setattr(self.experimentsModif[experiment_name], variable, value)

        return None

    def get_experiment_properties(self, variable):

        return [getattr(self.experimentsOri[experimentName], variable) for experimentName in self.experimentNames]

    def get_experiment_properties_modif(self, variable):

        return [getattr(self.experimentsModif[experimentName], variable) for experimentName in self.experimentNames]

    def filter_data_by_wavelength(self, min_wl, max_wl):

        self.initialize_experiment_modif()

        experiment_names = self.experimentNames

        for expName in experiment_names:

            # Use the original signal! - Allows going back
            wl = self.experimentsOri[expName].wavelength
            signal_desired = self.experimentsOri[expName].signalDesiredUnit
            signal_input = self.experimentsOri[expName].signalInput
            signal_abs = self.experimentsOri[expName].signalAbs
            signal_ht = self.experimentsOri[expName].signalHT

            if not np.isnan(signal_desired).all():
                self.experimentsModif[expName].signalDesiredUnit = filter_matrix_by_vector(signal_desired, wl, min_wl,
                                                                                           max_wl)
            else:
                self.experimentsModif[expName].signalDesiredUnit = np.nan

            if not np.isnan(signal_input).all():
                self.experimentsModif[expName].signalInput = filter_matrix_by_vector(signal_input, wl, min_wl, max_wl)
            else:
                self.experimentsModif[expName].signalInput = np.nan

            if not np.isnan(signal_abs).all():
                self.experimentsModif[expName].signalAbs = filter_matrix_by_vector(signal_abs, wl, min_wl, max_wl)
            else:
                self.experimentsModif[expName].signalAbs = np.nan

            self.experimentsModif[expName].signalHT = filter_matrix_by_vector(signal_ht, wl, min_wl, max_wl)
            self.experimentsModif[expName].wavelength = filter_vector_by_values(wl, min_wl, max_wl)

        return None



test4 = False

if test4:
    T = np.arange(20, 90, 2.5)
    D = np.arange(1, 6, 0.5)
    T_grid, D_grid = np.meshgrid(T, D)

    T, D = T_grid.flatten(), D_grid.flatten()

    yy = signal_two_state_tc_unfolding_monomer(T, D, 50, 50, 0, 2, 0.01, 10, 5, 0, 0, 0, 0)

    c = CdExperimentTcUnfolding()
    c.temperature = T
    c.chem_concentration = D
    c.signal_useful = np.array([yy]).reshape(1, -1)

    c.fit_signal()

    import matplotlib.pyplot as plt

    plt.scatter(T, y, c=D, cmap='viridis')
    #plt.plot(T,signal_two_state_tc_unfolding_monomer(T, D,*fit))
    plt.show()

    #plt.scatter(D,y,c=T, cmap='viridis')
    #plt.scatter(T,)
    #plt.show()
