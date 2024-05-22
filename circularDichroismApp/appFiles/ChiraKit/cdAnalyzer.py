import copy
import sys
import re
import warnings

from loadCDfilesHelpers  import *
from cdUnitsConverter    import *
from helpers             import *
from fitting_helpers     import *

from spectra_comparison_helpers     import *

from scipy.stats         import t
from scipy.stats         import chi2
from scipy.signal        import savgol_filter

import numpy  as np 
import pandas as pd 

import numpy.polynomial.polynomial as poly

sys.path.append('./secondary_structure_estimation_files')

from SelconsFunction                 import *

"""
Classes for the analysis of circular dichorism data
Code written by Osvaldo Burastero 

If you use this script please use a reference similar to the following:

    
Osvaldo Burastero, Nykola C. Jones, Søren Vrønning Hoffmann, & Maria M. Garcia-Alai (2024). 
ChiraKit (Version 1.0). Manuscript in preparation. https://spc.embl-hamburg.de/app/chirakit

No warranty whatsoever
If you have questions please contact me:
    oburastero@gmail.com

"""
 
# Matrix of secondary structure components (six components)
F1  = np.loadtxt("./secondary_structure_estimation_files/AU-F128_T-Nov22.txt",          dtype='f',   delimiter='\t') 
# Matrix of known (reference) CD spectra     
A1  = np.loadtxt("./secondary_structure_estimation_files/AU-A128_PCDDB-Nov22.txt",      dtype='f',   delimiter='\t')

# F1 and A1 contain the AU-SP175 and AU-SMP180 reference datasets

# Labels of the known CD spectra
#Lbl = np.loadtxt("./secondary_structure_estimation_files/Labels-SMP180_PCDDBOrder.txt", dtype='str', delimiter='\t')

class cd_experiment_general:

    '''
    Basic class for all CD spectra. It allows converting between different CD units
    '''

    def __init__(self):

        self.wavelength         = None  # 1D numpy array 
        self.signalInput        = None  # 2D dataframe with one column per measurement
        self.signalAbs          = None  # 2D dataframe with one column per measurement, absorbance units 
        self.signalDesiredUnit  = None  # 2D dataframe with one column per measurement, signal in user selected units 
        self.internalID         = None  # Internal ID based on the experiment name and spectra names 
        self.spectraNames       = None  # Names of the measured spectra, e.g., 

        self.signalHT           = None  # 2D dataframe with one column per measurement, unknown units

        # Metadata of the experiment, e.g, concentration, path length ... 
        self.metadata = None

        # Input units of the experiment 
        self.units              = 'millidegrees'  # Default Input units

        # Working units, selected by the user
        # String, one of: 'milliabsorbance', 'absorbance', molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity', 
        #                 'meanResidueMolarExtinction' or 'meanResidueMolarEllipticity'
        self.desiredUnits       = 'millidegrees'  

        # Start with an impossible temperature (in degree Celsius), so the user has to change it.
        self.temperature        = np.NaN

        self.numberOfResidues   = 0   # If protein sample, number of residues
        self.concentration      = 0   # Float, in mg/ml
        self.pathLength         = 0   # Float, in cm
        self.molecularWeight    = 0   # Float, in Dalton

        # Boolean to decide if the user needs to modify the parameters numberOfResidues, concentration, pathLength and molecularWeight
        # if 'isFakeExperiment' is set to True, then the matrix signalDesiredUnit won't be changed once it is created for the first time 
        self.isFakeExperiment       = False 

        # String, to know how the 'fake experiment' matrix 'signalDesiredUnit' was created
        # For example, it should be either 'molarEllipticity',  'molarExtinction', 'meanResidueMolarExtinction' or 'meanResidueMolarEllipticity'
        self.fakeExperimentSignal   = 'unknown' 

        self.isGenerated        = False # Boolean to help exporting only the generated data
        
        # List of dataframes containing the secondary structure content, one df per CD spectrum
        self.secondary_structure_content = None

        return None

    def experiment2absorbanceUnits(self,unitsBegin):

        """
        unitsBegin should be: 
            'absorbance' milliabsorbance', 'molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity',
            'meanResidueMolarExtinction' or 'meanResidueMolarEllipticity'
        """

        self.units     = unitsBegin

        signalInAbsUnits = convert2absorbance(self.signalInput,
            unitsBegin, self.concentration, self.pathLength, self.molecularWeight, self.numberOfResidues)

        self.signalAbs = signalInAbsUnits
          
        return None

    def experimentFromAbsorbanceUnits2otherUnits(self,unitsEND):

        """
        unitsEND should be: 
            'absorbance' milliabsorbance', 'molarExtinction', 'degrees', 'millidegrees', 'molarEllipticity',
            'meanResidueMolarExtinction' or meanResidueMolarEllipticity

        """

        self.desiredUnits      = unitsEND

        # Just copy the matrix in case the input and desired units are the same
        if self.desiredUnits == self.units:

            self.signalDesiredUnit = self.signalInput

        # Otherwise, compute the signal using the desired CD units
        else:

            signalNew = absorbance2desiredUnits(self.signalAbs,
                unitsEND, self.concentration, self.pathLength, self.molecularWeight, self.numberOfResidues)

            self.signalDesiredUnit = signalNew

        return None

    def load_data(self,file,name):

        fileType = detectFileType(file)

        # Define the dictionaries to call the data loading function 

        readDataFunctionDict = {

        "ChirascanFile" : readChirascanFileData,
        "PCCDBFile"     : readPCCDBFileData,
        "GenCDfile"     : readGenFileData,
        "DatFile"       : readDatFileData,
        "d0xFile"       : readD0xFileData,
        "plain_csv"     : read_custom_csv,
        "jasco_simple"  : read_jasco_single_sample_csv,
        "chirakit_txt_with_header"      : readChiraKitTxtData,
        'ChirascanFileTemperatureRamp'  : readChirascanFileDataThermalRamp
        }

        readMetadataFunctionDict = {

        "ChirascanFile" : readChirascanFileMetaData,
        "PCCDBFile"     : readPCCDBFileMetaData,
        "GenCDfile"     : readGenFileMetaData,
        "DatFile"       : readDatFileMetaData,
        "d0xFile"       : readD0xFileMetaData,
        "plain_csv"     : readChiraKitTxtMetaData,
        "jasco_simple"  : read_jasco_single_MetaData,
        "chirakit_txt_with_header"      : readChiraKitTxtMetaData,
        "ChirascanFileTemperatureRamp"  : readChirascanFileMetaData
        }

        # Read the metadata
        self.metadata   = readMetadataFunctionDict[fileType](file)

        # If the file is of type Chirascan temperature ramp, assign spectra names based on temperature
        if fileType == 'ChirascanFileTemperatureRamp':

            self.wavelength, self.signalInput, self.temperature, self.signalHT  = readDataFunctionDict[fileType](file)
        
            self.internalID   = [name + ' ' + str(round(t,1)) for t in self.temperature]
            self.spectraNames = [name + ' ' + str(round(t,1)) for t in self.temperature]

        else:

            self.wavelength, self.signalInput, self.spectraNames, self.signalHT = readDataFunctionDict[fileType](file)

            # If the file is of type d0x, load the temperature data directly from the file
            if fileType == 'd0xFile':

                self.temperature = readTemperatureD0x(file)

            # If the file is not of type d0x and not of type Chirascan temperature ramp, 
            # try to load the temperature data from the metadata 
            else:

                try:

                    for key, value in self.metadata.items():

                        if 'temp' in key.lower():

                            value = value.replace(",", " ")

                            digits = ''.join(c for c in value if c.isdigit() or c == '.')
 
                            # In case we have a numeric like string o less than 6 characters, we assume we have only one temperature value
                            if is_float(digits) and len(digits) < 6:
                                
                                self.temperature = float(digits)
                                break

                            # Try to split based on spaces
                            if are_all_strings_numeric(value.split()):

                                self.temperature = np.array([float(x) for x in value.split()])
                                break
                            
                except:

                    pass # Keep default temperature 

            # Assign internalID if we only have one spectrum

            if self.signalInput.shape[1] == 1:
                self.internalID   = [name]
    
            else:

                # Expand the temperature float into a vector, if required. 
                if isinstance(self.temperature, float):

                    self.temperature = np.array([self.temperature for _ in self.spectraNames])
                
                # Assign internalID if we have many spectra
                self.internalID = [name + ' ' + x for x in self.spectraNames]

            self.spectraNames = self.internalID

        self.units = guess_input_units_from_metadata_dictionary(self.metadata)

        self.numberOfResidues   = guess_parameter_from_metadata_dictionary(self.metadata,['number of residues','nres','aminoacids'])
        self.concentration      = guess_parameter_from_metadata_dictionary(self.metadata,['concentration','conc','mg'])
        self.pathLength         = guess_parameter_from_metadata_dictionary(self.metadata,['path length','path',' length'])
        self.molecularWeight    = guess_parameter_from_metadata_dictionary(self.metadata,['molecular weight','mw','dalton'])

        self.name = name

        # Sort in increasing order
        # Get the indices that would sort the wavelength vector
        sorted_indices = np.argsort(self.wavelength)

        # Use the sorted_indices to rearrange the rows of the matrix
        self.signalInput    = self.signalInput[sorted_indices]
        self.wavelength     = self.wavelength[sorted_indices]
        self.signalHT       = self.signalHT[sorted_indices]

        return None

    def init_and_check_secondary_str_method(self,lowerWL):

        # Verify that we have the protein concentration, path length, number of aminoacids, and weight (if needed)
        # Verify that the higher WL >= 240 nm and lower WL <= 190 nm

        self.secondary_structure_content = []
        self.fitted_spectra_sec_str      = {}
        self.query_spectra_sec_str       = {}
        self.lowerWL_sec_str             = lowerWL

        self.currentUnits = self.desiredUnits
        shouldChangeUnits = self.currentUnits != 'meanResidueMolarExtinction'

        shouldStop1 = any([
            self.numberOfResidues == 0, self.concentration   == 0,
            self.pathLength       == 0, self.molecularWeight == 0
            ])

        # Check that we actually need the parameters: protein concentration, path length, number of aminoacids, and weight
        shouldStop2 = self.units != 'meanResidueMolarExtinction'

        shouldStop3 = any([
            np.max(self.wavelength) < 240, 
            np.min(self.wavelength) > 190
            ])

        if (shouldStop1 and shouldStop2) or shouldStop3:

            return False

        return True

    def set_secondary_structure_method_references_default(self):

        wavelength_temp = self.wavelength[self.wavelength >= self.lowerWL_sec_str]
        
        # Obtain the lower wavelength limit (can't be lower than 175) 
        LimWL = int(np.max([np.min(wavelength_temp),175]))

        # Select the most adequate spectra reference database
        if LimWL >= 180:

            F          = F1[:,:] 
            self.A_sel = A1[:(240-LimWL+1),:] # At most, the first 61 WL (240-180), all 128 spectra (AU-SMP180)

            # Create a new F matrix with only four secondary structure components
            F_alpha        = F[0,:] + F[1,:]
            F_beta         = F[2,:] + F[3,:]
            self.F_sel     = np.array([F_alpha, F_beta, F[4,:], F[5,:]])
            self.SStruct_labels = ['Alpha', 'Beta', 'Turns', 'Unord']
            
        else:
            
            self.F_sel          = F1[:,:71] 
            self.A_sel          = A1[:(240-LimWL+1),:71] # At most, all 66 WL (240-175), first 71 spectra (AU-SP175)
            self.SStruct_labels = ['Alpha-r','Alpha-d', 'Beta-r', 'Beta-d', 'Turns', 'Unord']

        # Wavelength even sequence with 1 nm step
        self.wl_even_sequence = np.arange(240, LimWL-1, -1)

        return None

    def set_secondary_structure_method_references_user(self,matrixF,matrixA,max_wl_ref,wl_ref_step,SStruct_labels):

        '''
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

        '''

        wl_even_sequence = max_wl_ref - ( np.arange(0,matrixA.shape[0]) * wl_ref_step )

        cond1 = wl_even_sequence >= np.min(self.wavelength)
        cond2 = wl_even_sequence <= np.max(self.wavelength)
        cond3 = wl_even_sequence >= self.lowerWL_sec_str

        selected_idx = np.logical_and(cond1, np.logical_and(cond2, cond3))

        self.F_sel            = matrixF
        self.A_sel            = matrixA[selected_idx,:]
        self.wl_even_sequence = wl_even_sequence[selected_idx]
        self.SStruct_labels = SStruct_labels

        return None

    def estimate_secondary_structure(self):

        '''
        Run the Selcon3 algorithm to estimate the percentage of different secondary structure components
        '''

        self.experimentFromAbsorbanceUnits2otherUnits('meanResidueMolarExtinction')

        # Do secondary structure fitting
        signal_temp     = self.signalDesiredUnit[ self.wavelength >= self.lowerWL_sec_str  ]
        wavelength_temp = self.wavelength[        self.wavelength >= self.lowerWL_sec_str  ]

        # Sort in increasing order, in case that the data  was not loaded with the self.load_data() function
        # Get the indices that would sort the wavelength vector
        sorted_indices = np.argsort(wavelength_temp)

        # Use the sorted_indices to rearrange the rows of the matrix and vector
        signal_temp     = signal_temp[sorted_indices]
        wavelength_temp = wavelength_temp[sorted_indices]

        # Total CD spectra of this experiment
        nSpectra = self.signalDesiredUnit.shape[1]

        # Initalize list of empty lists
        # Each sublist will contain the results of the secondary structure fitting
        self.secondary_structure_content = [ None for _ in range(nSpectra)]

        # Iterate over the columns 
        for i in range(nSpectra):
        
            # Interpolate 
            delta_epsilon = np.interp(self.wl_even_sequence, wavelength_temp, signal_temp[:,i], left=None, right=None, period=None)

            try:

                return_List, sec_str_df, mean_refit_prot, query_spectrum =  SelconsPy(self.A_sel, self.F_sel, 
                    delta_epsilon, self.SStruct_labels)

                self.secondary_structure_content[i]               = sec_str_df
                self.fitted_spectra_sec_str[self.spectraNames[i]] = mean_refit_prot
                self.query_spectra_sec_str[self.spectraNames[i]]  = query_spectrum

            except:

                pass

        # End of - Do secondary structure fitting
       
        self.experimentFromAbsorbanceUnits2otherUnits(self.currentUnits)

        return None

class cd_experiment_fitting_model(cd_experiment_general):

    '''
    Advanced class for CD spectra with more measurements dimensions. For example, to load spectra at different temperatures or urea concentration
    This class is used as a starting point to build the class for thermal ramp analysis and for chemical denaturation analysis
    '''

    def __init__(self):

        super().__init__()
        self.fit_params            = None   # Values of the fitted parameters
        self.fit_rel_errors        = None   # Relative error of the fitted parameters
        self.signal_predicted      = None   # Predicted signal
        self.name                  = None   # String, experiment name
        self.decompositionDone     = False
        return None

    def decompose_spectra_pca(self):

        try:

            X         = self.signalDesiredUnit.T

            explained_variance, basis_spectra, coefficients = apply_pca(X) 

            self.coefficients_all     = coefficients
            self.basis_spectra_all    = basis_spectra
            
            # Cumulated explained variance of the components
            self.explained_variance  = explained_variance

            self.pca_based = True
            self.decompositionDone = True

        except:

            self.decompositionDone = False
            pass

        return None

    def decompose_spectra_svd(self):

        try:

            explained_variance, basis_spectra, coefficients = apply_svd(self.signalDesiredUnit) 

            # basis spectra and associated coefficients   
            self.basis_spectra_all    = basis_spectra       
            self.coefficients_all     = coefficients

            # Cumulated explained variance of the components
            self.explained_variance  = explained_variance

            self.pca_based = False
            self.decompositionDone = True

        except:

            self.decompositionDone = False
            pass

        return None

    def filter_basis_spectra(self,explained_variance_threshold=99):

        """
        Should be always run after decompose_spectra_svd or decompose_spectra_pca 
        so we create the matrices self.basis_spectra and self.coefficients
        based on the selected variance threshold
        """

        # Find the number of components (k) that capture at least threshold*100 percent of the variance or correlation
        k = np.sum(self.explained_variance < explained_variance_threshold) + 1

        self.k                = k
        
        self.basis_spectra    = self.basis_spectra_all[:,:k]
        self.coefficients     = self.coefficients_all[:k,:]

        return None

    def align_basis_spectra_and_coefficients(self):

        """
        Try to align the n selected basis spectra and the associated coefficients to the observed peak
        """

        self.basis_spectra, self.coefficients = align_basis_spectra_and_coefficients(self.signalDesiredUnit,self.basis_spectra,self.coefficients)

        return None

    def invert_selected_spectrum(self,n):

        n = int(n)

        if n <= self.basis_spectra.shape[1]:

            self.basis_spectra[:,n]  = - self.basis_spectra[:,n]
            self.coefficients[n,:]   = - self.coefficients[n,:]

        return None

    def rotate_basis_spectra(self):

        n_basis_spectra = self.basis_spectra.shape[1]

        if  n_basis_spectra == 2:

            self.basis_spectra, self.coefficients = rotate_two_basis_spectra(self.signalDesiredUnit,self.basis_spectra_all,self.pca_based)

        if n_basis_spectra == 3:

            self.basis_spectra, self.coefficients = rotate_three_basis_spectra(self.signalDesiredUnit,self.basis_spectra_all,self.pca_based)

        # Compute again the explained variance after rotation!

        data = self.signalDesiredUnit

        if self.pca_based:

            data_mean    = np.mean(data, axis=1,keepdims=True)
            data         = data - data_mean

        total_variance  = np.linalg.norm(data)**2  # Total variance in the data

        self.explained_variance = explained_variance_from_orthogonal_vectors(self.basis_spectra,self.coefficients,total_variance)

        return None

    def reconstruct_spectra(self):

        self.fitted_spectra = reconstruct_spectra(self.basis_spectra,self.coefficients,self.signalDesiredUnit,self.pca_based)

        return None

    def assign_useful_signal_svd(self,relevant_k=1):
        '''
        Trick to be able to fit the svd/pca coefficients
        This method will add the properties
            signal_useful and wavelength_useful

            signal_useful       contains the change in the coefficients
            wavelength_useful   contains the number of relevant coefficients
        '''

        self.signal_useful     = self.coefficients[0:relevant_k,:]
        self.wavelength_useful = np.arange(1,relevant_k+1)

        return None

    def assign_useful_signal(self,user_selected_wavelength = None):

        # Default method, use the filtered wavelength
        if user_selected_wavelength is None:

            self.signal_useful     = self.signalDesiredUnit[np.isin(self.wavelength, self.wavelength_filtered),:]
            self.wavelength_useful = self.wavelength_filtered

        # Alternative, subset the CD data based on user given wavelengths    
        else:

            # Conver to list if required
            if not isinstance(user_selected_wavelength,list):

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

            self.signal_useful     = self.signalDesiredUnit[np.isin(self.wavelength, available_wavelength),:]
            # Remove wavelengths that are not in the data
            self.wavelength_useful = available_wavelength[np.isin(available_wavelength, self.wavelength)]

        return None

    def estimate_useful_signal_based_on_snr_and_amplitude(self,measurement_factor):

        # 1) Apply the noise filtering algorithm
        wavelength_filtered = estimate_useful_wavelength_based_on_snr(self.wavelength,self.signalDesiredUnit)

        signal_useful     = self.signalDesiredUnit[np.isin(self.wavelength, wavelength_filtered),:]
    
        # 2) Apply the amplitude filtering algorithm (on the already filtered data!)
        wavelength_filtered = estimate_useful_wavelength_based_on_amplitude(wavelength_filtered,signal_useful,measurement_factor)

        self.wavelength_filtered = wavelength_filtered

        return None

    def estimate_baselines_parameters(self,measurement_factor,baseline_degree_window=8,bootstrap_intervals=False):

        """

        Estimate the temperature / chemical dependence of the native/unfolded state signal.

        Input - an integer 'baseline_degree_window' that is the temperature / chemical window to fit the equation of a line

        Result - we assign self.bNs, self.bUs, self.kNs, self.kUs
    
        kN, bN: slope and intercept of the pre-transition baseline Native State 
        kU, bU: slope and intercept of the post-transition baseline, Unfolded State 

        """

        max_native   = min(measurement_factor) + baseline_degree_window
        min_unfolded = max(measurement_factor) - baseline_degree_window

        native_signal       = filter_matrix_by_vector((self.signal_useful).T,measurement_factor,
            min(measurement_factor),max_native)

        native_measurement_factor      = filter_vector_by_values(measurement_factor,min(measurement_factor),max_native)

        unfolded_signal     = filter_matrix_by_vector((self.signal_useful).T,measurement_factor,
            min_unfolded,max(measurement_factor))

        unfolded_measurement_factor    = filter_vector_by_values(measurement_factor,min_unfolded,max(measurement_factor))

        fitNative   = [poly.Polynomial.fit(native_measurement_factor,  native_signal[:,i],1,
            full=True) for i in range(native_signal.shape[1])]

        fitUnfolded = [poly.Polynomial.fit(unfolded_measurement_factor,  unfolded_signal[:,i],1,
            full=True) for i in range(unfolded_signal.shape[1])]
        
        coefNative   = [x[0].coef for x in fitNative]
        bN           = np.array( [c[0] for c in coefNative] )
        kN           = np.array( [c[1] for c in coefNative] )

        coefUnfolded = [x[0].coef for x in fitUnfolded]
        bU           = np.array( [c[0] for c in coefUnfolded] )
        kU           = np.array( [c[1] for c in coefUnfolded] )

        self.bNs, self.kNs  = bN, kN
        self.bUs, self.kUs  = bU, kU

        # Bootstrap method for errors - not used in production. The idea was to analyse if we need to fit or not the slopes
        if bootstrap_intervals:

            x_data         = native_measurement_factor
            conf_intervals = []

            for i in range(native_signal.shape[1]):

                y_data = native_signal[:,i]

                conf_interval = compute_bootstrap_regression_slope(x_data,y_data)
                conf_intervals.append(conf_interval)

                self.conf_interval = conf_interval
        
            self.baselineNativeIsZero = 0 > self.conf_interval[0] and 0 < self.conf_interval[1]

        return None 

class cd_experiment_comparison(cd_experiment_fitting_model):

    '''
    Advanced class to handle CD spectra with categorical labels. 
    For example, 'control' and 'treatment'.

    We inherit the class 'cd_experiment_fitting_model' to take advantage of the PCA functions (for future developments)

    Some lines are commented with the pattern 'HW_Method'. I implemented the method 
    proposed by Hristova, Kalina, and William C. Wimley. 
    "Determining the statistical significance of the difference between arbitrary curves: A spreadsheet method." 
    Plos one 18.10 (2023): e0289619.

    The code works. However, the results were not satisfactory so it is not used in production.
    To try it, uncomment the lines with that pattern
    '''

    def __init__(self):

        self.signalDesiredUnit     = None  # 2D (size n x m) np array, CD signal matrix, needs to be called 'signalDesiredUnit' to use the PCA functions 
        self.wavelength            = None  # 1D (length n)   np array, wavelength
        self.labels                = None  # 1D (length m)   np array, labels of the CD signal matrix - one label per column

        self.workingUnits          = None  # string, selected working units to create the dataset for comparisons

        self.labels_unique         = None  # 1D (length z)   np array, unique labels
        self.labels_unique_N       = None  # 1D (length z)   np array, number of curves per label

        self.means                 = None  # 2D (size n x z) np array, CD signal average            (per wavelength) of each category
        self.sds                   = None  # 2D (size n x z) np array, CD signal standard deviation (per wavelength) of each category
        #HW_Method self.standard_err          = None  # 2D (size n x z) np array, CD signal standard error     (per wavelength) of each category


        self.distance_matrix       = None  # 2D (size m x m) np array, containing the all versus all comparisons
        self.distances             = None  # list of lists containing the all versus all comparisons. One sublist per comparison
        self.comparison_labels     = None  # list containing the labels of the comparisons.   One element per comparison

                                           # Let 'h' be the total number of comparisons to be done,
        self.difference_spectra     = None # 2D (size n x h) np array,
        self.difference_spectra_sd  = None # 2D (size n x h) np array,
        self.difference_spectra_lbl = None # length h, np array, labels of the difference spectra

        return None

    def summarise_signal_per_label(self):

        '''
        Obtain the average and standard deviation per unique label
        '''

        self.labels_unique   = np.unique(self.labels)
        
        labels_unique_N = []
        means_by_label  = []
        sds_by_label    = []
        #HW_Method ses_by_label    = []

        for label in self.labels_unique:

            sel_idx = self.labels == label

            mean_by_wl = np.mean( self.signalDesiredUnit[:,sel_idx], axis=1 )

            labels_unique_N.append(sum(sel_idx))

            if np.sum(sel_idx) > 1:
                
                sd_by_wl   = np.std(  self.signalDesiredUnit[:,sel_idx], axis=1 , ddof=1) # ddof=1 to use (N - 1) in the denominator
                #HW_Method se_by_wl   = sd_by_wl / np.sqrt(len(sel_idx))

            else:

                sd_by_wl   = np.full(self.signalDesiredUnit.shape[0], np.nan) # maybe np.nan would be better ...
                #HW_Method se_by_wl   = np.full(self.signalDesiredUnit.shape[0], np.nan) # maybe np.nan would be better ...

            means_by_label.append(mean_by_wl)
            sds_by_label.append(sd_by_wl)
            #ses_by_label.append(se_by_wl)

        self.means           = np.array(means_by_label).T
        self.sds             = np.array(sds_by_label).T
        #HW_Method self.standard_err    = np.array(ses_by_label).T
        self.labels_unique_N = np.array(labels_unique_N)

        return None

    def generate_comparison_labels(self):

        '''
        Create the comparison labels, for example, if we have the unique labels 'Treatment' and 'Control'

        Then, the 'comparison_labels' will be 'Treatment versus Treatment', 'Treatment versus Control' and 'Control versus Control'
        '''

        comparison_labels = []

        for i,label1 in enumerate(self.labels_unique):
            
            for label2 in self.labels_unique[i:]:

                comparison_labels.append(label1 + ' versus ' + label2)

        self.comparison_labels = np.array(comparison_labels)

        return None

    def generate_difference_spectra(self):

        '''
        Generate all possible 'difference' spectra. For instance, if we have the unique labels
        'Control' and 'Treatment', we will create the spectrum 'Control - Treatment'.   
        '''

        difference_spectra        = []
        difference_spectra_sd     = []
        difference_spectra_lbl    = []

        # Uncomment the next lines to use test the method proposed by Hristova, Kalina, and William C. Wimley (2023). Plos one
        # Uncomment 
        #HW_Method difference_spectra_se     = []
        #HW_Method difference_spectra_nsmall = []
        #HW_Method difference_dfreedom       = []
        
        for i,label1 in enumerate(self.labels_unique[:-1]):

            for ii,label2 in enumerate(self.labels_unique[(i+1):]):

                difference_spectra.append(self.means[:,i] - self.means[:,(ii+i+1)])
                difference_spectra_sd.append(np.sqrt(self.sds[:,i]**2          + self.sds[:,(ii+i+1)]**2))
                
                difference_spectra_lbl.append(label1 + ' - ' + label2)

                # Uncomment the next lines to use test the method proposed by Hristova, Kalina, and William C. Wimley (2023). Plos one
                #HW_Method difference_spectra_se.append(np.sqrt(self.standard_err[:,i]**2 + self.standard_err[:,(ii+i+1)]**2))
                #HW_Method nsmall = np.min([self.labels_unique_N[i],self.labels_unique_N[i+ii+1]])
                #HW_Method df     = nsmall*2 - 2

                #HW_Method difference_spectra_nsmall.append(nsmall)
                #HW_Method difference_dfreedom.append(df)

        self.difference_spectra        =  np.array(difference_spectra).T 
        self.difference_spectra_sd     =  np.array(difference_spectra_sd).T
        self.difference_spectra_lbl    =  np.array(difference_spectra_lbl)
        
        ''' HW_Method lines

        self.difference_spectra_se     =  np.array(difference_spectra_se).T
        self.difference_spectra_abs    =  np.abs(self.difference_spectra)
        self.difference_spectra_nsmall =  np.array(difference_spectra_nsmall)
        self.difference_dfreedom       =  np.array(difference_dfreedom)

        self.difference_t_value        = self.difference_spectra_abs / self.difference_spectra_se

        difference_probabilities  = []
        difference_chiSquares     = []

        for i,df in enumerate(difference_dfreedom):

            probabilities = 2 * (1 - t.cdf(abs(self.difference_t_value[:,i]), df))
            chiSquare     = chi2.isf(probabilities, 1)

            difference_probabilities.append(probabilities)
            difference_chiSquares.append(chiSquare)

        self.difference_probabilities = np.array(difference_probabilities).T
        self.difference_chiSquares    = np.array(difference_chiSquares).T

        difference_chiSquares_sum      = np.sum(self.difference_chiSquares,axis=0)
        df_for_chiSq                   = self.difference_chiSquares.shape[0]

        right_tailed_probability       = 1 - chi2.cdf(difference_chiSquares_sum, df_for_chiSq)
        ''' 

        return None

    def find_distances(self):

        '''
        For each spectrum, compute the normalised euclidean distance against all other spectra
        
        As a result of running this method, we'll add to the class 

            self.distances  -   list of numpy arrays. Each numpy array contains the computed distances of one category. 
                                For instance, if the unique labels are 'Treatment' and 'Control', 
                                self.distances will have length 3, 
                                one array will have the distances of all the 'Treatment' spectra against the 'Treatment' spectra,   
                                one array will have the distances of all the 'Treatment' spectra against the 'Control'   spectra, and
                                one array will have the distances of all the 'Control'   spectra against the 'Control'   spectra.

            self.distance_matrix - matrix of the normalised euclidean distance for all the spectra
        '''

        n               = len(self.labels)
        distance_matrix = np.full((n, n), np.nan)

        comparison_labels_lst = self.comparison_labels.tolist()

        distances        = [ [] for x in self.comparison_labels]

        for i in range(self.signalDesiredUnit.shape[1] - 1):

            distance_matrix[i,i] = 0 # Fill the diagonal

            label1 = self.labels[i]

            for ii in range(i+1,self.signalDesiredUnit.shape[1]):

                label2   = self.labels[ii]
                match_id = comparison_labels_lst.index(label1 + ' versus ' + label2)

                distance = normalised_euclidean_distance(self.signalDesiredUnit[:,i],self.signalDesiredUnit[:,ii])
                distances[match_id].append(distance)
                
                distance_matrix[i,ii] = distance

        distances            = [ np.array(x) for x in distances]
        self.distances       = distances
        self.distance_matrix = distance_matrix

        return None

class cd_experiment_chemical_unfolding(cd_experiment_fitting_model):

    '''
    Advanced class to analyse chemical unfolding
    '''

    def __init__(self):

        super().__init__()

        self.chem_concentration = None # 1D numpy array. It could be for example, the concentration of Urea, Guanidim chloride or pH

        return None

    def fit_signal(self,fitSlopeNative=True,fitSlopeUnfolded=True): 

        wavelengths = self.wavelength_useful
        signal      = self.signal_useful

        # Initial parameters have to be in order: 
            # Global M, Global D50
            # Single intercepts folded, Single intercepts unfolded
            # Single slopes folded    , Single slopes unfolded

        # Set initial parameters
        p0          = np.concatenate(((1,np.median(self.chem_concentration)),self.bNs,self.bUs,self.kNs,self.kUs))

        low_bound    =  np.array([0  , np.min(self.chem_concentration) + 1 ] + [x/100-10   if x>0 else 100*x-10  for x in p0[2:]])
        high_bound   =  np.array([1e5, np.max(self.chem_concentration) - 1 ] + [100*x+10   if x>0 else x/100+10  for x in p0[2:]]) 

        listOfChemConcentration = [self.chem_concentration for _ in wavelengths]
        listOfSignals           = signal.tolist()
        totalDataSets           = len(listOfChemConcentration)

        if not fitSlopeUnfolded:

            p0, low_bound, high_bound = p0[:-totalDataSets], low_bound[:-totalDataSets], high_bound[:-totalDataSets]

        if not fitSlopeNative:

            l1 = (2+totalDataSets*2)
            l2 = (2+totalDataSets*3)

            indices_to_delete = [x for x in range(l1,l2)]

            p0           = np.delete(p0, indices_to_delete) 
            low_bound    = np.delete(low_bound, indices_to_delete) 
            high_bound   = np.delete(high_bound, indices_to_delete)

        global_fit_params, cov = fit_chemical_unfolding(listOfChemConcentration,listOfSignals,self.temperature,
            p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        ## Try to re fit the slopes and/or baselines if the fitted parameters are close to the boundaries

        diffMin = (global_fit_params - low_bound)[2:]
        diffMax = (high_bound        - global_fit_params)[2:]

        cond1 = np.any(diffMin < 1)
        cond2 = np.any(diffMax < 1)

        if cond1:

            for i,value in enumerate(diffMin):

                if value < 1:

                    low_bound[i+2] = low_bound[i+2] - 500

        if cond2:

            for i,value in enumerate(diffMax):

                if value < 1:

                    high_bound[i+2] = high_bound[i+2] + 500

        if cond1 or cond2:

            p0 = global_fit_params

            global_fit_params, cov = fit_chemical_unfolding(listOfChemConcentration,listOfSignals,self.temperature,
                p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        errors                 = np.sqrt(np.diag(cov))

        M             = global_fit_params[0] # M
        D50           = global_fit_params[1] # Concentration at which the sample is 50 % unfolded

        interceptsFolded    = global_fit_params[2:(2+totalDataSets)]
        interceptsUnfolded  = global_fit_params[(2+totalDataSets):(2+totalDataSets*2)]

        slopesFolded   = global_fit_params[(2+totalDataSets*2):(2+totalDataSets*3)] if fitSlopeNative   else np.full(len(interceptsFolded), 0)
        slopesUnfolded = global_fit_params[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), 0)

        M_error            = errors[0] # m parameter
        d50_error          = errors[1] # concentration of denaturant at which the sample is 50 % unfolded

        interceptsFoldedError    = errors[2:(2+totalDataSets)]
        interceptsUnfoldedError  = errors[(2+totalDataSets):(2+totalDataSets*2)]

        slopesFoldedError   = errors[(2+totalDataSets*2):(2+totalDataSets*3)] if fitSlopeNative   else np.full(len(interceptsFolded), np.nan)
        slopesUnfoldedError = errors[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), np.nan)

        predicted = []

        # To be filled in this order 
        fit_params = []
        fit_errors = []

        # Fill the parameters dataframe
        i = 0

        for bN, kN, bU, kU in zip(interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded):

            Y = chemical_unfolding_with_linear_dependence_one_curve(self.chem_concentration,self.temperature,D50,M,bN,kN,bU,kU)

            predicted.append(Y)

            fit_params.append([kN,bN,kU,bU,D50,M, str(wavelengths[i]) + ' ' + self.name])

            i +=1

        # Fill the errors dataframe
        i = 0

        for bNe, kNe, bUe, kUe in zip(interceptsFoldedError,slopesFoldedError,interceptsFoldedError,slopesUnfoldedError):

            fit_errors.append([kNe,bNe,kUe,bUe,d50_error,M_error, str(wavelengths[i]) + ' ' + self.name])

            i +=1

        fit_params = np.array(fit_params)
        fit_errors = np.array(fit_errors)

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'D50', 'M', 'Condition']

        # Convert the NumPy array to a Pandas DataFrame with column names
        df_fit_params = pd.DataFrame(fit_params, columns=column_names)
        df_fit_errors = pd.DataFrame(fit_errors, columns=column_names)

        # Remove non necessary columns
        if not fitSlopeNative:

            df_fit_params = df_fit_params.drop('kN', axis=1)
            df_fit_errors = df_fit_errors.drop('kN', axis=1)
            column_names.remove('kN')

        if not fitSlopeUnfolded:

            df_fit_params = df_fit_params.drop('kU', axis=1)
            df_fit_errors = df_fit_errors.drop('kU', axis=1)
            column_names.remove('kU')

        # Convert to numeric (all columns except 'Condition')
        df_fit_params[ column_names[:-1] ] = df_fit_params[ column_names[:-1] ].apply(pd.to_numeric)
        df_fit_errors[ column_names[:-1] ] = df_fit_errors[ column_names[:-1] ].apply(pd.to_numeric)

        # Compute relative errors (in percentage)
        for col in column_names[:-1]:

            df_fit_errors[col] = ( df_fit_errors[col]  / df_fit_params[col] * 100 ).abs()

        predicted = np.array(predicted)
        
        self.signal_predicted      = np.array(predicted)
        self.fit_params            = df_fit_params
        self.fit_rel_errors        = df_fit_errors

        self.fractions  = chem_two_state_rev_unfolding_fractions(self.temperature,self.chem_concentration,M,D50)

        return None

    def fit_signal_three_state(self,fitSlopeNative=True,fitSlopeUnfolded=True,D50v1_init=0,D50v2_init=0): 

        wavelengths = self.wavelength_useful
        signal      = self.signal_useful

        # Initial parameters have to be in order: 
            # Global M, Global D50
            # Single intercepts folded, Single intercepts unfolded
            # Single slopes folded    , Single slopes unfolded

        listOfChemConcentration = [self.chem_concentration for _ in wavelengths]
        listOfSignals           = signal.tolist()
        totalDataSets           = len(listOfChemConcentration)

        bIs   = (self.bNs + self.kUs) * 0.5

        minC  = np.min(self.chem_concentration)
        maxC  = np.max(self.chem_concentration)
        medC  = np.median(self.chem_concentration)

        # Set initial parameters
        p0          = np.concatenate(((0.25,minC+1,0.25,maxC-1),self.bNs,self.bUs,bIs,self.kNs,self.kUs))

        low_bound    =  np.array([0.2  , minC + 0.5,0.2  , minC + 0.5 ] + [x/100-100   if x>0 else 100*x-100  for x in p0[4:]])
        high_bound   =  np.array([20, maxC - 2,20, maxC - 0.5 ]         + [100*x+100   if x>0 else x/100+100  for x in p0[4:]]) 

        # If required, change the initial guess and fitting limits
        if D50v1_init != 0:

            p0[1], low_bound[1], high_bound[1] = D50v1_init, D50v1_init - 1.75, D50v1_init + 1.75

        if D50v2_init != 0:

            p0[3], low_bound[3], high_bound[3] = D50v2_init, D50v2_init - 1.75, D50v2_init + 1.75

        if not fitSlopeUnfolded:

            p0, low_bound, high_bound = p0[:-totalDataSets], low_bound[:-totalDataSets], high_bound[:-totalDataSets]

        if not fitSlopeNative:

            l1 = (4+totalDataSets*3)
            l2 = (4+totalDataSets*4)

            indices_to_delete = [x for x in range(l1,l2)]

            p0           = np.delete(p0, indices_to_delete) 
            low_bound    = np.delete(low_bound, indices_to_delete) 
            high_bound   = np.delete(high_bound, indices_to_delete)

        global_fit_params, cov = fit_chemical_unfolding_three_species(listOfChemConcentration,listOfSignals,self.temperature,
            p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        ## Try to re fit the slopes and/or baselines if the fitted parameters are close to the boundaries

        diffMin = (global_fit_params - low_bound)[4:]
        diffMax = (high_bound        - global_fit_params)[4:]

        cond1 = np.any(diffMin < 1)
        cond2 = np.any(diffMax < 1)

        if cond1:

            for i,value in enumerate(diffMin):

                if value < 1:

                    low_bound[i+4] = low_bound[i+4] - 1e3

        if cond2:

            for i,value in enumerate(diffMax):

                if value < 1:

                    high_bound[i+4] = high_bound[i+4] + 1e3

        if cond1 or cond2:

            p0 = global_fit_params

            global_fit_params, cov = fit_chemical_unfolding_three_species(listOfChemConcentration,listOfSignals,self.temperature,
                p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

            ## Try a second re fit the slopes and/or baselines if the fitted parameters are still close to the boundaries

            diffMin = (global_fit_params - low_bound)[4:]
            diffMax = (high_bound        - global_fit_params)[4:]

            cond1 = np.any(diffMin < 10)
            cond2 = np.any(diffMax < 10)

            if cond1:

                for i,value in enumerate(diffMin):

                    if value < 10:

                        low_bound[i+4] = low_bound[i+4] - 1e5

            if cond2:

                for i,value in enumerate(diffMax):

                    if value < 10:

                        high_bound[i+4] = high_bound[i+4] + 1e5

            if cond1 or cond2:

                p0 = global_fit_params

                global_fit_params, cov = fit_chemical_unfolding_three_species(listOfChemConcentration,listOfSignals,self.temperature,
                    p0,low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        errors                 = np.sqrt(np.diag(cov))

        M1             = global_fit_params[0] # M
        D50v1          = global_fit_params[1] # Concentration at which the sample is 50 % unfolded

        M2             = global_fit_params[2] # M
        D50v2          = global_fit_params[3] # Concentration at which the sample is 50 % unfolded

        interceptsFolded       = global_fit_params[(4+totalDataSets*0):(4+totalDataSets*1)]
        interceptsUnfolded     = global_fit_params[(4+totalDataSets*1):(4+totalDataSets*2)]
        interceptsIntermediate = global_fit_params[(4+totalDataSets*2):(4+totalDataSets*3)]

        slopesFolded   = global_fit_params[(4+totalDataSets*3):(4+totalDataSets*4)] if fitSlopeNative   else np.full(len(interceptsFolded), 0)
        slopesUnfolded = global_fit_params[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), 0)

        M1_error            = errors[0] # m parameter
        d50v1_error         = errors[1] # concentration of denaturant at which the sample is 50 % unfolded

        M2_error            = errors[2] # m parameter
        d50v2_error         = errors[3] # concentration of denaturant at which the sample is 50 % unfolded

        interceptsFoldedError    = errors[(4+totalDataSets*0):(4+totalDataSets*1)]
        interceptsUnfoldedError  = errors[(4+totalDataSets*1):(4+totalDataSets*2)]
        interceptsIntermError    = errors[(4+totalDataSets*2):(4+totalDataSets*3)]

        slopesFoldedError   = errors[(4+totalDataSets*3):(4+totalDataSets*4)] if fitSlopeNative   else np.full(len(interceptsFolded), np.nan)
        slopesUnfoldedError = errors[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), np.nan)

        predicted = []

        # To be filled in this order 
        fit_params = []
        fit_errors = []

        # Fill the parameters dataframe
        i = 0

        for bN, kN, bU, kU, bI in zip(interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded,interceptsIntermediate):

            Y = chemical_unfolding_with_linear_dependence_one_curve_three_species(self.chem_concentration,self.temperature,D50v1,M1,D50v2,M2,bN,kN,bU,kU,bI)

            predicted.append(Y)

            fit_params.append([kN,bN,kU,bU,bI,D50v1,M1,D50v2,M2, str(wavelengths[i]) + ' Dataset ' + self.name])

            i +=1

        # Fill the errors dataframe
        i = 0

        for bNe, kNe, bUe, kUe, bIe in zip(interceptsFoldedError,slopesFoldedError,interceptsFoldedError,slopesUnfoldedError,interceptsIntermError):

            fit_errors.append([kNe,bNe,kUe,bUe,bIe,d50v1_error,M1_error,d50v2_error,M2_error, str(wavelengths[i]) + ' Dataset ' + self.name])

            i +=1

        fit_params = np.array(fit_params)
        fit_errors = np.array(fit_errors)

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU','bI', 'D50_step1', 'M1','D50_step2', 'M2', 'Condition']

        # Convert the NumPy array to a Pandas DataFrame with column names
        df_fit_params = pd.DataFrame(fit_params, columns=column_names)
        df_fit_errors = pd.DataFrame(fit_errors, columns=column_names)

        # Remove non necessary columns
        if not fitSlopeNative:

            df_fit_params = df_fit_params.drop('kN', axis=1)
            df_fit_errors = df_fit_errors.drop('kN', axis=1)
            column_names.remove('kN')

        if not fitSlopeUnfolded:

            df_fit_params = df_fit_params.drop('kU', axis=1)
            df_fit_errors = df_fit_errors.drop('kU', axis=1)
            column_names.remove('kU')

        # Convert to numeric (all columns except 'Condition')
        df_fit_params[ column_names[:-1] ] = df_fit_params[ column_names[:-1] ].apply(pd.to_numeric)
        df_fit_errors[ column_names[:-1] ] = df_fit_errors[ column_names[:-1] ].apply(pd.to_numeric)

        # Compute relative errors (in percentage)
        for col in column_names[:-1]:

            df_fit_errors[col] = ( df_fit_errors[col]  / df_fit_params[col] * 100 ).abs()

        predicted = np.array(predicted)
        
        self.signal_predicted      = np.array(predicted)
        self.fit_params            = df_fit_params
        self.fit_rel_errors        = df_fit_errors

        self.fractions  = chem_three_state_rev_unfolding_fractions(self.temperature,self.chem_concentration,M1,D50v1,M2,D50v2)

        return None

class cd_experiment_thermal_ramp(cd_experiment_fitting_model):

    def __init__(self):

        super().__init__()       # Load attributes from the parent class
        
        return None
        
    def estimate_signal_derivates(self,temp_window_length=5,dt=0.2):

        """
        Compute the 1st and 2nd derivative of self.signal using the Savitzky-Golay-Filter
        Estimate the melting temperature based on the derivative peaks (maximums or minimums)

        Input:
                - An integer 'temp_window_length' that is the window_length for the Savitzky-Golay-Filter
                - A float 'dt' to perform a linear interpolation along the temperature axis

        Result:
                - We assign self.derivative and self.tms_from_deriv

        """

        temperature_even_sequence = np.arange(min(self.temperature), max(self.temperature), dt)
        
        # Calculate the number of rows in self.signal_useful
        num_rows = self.signal_useful.shape[0]

        # Initialize an empty array to store the interpolated values for each row
        # interpolated_sequences will have dimensions x*m where x is the length of the even_sequence and m
        #   is the number of wavelengths data were measured

        # self.temperature (x-coordinates) should be in increasing order!!

        interpolated_sequences = np.empty((len(temperature_even_sequence),self.signal_useful.shape[0]))

        # Loop through each row and apply interpolation
        for row in range(num_rows):
            interpolated_sequence = np.interp(temperature_even_sequence, self.temperature, self.signal_useful[row, :], left=None, right=None, period=None)
            # Transpose results
            interpolated_sequences[:,row] = interpolated_sequence

        odd_n_data_points_window_len         = np.ceil(temp_window_length / dt) // 2 * 2 + 1

        self.derivative = savgol_filter(interpolated_sequences,axis=0,
            window_length=odd_n_data_points_window_len,polyorder=4,
            deriv=1,mode="nearest")

        """
        Estimate Tm from the derivative curve
        To do this, we first shift the first derivative by using the mean of the 
            median of the first and last 6 degrees of each melting curve.
        """

        der_temp_init = filter_matrix_by_vector(self.derivative,
            temperature_even_sequence,min(temperature_even_sequence)+6,min(temperature_even_sequence)+11)

        der_temp_end  = filter_matrix_by_vector(self.derivative,
            temperature_even_sequence,max(temperature_even_sequence)-11,max(temperature_even_sequence)-6)

        med_init  = np.median(der_temp_init,axis=0)
        med_end   = np.median(der_temp_end ,axis=0)   
        mid_value = np.array([(x+y)/2 for x,y in zip(med_init,med_end)])

        der_temp = filter_matrix_by_vector(self.derivative,
            temperature_even_sequence,min(temperature_even_sequence)+6,max(temperature_even_sequence)-6)

        mid_value = mid_value * np.where(mid_value > 0,1,-1)

        der_temp = np.add(der_temp,mid_value) 

        temp_temp = filter_vector_by_values(temperature_even_sequence,min(temperature_even_sequence)+6,max(temperature_even_sequence)-6)

        max_der = np.amax(der_temp,axis=0)
        min_der = np.amin(der_temp,axis=0)

        der_direction_temp = [abs(maxd) > abs(mind) for maxd,mind in zip(max_der,min_der)]

        if  sum(der_direction_temp) > (len(der_direction_temp) / 2):

            self.tms_from_deriv = get_temp_at_maximum_of_derivative(temp_temp,der_temp)

        else:

            self.tms_from_deriv = get_temp_at_minimum_of_derivative(temp_temp,der_temp)
        
        self.temperature_even_sequence = temperature_even_sequence

        self.global_Tm = np.mean(self.tms_from_deriv)

        return None

    def fit_signal_three_state(self,fitSlopeNative=False,fitSlopeUnfolded=False,T1_init=0,T2_init=0): 

        wavelengths = self.wavelength_useful
        signal      = self.signal_useful
    
        # Initial parameters have to be in order: 
            # Global melting temperature 1, Global enthalpy of unfolding 1, 
            # Global melting temperature 2, Global enthalpy of unfolding 2, 
            # Single intercepts folded, Single slopes folded, 
            # Single intercepts unfolded, Single slopes unfolded

        # Set initial parameters

        listOfTemperatures = [self.temperature for _ in wavelengths]
        listOfSignals      = signal.tolist()
        totalDataSets      = len(listOfTemperatures)

        bIs   = (self.bNs + self.kUs) * 0.5

        minT  = np.min(self.temperature)
        maxT  = np.max(self.temperature)
        meanT = np.median(self.temperature)

        p0           = np.concatenate( ( ( minT+10,10,maxT-20,10),self.bNs,self.bUs,bIs,self.kNs,self.kUs))

        low_bound    =  [minT+4,1,minT+4,1 ]     + [x/100-100 if x>0 else 100*x-100  for x in p0[4:]]
        high_bound   =  [maxT-7,250,maxT,250]    + [100*x+100 if x>0 else x/100+100  for x in p0[4:]] 

        if T1_init != 0:

            p0[0], low_bound[0], high_bound[0] = T1_init, T1_init - 15, T1_init + 15

        if T2_init != 0:

            p0[2], low_bound[2], high_bound[2] = T2_init, T2_init - 15, T2_init + 15

        if not fitSlopeUnfolded:

            p0, low_bound, high_bound = p0[:-totalDataSets], low_bound[:-totalDataSets], high_bound[:-totalDataSets]

        if not fitSlopeNative:

            l1 = (4+totalDataSets*3)
            l2 = (4+totalDataSets*4)

            indices_to_delete = [x for x in range(l1,l2)]

            p0           = np.delete(p0, indices_to_delete) 
            low_bound    = np.delete(low_bound, indices_to_delete) 
            high_bound   = np.delete(high_bound, indices_to_delete)

        global_fit_params, cov = fit_thermal_unfolding_three_species(listOfTemperatures,listOfSignals,p0,
            low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        ## Re fit the slopes and/or baselines if the fitted parameters are close to the boundaries

        diffMin = (global_fit_params - low_bound)[4:]
        diffMax = (high_bound        - global_fit_params)[4:]

        cond1 = np.any(diffMin < 1)
        cond2 = np.any(diffMax < 1)

        if cond1:

            for i,value in enumerate(diffMin):

                if value < 1:

                    low_bound[i+4] = low_bound[i+4] - 500

        if cond2:

            for i,value in enumerate(diffMax):

                if value < 1:

                    high_bound[i+4] = high_bound[i+4] + 500

        if cond1 or cond2:

            p0 = global_fit_params

            global_fit_params, cov = fit_thermal_unfolding_three_species(listOfTemperatures,listOfSignals,p0,
                low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        errors                 = np.sqrt(np.diag(cov))

        T1            = global_fit_params[0] # Temperature of melting
        DH1           = global_fit_params[1] # Enthalpy of unfolding

        T2            = global_fit_params[2] # Temperature of melting
        DH2           = global_fit_params[3] # Enthalpy of unfolding

        interceptsFolded       = global_fit_params[(4+totalDataSets*0):(4+totalDataSets*1)]
        interceptsUnfolded     = global_fit_params[(4+totalDataSets*1):(4+totalDataSets*2)]
        interceptsIntermediate = global_fit_params[(4+totalDataSets*2):(4+totalDataSets*3)]

        slopesFolded   = global_fit_params[(4+totalDataSets*3):(4+totalDataSets*4)] if fitSlopeNative   else np.full(len(interceptsFolded), 0)
        slopesUnfolded = global_fit_params[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), 0)

        T1_error            = errors[0] # Temperature of melting
        DH1_error           = errors[1] # Enthalpy of unfolding

        T2_error            = errors[2] # Temperature of melting
        DH2_error           = errors[3] # Enthalpy of unfolding

        interceptsFoldedError    = errors[(4+totalDataSets*0):(4+totalDataSets*1)]
        interceptsUnfoldedError  = errors[(4+totalDataSets*1):(4+totalDataSets*2)]
        interceptsIntermError    = errors[(4+totalDataSets*2):(4+totalDataSets*3)]

        slopesFoldedError   = errors[(4+totalDataSets*3):(4+totalDataSets*4)] if fitSlopeNative   else np.full(len(interceptsFolded), np.nan)
        slopesUnfoldedError = errors[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), np.nan)

        predicted = []

        # To be filled in this order - kN bN kU bU dHm Tm Condition
        fit_params = []
        fit_errors = []

        # Fill the parameters dataframe
        i = 0

        for bN, kN, bU, kU, bI in zip(interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded,interceptsIntermediate):

            Y = thermal_unfolding_one_curve_three_species(self.temperature,T1,DH1,T2,DH2,bN,kN,bU,kU,bI)
            predicted.append(Y)

            fit_params.append([kN,bN,kU,bU,bI,DH1,T1,DH2,T2, str(wavelengths[i]) + '. Dataset: ' + self.name])

            i +=1

        # Fill the errors dataframe
        i = 0

        for bNe, kNe, bUe, kUe, bIe in zip(interceptsFoldedError,slopesFoldedError,interceptsUnfoldedError,slopesUnfoldedError,interceptsIntermError):

            fit_errors.append([kNe,bNe,kUe,bUe,bIe,DH1_error,T1_error,DH2_error,T2_error, str(wavelengths[i]) + '. Dataset: ' + self.name])

            i +=1

        fit_params = np.array(fit_params)
        fit_errors = np.array(fit_errors)

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU','bI','DH1', 'T1','DH2', 'T2', 'Condition']

        # Convert the NumPy array to a Pandas DataFrame with column names
        df_fit_params = pd.DataFrame(fit_params, columns=column_names)
        df_fit_errors = pd.DataFrame(fit_errors, columns=column_names)

        # Remove non necessary columns
        if not fitSlopeNative:

            df_fit_params = df_fit_params.drop('kN', axis=1)
            df_fit_errors = df_fit_errors.drop('kN', axis=1)
            column_names.remove('kN')

        if not fitSlopeUnfolded:

            df_fit_params = df_fit_params.drop('kU', axis=1)
            df_fit_errors = df_fit_errors.drop('kU', axis=1)
            column_names.remove('kU')

        # Convert to numeric (all columns except 'Condition')
        df_fit_params[ column_names[:-1] ] = df_fit_params[ column_names[:-1] ].apply(pd.to_numeric)
        df_fit_errors[ column_names[:-1] ] = df_fit_errors[ column_names[:-1] ].apply(pd.to_numeric)

        # Compute relative errors (in percentage)
        for col in column_names[:-1]:

            df_fit_errors[col] = ( df_fit_errors[col]  / df_fit_params[col] * 100 ).abs()

        self.signal_predicted      = np.array(predicted)
        self.fit_params            = df_fit_params
        self.fit_rel_errors        = df_fit_errors

        self.fractions             = three_state_rev_unfolding_fractions(self.temperature,DH1,DH2,T1,T2)

        return None

    def fit_signal(self,fitSlopeNative=True,fitSlopeUnfolded=True): 

        wavelengths = self.wavelength_useful
        signal      = self.signal_useful
    
        # Initial parameters have to be in order: 
            # Global melting temperature, Global enthalpy of unfolding, 
            # Single intercepts folded, Single slopes folded, 
            # Single intercepts unfolded, Single slopes unfolded

        # Set initial parameters
        p0          = np.concatenate(((self.global_Tm,100),self.bNs,self.bUs,self.kNs,self.kUs))

        low_bound    =  [min(self.temperature)+7,10 ] + [x/100-10 if x>0 else 100*x-10  for x in p0[2:]]
        high_bound   =  [max(self.temperature)-6,250] + [100*x+10 if x>0 else x/100+10  for x in p0[2:]] 

        listOfTemperatures = [self.temperature for _ in wavelengths]
        listOfSignals      = signal.tolist()
        totalDataSets = len(listOfTemperatures)

        if not fitSlopeUnfolded:

            p0, low_bound, high_bound = p0[:-totalDataSets], low_bound[:-totalDataSets], high_bound[:-totalDataSets]

        if not fitSlopeNative:

            l1 = (2+totalDataSets*2)
            l2 = (2+totalDataSets*3)

            indices_to_delete = [x for x in range(l1,l2)]

            p0           = np.delete(p0, indices_to_delete) 
            low_bound    = np.delete(low_bound, indices_to_delete) 
            high_bound   = np.delete(high_bound, indices_to_delete)

        global_fit_params, cov = fit_thermal_unfolding(listOfTemperatures,listOfSignals,p0,
            low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        ## Re fit the slopes and/or baselines if the fitted parameters are close to the boundaries

        diffMin = (global_fit_params - low_bound)[2:]
        diffMax = (high_bound        - global_fit_params)[2:]

        cond1 = np.any(diffMin < 1)
        cond2 = np.any(diffMax < 1)

        if cond1:

            for i,value in enumerate(diffMin):

                if value < 1:

                    low_bound[i+2] = low_bound[i+2] - 500

        if cond2:

            for i,value in enumerate(diffMax):

                if value < 1:

                    high_bound[i+2] = high_bound[i+2] + 500

        if cond1 or cond2:

            p0 = global_fit_params

            global_fit_params, cov = fit_thermal_unfolding(listOfTemperatures,listOfSignals,p0,
                low_bound,high_bound,fitSlopeNative,fitSlopeUnfolded)

        errors                 = np.sqrt(np.diag(cov))

        Tm            = global_fit_params[0] # Temperature of melting
        dh            = global_fit_params[1] # Enthalpy of unfolding

        interceptsFolded    = global_fit_params[2:(2+totalDataSets)]
        interceptsUnfolded  = global_fit_params[(2+totalDataSets):(2+totalDataSets*2)]

        slopesFolded   = global_fit_params[(2+totalDataSets*2):(2+totalDataSets*3)] if fitSlopeNative   else np.full(len(interceptsFolded), 0)
        slopesUnfolded = global_fit_params[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), 0)

        Tm_error            = errors[0] # Temperature of melting
        dh_error            = errors[1] # Enthalpy of unfolding

        interceptsFoldedError    = errors[2:(2+totalDataSets)]
        interceptsUnfoldedError  = errors[(2+totalDataSets):(2+totalDataSets*2)]

        slopesFoldedError   = errors[(2+totalDataSets*2):(2+totalDataSets*3)] if fitSlopeNative   else np.full(len(interceptsFolded), np.nan)
        slopesUnfoldedError = errors[(len(global_fit_params)-totalDataSets):] if fitSlopeUnfolded else np.full(len(interceptsFolded), np.nan)

        predicted = []

        # To be filled in this order - kN bN kU bU dHm Tm Condition
        fit_params = []
        fit_errors = []

        # Fill the parameters dataframe
        i = 0

        for bN, kN, bU, kU in zip(interceptsFolded,slopesFolded,interceptsUnfolded,slopesUnfolded):

            Y = two_state_thermal_unfolding_one_curve(self.temperature,Tm,dh,bN,kN,bU,kU,0)
            predicted.append(Y)

            fit_params.append([kN,bN,kU,bU,dh,Tm, str(wavelengths[i]) + ' ' + self.name])

            i +=1

        # Fill the errors dataframe
        i = 0

        for bNe, kNe, bUe, kUe in zip(interceptsFoldedError,slopesFoldedError,interceptsUnfoldedError,slopesUnfoldedError):

            fit_errors.append([kNe,bNe,kUe,bUe,dh_error,Tm_error, str(wavelengths[i]) + ' ' + self.name])

            i +=1

        fit_params = np.array(fit_params)
        fit_errors = np.array(fit_errors)

        # Column names
        column_names = ['kN', 'bN', 'kU', 'bU', 'dH', 'Tm', 'Condition']

        # Convert the NumPy array to a Pandas DataFrame with column names
        df_fit_params = pd.DataFrame(fit_params, columns=column_names)
        df_fit_errors = pd.DataFrame(fit_errors, columns=column_names)

        # Remove non necessary columns
        if not fitSlopeNative:

            df_fit_params = df_fit_params.drop('kN', axis=1)
            df_fit_errors = df_fit_errors.drop('kN', axis=1)
            column_names.remove('kN')

        if not fitSlopeUnfolded:

            df_fit_params = df_fit_params.drop('kU', axis=1)
            df_fit_errors = df_fit_errors.drop('kU', axis=1)
            column_names.remove('kU')

        # Convert to numeric (all columns except 'Condition')
        df_fit_params[ column_names[:-1] ] = df_fit_params[ column_names[:-1] ].apply(pd.to_numeric)
        df_fit_errors[ column_names[:-1] ] = df_fit_errors[ column_names[:-1] ].apply(pd.to_numeric)

        # Compute relative errors (in percentage)
        for col in column_names[:-1]:

            df_fit_errors[col] = ( df_fit_errors[col]  / df_fit_params[col] * 100 ).abs()

        self.signal_predicted      = np.array(predicted)
        self.fit_params            = df_fit_params
        self.fit_rel_errors        = df_fit_errors

        self.fractions             = two_state_rev_unfolding_fractions(self.temperature,dh,Tm)

        return None

class cd_experiment_custom_analysis(cd_experiment_fitting_model):

    '''
    Advanced class to analyse the CD data using custom models. 
    For example, binding, or thermochemical unfolding
    '''

    def __init__(self):

        super().__init__()       # Load attributes from the parent class
        self.model_function   = None

        # Numpy array containing a certain measurement dimension
        self.first_measurement_dimension  = None # E.g., temperature or ligand concentration
        self.second_measurement_dimension = None 

        # experimental parameters names
        self.first_exp_param_name         = None
        self.second_exp_param_name        = None

        return None

    def generate_model_function(self,text):

        cleaned_text = clean_function_text(text)

        self.cleaned_text_function = cleaned_text

        fit_params_names = extract_words(cleaned_text)

        fit_params_names = [ x for x in fit_params_names if x not in ['sqrt','log','np','exp',self.first_exp_param_name,self.second_exp_param_name] ]

        # remove duplicates and sort it
        fit_params_names = sorted(list(set(fit_params_names)))

        # find global and local parameters
        global_params_names  = [ x for x in fit_params_names if 'global'     in x.lower() ]
        local_params_names   = [ x for x in fit_params_names if 'global' not in x.lower() ]
        
        local_params_names_extended   = [ ]

        for wl in self.wavelength_useful:

            local_params_names_extended += [str(int(wl)) + ' ' + x for x in local_params_names]

        paramsSign           = []

        for param in global_params_names + local_params_names:

            if 'Pos' in param:

                paramsSign.append('positive')

            elif 'Neg' in param:

                paramsSign.append('negative')

            else:

                paramsSign.append('unknown')

        self.paramsSign = paramsSign

        wavelengths  = self.wavelength_useful

        self.global_params_names         = global_params_names
        self.local_params_names          = local_params_names
        self.local_params_names_extended = local_params_names_extended

        self.n_global_params = len(global_params_names)
        self.n_local_params  = len(local_params_names)

        # Create a Python function using eval
        def custom_function(dummyVariable,*args):

            """
            args should be in this order: 
                global param 1, ... , global param N, 
                local param 1 (first wavelength), local param 1 (second wavelength), 
                ...  local param Z (last wavelength).

            The total number of arguments should equal 

                #global_params + (#local_params * #wavelengths)

            """

            signal = []

            dInit = {'np'              : np,
            self.first_exp_param_name  : self.first_measurement_dimension,
            self.second_exp_param_name : self.second_measurement_dimension
            }

            # initialize the i counter, just in case there are no global params
            i = -1
            for i,fit_param in enumerate(global_params_names):

                dInit[fit_param] = args[i]

            counter = 0

            for _ in wavelengths:

                d = dInit.copy()

                for fit_param in local_params_names:
                
                    d[fit_param] = args[i + counter + 1]

                    counter += 1

                signal.append(eval(cleaned_text, d))

            signal = np.array(signal).flatten()

            return signal

        self.model_function = custom_function

        return None

    def get_initial_values(self,leftLimitLogSearch=-3,rightLimitLogSearch=3):

        signal = self.signal_useful
        wls    = self.wavelength_useful

        nParamsPerCurve = len(self.local_params_names) + len(self.global_params_names)

        # Generate a list of list containing a
        # combination of log-spaced values
        combis = generate_initial_params_combinations(nParamsPerCurve,self.paramsSign,
            leftLimitLogSearch,rightLimitLogSearch)
    
        # Try to find acceptable initial values by running the user defined function with different parameter values

        dCombInit = {'np'          : np,
        self.first_exp_param_name  : self.first_measurement_dimension,
        self.second_exp_param_name : self.second_measurement_dimension
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

                predicted.append(eval(self.cleaned_text_function,dComb))
        
        all_wl_rss = []
        sel_combi  = []    

        for i,wl in enumerate(wls):

            wl_rss = []

            signalSlice = signal[i,:]

            for prediction in predicted:

                if np.isnan(prediction).any():
                    wl_rss.append(np.inf)
                
                else:

                    residuals = signalSlice - prediction
                    wl_rss.append(np.sum(residuals**2))

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

        tol_factor   = 500

        low_bound    =  [x/tol_factor  if x>0 else tol_factor*x  for x in p0]
        high_bound   =  [tol_factor*x  if x>0 else x/tol_factor  for x in p0]

        self.p0         = round_to_significant_digits(p0,2)
        self.low_bound  = round_to_significant_digits(low_bound,2)
        self.high_bound = round_to_significant_digits(high_bound,2)

        return None

    def fit_signal(self):

        """
        Find the fitting parameters by using the function generated through self.generate_model_function
        """

        # Reset the fitting params and the predicted signal
        self.fit_params         = None
        self.fit_rel_errors     = None
        self.signal_predicted   = None
        self.boundaries_updated = False

        if self.p0 is None:

            self.get_initial_values()

        signal = self.signal_useful
        wls    = self.wavelength_useful

        c1 = np.all(self.p0 <= self.high_bound)
        c2 = np.all(self.p0 >= self.low_bound)

        if not c1 or not c2:

            return None

        try:

            all_fitted_params, cov = curve_fit(self.model_function,1,signal.flatten(),
                p0 = self.p0,
                bounds=(tuple(self.low_bound), tuple(self.high_bound)))

            # Try to re-fit if the values are close to the boundaries +/- 500

            diffMin = (all_fitted_params - self.low_bound)
            diffMax = (self.high_bound   - all_fitted_params)

            cond1 = np.any(diffMin < 1)
            cond2 = np.any(diffMax < 1)

            if cond1:

                for i,value in enumerate(diffMin):

                    if value < 1:

                        self.low_bound[i] = self.low_bound[i] - 500

            if cond2:

                for i,value in enumerate(diffMax):

                    if value < 1:

                        self.high_bound[i] = self.high_bound[i] + 500

            if cond1 or cond2:

                self.boundaries_updated = True

                self.p0 = all_fitted_params

                all_fitted_params, cov = curve_fit(self.model_function,1,signal.flatten(),
                    p0 = self.p0,
                    bounds=(tuple(self.low_bound), tuple(self.high_bound)))

            self.all_fitted_params = all_fitted_params
            errors                 = np.sqrt(np.diag(cov))

            predicted      = self.model_function(None,*all_fitted_params)
            predicted      = predicted.reshape((self.signal_useful.shape[0], -1))
            
            self.signal_predicted = predicted

            # initialize the dataframe that will have the fitted parameters values
            local_params_df  = pd.DataFrame({'Condition': wls})
            global_params_df = pd.DataFrame({'Condition': wls})

            local_params_df_error  = pd.DataFrame({'Condition': wls})
            global_params_df_error = pd.DataFrame({'Condition': wls})

            if self.n_local_params > 0:

                column_names    = self.local_params_names
                local_params1    = all_fitted_params[self.n_global_params:]
                local_params2    = local_params1.reshape((-1, self.n_local_params)) 
                local_params3    = pd.DataFrame(local_params2, columns=column_names)
                local_params_df  = pd.concat([local_params3, local_params_df], axis=1)

                local_params_err1     = errors[self.n_global_params:]
                local_params_err2     = local_params_err1.reshape((-1, self.n_local_params)) 
                local_params_err3     = np.abs(local_params_err2 / local_params2 * 100)
                local_params_err4     = pd.DataFrame(local_params_err3, columns=column_names)
                local_params_df_error = pd.concat([local_params_err4, local_params_df_error], axis=1)

            if self.n_global_params > 0:
                 
                for i,param in enumerate(all_fitted_params[:self.n_global_params]):
                                               
                    global_params_df[self.global_params_names[i]] = param 

                    rel_err = np.abs(errors[i] / param * 100)

                    global_params_df_error[self.global_params_names[i]] = rel_err

            all_params_df     = pd.merge(global_params_df,local_params_df)
            all_params_df_err = pd.merge(global_params_df_error,local_params_df_error)

            self.fit_params     = all_params_df
            self.fit_rel_errors = all_params_df_err

        except:

            pass
         
        return None

class cdAnalyzer:

    """
    Useful to work with many different cd experiments
    """

    def __init__(self):

        """
        Create 
        """

        self.experimentsOri           = {}  # Dictionary where each key value pair corresponds to one CD experiment 
        self.experimentsModif         = {}  # To subset the wavelength range of the CD experiments inside experimentsOri  
        self.experimentNames          = []
        self.sharedParameters         = False # To modify all the experiment at the same time with shared parameters

        self.experimentsThermal       = {} # Dictionary where each key value pair corresponds to one CD based thermal unfolding experiment
        self.experimentNamesThermal   = []

        self.experimentsChemical      = {} # Dictionary where each key value pair corresponds to one CD based chemical unfolding experiment
        self.experimentNamesChemical  = []

        self.experimentsCustom      = {} # Dictionary where each key value pair corresponds to one CD experiment with a custom formula for analysis
        self.experimentNamesCustom  = []

        return None

    def loadExperiment(self,file,name):

        """
        Append one experiment to experimentsOri 
        """

        if name in self.experimentNames:

            return "Experiment name already selected!"

        try:

            self.experimentsOri[name] = cd_experiment_general()         
            self.experimentsOri[name].load_data(file,name)

            self.experimentNames.append(name)

            # Adapt the number of significant digits

            return "Data loaded successfully!!!"

        except:

            pass
        
        return "Data could not be loaded"

    def deleteExperiment(self,names):

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

    def clean_experiments(self,experiment_type):

        if experiment_type == 'custom':

            self.experimentsCustom      = {} 
            self.experimentNamesCustom  = []

        if experiment_type == 'thermal':

            self.experimentsThermal       = {} 
            self.experimentNamesThermal   = []

        if experiment_type == 'chemical':

            self.experimentsChemical      = {} 
            self.experimentNamesChemical  = []

        return None

    def experiments2absorbanceUnits(self,unitsDic):

        """

        Sets the proper units for each experiment in experimentsOri
            e.g. of unitsList - {'exp' : 'milliabsorbance','blank' : 'degrees', ... }

        """

        lst_of_experiments_with_new_units = []

        for experimentName , experimentUnits in unitsDic.items():

            if not self.experimentsOri[experimentName].isFakeExperiment:

                if experimentUnits != self.experimentsOri[experimentName].units:
                    lst_of_experiments_with_new_units.append(experimentName + ' : ' + experimentUnits)

                self.experimentsOri[experimentName].experiment2absorbanceUnits(experimentUnits)

        # Return which experiments were modified
        return lst_of_experiments_with_new_units

    def allExperiments2absorbanceUnits(self,experimentUnits):

        for experimentName in self.experimentNames:

            if not self.experimentsOri[experimentName].isFakeExperiment:

                self.experimentsOri[experimentName].experiment2absorbanceUnits(experimentUnits)

        return None

    def experimentsAbsorbanceUnits2otherUnits(self,unitsDic):

        """

        Transform absorbance units to desired units
            e.g. of unitsList - {'exp' : 'milliabsorbance','blank' : 'degrees', ... }

        This function must be always called after running setExperimentUnits !!!

        """

        for experimentName , experimentUnits in unitsDic.items():

            if not self.experimentsOri[experimentName].isFakeExperiment:

                self.experimentsOri[experimentName].experimentFromAbsorbanceUnits2otherUnits(experimentUnits)

            else:

                self.experimentsOri[experimentName].signalDesiredUnit = self.experimentsOri[experimentName].signalInput

        return None

    def allExperimentsAbsorbanceUnits2otherUnits(self,experimentUnits):

        for experimentName in self.experimentNames:

            if not self.experimentsOri[experimentName].isFakeExperiment:

                self.experimentsOri[experimentName].experimentFromAbsorbanceUnits2otherUnits(experimentUnits)

            else:

                self.experimentsOri[experimentName].signalDesiredUnit = self.experimentsOri[experimentName].signalInput

        return None

    def initializeExperimentModif(self):

        """

        Required to store the original data and allow post-processing

        """

        self.experimentsModif = copy.deepcopy(self.experimentsOri)

        return None

    def setExperimentProperties(self,experimentName,variable,value):

        """
        experimentName must be in self.experimentNames
        variable can be 'replicates', 'reads', or 'scans'
        value is a number
        """

        setattr(self.experimentsOri[experimentName], variable, value)
        setattr(self.experimentsModif[experimentName], variable, value)

        return None

    def getExperimentProperties(self,variable):

        """
        variable can be 'replicates', 'reads', or 'scans'
        """

        return [getattr(self.experimentsOri[experimentName], variable) for experimentName in self.experimentNames]

    def getExperimentPropertiesModif(self,variable):

        """
        variable can be 'replicates', 'reads', or 'scans'
        """

        return [getattr(self.experimentsModif[experimentName], variable) for experimentName in self.experimentNames]

    def filterDataByWavelength(self,minWL,maxWL):

        self.initializeExperimentModif()

        experimentNames = self.experimentNames

        for expName in experimentNames:

            # Use the original signal! - Allows going back
            wl            = self.experimentsOri[expName].wavelength
            signalDesired = self.experimentsOri[expName].signalDesiredUnit
            signalInput   = self.experimentsOri[expName].signalInput
            signalAbs     = self.experimentsOri[expName].signalAbs
            signalHT      = self.experimentsOri[expName].signalHT

            if not np.isnan(signalDesired).all():
                self.experimentsModif[expName].signalDesiredUnit = filter_matrix_by_vector(signalDesired,wl,minWL,maxWL)
            else:
                self.experimentsModif[expName].signalDesiredUnit = np.nan

            if not np.isnan(signalInput).all():
                self.experimentsModif[expName].signalInput = filter_matrix_by_vector(signalInput,wl,minWL,maxWL)
            else:
                self.experimentsModif[expName].signalInput = np.nan

            if not np.isnan(signalAbs).all():
                self.experimentsModif[expName].signalAbs = filter_matrix_by_vector(signalAbs,wl,minWL,maxWL)
            else:
                self.experimentsModif[expName].signalAbs = np.nan

            self.experimentsModif[expName].signalHT          = filter_matrix_by_vector(signalHT,wl,minWL,maxWL)
            self.experimentsModif[expName].wavelength        = filter_vector_by_values(wl,minWL,maxWL)

        return None


if False:
    
    t1 = cd_experiment_thermal_ramp()
    t1.load_data('/home/os/Downloads/zero_new.csv','t')
    t1.signalDesiredUnit = t1.signalInput
    t1.assign_useful_signal([195,222,200])
    t1.estimate_signal_derivates()
    t1.estimate_baselines_parameters(t1.temperature)

    #t1.fit_signal_three_state(True,True)
    #t1.fit_signal_three_state(True,False)
    t1.fit_signal_three_state(False,True)
    t1.fit_signal_three_state(False,False)

if False:
    
    t1 = cd_experiment_custom_analysis()
    t1.load_data('/home/os/Downloads/zero_new.csv','t') 

    t1.temperature = np.array([5.7, 10.3, 15,19.6,24.2,28.9,33.5,38.2,42.8,47.4,52.1,56.7,61.4,66,70.6,75.3,79.9,84.5])  + 273.15

    t1.first_measurement_dimension = t1.temperature 
    t1.first_exp_param_name         = 'T'

    t1.signalDesiredUnit = t1.signalInput
    t1.assign_useful_signal([183,193,203,213,223,233])
    t1.generate_model_function('(bn+kOne*e^(-deltaHGlobalAPos/(0.001987*T)*(1-T/TGlobalAPos))+(bu)*e^(-deltaHGlobalAPos/(0.001987*T)*(1-T/TGlobalAPos))*e^(-deltaHGlobalBPos/(0.001987*T)*(1-T/TGlobalBPos)))/(1+e^(-deltaHGlobalAPos/(0.001987*T)*(1-T/TGlobalAPos))*e^(-deltaHGlobalBPos/(0.001987*T)*(1-T/TGlobalBPos)))')

    p0  = np.repeat([1],18)
    lb  = np.repeat([-500],18)
    hb  = np.repeat([500],18)

    t1.p0         = np.concatenate((np.array([305,330,50,50]), p0), axis=0)
    t1.low_bound  = np.concatenate((np.array([295,320,2,2]), lb), axis=0)
    t1.high_bound = np.concatenate((np.array([320,380,100,100]), hb), axis=0)

    t1.fit_signal()

    tm1 = t1.fit_params.iloc[0,1]
    tm2 = t1.fit_params.iloc[0,2]
    dh1 = t1.fit_params.iloc[0,3]
    dh2 = t1.fit_params.iloc[0,4]

    R = 1.987 / 1000 # kcal/mol

    A = np.exp(-dh1*(1-t1.temperature/tm1)/(R*t1.temperature))
    B = np.exp(-dh2*(1-t1.temperature/tm2)/(R*t1.temperature))

    xN = 1 / (1+A+A*B)
    xi = A / (1+A+A*B)
    xD = A*B / (1+A+A*B)

    rss = []

    for t in np.arange(305,320):

            t1.p0         = np.concatenate((np.array([t,330,20,20]), p0), axis=0)
            t1.low_bound  = np.concatenate((np.array([t,320,1,1]), lb), axis=0)
            t1.high_bound = np.concatenate((np.array([t,380,100,100]), hb), axis=0)

            rss.append(np.sum(np.square(t1.signal_predicted - t1.signal_useful)))

    print(t1.fit_params)
    print(t1.fit_rel_errors)
    

    matrix = t1.signal_predicted
    import matplotlib.pyplot as plt

    num_columns = matrix.shape[0]

    # Plot each column of the matrix against the vector
    for i in range(num_columns):
        # Extract the i-th column from the matrix
        column_data1 = matrix[i, :]
        column_data2 = t1.signal_useful[i, :]

        # Plot the column against the vector
        plt.plot(t1.temperature, column_data1, label=f'Column {i+1}')
        plt.scatter(t1.temperature, column_data2, label=f'Column {i+1}')

    # Add labels and title
    plt.xlabel('Vector')
    plt.ylabel('Column Values')
    plt.title('Plot of Each Column of a Matrix Against a Vector')

    # Add a legend
    plt.legend()

    # Display the plot
    plt.show()

    df = pd.DataFrame({'xN':xN,'xI':xi,'xD':xD,'T':t1.temperature})

        # Plot each column against the specified x column
    ax = df.plot(x='T', y=['xN', 'xI', 'xD'], kind='line')

    # Set labels and title
    ax.set_xlabel('X Column')
    ax.set_ylabel('Values')

    # Add a legend to display the labels for the vertical lines
    ax.legend()

    # Show the plot
    plt.show()