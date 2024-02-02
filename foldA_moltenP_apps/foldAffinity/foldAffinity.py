#!/usr/bin/python3
import copy
import json

import numpy  as np
import pandas as pd

import scipy.stats
import scipy.optimize

from scipy.signal     import savgol_filter
from scipy.optimize   import curve_fit

from helpers                                  import *
from fitting_helpers_unfolded_fraction        import *
from fitting_helpers_thermal_unfolding        import *

class DSF_binding:

    '''
    Class for binding analysis from fluorescence melting curves
    This class was first written by Stephan Niebling and 
    	later modified by Osvaldo Burastero to adapt it to a shiny server 
    It is based on isothermal analysis by Bai et al., Sci. Rep. 2019.

    No warranty whatsoever
    If you have questions please contact 
    stephan.niebling@embl-hamburg.de or oburastero@gmail.com

    References
    Bai, N.; Roder, H.; Dickson, A. & Karanicolas, J. Scientific Reports, 9, 2019.
    Niebling, Stephan, et al.  Scientific Reports 11.1, 2021.
    '''

    def __init__(self):

        # Heat capacity of unfolding
        self.cp    = 0

        # Protein concentration
        self.pconc = None

        # List of the possible signals to be analysed, e.g., [330 nm, 350nm]
        self.signals = None

        # Selected signal 
        self.signal_type = None

        # Array of the conditions names
        self.conditions_original    = None

        # Boolean array to decide which positions should be used
        self.include_vector = None

        # Dictionaries where the keys are the different signals and the values have the signal and temperature data
        self.signal_data_dictionary = {}
        self.temp_data_dictionary   = {}

        # Values of the selected signal
        self.fluo                   = None
        # Temperature values 
        self.temps                  = None

        # Derivative values of the selected signal
        self.derivative             = None

        # Maximum of the first derivative
        self.max_derivative         = None

        # Tm from the first derivative
        self.tmsFromDerivative = None

        # Fitted parameters of the melting curves - local, global or global_cp options
        self.local_fit_params     = None
        self.local_fit_errors     = None
        self.global_fit_params    = None
        self.global_fit_errors    = None
        self.global_fit_cp_params = None
        self.global_fit_cp_errors = None

        #  Fitted parameters of the melting curves - one from the local, global or global_cp options   
        self.fit_fluo_params =  None
        self.fit_fluo_errs   =  None

        # Fine conc grid for plotting
        self.kd_model_conc   = None
        # Fitted binding curves
        self.kd_models       = None
        
        # Lower limits based on fit errors (covariance matrix)
        self.kd_models_lower = None

        # Upper limits based on fit errors (covariance matrix)
        self.kd_models_upper = None

        # Fitted signal
        self.fit_fluo_pred   =  None  

        # Limits of the asymmetric confidence intervals
        self.bind_ci95_asymmetric_low = None
        self.bind_ci95_asymmetric_up  = None

        # Fitted parameters and errors of the isothermal curves 
        self.bind_params          = None
        self.bind_errors          = None
        
        self.isothermal_data = None
        self.unique_concs    = None
        self.isothermal_fits = None
        self.concentrations  = None
        self.isothermal_ts   = None

        # Fitted data and parameters from an alternative model - Tm shift fitting - done in the R shiny app only!!!
        self.tms_fit_pred             = None
        self.tms_fit_info             = None
        self.tms_fit_asymmetricCI95   = None
        self.tms_fit_model            = None

        return None

    def load_JSON_file(self,JSON_file):

        # Read JSON data from a file
        with open(JSON_file, "r") as file:
            json_data = file.read()

        # Parse JSON data into a dictionary
        data_dict     = json.loads(json_data)
        data_dict_new = {}

        # Decode file
        for attribute, value in data_dict.items():
            
            if 'comment' in attribute:
                continue

            if isinstance(value, dict) :

                data_dict_new[attribute] = {}

                for subAttribute, subValue in value.items():
                        
                    if isinstance(subValue,list) and isinstance(subValue[0], str) and ';' in subValue[0]:
    
                        subValue = [x.split(';') for x in subValue]

                        data_dict_new[attribute][subAttribute] = np.array(subValue, dtype=float)

                    else:

                        data_dict_new[attribute][subAttribute] = np.array(subValue)

            elif isinstance(value, list) and isinstance(value[0], str) and ';' in value[0]:

                value      = [x.split(';') for x in value]

                data_dict_new[attribute] = np.array(value, dtype=float)

            elif isinstance(value, list):

                data_dict_new[attribute] = np.array(value)

            else:

                # Leave as it is
                data_dict_new[attribute] = value 

        for key, value in data_dict_new.items():
            setattr(self, key, value)

        return None

    def export_JSON_file(self,JSON_file):

        my_object_dict = vars(self)
        commentDictionary = {'_comment' : 'Use this JSON file to load the session in FoldAffinity. \
To parse the data, consider that \
all 2D arrays where converted into lists with sublists, and each sublist is now a string with \
elements separated by semicolons.'}

        my_object_dict = {**commentDictionary,**my_object_dict}
        
        json_data = json.dumps(my_object_dict,cls=NumpyEncoder,indent=2)

        # Write JSON data to a file
        with open(JSON_file, "w") as file:
            file.write(json_data)

        return None

    def load_example_data(self):

        fluo  = np.load('testFluo.npy')    
        temps = np.load('testTemp.npy')
    
        self.signals        = np.array(["350nm"])

        self.signal_data_dictionary["350nm"]   = fluo
        self.temp_data_dictionary["350nm"]     = temps

        self.conditions_original = np.loadtxt("testConcs.txt")

        return None

    def load_nanoDSF_xlsx(self, processed_dsf_file,sheet_names=["350nm","330nm","Scattering","Ratio"]):

        """
        Load nanotemper processed xlsx file and extract data
            
        """

        # Load excel file
        # this needs to be the processed file!
        xls = pd.ExcelFile(processed_dsf_file)
        
        conditions_df = pd.read_excel(xls, "Overview")
        conditions    = conditions_df[['Sample ID']].values.flatten()

        self.conditions_original = conditions
        # Add position index to avoid problems with empty names

        possible_signals    = ["350nm","330nm","Scattering","Ratio"]

        include             = []

        # Change signal name if unfolding curve is present
        for sn in sheet_names:

            include_value = any([ps in sn and "deriv" not in sn.lower() and "fold" not in sn.lower() for ps in possible_signals])

            include.append(include_value)

        sheet_names_to_load = np.array([s for (i, s) in zip(include, sheet_names) if i])

        self.signals        = np.array([" ".join(sn.split()) for sn in sheet_names_to_load])

        for sn, signal in zip(sheet_names_to_load, self.signals):

            dat = pd.read_excel(xls, sn, index_col=None, header=None)

            # Find first index with numbers (above is header)
            first_row = int(np.argwhere(list(dat.iloc[:, 0] == 'Time [s]'))) + 1

            fluo   = np.array(dat.iloc[first_row:, 2:]).astype('float')
            temp   = np.array(dat.iloc[first_row:, 1]).astype('float') 

            if (len(temp)) > 900:
                fluo = fluo[::4]
                temp = temp[::4]

            while (len(temp)) > 400:

                fluo = fluo[::2]
                temp = temp[::2]

            self.signal_data_dictionary[signal]   = fluo
            self.temp_data_dictionary[signal]     = temp

        return None

    def load_tycho_xlsx(self,file):
        
        xls = pd.ExcelFile(file)
        
        # Retrieve the conditions names - in sheet "Results"
        df = pd.read_excel(xls, "Results")

        for colIndex, column in enumerate(df):

            columnValues      = df[column]
            has_desired_value = columnValues.isin(["Capillary label"])

            if any(has_desired_value):

                rowIndexBegin   = np.flatnonzero(has_desired_value)[0]+1
                initialRatioCol = df.iloc[rowIndexBegin:,colIndex+4]
                try:
                    rowIndexEnd     = np.flatnonzero(initialRatioCol.isnull())[0] + rowIndexBegin
                except:
                    rowIndexEnd     = rowIndexBegin+6
                break

        conditions         = np.array(df.iloc[rowIndexBegin:rowIndexEnd,colIndex])
        numberOfConditions = len(conditions)

        self.conditions_original = conditions

        # Retrieve the fluorescence signal - in sheet "Profiles_raw"
        df = pd.read_excel(xls, "Profiles_raw")

        for colIndex, column in enumerate(df):

            columnValues      = df[column]
            has_desired_value = columnValues.isin(["Temperature [°C]"])

            if any(has_desired_value):

                rowIndexBegin   = np.flatnonzero(has_desired_value)[0]+1
                rowIndexEnd     = np.flatnonzero(columnValues[(rowIndexBegin):].isnull())[0]+rowIndexBegin

                temperature     = np.array(df.iloc[rowIndexBegin:rowIndexEnd,colIndex]).astype('float')

                break

        columnNames  = df.iloc[rowIndexBegin-1,:]

        startIndexes = [i for i, j in enumerate(columnNames) if j == conditions[0]]

        signals      = np.array(df.iloc[rowIndexBegin-3,startIndexes])
        signalsClean = []
        
        for s in signals:
            if "ratio" in s.lower():
                signalsClean.append("Ratio")
            else:  
                if "330" in s.lower() :
                    signalsClean.append("330nm")
                if "350" in s.lower():
                    signalsClean.append("350nm")

        for i, signal in enumerate(signalsClean):

            colStartIndex   = startIndexes[i]
            colEndIndex     = colStartIndex + numberOfConditions

            fluo            = np.array(df.iloc[rowIndexBegin:rowIndexEnd,colStartIndex:colEndIndex]).astype('float')
            temp            = temperature

            while (len(temp)) > 400:

                fluo = fluo[::2]
                temp = temp[::2]

            self.signal_data_dictionary[signal]   = fluo
            self.temp_data_dictionary[signal]     = temp

        self.signals        = np.array(signalsClean)

        return None

    def load_Thermofluor_xlsx(self,thermofluor_file):

        """

        Load DSF Thermofluor file and extract data

        """

        xls            = pd.ExcelFile(thermofluor_file)
        dat            = pd.read_excel(xls, "RFU",header=None)
        conditions     = np.array(dat.iloc[0, 1:])

        self.conditions_original   = conditions

        signal = "DSF_RFU"

        self.signal_data_dictionary[signal]   = np.array(dat.iloc[1:,1:]).astype('float')
        self.temp_data_dictionary[signal]     = np.array(dat.iloc[1:, 0]).astype('float') 

        self.signals = np.array([signal])

        return None

    def load_panta_xlsx(self,pantaFile):

        sheet_names   = get_sheet_names_of_xlsx(pantaFile)

        try:

            data          = pd.read_excel(pantaFile, "Data Export")

        except:

            data          = pd.read_excel(pantaFile, "melting-scan") # Alternative format of PANTA

        column_names  = [str.lower(c) for c in data.columns]

        pos_350       = [i for i,x in enumerate(column_names) if "350"   in x and "deriv" not in x and "330" not in x]
        pos_330       = [i for i,x in enumerate(column_names) if "330"   in x and "deriv" not in x and "350" not in x]
        scattering    = [i for i,x in enumerate(column_names) if "scattering"   in x and "deriv" not in x]
        pos_ratio     = [i for i,x in enumerate(column_names) if "ratio" in x and "deriv" not in x]

        possible_signals    = ["350nm","330nm","Scattering","Ratio"]
        signals             = []

        all_positions = [pos_350,pos_330,scattering,pos_ratio]

        for positions,signal in zip(all_positions,possible_signals):

            if len(positions) > 0:

                fluo, temp = subset_panta_data(data,positions)

                self.signal_data_dictionary[signal]   = fluo
                self.temp_data_dictionary[signal]     = temp

                signals.append(signal)

        self.signals        = np.array(signals)

        try:

            conditions_df = pd.read_excel(pantaFile, "Overview")
            conditions    = conditions_df[['Sample ID']].values.flatten().astype(str)

        except:

            conditions = np.repeat(0,fluo.shape[1])


        self.conditions_original = conditions

        
        return None

    def load_QuantStudio_txt(self,QSfile):

        """

        Input: A txt file ('QSfile') where column 2 has the well position, 
        column 3 the temperature and column 4 the fluorescence signal. Index starts at 1!!!

        --- Caution ---
        The first rows of the file are comments that start with '*' and are not readed
        The temperature and signal column have commas that need to be deleted

        Columns are separated by spaces
        """

        start_row = getStartLineQuantStudioFile(QSfile)
        data      = pd.read_csv(QSfile,skiprows=start_row,sep="\s+",header=None)

        u, ind     = np.unique(data.iloc[:,1], return_index=True)
        conditions = u[np.argsort(ind)]

        self.conditions_original = conditions

        fluo , temp = generateMergedQuantStudioDataFrame(data)

        signal = "Fluorescence"

        self.signal_data_dictionary[signal]   = fluo
        self.temp_data_dictionary[signal]     = temp
        self.signals = np.array([signal])

        return None

    def load_Agilents_MX3005P_qPCR_txt(self,filename):

        """
        Input: A txt file where the 2nd column has the fluorescence data, and 
        the 3rd column the temperature. 
        Wells are separated by rows containing a sentence like this one: 'Segment  2 Plateau  1 Well  1' 
        """

        dfs      = []
        well_num = []

        with open(filename, 'r') as f:

            ls  = f.read().splitlines()

            for i,line in enumerate(ls):
                if line.startswith('Segment') and 'Well' in line:
                    # Get the well number
                    well_num .append(line.split()[-1])

                    fluorescence = []
                    temperature  = []

                    for line2 in ls[i+2:]:
                        if line2.startswith('Segment'):
                            break
                        else:
                            data = line2.split()
                            fluorescence.append(float(data[1]))
                            temperature.append(float(data[2]))

                    df = pd.DataFrame({'temperature':temperature,'signal'+str(well_num):fluorescence})
                    df.sort_values('temperature',inplace=True)
                    dfs.append(df)

        # Combine dataframes so we can obtain a vector of temperatures and a matrix of fluorescence signal
        merged = pd.merge_asof(dfs[0].dropna(),dfs[1].dropna(),on="temperature", direction='nearest',allow_exact_matches=True)

        for df in dfs[2:]:
            
            merged = pd.merge_asof(merged,df.dropna(),on="temperature", direction='nearest',allow_exact_matches=True)       

        fluo   = np.array(merged.iloc[:, 1:]).astype('float')
        temp   = np.array(merged.iloc[:, 0]).astype('float') 

        # Reduce data so we can plot and fit the data faster.
        while len(temp) > 700:
            fluo = fluo[::2]
            temp = temp[::2]

        signal = "Fluorescence"

        self.conditions_original            = well_num
        self.signal_data_dictionary[signal] = fluo
        self.temp_data_dictionary[signal]   = temp
        self.signals = np.array([signal])

        return None

    def load_csv_file(self,file):

        """

        Input: A csv file where the first column has the temperature and all the next columns the 
        fluorescence data, header is required

        """

        dat            = pd.read_csv(file)

        conditions  = [str(c) for c in dat.columns[1:]]
        
        self.conditions_original   = conditions

        signal = "Fluorescence"

        self.signal_data_dictionary[signal]   = np.array(dat.iloc[:,1:]).astype('float')
        self.temp_data_dictionary[signal]     = np.array(dat.iloc[:, 0]).astype('float')  
        self.signals = np.array([signal])

        return None

    def set_signal(self,which):

        """

        Assign self.fluo, self.temps according to the desired signal

        """

        self.fluo   = self.signal_data_dictionary[which]
        self.temps  = self.temp_data_dictionary[which]
        
        self.signal_type = which
        self.set_dt()

        # Move the selectes signal to the first position - useful to export the JSON
        self.signals = self.signals[self.signals != which]
        self.signals = np.insert(self.signals, 0, which)

        return None

    def set_dt(self):

        self.dt         = ( max(self.temps) - min(self.temps) ) / (len(self.temps) - 1)

        return None

    def median_filter(self,n_degree_window):

        """

        Use this function if the fluorescence curves present spikes. It applies a
        rolling median window filter 

        Requires:

        An integer 'n_degree_window'


        """

        self.fluo = np.apply_along_axis(median_filter_from_fluo_and_temp_vectors,0,
            self.fluo,self.temps,n_degree_window)

        return None

    def estimate_fluo_derivates(self,temp_window_length):

        odd_n_data_points_window_len = np.ceil(temp_window_length / self.dt) // 2 * 2 + 1

        self.derivative = savgol_filter(self.fluo,axis=0,
            window_length=odd_n_data_points_window_len,polyorder=4,
            deriv=1,mode="nearest")

        """
        
        Estimate Tm from the derivative curve
        To do this, we first shift the first derivative by using the mean of the 
            median of the first and last 5 degrees of each melting curve.

        """

        der_temp_init = filter_fluo_by_temp(self.derivative,
            self.temps,min(self.temps)+6,min(self.temps)+11)

        der_temp_end  = filter_fluo_by_temp(self.derivative,
            self.temps,max(self.temps)-11,max(self.temps)-6)

        med_init  = np.median(der_temp_init,axis=0)
        med_end   = np.median(der_temp_end ,axis=0)   
        mid_value = np.array([(x+y)/2 for x,y in zip(med_init,med_end)])

        der_temp = filter_fluo_by_temp(self.derivative,
            self.temps,min(self.temps)+1,max(self.temps)-1)

        mid_value = mid_value * np.where(mid_value>0,1,-1)

        der_temp = np.add(der_temp,mid_value) 

        temp_temp = filter_temp_by_temp(self.temps,min(self.temps)+1,max(self.temps)-1)

        max_der = np.amax(der_temp,axis=0)
        min_der = np.amin(der_temp,axis=0)

        der_direction_temp = [abs(maxd) > abs(mind) for maxd,mind in zip(max_der,min_der)]

        if  sum(der_direction_temp) > (len(der_direction_temp) / 2):

            self.tms = get_temp_at_maximum_of_derivative(temp_temp,der_temp)

        else:

            self.tms = get_temp_at_minimum_of_derivative(temp_temp,der_temp)

        self.tmsFromDerivative = self.tms

        return None

    def fit_fluo_local(self):
        '''
        This is the first step of fitting the melting curves
        It uses the deltaCp specified before or otherwise 
        assumes a deltaCp of 0
        
        '''

        num_datasets     = len(self.concentrations)
        local_fit_params = np.empty([7, num_datasets])
        pcovs = []
        errors = np.empty([7, num_datasets])
           
        for index in range(num_datasets):
            fit_params, err = fit_single_thermal_curve(self.temps, self.fluo[:, index], 
                self.cp)

            for index2 in range(6):
            # Fill parameters
                local_fit_params[index2, index] = fit_params[index2]
                errors[index2, index]           = err[index2]

            local_fit_params[6, index] = self.cp
            errors[6, index]           = np.nan # Dummy entry

        self.local_fit_params = local_fit_params
        self.local_fit_errors = errors
        self.fitting_option   = "Local"

        return None

    def fit_fluo_global(self):
        '''
        This is the second step of fitting the melting curves:
        Now the thermal curves are fitted with the global slopes
        All other parameters are local
        This function requires the local fit as input (starting values)
        '''
        # Get existing parameters
        old_params   =  self.local_fit_params
        num_datasets = len(self.concentrations)

        # Initial parameters have to be in order: Tms, dHs, unfolded_intercepts, folded_intercepts, unfolded slope, folded_slope
        Tms = old_params[0, :]
        dHs = old_params[1, :]
        unfolded_intercepts = old_params[2, :]
        folded_intercepts = old_params[3, :]
        # For the shared slopes we could also use the average
        # For now will stick to what was done in the original code
        folded_slope   = np.mean(old_params[4, :])
        unfolded_slope = np.mean(old_params[5, :])
        last_temp = len(self.temps) - 1
        # Instead of looping through the datasets, we will concatenate all data
        temp_concat, fluo_concat = concatenate_fluorescence(self.temps, self.fluo)

        # Starting values for fit
        p0 = (*Tms, *dHs, *unfolded_intercepts, *folded_intercepts, folded_slope, unfolded_slope)

        # Specify limits for fit
        low_bound = [-np.inf] * (4 * num_datasets + 2)
        high_bound = [np.inf] * (4 * num_datasets + 2)
        # Specify lower limits for Tms
        low_bound[:num_datasets]  = [np.min(self.temps)] * num_datasets
        high_bound[:num_datasets] = [np.max(self.temps)] * num_datasets
        # Specify lower limits for dH
        low_bound[num_datasets:(2 * num_datasets)] = [0] * num_datasets

        # Do global fit
        params, cov = curve_fit(global_thermal_curves_woCp(self.cp), temp_concat, fluo_concat, p0=p0, bounds=(low_bound, high_bound), max_nfev=1E3)
        # Now bring parameter matrix in shape needed by program
        fit_params = np.copy(old_params)
        fit_params[0, :] = params[:num_datasets]  # Tm
        fit_params[1, :] = params[num_datasets:2 * num_datasets]  # dH
        fit_params[2, :] = params[2 * num_datasets:3 * num_datasets]  # unfolded intercepts
        fit_params[3, :] = params[3 * num_datasets:4 * num_datasets]  # folded intercepts
        fit_params[4, :] = np.tile(params[4 * num_datasets], num_datasets)
        fit_params[5, :] = np.tile(params[4 * num_datasets + 1], num_datasets)
        fit_params[6, :] = np.tile(self.cp, num_datasets) # Cp
        # Same for errors
        errors = np.sqrt(np.diag(cov))
        fit_errors = np.copy(old_params) * np.nan
        fit_errors[0, :] = errors[:num_datasets]  # Tm
        fit_errors[1, :] = errors[num_datasets:2 * num_datasets]  # dH
        fit_errors[2, :] = errors[2 * num_datasets:3 * num_datasets]  # unfolded intercepts
        fit_errors[3, :] = errors[3 * num_datasets:4 * num_datasets]  # folded intercepts
        fit_errors[4, :] = np.tile(errors[4 * num_datasets], num_datasets)
        fit_errors[5, :] = np.tile(errors[4 * num_datasets + 1], num_datasets)
        fit_errors[6, :] = np.tile(np.nan, num_datasets) # Dummy value for Cp
        self.global_fit_params = fit_params
        self.global_fit_errors = fit_errors
        # Update melting temperatures
        self.tms = fit_params[0]

        self.fitting_option = "Global"

        return None

    def fit_fluo_global_cp(self): # , correct_local_fit=True
        '''

        '''

        # Get existing parameters
        old_params = self.global_fit_params

        num_datasets = len(self.concentrations)

        # Initial parameters have to be in order: Tms, dHs, unfolded_intercepts, folded_intercepts, unfolded slope, folded_slope
        Tms = old_params[0, :]
        dHs = old_params[1, :]
        unfolded_intercepts = old_params[2, :]
        folded_intercepts = old_params[3, :]
        # For the shared slopes we use the average
        folded_slope = np.mean(old_params[4, :])
        unfolded_slope = np.mean(old_params[5, :])
        last_temp = len(self.temps) - 1
        # Instead of looping through the datasets, we will concatenate all data
        temp_concat, fluo_concat = concatenate_fluorescence(self.temps, self.fluo)

        # Starting values for fit
        p0 = (*Tms, *dHs, *unfolded_intercepts, *folded_intercepts, folded_slope, unfolded_slope, self.cp)

        # Specify limits for fit
        low_bound = [-np.inf] * (4 * num_datasets + 3)
        high_bound = [np.inf] * (4 * num_datasets + 3)
        # Specify lower limits for Tms
        low_bound[:num_datasets] = [np.min(self.temps)] * num_datasets
        high_bound[:num_datasets] = [np.max(self.temps)] * num_datasets
        # Specify lower limits for dH
        low_bound[num_datasets:(2 * num_datasets)] = [0] * num_datasets
        # Specify lower limit for deltaCP

        try:
            global_fit_params, cov = curve_fit(global_thermal_curves, temp_concat, fluo_concat, 
                p0=p0,bounds=(low_bound, high_bound),max_nfev=500)

            # Do global fit
            # Now bring parameter matrix in shape needed by program
            fit_params = np.zeros((7, num_datasets)) * np.nan  # np.copy(old_params)
            fit_params[0, :] = global_fit_params[:num_datasets]  # Tm
            fit_params[1, :] = global_fit_params[num_datasets:2 * num_datasets]  # dH
            fit_params[2, :] = global_fit_params[2 * num_datasets:3 * num_datasets]  # unfolded intercepts
            fit_params[3, :] = global_fit_params[3 * num_datasets:4 * num_datasets]  # folded intercepts
            fit_params[4, :] = np.tile(global_fit_params[4 * num_datasets], num_datasets)
            fit_params[5, :] = np.tile(global_fit_params[4 * num_datasets + 1], num_datasets)
            fit_params[6, :] = np.tile(global_fit_params[4 * num_datasets + 2], num_datasets)
            # Same for errors
            errors = np.sqrt(np.diag(cov))
            fit_errors = np.copy(fit_params) * np.nan
            fit_errors[0, :] = errors[:num_datasets]  # Tm
            fit_errors[1, :] = errors[num_datasets:2 * num_datasets]  # dH
            fit_errors[2, :] = errors[2 * num_datasets:3 * num_datasets]  # unfolded intercepts
            fit_errors[3, :] = errors[3 * num_datasets:4 * num_datasets]  # folded intercepts
            fit_errors[4, :] = np.tile(errors[4 * num_datasets], num_datasets)
            fit_errors[5, :] = np.tile(errors[4 * num_datasets + 1], num_datasets)
            fit_errors[6, :] = np.tile(errors[4 * num_datasets + 2], num_datasets)
            # Send to instance
            self.global_fit_cp_params = fit_params
            self.global_fit_cp_errors = fit_errors
            # Update deltaCp
            self.cp = fit_params[6, 0]
            # Update melting temperatures
            self.tms = fit_params[0]

            self.fitting_option = "Global_CP"

        except RuntimeError:
            self.self.fitting_option = "Global"

        return None

    def fit_isothermal(self,model_type):  
        '''
        This is an alternative Kd fit directly using
        to fit the fus.


        Output
        It outputs a fit with error from the covariance matrix
        self.bind_params: Ku and Kd (latter in M)
        self.bind_errors: Ku and Kd errors
        '''

        # Initialize isothermal data and fits
        self.isothermal_data = np.empty([len(self.concentrations), len(self.isothermal_ts)])
        self.unique_concs = np.unique(self.concentrations)
        self.isothermal_fits = np.empty([len(self.unique_concs), len(self.isothermal_ts)])
        bind_params, bind_errors, bind_ci95_asymmetric_low, bind_ci95_asymmetric_up = [], [], [], []
        # Select fitting parameters
        fitting_params, fitting_errors, fit_melting = self.select_fitting_params()

        # Loop through temperatures and fit Ku and Kd
        for i, binding_temp in enumerate(self.isothermal_ts):
            # Select fu
            self.isothermal_data[:, i] = calculate_fraction_unfolded(binding_temp,
                                                                          self.concentrations,
                                                                          fitting_params,
                                                                          self.cp).transpose()
            fu = self.isothermal_data[:, i]
            # In case there are nans in fu: remove them
            if np.sum(np.isnan(fu)) > 0:
                inds = ~np.isnan(fu)
                concs_temp = copy.deepcopy(self.concentrations[inds])
                fu = fu[inds]
            else:
                concs_temp = self.concentrations
            # Determine start values for Ku and Kd
            start_Ku = np.max(fu) / (1.01 - np.max(fu))  # Based on maximum fu since ligand stabilizes folded, 1.01 to avoid division by zero
            half_pos = np.argmin(np.abs(fu - .5 * (np.max(fu) - np.min(fu))))
            start_Kd = concs_temp[half_pos]
            
            # Do curve fitting
            if model_type == "One_Site":
                fit_func = calculate_fitted_isothermal_simple(self.pconc)
                params, pcov = curve_fit(fit_func, concs_temp, fu, max_nfev=5E3, method='trf',
                    p0=[start_Ku, start_Kd], bounds=((0, 0), (np.inf, 0.1)))    

            # Protein dimerization, Equilibria are U <-> F & F + F <-> FF
            if model_type == "Dimerization":
                fit_func = calculate_fitted_isothermal_simple_dim()
                params, pcov = curve_fit(fit_func, concs_temp, fu, max_nfev=5E3, method='trf',
                    p0=[start_Ku, start_Kdim], bounds=((0, 0), (np.inf, np.inf)))    

            if model_type == "Two_Sites_One_Kd":
                fit_func = calculate_fitted_isothermal_2kds_simple2(self.pconc)
                params, pcov = curve_fit(fit_func, concs_temp, fu, max_nfev=5E3, method='trf',
                    p0=[start_Ku, start_Kd], bounds=((0, 0), (np.inf, 1)))  

            if model_type == "Two_Sites_Two_Kd":
                fit_func = calculate_fitted_isothermal_2different_kds_simple2(self.pconc)
                params, pcov = curve_fit(fit_func, concs_temp, fu, max_nfev=5E3, method='trf',
                    p0=[start_Ku, start_Kd,start_Kd], bounds=((0, 0, 0), (np.inf, 2, 2)))  #2 M upper bound

            self.isothermal_fits[:, i] = fit_func(self.unique_concs, *params)

            # Only working for one Kd models for now !!!
            fuPred          = fit_func(concs_temp, *params)
            desired_rss     = get_desired_rss(fu,fuPred,len(concs_temp),2)
            kd_estimated    = params[1]
            asymmetric_ci95 = get_asymmetric_ci95(kd_estimated,fu,concs_temp,start_Ku,desired_rss,self.pconc,model_type)

            bind_ci95_asymmetric_low.append(asymmetric_ci95[0])
            bind_ci95_asymmetric_up.append(asymmetric_ci95[1])

            bind_params.append(params)
            bind_errors.append(np.sqrt(np.diag(pcov)))
        # Write to instance
        self.bind_ci95_asymmetric_low = np.array(bind_ci95_asymmetric_low)
        self.bind_ci95_asymmetric_up = np.array(bind_ci95_asymmetric_up)

        self.bind_params          = np.array(bind_params)
        self.bind_errors          = np.array(bind_errors)
        return None

    ## Start of Functions for shiny server

    def shiny_export_fit_fluo(self):
        '''
        This is a helper function to export fluorescence fits
        that can then be plotted in the web app

        It writes the following arrays to the instance:
        self.fit_fluo_pred   
        self.fit_fluo_params 
        self.fit_fluo_errs   
        '''    

        if not hasattr(self, 'select_fitting_params'):
            return None, None

        fitting_params, fitting_errors, fit_melting = self.select_fitting_params()
        num_datasets = len(self.concentrations)

        fit_fluo = np.empty((self.fluo).shape)

        for c_index in range(num_datasets):

            Tm = fitting_params[0, c_index]  # Tm
            dH = fitting_params[1, c_index]  # dH
            ti = fitting_params[2, c_index]  # unfolded intercept
            bi = fitting_params[3, c_index]  # folded intercept
            ts = fitting_params[4, c_index]  # unfolded slope
            bs = fitting_params[5, c_index]  # folded slope

            for t_index in range(len(self.temps)):
                T = self.temps[t_index] + 273.15
                R = 1.987 / 1000
                dG = dH * (1 - T / (Tm + 273.15)) - self.cp * (Tm + 273.15 - T + T * np.log(T / (Tm + 273.15)))
                try:
                    Ku = np.exp(-dG / (R * T))
                except RuntimeWarning:
                    print(dG, -dG / (R*T))

                fit_fluo[t_index,c_index] = (Ku / (1 + Ku)) * (ts * T + ti) + (1 / (1 + Ku)) * (bs * T + bi)
        
        self.fit_fluo_pred   =  fit_fluo     
        self.fit_fluo_params =  fitting_params
        self.fit_fluo_errs   =  fitting_errors

    def pre_fit_isothermal(self):
        """
        Fill self.isothermal_data with the predicted fraction unfolded
        """

        # Initialize isothermal data and fits
        self.isothermal_data = np.empty([len(self.concentrations), len(self.isothermal_ts)])
        fitting_params, fitting_errors, fit_melting = self.select_fitting_params()
        # Loop through temperatures and fill self.isothermal_data
        for i, binding_temp in enumerate(self.isothermal_ts):
            # Select fu
            self.isothermal_data[:, i] = calculate_fraction_unfolded(binding_temp,
                self.concentrations,fitting_params,self.cp).transpose()

    def shiny_export_isothermal(self):
        '''
        This is a helper function to export isothermal fits
        that can then be plotted in the web app

        It write the following arrays to the instance:
        self.kd_model_conc: fine conc grid for plotting
        self.kd_models: Fitted binding curves
        self.kd_models_lower: Lower limits based on fit errors
        self.kd_models_upper: Upper limits based on fit errors
        '''
        if not hasattr(self, 'bind_params'):
            return None, None
        # Determine fit_melting
        fitting_params, fitting_errors, fit_melting = self.select_fitting_params()
        if len(fitting_params) == 0:
            return None, None
        # Initialize arrays
        kd_fit_conc = np.unique(self.concentrations)
        # Finer grid for lines
        kd_model_conc = 10 ** (np.linspace(np.log10(kd_fit_conc[1] * .5), np.log10(kd_fit_conc[-1] * 2), 100))
        kd_models = np.ones((len(kd_model_conc), len(self.isothermal_ts))) * np.nan
        kd_models_lower = copy.deepcopy(kd_models)
        kd_models_upper = copy.deepcopy(kd_models)
        # Loop through isothermal ts and fill arrays
        for t_index in range(len(self.isothermal_ts)):
            # Choose right parameters and errors
            params = self.bind_params[t_index]
            err    = self.bind_errors[t_index]
            fu     = self.isothermal_data[:, t_index]
            ku, kd = params[0], params[1]
            ku_err, kd_err = err[0], err[1]
            # Model
            model = calculate_fitted_isothermal_simple(self.pconc)
            kd_models[:, t_index] = model(kd_model_conc, *params)
            # Errors
            params_lower = params - err
            params_lower[params_lower< 0] = 1E-20

            #   This are not real confidence bands from non linear regression !!!!!!!
            #   We are not using them in the shiny server

            kd_models_upper[:, t_index] = model(kd_model_conc, *(params + err))
            kd_models_lower[:, t_index] = model(kd_model_conc, *(params_lower))
        # Write to class instance
        self.kd_model_conc = kd_model_conc
        self.kd_models = kd_models
        self.kd_models_lower, self.kd_models_upper = kd_models_lower, kd_models_upper


    ## End of Functions for shiny server

    def select_fitting_params(self):
        '''
        Checks if temperature fit was done and selects respective fitting parameters
        If global fit exists, it will choose these parameters
        Otherwise it selects local fitting parameters
        Args:
        Returns:
        '''
        # Check if global fit with variable Cp was done

        if self.fitting_option == "Global_CP":

            fit_params = self.global_fit_cp_params
            fit_errors = self.global_fit_cp_errors
            fit_melting = 'global_cp'

        elif self.fitting_option == "Global":

            fit_params = self.global_fit_params
            fit_errors = self.global_fit_errors
            fit_melting = 'global'

        elif self.fitting_option == "Local":

            fit_params = self.local_fit_params
            fit_errors = self.local_fit_errors
            fit_melting = 'local'

        else:

            return [], [], []
        # Update melting temperatures
        self.tms = fit_params[0]
        
        return fit_params, fit_errors, fit_melting

debug_run = False

xls_file1 = "./www/nDSFdemoFile.xlsx"

if debug_run:

    concentrations = [5,5,2.5,2.5,1.3,1.3,0.63,0.63,0.31,0.31,0.16,0.16,0.078,0.078,0.039,0.039,0.02,0.02,
    0.0098,0.0098,0.0049,0.0049,0.0024,0.0024,0.0012,0.0012,0.00061,0.00061,0.00031,0.00031,0.00015,0.00015]

    concentrations = [x / 1000 for x in concentrations]

    t = DSF_binding()

    t.load_nanoDSF_xlsx(xls_file1)
    t.set_signal('Ratio')

    t.concentrations = np.array(concentrations)

    #t.median_filter(2)

    t.fluo  = filter_fluo_by_temp(t.fluo,t.temps,40,90)
    t.temps = filter_temp_by_temp(t.temps,40,90)

    t.fluo  = t.fluo[:,:32]

    t.cp = 0

    t.fit_fluo_local()
    t.fit_fluo_global()
    #t.fit_fluo_global_cp()

    #print(t.cp)

    t.pconc = 6.5/(1E6)
    t.isothermal_ts = np.array([65,66])
    t.fit_isothermal("One_Site")
    t.shiny_export_fit_fluo()
    t.shiny_export_isothermal()

    t.export_JSON_file('person.json')

if False:

    t2 = DSF_binding()
    t2.load_JSON_file("/home/osvaldo/Downloads/fAffinitySession_2023-05-21.json")

    t2.export_JSON_file("/home/osvaldo/Downloads/test.json")

    #kds = [x[1]*1E6 for x in t.bind_params ]

    #errors = [x[1]*1E6 for x in t.bind_errors ]

    #errors = (np.array(errors) / np.array(kds)) * 100

    #print(t.bind_ci95_asymmetric)
    #print(errors)


#t = DSF_binding()
#folder="/home/osvaldo/pkng_inhibitors/SupplementaryMaterial/melting_curves_analysis"
#file=folder+"/"+"2.3.20-lig04-m9-zn-dtrp-batch2.xlsx"
#t.load_nanoDSF_xlsx(file)
#t.set_signal('350nm')

#t.fluo  = filter_fluo_by_temp(t.fluo,t.temps,30,50)
#t.temps = filter_temp_by_temp(t.temps,30,50)

#np.save('testFluo.npy', t.fluo)    # .npy extension is added if not given
#np.save('testTemp.npy', t.temps) 

#file = "/home/osvaldo/Downloads/2022-01-19_1616CET_T6-039.xlsx"
#snames = (get_sheet_names_of_xlsx(file))
#print(snames)
#t.load_tycho_xlsx(file)

#print(t.fluo)
#t.set_signal('Ratio')
#t.median_filter(5)
#t.estimate_fluo_derivates(10)



#print(t.min_derivative)
