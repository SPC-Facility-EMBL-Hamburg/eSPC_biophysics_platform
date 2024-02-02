#!/usr/bin/python3
import pandas as pd
import numpy as np
import numpy.polynomial.polynomial as poly

from helpers         import *


class MST_fit:

    """
    Class for loading the data from MST protein-ligand binding studies
    This class was written by Osvaldo Burastero based on code done by Stephan Niebling

    No warranty whatsoever
    If you have questions please contact me:
    	oburastero@gmail.com
    """

    def __init__(self):

        self.expID_vector = np.array("A") # Default value, changed later by the user

        return None

    """
            self.signal_data_dictionary = None
            self.time_data_dictionary   = None
    """

    def init_dictionary_to_store_original_fluo_time_data(self):

        self.signal_data_dictionary = {}
        self.time_data_dictionary   = {}

    def load_example_data(self):

        """ Requires signal.npy and time.npy !!!

        """

        self.init_dictionary_to_store_original_fluo_time_data()
        self.signal_data_dictionary["Raw Fluorescence"] = np.load("signal.npy")
        self.time_data_dictionary["Raw Fluorescence"]   = np.load("time.npy")
        self.concs  = 150 / np.power(2,np.arange(0,16))

        return None

    def load_MST_xlsx(self, mst_file_xlsx):

        """
        Load nanotemper xlsx file and extract data
        """

        self.init_dictionary_to_store_original_fluo_time_data()

        xls = pd.ExcelFile(mst_file_xlsx)
        
        dat = pd.read_excel(xls, 'RawData', index_col=None, header=None)

        first_row_signal = int(np.argwhere(list(dat.iloc[:,0] == 'Time [s]'))) + 1
        row_ligand_conc  = int(np.argwhere(list(dat.iloc[:,0] == 'Ligand Concentration:'))) 
        row_cap_pos      = int(np.argwhere(list(dat.iloc[:,0] == 'Capillary Position:'))) 
        row_ligand       = int(np.argwhere(list(dat.iloc[:,0] == 'Ligand:'))) 

        total_measurements = int((len(dat.iloc[first_row_signal,:]) + 1) / 3)
        ligand_columns     = [1+x*3 for x in range(total_measurements)]

        ligand_concentrations = np.array(dat.iloc[row_ligand_conc,ligand_columns]).astype('float32')
        ligand_capillaries    = np.array(dat.iloc[row_cap_pos,ligand_columns])
        ligand_name           = np.array(dat.iloc[row_ligand,1])

        signal = np.array(dat.iloc[first_row_signal:,ligand_columns]).astype('float32')
        time   = np.array(dat.iloc[first_row_signal:,0]).astype('float32')

        #Remove rows with empty values
        no_nan_rows = ~np.isnan(signal).any(axis=1)

        signal = signal[no_nan_rows]
        time   = time[no_nan_rows]

        # Default value in case we can't read it from the file
        protConc = np.ones(len(ligand_concentrations)) 

        try: 
            row_prot_conc    =  int(np.argwhere(list(dat.iloc[:,0] == 'TargetConcentration:')))
            possibleProtConc =  np.array(dat.iloc[row_prot_conc,ligand_columns]).astype('float32')
            if len(possibleProtConc) == len(ligand_concentrations) and len(possibleProtConc) > 5:
                protConc = possibleProtConc

        except:
            pass
        
        self.protConc      = protConc
        self.experimentID  = np.array(["A" for _ in ligand_concentrations])
        self.concs         = ligand_concentrations
        
        self.signal_data_dictionary["Raw Fluorescence"] = signal
        self.time_data_dictionary["Raw Fluorescence"]   = time

        return None

    def load_MST_csv(self, file, sep,header):

        """
        Load CSV with two columns, ligand concentration and signal
        """

        self.init_dictionary_to_store_original_fluo_time_data()

        ligand_concentrations, signal, proteinConc, experimentID = readCSV(file, sep,header)

        self.concs                                      = ligand_concentrations
        self.protConc                                   = proteinConc
        self.experimentID                               = experimentID

        self.signal_data_dictionary["Raw Fluorescence"] = signal
        self.time_data_dictionary["Raw Fluorescence"]   = np.arange(-1,2) # placeholder only

        return None

    def set_signal(self,which):

        """
        Assign self.signal, self.times according to the desired signal
        """
        self.signal    = self.signal_data_dictionary[which]
        self.times     = self.time_data_dictionary[which]

    def sort(self):
        '''
        Sort ligand concentrations and fluorescence
        with ascending ligand concentration
        '''
        sort_ind = np.argsort(self.concs)
        self.concs = self.concs[sort_ind]
        self.signal = self.signal[:, sort_ind]
        return None

    def median_filter(self,n_degree_window):

        """

        Use this function if the fluorescence curves present spikes. It applies a
        rolling median window filter 

        Requires:

        An integer 'n_degree_window'


        """

        self.signal = np.apply_along_axis(median_filter_from_fluo_and_time_vectors,0,
            self.signal,self.times,n_degree_window)

        return None

    def get_cold_fluo(self,cold_min,cold_max):


        ind_cold    = (self.times >= cold_min) * (self.times <= cold_max)

        if not any(ind_cold):

            self.F_cold = self.signal[np.abs(self.times - cold_max).argmin()]

        else:

            self.F_cold = np.mean(self.signal[ind_cold,:], axis=0)

        return None

    def get_fnorm(self,hot_min,hot_max):

        ind_hot     = (self.times >= hot_min) * (self.times <= hot_max)

        if not any(ind_hot):

            F_hot = self.signal[np.abs(self.times - hot_max).argmin()]

        else:

            F_hot       = np.mean(self.signal[ind_hot,:], axis=0)

        self.F_norm = F_hot / self.F_cold

        return None

run_debug = False

if run_debug:

    # Hot region
    hot = 25
    # File names, adjust these!
    mst = MST_fit()
    mst.load_MST_xlsx("./www/demo.xlsx")
    mst.set_signal("Raw Fluorescence")
    print(mst.times)
    #mst.sort()
    #mst.normalize()
    # Calculate fnorm and chose the "hot" region
    #mst.calc_fnorm(hot=hot)

#mst = MST_fit()
#mst.load_MST_csv("./www/test3.csv")

#for file in ["./www/test3.csv","./www/test.csv","./www/test2.csv","./www/test4.csv"]:

#	s = getSepCharacter(file)

#	for k in s:
#		h  = csvHasHeader(file,k)
#		res = readCSV(file,k,h)
#		print(res)

if False:

    print(csvHasHeader("./www/test4.csv",","))
    print(csvHasHeader("./www/test3.csv",","))
    print(csvHasHeader("./www/test2.csv",","))
    print(csvHasHeader("./www/test.csv",","))

    readCSV("./www/test4.csv",",",True)
    readCSV("./www/test3.csv",",",False)
    readCSV("./www/test2.csv",",",False)
    readCSV("./www/test.csv",",",False)
