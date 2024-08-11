#!/usr/bin/python3
import os
import pandas as pd
import numpy  as np
from pandas.api.types import is_string_dtype
from pandas.api.types import is_numeric_dtype

# Function to find the first occurrence of a value in an array
def find_first_occurrence(array, value):

    try:

        return int(np.where(array == value)[0][0])

    except:

        return None

def fileIsNanoTemperCSV(file):

    with open(file) as f:

        l0 = f.read().splitlines()[0].lower()

    return 'capillary index' in l0 and 'response' in l0 

#  Load CSV with five columns, Capillary Index;Response;Ligand Concentration;Target Concentration;Buffer

def readNanoTemperCSV(file,sep,expid='A'):

    csv                    = pd.read_csv(file,sep=sep)
    ligand_concentrations  = np.array(csv.iloc[:,2]).astype('float32')
    proteinConc            = np.array(csv.iloc[:,3]).astype('float32')

    signal                 =  np.tile(np.array(csv.iloc[:,1]).astype('float32'),(3, 1))

    experimentID  = np.array([expid for _ in ligand_concentrations])

    return ligand_concentrations, signal, proteinConc, experimentID

def readManyNanoTemperCSV(files,sep):

    ligand_concentrations, signal, proteinConc, experimentID = [], [], [], []

    for file in files:

        file_name_without_extension = os.path.splitext(os.path.basename(file))[0]

        lt, st, pt, et = readNanoTemperCSV(file,sep,file_name_without_extension)

        ligand_concentrations.append(lt)
        signal.append(st)
        proteinConc.append(pt)
        experimentID.append(et)

    ligand_concentrations = np.concatenate(ligand_concentrations)  
    proteinConc           = np.concatenate(proteinConc)  
    experimentID          = np.concatenate(experimentID)
    signal                = np.hstack(signal)

    return ligand_concentrations, signal, proteinConc, experimentID

def getSepCharacter(file):

    """
    
    Get the character that splits the csv into 2-4 columns

    Input 
        - file: path of the input file to be read

    Output
        - string, one of this: [" ",",",";"])
          
    """

    with open(file) as f:

        ls = f.read().splitlines()
        
        for sep in [",",";"," "]:

            lens = [ len(l.split(sep)) for l in ls ]
            
            c1 = (lens.count(lens[0]) == len(lens)) # Check all columns have equal length
            c2 = lens[0] >=2 and lens[0] <= 5

            if c1 and c2:
                return sep

    return None

def csvHasHeader(file,sep):

    """
    Split the first line using a custom separator and 
    check that it contains, at most, one string that represents a number (this could be the ligand name)

    Input 
        - file: path of the input file to be read
        - sep:  string, one of this: [" ",",",";"], used to split the first line into a list

    Output
        - boolean

    """

    with open(file) as f:

        ls = f.read().splitlines()
        firstLine = ls[0].split(sep)

        # Dummy replace to make the method isnumeric work
        x = [ l.replace("-", "9").replace(".", "9").replace("e", "9").replace("+", "9").replace("E", "9") for l in firstLine ]
        x = [ l.isnumeric() for l in x ]

        return sum(x) <= 1 

def readCSV(file,sep,header):

    """
    Try to read a csv with the ligand concentration, signal, and protein concentration and/or experiment ID

    Input
        - file: path of the input file to be read
        - sep: dictionary where the keys are possible strings (one of this: [" ",",",";"])
          and the values are the number of detected columns (minimum 2) 
        - header: boolean that describes if the file has a header

    Output
        - numpy arrays of the ligand concentration, signal, 
          protein concentration and experiment ID
    """

    # I tried if pandas could infer the header - it didn't work!
    if header:
        csv = pd.read_csv(file,sep=sep)
    else:
        csv = pd.read_csv(file,sep=sep,header=None)

    # Read Nanotemper csvs

    # By default, the first column has the ligand concentration
    # and the second column, the signal value

    ligand_concentrations  =  np.array(csv.iloc[:,0]).astype('float32')
    signal                 =  np.tile(np.array(csv.iloc[:,1]).astype('float32'),(3, 1))

    # default values
    proteinConc   = np.ones(len(ligand_concentrations))
    experimentID  = np.array(["A" for _ in ligand_concentrations])

    # Let's check the 3rd and 4th columns 
    nColumns = csv.shape[1]

    if nColumns >= 3:

        # Assign the protein concentration to the 3rd column if this column has numeric values
        if is_numeric_dtype(csv.iloc[:, 2]):

            proteinConc = np.array(csv.iloc[:,2]).astype('float')

            # Assign the experiment ID (4th column, if present)
            if nColumns == 4:

                experimentID   = np.array(csv.iloc[:,3]).astype('str')

        else:

            # Assign the experiment ID (3rd column if it isn't numeric)
            experimentID   = np.array(csv.iloc[:,2]).astype('str')
    
            if nColumns == 4:

                # Catch exception in case the fourth column is not numeric! (error in the input file)
                try:
                    proteinConc   = np.array(csv.iloc[:,3]).astype('float')
                except:
                    pass

    return ligand_concentrations, signal, proteinConc, experimentID

def median_filter_from_fluo_and_time_vectors(fluo_vec,time_vec,rolling_window):

    """

    Compute the median filter of the fluorescence vector using a time rolling window


	Requires: 

		1) The fluorescence 1D vector 'fluo_vec'
		2) The temperature  1D vector 'time_vec'
		3) The size of the rolling window in celsius 'rolling_window'

	Returns the fluorescence vector passed through the median filter


    """

    scaling_factor = 10000 

    time_vec     =  np.multiply(time_vec,scaling_factor).astype(int) 
    series       =  pd.Series(fluo_vec,index=time_vec,dtype=float)
    series.index =  pd.to_datetime(series.index,unit='s')

    roll_window  = str(int(rolling_window*scaling_factor))+"s"

    fluo_df_median_filt = series.rolling(roll_window).median().to_numpy()

    return fluo_df_median_filt

def normalize_with_cold_region_mean(signal,times):
    '''
    Normalize with values at t <= 0
    '''

    normalized_signal = signal / np.mean(signal[times <=0,:], axis=0)

    return normalized_signal