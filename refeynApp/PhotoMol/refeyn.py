from scipy.optimize import curve_fit
import pandas as pd
import numpy  as np
import h5py

from helpers import * 

class Refeyn:
    '''
    Simple class to load refeyn eventsFitted.h5 files from mass photometer

    If no information about the masses is given, then calibration parameters can be used 
    to convert contrasts to masses.
    '''

    def __init__(self):

        self.massesLoaded = False

        return None
    
    def load_data_h5(self,filename):
        self.fn           = filename
        self.massesLoaded = False

        data = h5py.File(self.fn, 'r')
        contrasts       = np.array(data['contrasts']).squeeze()
        self.contrasts  = contrasts[~np.isnan(contrasts)]

        if 'masses_kDa' in data.keys():

            masses_kDa      = np.array(data['masses_kDa']).squeeze()
            self.masses_kDa = masses_kDa[~np.isnan(masses_kDa)]

            # Check we have some data
            if len(self.masses_kDa) > 50:

                self.n_binding    = np.sum(self.masses_kDa >= 0)
                self.n_unbinding  = np.sum(self.masses_kDa <  0)

                ## Find peaks for big masses

                self.create_histo([0,max(self.masses_kDa)],bin_width=10)
                self.findInitialPeaks()

                self.massesLoaded = True

        return None

    def load_data_csv(self,filename):
        self.fn           = filename
        self.massesLoaded = False

        data = pd.read_csv(filename)
        
        contrasts       = np.array(data['contrasts']).squeeze()
        self.contrasts  = contrasts[~np.isnan(contrasts)]

        if 'masses_kDa' in data.columns:

            self.masses_kDa   = np.array(data['masses_kDa']).squeeze()

            try:

                self.n_binding    = np.sum(self.masses_kDa >= 0)
                self.n_unbinding  = np.sum(self.masses_kDa <  0)
                self.create_histo([0,max(self.masses_kDa)],bin_width=10)
                self.findInitialPeaks()
                self.massesLoaded = True

            except:
                
                pass

        return None

    def findInitialPeaks(self):
        
        pks_initial       = findPeaks(self.hist_counts,self.hist_mass)
        self.pks_initial  = [int(p) for p in pks_initial if p >= 0]             

        return None

    def contrastsToMass(self,slope,intercept):
        '''
        Function to convert masses from
        contrasts using known calibration parameters 

        Caution! slope and intercept are based on f(mass) = contrast !!!! 
        In other words, contrast = slope*mass + intercept

        '''
        
        interceptInverse = -intercept / slope 
        slopeInverse     = 1 / slope          

        self.masses_kDa   = np.polyval(np.array([slopeInverse,interceptInverse]), self.contrasts)
        self.n_binding    = np.sum(self.masses_kDa >= 0)
        self.n_unbinding  = np.sum(self.masses_kDa <  0)

        self.create_histo([0,max(self.masses_kDa)],bin_width=10)
        self.findInitialPeaks()
        self.massesLoaded = True
        
        return None

    def export_h5_dataset(self,filename):

        '''
        Creates a h5 file with a dataset called 'masses_kDa' 
        '''

        assert self.massesLoaded

        hf = h5py.File(filename, 'w')
        hf.create_dataset('masses_kDa', data=self.masses_kDa)
        hf.create_dataset('contrasts', data=self.contrasts)
        hf.close()

        return None

    def create_histo(self, window=[0,2000], bin_width=10):
        '''
        Creates histogram of masses
        '''
        # Determine number of bins based on bin_width
        nbins = (window[1] - window[0]) // bin_width
        nbins = int(nbins)
        # Create histogram
        self.hist_counts, self.hist_bins = np.histogram(self.masses_kDa, range=window, bins=nbins)
        self.hist_mass = (self.hist_bins[1:] + self.hist_bins[:-1]) / 2.0
        # Write parameters to instance
        self.bin_width    = bin_width
        self.hist_window  = window
        self.hist_nbins   = nbins

        return None

    def create_fit_table(self):
        '''
        Uses info in self.fit to generate a 
        pandas DataFrame that summarizes fit results
        '''
        # Create lists with fitting parameters
        # These are later used to create a pandas DataFrame
        list_pos, list_sigma, list_ampl, list_counts = [], [], [], []
        # Loop over entries in optimized parameters
        for i in range(int(len(self.popt)/3)):
            list_pos.append(round(self.popt[3*i]))
            list_ampl.append(round(self.popt[3*i+1]))
            list_sigma.append(round(self.popt[3*i+2]))
            list_counts.append(round(np.trapz(self.fit[:,i+1], x=self.fit[:,0]) / np.diff(self.hist_mass)[0]))
        # Create Pandas Dataframe

        total_counts = self.n_binding

        ## Add the counts from the unbinding data, if necessary
        has_negative_elements = any(element < 0 for element in list_pos)

        if has_negative_elements:
            total_counts = total_counts + self.n_unbinding

        self.fit_table = pd.DataFrame(data={'Position / kDa': list_pos,
                                            'Sigma / kDa': list_sigma,
                                            'Counts' : list_counts,
                                            'Counts / %': np.round(np.array(list_counts)/total_counts*100),
                                            'Amplitudes' : list_ampl}
                                      ) 
        return None
    
    def fit_histo(self, guess_pos=[66,148,480], tol=100, max_std=200,min_observed_mass=40,baseline=0):
        '''
        Fit gaussians to histogram
        guess: list with guessed centers, defines the number of gaussians to be used, 
        '''
        # If no guess are taken, only fit one gaussian and use maximum in histogram as guess

        if len(guess_pos) == 0:
            return None
        else:
            # Get amplitude for each guess position
            guess_amp = []
            for pos in guess_pos:
                ind = np.argmin(np.abs(self.hist_mass - pos))
                guess_amp.append(self.hist_counts[ind])
            fit_guess    = np.column_stack((np.array(guess_pos), np.array(guess_amp), np.array([5]*len(guess_pos)))).flatten()
            lower_bounds = np.column_stack((np.array(guess_pos) - tol , np.array([0]*len(guess_pos)), np.array([0]*len(guess_pos)))).flatten()
            upper_bounds = np.column_stack((np.array(guess_pos) + tol , np.array([np.max(self.hist_counts)*1.2]*len(guess_pos)), np.array([max_std]*len(guess_pos)))).flatten()
            bounds = (tuple(lower_bounds), tuple(upper_bounds))

        # Use truncated gaussian function
        def fitting_helper_func(x,*params):
            return truncated_multi_gauss_with_baseline(x,min_observed_mass,baseline,*params)

        func = fitting_helper_func
        
        values_have_sense = np.greater(np.abs(self.hist_mass),min_observed_mass)

        # Restrict data range for the fitting
        hist_counts = self.hist_counts * values_have_sense
        hist_mass   = self.hist_mass   

        # Do fit
        self.popt, self.pcov = curve_fit(func, hist_mass, hist_counts, p0=fit_guess, bounds=bounds)  #, method='dogbox', maxfev=1E5)
        # Create fit and individual gaussians for plotting
        # Finer grid
        x = np.linspace(np.min(self.hist_mass), np.max(self.hist_mass), 800)
        single_gauss = []
        for i in range(0, len(self.popt), 3):
            ctr = self.popt[i]
            amp = self.popt[i+1]
            wid = self.popt[i+2]
            single_gauss.append(func(x, ctr, amp, wid))
        
        # Sum of all
        fit_sum = func(x, *self.popt)

        # Create one array for all
        self.fit = np.column_stack((x, np.array(single_gauss).T, fit_sum))
        # Errors
        self.fit_error = np.sqrt(np.diag(self.pcov))
        # Create fit table
        self.create_fit_table()

        return None

if False: 
    
    refeyn = Refeyn()

    #refeyn.load_data_h5("demoTestFile")

    refeyn.load_data_h5("/home/osvaldo/Downloads/test.h5")

    print(refeyn.pks_initial)

    #massActual = refeyn.masses_kDa

    #refeyn.contrastsToMass(-5.72173609e-05,4.73634440e-05)

    #massNew    = refeyn.masses_kDa

    #refeyn.create_histo(window=np.array([0,1000]),bin_width=4)

    #refeyn.load_data("/home/osvaldo/Downloads/eventsFitted.h5")
    #print(refeyn.hist_mass)
    #print(refeyn.pks_initial)

    #pks = findPeaks(refeyn.hist_counts,refeyn.hist_mass)
    #print(pks)
    #refeyn.fit_histo(guess_pos=np.array([66,148,480]),tol=20,max_std= 100) 

    #print(refeyn.fit_table)

    #data = h5py.File("/home/osvaldo/Downloads/eventsFitted.h5", 'r')
    #print(data.keys())

    #df = np.array(data['mass_model']).squeeze()
    #print(df)
