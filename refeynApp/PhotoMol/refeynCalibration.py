from scipy.optimize import curve_fit
import pandas as pd
import numpy  as np
import h5py
from helpers import * 

class RefeynCalib:
    '''
    class to load refeyn files from mass photometer where only the contrasts are present 
    Used with known masses and to fit a line to the contrast versus mass.
    '''

    def __init__(self):
        
        return None

    def load_data_h5(self,filename):
        self.fn = filename
        # Load data
        data = h5py.File(self.fn, 'r')
        
        self.contrasts  = np.array(data['contrasts']).squeeze()
        self.n_binding  = np.sum(self.contrasts<0)
        self.create_histo()
        self.findInitialPeaks()

        return None

    def load_data_csv(self,filename):
        self.fn = filename
        
        data = pd.read_csv(filename)
        
        contrasts  = np.array(data['contrasts']).squeeze()
        self.contrasts = contrasts[~np.isnan(contrasts)]
        self.n_binding  = np.sum(self.contrasts<0)
        self.create_histo()
        self.findInitialPeaks()

        return None

    def create_histo(self, window=[-1,0], bin_width=0.0004):
        '''
        Creates histogram 
        '''
        # Determine number of bins based on bin_width
        nbins = (window[1] - window[0]) // bin_width
        nbins = int(nbins)
        # Create histogram
        self.hist_counts_contrasts, self.hist_bins = np.histogram(self.contrasts, range=window, bins=nbins)
        self.hist_centers_contrasts = (self.hist_bins[1:] + self.hist_bins[:-1]) / 2.0
        # Write parameters to instance
        self.bin_width    = bin_width
        self.hist_window  = window
        self.hist_nbins   = nbins

        return None

    def findInitialPeaks(self):
        
        pks_initial       = findPeaks(self.hist_counts_contrasts,self.hist_centers_contrasts,
            height=10, distance=4,prominence=4,masses=False)
        self.pks_initial  = [p for p in pks_initial if p <= 0]             

        return None

    def create_fit_table(self):
        '''
        Uses info in self.fit to generate a 
        pandas DataFrame that summarizes fit results
        '''

        list_pos, list_sigma, list_ampl, list_counts = [], [], [], []
        # Loop over entries in optimized parameters
        for i in range(int(len(self.popt)/3)):
            list_pos.append((self.popt[3*i]))
            list_ampl.append((self.popt[3*i+1]))
            list_sigma.append((self.popt[3*i+2]))
            list_counts.append(round(np.trapz(self.fit[:,i+1], x=self.fit[:,0]) / np.diff(self.hist_centers_contrasts)[0]))
        # Create Pandas Dataframe
        self.fit_table = pd.DataFrame(data={'Position / contrast': list_pos,
                                            'Sigma / contrast': list_sigma,
                                            'Counts' : list_counts,
                                            'Counts / %': np.round(np.array(list_counts)/self.n_binding*100),
                                            'Amplitudes' : list_ampl}
                                      ) 

        return None

    def fit_histo(self, guess_pos=[0], max_std=0.1,tol=0.05):
        '''
        Fit gaussians to histogram
        guess: list with guessed centers, defines the number of gaussians to be used, 
        '''

        # Get amplitude for each guess position
        guess_amp = []
        for pos in guess_pos:
            ind = np.argmin(np.abs(self.hist_centers_contrasts - pos))
            guess_amp.append(self.hist_counts_contrasts[ind])

        fit_guess    = np.column_stack((np.array(guess_pos), np.array(guess_amp), np.array([0.001]*len(guess_pos)))).flatten()
        lower_bounds = np.column_stack((np.array(guess_pos) - tol , np.array([0]*len(guess_pos)), np.array([0]*len(guess_pos)))).flatten()
        upper_bounds = np.column_stack((np.array(guess_pos) + tol , np.array([np.max(self.hist_counts_contrasts)*1.2]*len(guess_pos)), np.array([max_std]*len(guess_pos)))).flatten()
        bounds = (tuple(lower_bounds), tuple(upper_bounds))

        # Use truncated gaussian function
        def fitting_helper_func(x,*params):
            return multi_gauss(x,*params)

        func = fitting_helper_func

        self.popt, self.pcov = curve_fit(func, self.hist_centers_contrasts, self.hist_counts_contrasts, p0=fit_guess, bounds=bounds)  #, method='dogbox', maxfev=1E5)
        
        # Create fit and individual gaussians for plotting
        # Finer grid
        x = np.linspace(np.min(self.hist_centers_contrasts), np.max(self.hist_centers_contrasts), 1000)
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

    def calibrate(self,calib_floats):
        ''' 
        Calibration based on contrasts histogram
        You need to have run create_histo and fit_histo(contrasts=True) before
        calib_floats: List of calibration standards defined by floats: e.g. 66, 146, 480
        '''
        # Dictionary with masses

        calib_stand = np.array(calib_floats)
        # Sort in reversed order
        #calib_stand.sort()

        # Calibration points from fitting
        calib_points = self.fit_table['Position / contrast']
        # Check if number of fitted gaussians fits to number of calibrants
        #if not (len(calib_points) == len(calib_stand)):
        #    print("%i calibration standards given, but %i gaussians fitted to contrasts!" % (len(calib_stand), len(calib_points)))
        # Now use these to calibrate system
        # First order polynomial fit

        params = np.polyfit(calib_stand, calib_points, 1) # x-axis are the masses (calib_stand), y-axis are the contrasts
        # returns slope, intercept
        # Calculate R2
        calc = np.polyval(params, calib_stand)
        r2   = r_squared(calib_points, calc)

        self.calib = {'standards':  calib_stand,
                      'exp_points': calib_points,
                      'fit_params': params,
                      'fit_r2': r2}

        self.calib_stand    = calib_stand
        self.calib_points   = calib_points
        self.calib_params   = params
        self.calib_r2       = r2

        return None

#refeyn = RefeynCalib()
#refeyn.load_data_csv("/home/osvaldo/Downloads/eventsFound.csv")
#print([round(p,3) for p in refeyn.pks_initial])
#print(min(refeyn.contrasts)*1e3)

#refeyn.fit_histo(guess_pos=[-.004, -.008, -.027])
#refeyn.calibrate([66, 146, 480])
#print(refeyn.calib_params)
#print(refeyn.fit_table)
