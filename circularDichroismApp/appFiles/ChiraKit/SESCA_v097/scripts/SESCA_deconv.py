#!/usr/bin/env python
#

import sys
import os
import math
import time
import random

workdir = os.getcwd()
stime = time.time()

Usage =  "*******************************************\n"
Usage += "SESCA spectrum deconvolution module to determine ideal secondary structure (SS) compositions:\n"
Usage += "Usage: SESCA_deconv.py <spectrum_file> <basis_set_file> @flag <argument>\n"
Usage += "Possible command line flags are:\n"
Usage += "   @spect <filename> specify spectrum file (default: None)\n"
Usage += "   @lib   <filename> specify basis set file (default: None)\n"
Usage += "   @corr  <filename> specify baseline/sidechain correction file (default: None)\n"
Usage += "   @ref   <filename> specify refernce SS file (default: None)\n"
Usage += "   @write <0 - off or filename> specify output file name (default: 'BS_deconv.out')\n"
Usage += "   @W0 <float,float,...,float or 'random' or '0' for test mode >, set initial coefficients (default: 'random')\n"
Usage += "   @SF0 <float>  set initial scaling factor (default: 1.0)\n"
Usage += "   @test <0,1,2> run test mode: 0 - off, 1 - test evaluation, 2 - use initial coefficient (no fitting), (default: 0) '\n"
Usage += "   @range <float,float> set minimum and maximum wavelength for the spectrum (default: None)\n"
Usage += "   @rep   <int> set the number of repeats for deconvolutions with random starting weights (default: 10)\n"
Usage += "   @err   <int> select calibration curve for error estimation (default: 2 for standard rescaled error estimates)\n"
Usage += "   @iter <int> set the maximum number of iterations during the SS-search (default: 5000)\n"
Usage += "   @mult1 <int> set Lagrange multiplier for forcing the sum of coefficients to be 1.0 (default: 0)\n"
Usage += "   @mult2 <int> set Lagrange multiplier for forcing coefficients to be positive (default: 1000)\n"
Usage += "   @nrmsd <0,1> set normalization for fit RMSD by: 0 - None, 1 - spectrum amplitude (default: 0) \n"
Usage +=  "*******************************************\n"

#Specify SESCA files and directories here:
#################################################
SESCA_Dir =     "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097"
SESCA_scripts =  os.path.join(SESCA_Dir,"scripts")
##################################################
#the code from here should not be changed

#Default variables:
#Basic I/O paramerters:
infile = ""
libfile = ""
reffile = ""
corrfile = ""
outfile = "BS_deconv.out"
Inp_files = [infile, libfile, reffile, corrfile]
Out_files = [outfile]
#deconvolution parameters:
W_custom = []
SF0 = 1.0
L_range = ["",""]
random_seed = 10
test = 0
error = 0
Dec_param = [W_custom, L_range, random_seed, test, error, SF0] 
#minimization parameters:
Lambda1 = 0
Lambda2 = 1000
i_steps = 5000
nrmsd = 0
Opt_param = [Lambda1, Lambda2, i_steps, nrmsd]
#misc. control flags:
failmark = 0
verbosity = 1

#Def_Args = [infile, libfile, outfile, W_custom, L_range, random_seed, test, Lambda, i_steps,failmark, verbosity]
Def_Args = [Inp_files,  Out_files, Dec_param, Opt_param, failmark, verbosity]

#function definitions:
def Pass_Defaults():
	return Def_Args

def Pass_Imports():
	return IMPORTS

#function to control verbosity:
def Vprint(level,*Messages):
	if level <= verbosity:
		string = ""
		for message in Messages:
			string += str(message)+" "
		print(string)
	else:
		pass

#function to set verbosity levels:
def Set_verb(int):
	global verbosity
	verbosity = int

#read arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(3, "Reading in %1d arguments:\n" % argnum, Args)
#	getting default arguments:
	Inp_files,  Out_files, Dec_param, Opt_param, failmark, verbosity = Pass_Defaults()
	New_Args = [Inp_files,  Out_files, Dec_param, Opt_param, failmark, verbosity]
	infile, libfile, reffile,corrfile = Inp_files
	outfile = Out_files[0]
	W_custom, L_range, random_seed, test, error, SF0 = Dec_param
	Lambda1, Lambda2, i_steps, nrmsd = Opt_param

#	define possible command flags:
	FLAGS = ["spect","lib","corr","write","W0","SF0","range","rep","ref","err","iter","mult1","nrmsd","mult2","test","verb"]
	flag = ""
	acnt = 0
#	modify defaults based on read arguments:
	for arg in Args:
		if arg.startswith("@"):
			flag = arg.strip("@")
			if not flag in FLAGS:
				Vprint(1, "Unknown flag:",flag)
				New_Args[4] = 1
 
		elif flag == "spect" or (flag == "" and acnt == 0):
			New_Args[0][0] = arg
			flag = ""
		elif flag == "lib" or (flag == "" and acnt == 1):
			New_Args[0][1] = arg
			flag = ""
		elif flag == "ref":
			New_Args[0][2] = arg
			flag = ""
		elif flag == "corr":
			New_Args[0][3] = arg
			flag = ""
		elif flag == "write":
			if arg == "0":
				New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "W0" or (flag == "" and acnt == 2):
			if arg == "random":
				if New_Args[2][2] == 0:
					New_Args[2][2] = 10
			elif arg == "0":
				New_Args[2][2] = 0
			else:
#				try:
					parts = arg.split(",")
					W_new = []
					for number in parts:
						value = float(number)	
						W_new.append(value)
					New_Args[2][0] = W_new
#				except Exception:
#					Vprint(1, "@W0 takes only '0','random', or comma separated floats as an argument!")
#					New_Args[4] = 1
			flag = ""
		elif flag == "rep":
			try:
				New_Args[2][2] = int(arg)
			except Exception:
				Vprint(1, "@rep takes only integers as an argument!")
				New_Args[4] = 1
			flag = ""
		elif flag == "err":
			try:
				New_Args[2][4] = int(arg)
			except Exception:
				Vprint(1, "@rep takes only integers as an argument!")
				New_Args[4] = 1
			flag = ""
		elif flag == "SF0":
			try:
				New_Args[2][5] = float(arg)
			except Exception:
				Vprint(1, "@SF0 takes only floats as an argument!")
				New_Args[4] = 1
			flag = ""
		elif flag == "range":
			try:
				parts = arg.split(",")
				lmin = float(parts[0])
				lmax = float(parts[1])
				if lmax >= lmin:
					New_Args[2][1] = [lmin,lmax]
				else:
					New_Args[2][1] = [lmax,lmin]
			except Exception:
				Vprint(1, "@range takes only two comma separated floats as an argument!")
				New_Args[4] = 1
			flag = ""
		elif flag == "test":
			if arg in ["0","1","2"]:
				New_Args[2][3] = int(arg)
			else:
				Vprint(1,"@test takes only <0,1,2> as arguments!")
				New_Args[4] = 1	
			flag = "" 
		elif flag == "mult1":
			try:
				New_Args[3][0] = float(arg)	
			except Exception:
				Vprint(1, "@mult1 only takes float arguments")
				New_Args[4] = 1	
			flag = ""
		elif flag == "mult2":
			try:
				New_Args[3][1] = float(arg)	
			except Exception:
				Vprint(1, "@mult2 only takes floats arguments")
				New_Args[4] = 1	
			flag = ""
		elif flag == "iter":
			try:
				New_Args[3][2] = int(arg)	
			except Exception:
				Vprint(1, "@iter only takes integer arguments")
				New_Args[4] = 1	
			flag = ""
		elif flag == "nrmsd":
			if arg in ["0","1","2"]:
				New_Args[3][3] = int(arg)	
			else:
				Vprint(1, "@iter only takes integer arguments")
				New_Args[4] = 1	
			flag = ""
		elif flag == "verb":
			try:
				New_Args[5] = int(arg)	
			except Exception:
				Vprint(1, "@verb only takes integer arguments")
				New_Args[4] = 1	
			flag = ""
		else:
			Vprint(1, "Unknown argument:", arg)
			New_Args[4] = 1
		acnt += 1

	return New_Args

#function to load non-standard modules:
def Import_Custom():
	IMPORTS = []
	Vprint(2, "Loading modules:")
	try:
		globals()["np"] = __import__("numpy")
		IMPORTS.append("numpy")
	except ImportError:
		Vprint(2,"\nWarning, Could not import numpy module!\n")
	try:
		globals()["spo"] = __import__("scipy.optimize", globals(),locals(),["minimize"])
		IMPORTS.append("minimize")
	except ImportError:
		Vprint(2,"\nWarning, minimizer from scipy.optimize not found,importing SESCA simplex module... \n")
		try:
			if not SESCA_scripts in sys.path:
				sys.path.append(SESCA_scripts)
			globals()["min2"] = __import__("SESCA_min")	
			IMPORTS.append("minimize2")
		except ImportError:
			Vprint(2,"\nWarning, Could not import simplex module!")
	Vprint(2, IMPORTS)
	return IMPORTS


#read in spectrum
def Read_Spectrum(Infile,mode):
	Vprint(1, "\nReading in reference CD file:",)    
	Spectrum = []
	if os.path.isfile(Infile) == False:
		Vprint(1, "\nFile not found")
		return "None"
	Vprint(1, Infile)
	r = open(Infile,"rb")
	ref_num = 0
	for line in r:
		line2 = line.decode("ascii").strip("\n")
		if not line2.startswith('#') and line2 != '':
#		getting spectrum information
			try:		
				parts = line2.split()
				wavelength = float(parts[0])
				value = float(parts[1])
				ref_num = ref_num + 1
				entry = [wavelength,value]
				Spectrum.append(entry)
				Vprint(4, entry)
			except Exception:
				Vprint(2, "cannot read line:",line2)
	r.close()	
	return Spectrum

# read in basis set:
def Read_BS(Infile,L_range):
	Vprint(1, "\nReading in CD basis spectrum library:")    
	if os.path.isfile(Infile) == False:
		Vprint(1, "\nFile not found",Infile)
		return ("None")
	Vprint(1, Infile)
	BS_DATA = []
	lib_num = 0
	l = open(Infile,"rb")
	for line in l:
		line2 = line.decode("ascii").strip("\n")
		Vprint(5,line2)
#		reading spectrum information information		
		if line2.startswith("#|"):
			Labels = line2.strip("#|").split()
			for label in Labels:
				entry = [str(label), []]
				Vprint(4,entry)
				BS_DATA.append(entry)
		elif line2.strip() == '':
			valid = 0
		elif not line2.startswith('#') and line2.strip() != '':
			parts = line2.split()
			Vprint(5, parts)
			valid = 1
			if L_range != ["",""]:
				if L_range[0] != "" and L_range[0] > float(parts[0]):
					valid = 0
				if L_range[1] != "" and L_range[1] < float(parts[0]):
                                        valid = 0
			if valid == 1:
				try:		 
					cnt = 0
					for item in parts:
						value = float(item)
						BS_DATA[cnt][1].append(value)
						cnt += 1
					lib_num = lib_num + 1
					Vprint(4, "Library entry read:",line2)
				except Exception:
					Vprint(2, "Library, unable to process:",line2)

# getting wavelength boundaries for the basis spectra:		
	spect_num = len(BS_DATA)-1
	if BS_DATA[0][1] != []:
		lmin = BS_DATA[0][1][0]
		lmax = BS_DATA[0][1][-1]
	else:
		lmin = 0.0
		lmax = 0.0
	Vprint(3, "\n%1d Basis spectra read\nWvavelength range: %3.1f - %3.1f nm" % (spect_num,lmin,lmax))
	DATA = [BS_DATA,spect_num,lmin,lmax]	
# returning info:
	l.close()
	return DATA

#function to read error parameters from the basis set:
def Read_Error_Par(Infile):
	CALIB = []
	Vprint(2, "Reading error parameters")

	f = open(Infile, "rb")
	flag = 0
	for line in f:
		line2 = str(line.decode("ascii")).strip("\n")
#		reading Error Calibration ers (linear error model, outdated):
		if line2.startswith("#Calibration parameters"):
			flag = 3
		elif line2.startswith("#") and flag == 3:
			parts = line2.split()
			try:
				param = (parts[1], float(parts[3]),float(parts[5]))
				Vprint(4,"Library calibration parameters:",param)
				CALIB.append(param)	
			except Exception:
				Vprint(2,"Library, unable to process:",line2)
	
#		reading Error Calibration parameters (nonlinear error model):
		elif line2.startswith("#Calibration curves"):
			flag = 4
		elif line2.startswith("#RMSDj/") and flag == 4:
#			mark calibration curve, set up minimum and maximum values:
			header = ["Curve", "", "", 0]
			CALIB = [header]
			flag = 5
		elif line2.startswith("#") and flag == 5:
#			record curve data:
			try:
				parts = line2.split()
				Curve_entry = []
				e_num =len(parts)
				for e in range(1,e_num):
					Curve_entry.append(float(parts[e]))
				Vprint(4,"Library calibration parameters:",Curve_entry)
				CALIB.append(Curve_entry)

#				record minimum and maximum RMSD range for the calibration:			
				R0 = Curve_entry[0]
				if CALIB[0][1] == "" or CALIB[0][1] > R0:
					CALIB[0][1] = R0
				if CALIB[0][2] == "" or CALIB[0][2] < R0:
					CALIB[0][2] = R0
#				record the number of calibration curves:
				curve_num = int((len(CALIB[1])-1)/2)
				CALIB[0][3] = curve_num
			except Exception:
				Vprint(2,"Library, unable to process:",line2)	
	f.close()

#	print parameters:
	if CALIB != []:
		Vprint(3,"\nError calibration parameters:")
		for entry in CALIB:
			Vprint(3, entry)

	return CALIB


#function to read in a reference SS (from SESCA_pred) 
def Read_Ref(ref_file):
	Vprint(1, "\nReading reference file:")
#	check if file exists:
	if os.path.isfile(ref_file) == False:
		Vprint(1, "\nFile not found",ref_file)
		return "None"
#	read file_
	Ref_SS = []
	flag = 0
	f = open(ref_file,"rb")
	for line in f:
		line2 = line.decode("ascii").strip("\n")
		Vprint(5,line2)
#		get reference coefficients:
		if line.startswith("#Weighting factors"):
			flag = 2
		elif flag == 2 and line.startswith("#"):
			parts = line.strip("\n").split()
			Ref_SS.append(float(parts[-1]))
		elif flag == 2 and not line.startswith("#"):
			flag = 0
	f.close()
	Vprint(4, "Reference structure:", Ref_SS)
	
	return Ref_SS

#function to calculate SS compositions:
def Compare_SS(SS1, SS2):
	dim1 = len(SS1)
	dim2 = len(SS2)
#	check dimensions:
	if dim1 != dim2:
		Vprint(1, "SS compositions have different dimensions!")
		return "None"

#	compute dSS:
	dSS = 0.0
	for i in range(dim1):
		dSS += math.fabs(SS1[i]-SS2[i])/2

	return dSS

#function to intrapolate spectrum intensity at wavelength X
def Intrapolate(Function, X):
			lower = ["",0.0,0.0]
			upper = ["",0.0,0.0]
			entry = [X,0.0,0]
			for point in Function:
				dist =  math.fabs(entry[0]-point[0])
				if entry[0] > point[0] and (lower[0] == "" or lower[0] > dist):
					lower[0] = dist
					lower[1] = point[0]
					lower[2] = point[1]
			
				if entry[0] < point[0] and (upper[0] == "" or upper[0] > dist):
					upper[0] = dist
					upper[1] = point[0]
					upper[2] = point[1]
			if lower[0] != "" and upper[0] != "" and upper[1] != lower[1]:
				div = math.fabs(upper[1]-lower[1])
				Wlow = upper[0]/div
				Wup = lower[0]/div
				value_int = upper[2]*Wup + lower[2]*Wlow
				entry[2] = 1
				Vprint(4, "value intrapolated from:\n  P1(%2.1f,%2.1f)\n  P2(%2.1f,%2.1f)" % (lower[1],lower[2],upper[1],upper[2]))
			else:
				Vprint(2, "Error, could not intrapolate at value:", X)
				value_int = 0.0
				entry[2] = 0
			entry[1] = value_int
			return entry

# function to filter wavelengths:
def Process_Data(Spectrum_data, BS_data, Corr_data):
	Vprint(1, "\nProcessing input data:\n")
	BS_num = BS_data[1]
	lmin = BS_data[2]
	lmax = BS_data[3]
	BS = BS_data[0]
#	extract spectrum range for exp.data:
	Lrange_spectrum = []
	for entry in Spectrum_data:
		Lrange_spectrum.append(entry[0])
#	print Lrange_spectrum
	
#	select BS wavelengths which match:
	Labels =  []
	Exp_Int = []
	BS_Ints = []
	Corr_Int = []
	Waves = []
	wcnt = 0
	for Wavelength in BS[0][1]:
		if Wavelength in Lrange_spectrum:
			Vprint(3, Wavelength, "added")
			Waves.append(Wavelength)
			TMP = []
#			if we have a match add spectrum intensity to exp. data
			for entry in Spectrum_data:
				if entry[0] == Wavelength:
					Exp_Int.append(entry[1])
			corr_found = 0
			for entry in Corr_data:
				if entry[0] == Wavelength:
					Corr_Int.append(entry[1])
					corr_found += 1
			if corr_found == 0:
				Corr_Int.append(0.0)
					
#			if we have a match add BS intensity values to Basis spectrum data
			for Column in BS:
				value = Column[1][wcnt]
				if value != Wavelength:
					TMP.append(value)
			Vprint(3, TMP)
			BS_Ints.append(TMP)	
		if wcnt == 0:
			for Column in BS:
				Labels.append(Column[0])	
		wcnt += 1
#	turning data into numpy arrays:
	Labels = np.array(Labels)
	BS_Ints = np.array(BS_Ints)
	Waves = np.array(Waves)
	Corr_Ints = np.array(Corr_Int)

	Vprint(3, Labels)
	Vprint(4, Exp_Int)
	Vprint(4, BS_Ints)
	Processed_data = [Exp_Int, Corr_Int, BS_Ints,Labels,Waves]
	return Processed_data

#function to produce calculated CD spectra:
def Compute_Spectra(Weights, BS_data, Corr_Spect):
	Spectrum = np.array([])
#	first compute the base spectrum 
	Base_Spect = np.dot(Weights, BS_data)

#	compute normalization term, and modify correction:
	if Corr_Spect != []:
		Norm = np.sum(Weights)
		Corr = np.multiply(Corr_Spect, Norm) 
#	note: correction term is exact if wieghts are normalized!

#	determine final spectrum
	Spectrum = np.add(Base_Spect, Corr)
	
	return Spectrum

#function to compute the estimation uncertainty:
def Compute_Error(Dev, Error_Par, curve):
	Header, e_num = Error_Par[0], len(Error_Par)
	Err_data = []
	check = 0
	Vprint(2, "Estimating uncertainty:")	
	Vprint(5, "Header:", Header)

#	Calculate error estimate with calibration curve
	if Header[0] == "Curve":
		Vprint(3,"Nonlinear error calibration")
		Error = []
		Bars = []
#		make sure we are in the calibration range
		Vprint(4,Header)
		Rmin, Rmax = Header[1],Header[2]
		if Dev < Rmin: 
			Dev = Rmin
			check = 2
		elif Dev > Rmax:
			Dev = Rmax
			check = 3
		else:
			check = 1
#		make sure a valid calibration curve is chosen
		if curve > Header[3]:
			k = 3
			Vprint(1,"Requested calibration curve not found, falling back to default")
		else:
			k = (curve-1)*2 + 1
#		determine estimates based on calibration:
		Vprint(3, "Using curve %1d"%curve)
		for j in range(1,e_num):
			E_data = [Error_Par[j][0], Error_Par[j][k]]
			B_data = [Error_Par[j][0], Error_Par[j][k+1]]
			Error.append(E_data)
			Bars.append(B_data)
		Err_est = Intrapolate(Error,Dev)[1]
		Bar_est = Intrapolate(Bars,Dev)[1]
		Err_data = [Err_est, Bar_est, curve, check]

#	Calculate error based on linear fits:
	if Header[0] == "mfit":
		Vprint(3,"Linear error estimation")
		mf,Sm = Error_Par[0][1], Error_Par[0][2]
		Rf,Sr = Error_Par[1][1], Error_Par[1][2]
		Sdev  = Error_Par[2][1]
		Fit = [mf,Sm, Rf, Sr, Sdev]
		Vprint(4, Fit)
		Error_est = float(Dev)/mf
		Bar_est = float(Sdev)
#		make sure the Error estimate makes sense:
		if Error_est < 0.0:
			check = -2
			Error_est = 0.0
		elif Error_est > 1.0:
			check = -3
			Error_est = 1.0
		else:
			check = -1
		Err_data = [Error_est, Bar_est, mf, check]

#	relay the data:					
	Vprint(3, "Estimated error:", Err_data)
	return Err_data

#helper function to multiply array elements by a scalar:
def Amult(Array, Scalar):
	Tmp = []
	for Element in Array:
		New_element = Scalar * Element
		Tmp.append(New_element)
	new_array = np.array(Tmp)
	return new_array

#function to calculate spectrum RMSD:
def Rmsd_Calc(Exp_data,Calc_data):
	data_num = len(Exp_data)
	if data_num != 0:
		Diff_data = np.subtract(Exp_data, Calc_data)
		Diff_square = np.multiply(Diff_data,Diff_data)
		Sum_sq = np.sum(Diff_square)
		Rmsd = math.sqrt(Sum_sq/data_num)
	else:
		Rmsd = 0.0
	return Rmsd

#function to calculate intensity normalized spectrum RMSD:
def Nrmsd_Calc(Exp_data,Calc_data):
	data_num = len(Exp_data)
	if data_num != 0:
#		Calclate RMSD as normal
		Diff_data = np.subtract(Exp_data, Calc_data)
		Diff_square = np.multiply(Diff_data,Diff_data)
		Sum_sq = np.sum(Diff_square)
		Rmsd = math.sqrt(Sum_sq/data_num)
#		divide by average amplitude
		Exp_sq = np.multiply(Exp_data,Exp_data)
		Sum_exp = np.sum(Exp_sq)
		Norm = math.sqrt(Sum_exp/data_num)
		Nrmsd = float(Rmsd)/Norm
	else:
		Nrmsd = 0.0
	return Nrmsd

#function to calculate spectrum mean unsigned error:
def Muse_Calc(Exp_data,Calc_data):
	data_num = len(Exp_data)
	if data_num != 0:
		Diff_data = np.subtract(Exp_data, Calc_data)
		Sum_abs = 0.0
		for value in Diff_data:
			Sum_abs += math.fabs(value)
		Muse = Sum_abs/data_num
	else:
		Muse = 0.0
	return Muse

#function to renormalize weights:
def Renorm_Weights(Weights):
	dim = len(Weights)
	W_new = []
#	first set any negative weights to 0:
	Sum = 0.0
	for Coeff in Weights:
		if Coeff < 0.0:
			W_new.append(0.0)
		else:
			W_new.append(Coeff)
			Sum += Coeff
#	normalize coeffs by the sum:
	if Sum != 0.0:
		for i in range(dim):
			W_new[i] = float(W_new[i])/Sum
			Scale = float(1.0)/Sum
	else:
		W_new = Weights
		Scale = 1.0

	Renorm = [W_new,Scale]
	
	return Renorm
		
#print fitted results:
def Print_fit(Waves,Exp_data,BS_data,Weights,Corr_data,Scaling):
	cnt = 0
	Sjl = Compute_Spectra(Weights,BS_data,Corr_data)
	Exp_scaled = np.multiply(Scaling, Exp_data)
	Sdiff = np.subtract(Exp_scaled, Sjl)
	All_lines =      "#Wl(nm)   Iexp     Icalc     Idev\n"
	for Wl in Waves:
		exp = Exp_scaled[cnt]
		calc = Sjl[cnt]
		diff = -1*Sdiff[cnt]
		data = (Wl,exp,calc,diff)
		string = "%1.1f   %1.4f    %1.4f    %1.4f\n" % data
		All_lines += string
		cnt += 1
	return All_lines

#function to print output data:
def Write_Outfile(Filenames, Aux_data, Exp_Int, Corr_Int, Bil, Results):
		workdir, infile, libfile, reffile, corrfile  = Filenames
		lnum_final, BS_num, Labels, Waves = Aux_data
		W_final, Rmsd_final, Muse_final, Nrmsd_final, Scaling, Error, dSS = Results
		Output = "#SESCA CD deconvolution results:\n#Workdir: %1s\n#Spectrum: %1s\n#Basis set: %1s\n" % (workdir,infile,libfile)
#		add reference file:
		if not reffile == "":
			Output += "#Reference: %1s\n" % reffile
#		add correction file:
		if not corrfile == "":
			Output += "#Correction: %1s\n" % corrfile

#		add scaling factor:
		if Scaling != 1.0:
			Output += "#Scaling factor from weight normalization: %1.3f\n" % Scaling
#		write SS coefficients:
		Output += "\n#Weighting factors for calc. CD spectrum:\n"
		for i in range(BS_num):
			data = (Labels[i+1],W_final[i])
			W_string = "# %12s  :   %1.4f\n" % data
			Output += W_string
		Output += "\n"
#		write original and calculated spectra:
		Output += Print_fit(Waves, Exp_Int, Bil, W_final, Corr_Int, Scaling)
		Rmsd_string = "\n#matches     RMSD     MUSE    NRMSD\n#    %1d     %1.3f    %1.3f    %1.3f\n" % (lnum_final, Rmsd_final, Muse_final, Nrmsd_final)
		Output += Rmsd_string
#		write uncertainty of the estimate:
		if Error != []:	
			Err_string = "\n#Typical uncertainty of the SS estimate based on residual error:\n"
			if Error[3] > 0:
				Err_string += "#Uncertianty determined from calibration curve %1d\n"%Error[2]
				if Error[3] == 2:
					Err_string = "#Warning, spectral deviation (RMSD: %1.3f) is smaller than the calibration range!\n" % Rmsd_final
				if Error[3] == 3:
					Err_string = "#Warning, spectral deviation (RMSD: %1.3f) is larger than the calibration range!\n" % Rmsd_final
			else:
				Err_string += "#Uncertianty determined from linear fit (mf: %2.1f kMRE)\n"%Error[2]
				if Error[3] == -2:
					Err_string = "#Warning, uncertainty is smaller than 0%, check error parameters\n"
				if Error[3] == -3:
					Err_string = "#Warning, uncertainty is larger than 100%, check error parameters\n"
			Err_string += "#dSS-est: %2.1f +/- %2.1f" % (Error[0]*100, Error[1]*100) + " %\n"
			
			Output += Err_string

#		write deviation from reference:
		if not dSS in ["","None"]:
			SS_string = "\n#Calculated deviation from reference structure:\n#dSS-ref: %2.1f "%(dSS*100) + "%\n"
			Output += SS_string

		return Output                
	
#function to generate random weights:
def Generate_Rweights(dimension,Wmin,Wmax):
	Weights = []
	Rand_limit = [Wmin,Wmax]
	for i in range(dimension-1):
		value = random.uniform(Rand_limit[0],Rand_limit[1])
		Weights.append(value)
		Rand_limit[1] -= value
	Weights.append(Rand_limit[1])
	Weights = np.array(Weights)
	np.random.shuffle(Weights)	
	return Weights

#function to Minimize RMSD as the function of Coeffs:
def Eval_function(Weights, Args):
#	import determined parameters:	
	Sexp, Bil, Scorr, Lambda1, Lambda2, Norm = Args
	Scalc = Compute_Spectra(Weights, Bil, Scorr)
#	Lagrange multipliers:
#	multiplier to keep the sum of all weights 1.0
	LG_mult = Lambda1 * math.fabs((1.0 - np.sum(Weights)))
#	extra terms to keep all weights positive:
	for Wi in Weights:
		if Wi >= 0.0:
			LG_mult += 0.0
		else:
			LG_mult += -Lambda2 * Wi
#	full function to return (including LG mulitplier):
	if Norm == 0:
		Full_function = Rmsd_Calc(Sexp, Scalc) + LG_mult
	if Norm == 1:
		Full_function = Nrmsd_Calc(Sexp, Scalc) + LG_mult
	return Full_function
	
def Run_function_tests(Winit,Exp_Int,Corr_Int,Bil,BS_num,Waves):
#	test evaluation function to make sure it works:
		Vprint(1, "\nEval. function test:")
		test_args = [Exp_Int, Bil, Corr_Int, 1000, 1000, 0]
		Value_init = Eval_function(Winit,test_args)
		Value_test = 2.593
		Scaling = 1.0
		Vprint(1, "Test weights:",Winit,"\nVexpected: %1.3f\nVcalculated: %1.3f" %(Value_test,Value_init))
		if math.fabs(Value_init-Value_test) <= 0.001:	
			Vprint(1, "Passed!")
		else:
			Vprint(1, "Failed!")
		Vprint(1,"\nTesting function with random weights:")
		for i in range(10):
			Weight_R = Generate_Rweights(BS_num,0.0,1.0)
			Value = Eval_function(Weight_R,test_args)
			string = "%1d [" % i + (BS_num*"%1.3f, ")%tuple(Weight_R) + "] Value: %1.3f" %Value
			Vprint(2, string)
#	test RMSD function:
		Sjl = np.dot(Winit,Bil)
		Sdiff = np.subtract(Exp_Int, Sjl)
		Vprint(1, "\nTesting CD/RMSD functions:")
		cnt = 0
		Sum = 0.0
		for Wl in Waves:
			exp = Exp_Int[cnt]
			calc = Sjl[cnt]
			diff = -1*Sdiff[cnt]
			data = (Wl,exp,calc,diff)
			string = "%1.1f %1.4f %1.4f %1.4f" % data
#			print string
			Sum += diff**2
			cnt += 1
		Vprint(1, Print_fit(Waves, Exp_Int, Bil, Winit, Corr_Int, Scaling))
		Rmsd = math.sqrt(Sum/cnt)
		Vprint(1, "#Final Rmsd:   %1.3f" % Rmsd) 
		Rmsd2 = Rmsd_Calc(Exp_Int,Sjl)
		Vprint(1, "#Rmsd2 check:  %1.3f" % Rmsd2)	
#	test random weights:
		Vprint(1,"\n#Random weights:")
		for i in range(10):
			Weight_R = Generate_Rweights(BS_num,0.0,1.0)
			W_sum = np.sum(Weight_R)
			Vprint(2, ("("+BS_num*"%1.3f "+")") % tuple(Weight_R), "sum: %1.3f" % W_sum)




		
#Main function for script execution:	
def Deconv_Main(Args):
# 	read default parameters
	Inp_files,  Out_files, Dec_param, Opt_param, failmark, verbosity = Args
	infile, libfile, reffile, corrfile = Inp_files
	outfile = Out_files[0]
	W_custom, L_range, random_seed, test, error, SF0 = Dec_param
	Lambda1, Lambda2, i_steps, nrmsd = Opt_param
	Set_verb(Args[5])


#	Import custom modules:
	IMPORTS = Import_Custom()


	#read in spectrum:
	Spect = []
	if failmark == 0 and test != 1:
		Spect = Read_Spectrum(infile,1)
		if Spect == "None": 
			failmark = 1

#	read in base line correction:
	Corr_Spect = []
	if corrfile != "" and failmark == 0 and test != 1:
		Corr_Spect = Read_Spectrum(corrfile, 1)
		if Corr_Spect == "None":
			failmark = 1
	 
	#read in basis spectra:
	BASIS = []
	EPAR =  []
	if libfile != "" and failmark == 0 and test != 1:
		BASIS = Read_BS(libfile,L_range)
		EPAR = Read_Error_Par(libfile)
		if BASIS == "None":
			failmark = 2

#	Parse Input information:
#	print BASIS
	if failmark != 0:
		print("\nError while reading input files, script stops!") 
		print(Usage) 
		sys.exit(failmark)


	if test == 1:
#		use spectrum CD_LYSM.out and basis set Map_DS-dT.dat
#		in wavelength range 200 - 205 nm
#		resultig RMSD should be 0.592 kMRE
#		Eval function should be 2.593 
		Winit = np.array([0.365, 0.04, 0.597])
		Labels = np.array(['Wavelength','Alpha','Beta','Coil'])
		Waves = np.array([  200.0,   201.0,   202.0,   203.0,   204.0,   205.0])
		Exp_Int_test = 	 [ -0.726,  -4.254,  -7.388,  -9.993, -11.972, -13.324]
		BS_test = [	 [ 27.445,  16.789,   7.005,  -1.956, -10.195, -17.508],
				 [ 48.892,  47.191,  44.444,  40.538,  35.539,  29.916],
				 [-21.712, -21.222, -20.406, -19.102, -17.323, -15.264]]
		Corr_test = np.array([ 0.0,  0.000,   0.000,  0.000,    0.000,   0.000])	
		Exp_Int = np.array(Exp_Int_test)
		Bil = np.array(BS_test)
		lnum_final, exp_check, BS_check, BS_num = [6,6,6,3]
		outfile = ""
		Run_function_tests(Winit, Exp_Int, Corr_test, Bil, BS_num, Waves)
		return "None"
	
	
	Parsed = Process_Data(Spect,BASIS, Corr_Spect)
#	declare internal variables:
	Exp_Int = Parsed[0]
	Corr_Int = Parsed[1]
	BS_Ints = Parsed[2]
	Labels = Parsed[3]
	Waves = Parsed[4]
	Bil = BS_Ints.transpose()
	lnum_final = len(Waves)
	exp_check = len(Exp_Int)
	BS_check = len(BS_Ints)
	BS_num = len(BS_Ints[0])

#	check dimensionality:
	if lnum_final != BS_check or lnum_final != exp_check:
		Vprint(1, "Error mismatch in array dimensions:")
		Vprint(1, "lnum: %1d, exp_int: %1d, BS_ints: %1d" % (lnum_final,exp_check,BS_check))
		failmark = 3
		print(Usage)
		sys.exit(failmark)
	else:
		Vprint(2, "number final wavelengths: %1d" % lnum_final)
			
# define initial spot:
	Seeds = []
	if W_custom != []:
		Seeds.append(np.array(W_custom))
	elif random_seed >= 1:
		for j in range(random_seed):
			W0 = Generate_Rweights(BS_num,0.0,1.0)
			Seeds.append(W0)
	else:
		Seeds.append(np.array([0.365, 0.04, 0.597]))

#	scale coefficients if SF0 was requested:
	if SF0 != 1.0:
		s_num = len(Seeds)
		c_num = len(Seeds[0])
		for s in range(s_num):
			for c in range(c_num):
				Seeds[s][c] = Seeds[s][c]/SF0
			

# execute deconvolution:	
	min_rep = False
	disp = 0
	if verbosity >= 3:
		min_rep = True
		disp = 1

	W_best = []
	W_final = []
	W_min = []
	cnt = 0
	for W0 in Seeds:
#		check starting point:
		Min_Args = [Exp_Int, Bil, Corr_Int, Lambda1, Lambda2,  nrmsd]
		Value = Eval_function(W0,Min_Args)
		Scaling = 1.0
		Vprint(3, "\nInitial weights:",("["+BS_num*"%1.3f "+"]")%tuple(W0)," Rmsd: %1.3f" % Value)
#		calculate best fit from each initial weight:
		if test != 2 and "minimize" in IMPORTS:
#	 		optimize weights with scipy simplex minimizer:
			Vprint(2, "\nFitting weights Sicpy.optimize:")
			Min_result = spo.minimize(Eval_function, W0, Min_Args, method='nelder-mead',
				 options={"ftol":1e-15,  "disp":min_rep, "maxiter":i_steps, "adaptive":True})
			W_min = Min_result.x
			Vprint(3, Min_result)
		elif test != 2 and "minimize2" in IMPORTS:
#	 		if scipy is not available use SESCA simplex minimizer instead:
			Vprint(2, "\nFitting weights with SESCA_min:")
			adaptive,tol = [1, 1e-15]	
			Min_param = [adaptive,tol, i_steps,disp]
			NM_param = min2.Pass_NMparam()
			Min_result = min2.Iterate_Min(Eval_function, W0, Min_Args, Min_param, NM_param)
			if Min_result != None:
				W_min = Min_result[2]
				Vprint(3, Min_result)
			else:
				W_min = W0
				Vprint(1, "Minimization failed, weights reset to:",W_final)				


		elif not "minimize" in IMPORTS and not "minimize2" in IMPORTS:
			Vprint(1, "\nMinimization was called but no minimizer was imported!")
			W_min = W0
			break
		else:
#			test input  weights without optimization
			W_min = W0

#		Calculate fit
		cnt += 1
		Calc_min = Compute_Spectra(W_min, Bil, Corr_Int)
		Rmsd_min = Rmsd_Calc(Exp_Int,Calc_min) 
		Muse_min = Muse_Calc(Exp_Int,Calc_min)
		Nrmsd_min = Nrmsd_Calc(Exp_Int,Calc_min)
		W_seed = [Rmsd_min, Nrmsd_min, Muse_min, W_min]
		Vprint(2, "Seed %1d: %1.3f "%(cnt,W_seed[nrmsd]),W_min) 
		if W_best == [] or W_best[nrmsd] > W_seed[nrmsd]:
			W_best = W_seed 

#	renormalize coefficients, in case the constraints are lifted:
	Rmsd_best,Nrmsd_best,Muse_best, W_uns = W_best
	W_final, scaling  = Renorm_Weights(W_uns)
	Exp_scaled = np.multiply(scaling, Exp_Int)
	Calc_final = Compute_Spectra(W_final, Bil, Corr_Int)

#	calculate deviation scores:
	Rmsd_final = Rmsd_Calc(Exp_scaled, Calc_final) 
	Muse_final = Muse_Calc(Exp_scaled, Calc_final)
	Nrmsd_final = Nrmsd_Calc(Exp_scaled, Calc_final)
	
#	print the final coefficients:
	Vprint(4, W_final)
	Vprint(2, "\nFinal weights:",("[ "+BS_num*"%1.3f "+"]")%tuple(W_final), "with scaling factor: %1.3f" % scaling)
	Vprint(2, " Rmsd: %1.3f,  Muse: %1.3f, Nrmsd: %1.3f" % (Rmsd_final, Muse_final, Nrmsd_final))
	Data = [W_final, Rmsd_final, Muse_final, Nrmsd_final, scaling]

#	compute typical uncertainty:
	Error_est = []
	if error != 0 and EPAR != []:
		Error_est = Compute_Error(Rmsd_final, EPAR, error)
	Data.append(Error_est)

#	Compare to reference structure if provided:
	dSS = ""
	if reffile != "":
#		read reference file:
		Ref_SS = Read_Ref(reffile)
#		compute SS deviation: 
		dSS = Compare_SS(Ref_SS, W_final)
		if not dSS in ["","None"]:
			Vprint(3, "deviation between estimated SS and reference: %1.4f"%dSS )
		else:
			dSS = ""
#	# append deviation data:
	Data.append(dSS)	

# 	Write output file:
	if outfile != "":
		Vprint(1, "Writing output data\n")
		Filenames = [workdir, infile, libfile, reffile, corrfile] 
		Aux_data = [lnum_final, BS_num, Labels, Waves] 
		Output_data = Write_Outfile(Filenames, Aux_data, Exp_Int, Corr_Int, Bil, Data)
		Output = Output_data.encode("ascii")
		o = open(outfile,"wb")
		o.write(Output)
		o.close()
		
	return Data

# executing standalone script:
if __name__ == '__main__':
#	handling command line arguments:	
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Inp_files,  Out_files, Dec_param, Opt_param, failmark, verbosity = Custom_Args
	outfile = Out_files[0]
	Set_verb(Custom_Args[5])
	
	Vprint(3,"\nDeconvolution module\nRun parameters:", Custom_Args)
	
#	executing main code
#	Data_main = Deconv_main(W_custom, infile, libfile, outfile, L_range, test)
	Data_main = Deconv_Main(Custom_Args)
#	print Data_main
	if test == 1:
		outfile = ""		

	ftime = time.time()
	runtime = ftime-stime
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",outfile)
else:
#	check dependencies, raise error if they are missing:
	IMPORTS = Import_Custom()
	if len(IMPORTS) >= 2:	
		Vprint(2,"SESCA CD deconvolution module (CD_deconv.py)")
	else:
		Vprint(1, IMPORTS)
		raise ImportError("Failed to load dependencies")
	
