#!/usr/bin/env python
#

import sys
import os
import math
import time

workdir = os.getcwd()
stime = time.time()


#script to calculate and compare theoretical CD based on secondary structure content
usage0  = "*******************************************************************\n"
usage   = "SESCA module for calculating theoretical CD spectra from secondary structure distributions\n"
usage  += "\nusage: CD_calc.py (<reference_file>) <target_file> @lib <basis set file> @flag <argument>\nPossible command flags are:\n"
usage  += "   @ref  <ref_file> specify reference file (CD spectrum)\n   @tar  <target_file> specify target file (structure composition file or CD spectrum)\n"
usage  += "   @lib  <BS_file> specify basis spectrum library (default is none, see libs subdir)\n"
usage  += "   @write <output_file> specify output file name (default: CD_comp.out)\n"
usage  += "   @mode <int> set work mode (defaul:t 1)\n      0 - generate CD spectrum only (no ref_file)\n"
usage  += "      1 - generate CD spectrum and compare to reference\n      2 - compare 2 CD spectra (no lib_file)\n"
usage  += "      3 - sum or scale CD spectra (no lib_file)\n      4 - recalculate spectra (from weights)\n"
usage  += "   @scale <0 / float> use scaling factor for the calculated CD spectrum, 0 means no scaling (default: 1.0)\n"
usage  += "   @range <float,float> limit wavelength range to work in (default: none)\n"
usage  += "   @norm <0,1 / float> normalize coefficients to 100%, if a float other than 1.0 is provided,\n"
usage  += "         the calculated CD spectrum is scaled by that amount (default: 0 - off)\n"
usage  += "   @err <int> select calibration curve for model error estimation (found in the basis set), 0 - no error estimation (default: 1)\n"
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 3)\n"

Usage = usage0 + usage + usage0


#default parameters:
ref_file = ""
ss_file = ""
lib_file = ""
out_file = "CD_comp.out"
mode = 1
L_range = ["",""]
sc_factor = 1.0
norm = 0
nrmsd = 0
L_range = ["",""]
error = 1

failmark = 0
verbosity = 3


Inp_files = [ss_file, lib_file, ref_file]
Out_files = [out_file]
Param = [mode, L_range, sc_factor, norm, nrmsd, error]
Def_Args = [Inp_files, Out_files, Param, failmark, verbosity]


#function definitions:
#function to pass on defaults:
def Pass_Defaults():
	return Def_Args

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

#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
	Inp_files, Out_files, Param, failmark, verbosity = Pass_Defaults()
#	New_Args = [ref_file, ss_file, out_file, lib_file, mode, L_range, scale, norm, sc_factor, write, failmark, verbosity]
	New_Args = [Inp_files, Out_files, Param, failmark, verbosity]
	ss_file, lib_file, ref_file = Inp_files
	out_file = Out_files[0]
	mode, L_range, sc_factor, norm, nrmsd, error = Param


	FLAGS = ["ref","tar","write","lib","mode","range","scale","norm","nrmsd","err","verb"]
	flag = ""
	acnt = 0
#	processing new passed arguments: 
	Vprint(2, "Recognized flags:")
	for arg in Args:
		if arg.startswith("@"):
			flag = arg.strip("@")
			if flag in FLAGS:
				Vprint(2, flag)
			else:
				Vprint(1,"Unknown flag:",flag)
				New_Args[3] = 1
#		standard I/O flags:
		elif flag == "tar":
			New_Args[0][0] = arg
			flag = ""
		elif flag == "lib":
			New_Args[0][1] = arg
			flag = ""
		elif flag == "ref":
			New_Args[0][2] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
				New_Args[1][0]= ""
			else:
				New_Args[1][0] = arg
			flag = ""
#		calculation parameters:
		elif flag == "mode":
			if arg in ["0","1","2","3","4"]:
				New_Args[2][0] = int(arg)
			else:
				Vprint(1, "invalid argument for @mode") 
				New_Args[3] = 1 
			flag = ""
		elif flag == "range":
			try:
				parts = arg.split(",")
				lower = float(parts[0])
				upper = float(parts[1])
				if upper >= lower:
					New_Args[2][1] = [lower,upper]
				else:
					New_Args[2][1] = [upper, lower]
			
			except Exception:
				Vprint(1, "@range only takes two comma-separated float arguments")  
				New_Args[3] = 1 
			flag = ""
		elif flag == "scale":
			if arg == "0":
				New_Args[2][2] = 0
			else:
				try:
					New_Args[2][2] = float(arg)
				except Exception:
					Vprint(1, "@scale only takes float arguments or '0'") 
					New_Args[3] = 1 
				flag = ""
		elif flag == "norm":
			if arg in ["0","1"]:
				New_Args[2][3] = int(arg)
			else:
				Vprint(1, "@norm only takes <1,0> as arguments")
				New_Args[3] = 1		
			flag = ""
		elif flag == "nrmsd":
			if arg in ["0","1","2"]:
				New_Args[2][4] = int(arg)
			else:
				Vprint(1, "@nrmsd only takes <1,0> as arguments")
				New_Args[3] = 1		
			flag = ""
		elif flag == "err":
			try:
				New_Args[2][5] = int(arg)
			except Exception:
				Vprint(1, "@err only takes integer arguments")
				New_Args[3] = 1		
			flag = ""
		elif flag == "verb":
			try:
				New_Args[4] = int(arg)
			except Exception:
				Vprint(1, "@verb only takes integer arguments") 
				New_Args[3] = 1	
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0][2] == "" and acnt == 0:
			New_Args[0][2] = arg
			flag = ""
		elif flag == "" and New_Args[0][0] == "" and acnt == 1:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[0][1] == "" and acnt == 2:
			New_Args[0][1] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args
		
#Function to read in structure data:
def Read_Struct_file(File):
	flag = 0
	Classes = []
	Vprint(1,"\nReading Structure File:",)
	if os.path.isfile(File) == False:
		Vprint(1, "\nFile not found")
		return "None"

	Vprint(1, File)
	i = open(File,"rb")
	for line0 in i:
		line = str(line0.decode("ascii")).strip("\n")
		if line.startswith("#Secondary structure"):
			flag = 1
		elif line.startswith("# Total residues counted"):
			flag = 1
		elif line.startswith("#Weighting factors"):
			flag = 2
		elif line.startswith('# ') and flag == 1:
			try:
				parts = line.strip("\n").split()
				ID = parts[1].strip(":")
				value = float(parts[2])
				data = (ID,value)
				Classes.append(data)
				Vprint(3, data)
			except Exception:
				Vprint(2, "Cannot to read line:",line)
		elif line.startswith('# ') and flag == 2:
			try:
				parts = line.strip("\n").split()
				ID = parts[1]
				value = float(parts[3])*100.0
				data = (ID,value)
				Classes.append(data)
				Vprint(3, data)
			except Exception:
				Vprint(2, "Cannot to read line:",line)
	i.close()
	return Classes


#function to read in spectrum data:
def Read_Spectrum_file(File,L_range):
	Vprint(1,"\nReading in spectrum file:",)
	Spectrum = []
	if os.path.isfile(File) == False:
		Vprint(1,"\nFile not found")
		return "None"

	Vprint(1, File)
	r = open(File,"rb")
	for line in r:
		line2 = str(line.decode("ascii")).strip("\n")
		if not (line2.startswith('#') or line2.startswith(";")) and line2 != '':
# 		getting spectrum information
			try:
				parts = line2.split()
				wavelength = float(parts[0])
				value = float(parts[1])
				entry = [wavelength,value]
				if L_range == ["",""]:
					Spectrum.append(entry)
					Vprint(4,entry)
				elif entry[0] >= L_range[0] and entry[0] <= L_range[1]:
					Spectrum.append(entry)
					Vprint(4, entry)
				else:
					Vprint(2, entry, "Out of range")
						
			except Exception:
				Vprint(2,"cannot read line:",line2)
	r.close()
	return Spectrum

# function to read in basis set parameters:
def Read_BS(File,L_range):
	Vprint(1,"\nReading in CD basis spectrum set:",)    
	if os.path.isfile(File) == False:
		Vprint(1,"\nFile not found: %1s"%File)
		return ("None")

	Vprint(1,File)
	BS_DATA = []
	lib_num = 0
	ASSIGN = []
	flag = 0
	COMBMAT = []
	class_num = 0
	BS_num = 0
	CALIB = []
	l = open(File,"rb")
	for line in l:
		line2 = str(line.decode("ascii")).strip("\n")
#		reading basis spectrum intensities:	
		if line2.startswith("#|"):
			flag = 0
			Labels = line2.strip("#|").split()
			for label in Labels:
#				create an array for each basis spectrum
				entry = [label, []]
				BS_DATA.append(entry)
			
		elif not line2.startswith('#') and line2.strip() != '':
			parts = line2.split()
			valid = 0
			try:
#				determine if this line is within the requested wavelength range:
				waveL = float(parts[0])	
				if L_range == ["",""]:
					valid = 1
				elif L_range[0] <= waveL and L_range[1] >= waveL:
					valid = 1
				if valid == 1:
#					add intensities to the basis spectra if the wavelength is OK
					cnt = 0
					for item in parts:
						value = float(item)
						BS_DATA[cnt][1].append(value)
						cnt += 1
					lib_num = lib_num + 1
					Vprint(4,"Library entry read:",line2)
			except Exception:
				Vprint(2,"Library, unable to process:",line2)
		elif line2 == "":
			valid = 0
			flag = 0

		#reading combination matrix (if available):
		elif line2.startswith("#Basis/") and COMBMAT == []:
#			record SS element codes:
			flag = 1
			parts = line2.split()
			Codes = ["Codes",parts[1:]]
			COMBMAT.append(Codes)
		elif line2.startswith("#") and flag == 1:
			try:
#				add assignent factors for each basis spectrum / SS element:
				first, second = line2.split(":")
				name = first.split()[1]
				values = second.split()
				for i  in range(len(values)):
					val = float(values[i])
					values[i] = val
				comb_entry = [name,values]
				COMBMAT.append(comb_entry)	
				Vprint(4,"Library combination map read:",comb_entry)
			except Exception:
				Vprint(2,"Library, unable to process:",line2)

#		reading Assignment information (if no combination matrix is available):
		elif line2.startswith("#Mapping") and COMBMAT == []:
			flag = 2
			Codes = ["Assign", []]
			ASSIGN.append(Codes)
		elif line2.startswith("#") and flag == 2:
			try:
				first, second = line2.split(":")
				name = first.split()[1]
				codes = second.split()
				map_entry = (name,codes)
				ASSIGN.append(map_entry)
				for Code in codes:
					if not Code in ASSIGN[0][1]:
						ASSIGN[0][1].append(Code)
				Vprint(4,"Library mapping read:",map_entry)
			except Exception:
				Vprint(2,"Library, unable to process:",line2)

#		reading Error Calibration parameters (linear error model, outdated):
		elif line2.startswith("#Calibration parameters"):
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
	l.close()

#	process basis set data:
	if ASSIGN != []:
		Vprint(3,"\nAssignment:")
		for entry in ASSIGN:
			Vprint(4,entry)

#	generate combination matrix from the assignment if necessary:
	if COMBMAT == [] and ASSIGN != []:
		COMBMAT.append(ASSIGN[0])
		COMBMAT[0][0] = "Codes"
		for i in range(1,len(ASSIGN)):
			Row = [ASSIGN[i][0],[]]
			for entry in COMBMAT[0][1]:
				if entry in ASSIGN[i][1]:
					Row[1].append(1.0)
				else:
					Row[1].append(0.0)
			COMBMAT.append(Row)	
		

	if COMBMAT != []:		
		Vprint(3,"\nCombination matrix:")		
		for entry in COMBMAT:
			Vprint(4,entry)
		BS_num = len(COMBMAT)-1
		class_num = len(COMBMAT[0][1])
	else:
		Vprint(1,"No assignment information was detected!")
		Vprint(1,"Please provide either an Assignment or Combination matrix block!")
		return ("None") 

	if CALIB != []:
		Vprint(3,"\nError calibration parameters:")
		for entry in CALIB:
			Vprint(4, entry)

# getting wavelength boundaries for the basis spectra:		
	spect_num = len(BS_DATA)-1
	if BS_DATA[0][1] != []:
		lmin = BS_DATA[0][1][0]
		lmax = BS_DATA[0][1][-1]
	else:
		lmin = 0.0
		lmax = 0.0

	if spect_num == BS_num:
		Vprint(3,"\n%1d Basis spectra read\nWavelength range: %3.1f - %3.1f nm" % (spect_num,lmin,lmax))
		DATA = [ASSIGN, COMBMAT, BS_DATA, CALIB, spect_num,lmin,lmax, class_num]
	
	else:
		Vprint(1,"Mismatch between the number of assigned classes (%1d) and the number of basis spectra (%1d)!" %(spect_num,BS_num))
		Vprint(1,"Please double check your basis set file!")
		DATA = ("None")

	
	return DATA


#function to generate theoretical CD spectrum:
def Generate_CD(SS_data, BS_data):
#		Check if SS classes are matching:
		Comb_Matrix = BS_data[1]
		Basis_Spect = BS_data[2]
		Labels = BS_data[1][0][1]
		BS_num = BS_data[4]
		SS_num = BS_data[7]
		SS_num2 = len(SS_data)
		
		SS_codes = []
		matches = 0
		Vprint(3, "Checking Class compatiblity with SS data:")
		for entry in SS_data:
			for k in range(SS_num2):
				Label = Labels[k]
				if entry[0] == Label:
					matches += 1
					SS_codes.append([entry[0],k])
					Vprint(4, "%1s Class in Basis set (%1d)" % (entry[0],k))
		if SS_num != matches:
			Vprint(1, "\nWarning, %1d of %1d matching classes were found!" %(matches,SS_num))
			Vprint(1, "Please check your target and basis set files!")
		else:
			Vprint(2, "\nAll %1d Class matched!" % SS_num)	
	
		CD_new = [[],[]]
#		Calculate Coefficients for the basis spectra:
		Vprint(2, "\nCalculating Coefficients:")
		for i in range(BS_num):
			Coeff = [Comb_Matrix[i+1][0],0.0]
			for k in range(SS_num2):
				code = SS_codes[k][1]
				Aki = Comb_Matrix[i+1][1][code]
				Wjk = SS_data[k][1]
				if Aki != 0.0:
					Coeff[1] += float(Wjk)/100 * Aki
#					print SS_data[k],Wjk, Labels[code],Aki
			CD_new[0].append(Coeff)	
			Vprint(3, "%1s : %1.3f" % tuple(Coeff))

#		Compute CD spectrum:
		CD_spect = Compute_CD(CD_new[0], Basis_Spect)
		CD_new[1] = CD_spect		
					
		return CD_new

#function to compute CD spectrum from Basis spectra:
def Compute_CD(Coeffs, Basis_Spect):	
		CD_spect = []
		Vprint(4, "\nTheoretical CD spectrum:")
		Wavelengths = Basis_Spect[0][1]
		BS_num = len(Basis_Spect)-1
		for l in range(len(Wavelengths)):
			Wave = [Wavelengths[l],0.0]
			for i in range(BS_num):
				Cji = Coeffs[i][1]
				Bil = Basis_Spect[i+1][1][l]
				Wave[1] += Cji * Bil
			CD_spect.append(Wave)	
			Vprint(4, "%6.1f   %3.4f" % tuple(Wave))		
					
		return CD_spect	

#helper function to divde coefficients (without normalization):
def Divide_Coeffs(SS_data, div= 100.0):
	New_Coeffs = []
	for Class in SS_data:
		value = float(Class[1])/div
		Temp = [Class[0],value]
		New_Coeffs.append(Temp)
	return New_Coeffs

#function to normalize weights:
def Normalize_Coeffs(SS_data):
	Old_Coeffs = [[],0.0]
#	summing coefficients
	for Class in SS_data:
		Old_Coeffs[0].append(Class[1])
		Old_Coeffs[1] += Class[1]
	Vprint(4, Old_Coeffs)
	New_Coeffs = []
	cnt = 0
#	make sure coefficients add up to 100%
	for Coeff in Old_Coeffs[0]:
		entry = ["",0.0]
		entry[0] = SS_data[cnt][0]
		entry[1] = Coeff*100/Old_Coeffs[1]
		Vprint(3, "%1s : %1.4f" % tuple(entry))
		New_Coeffs.append(entry)
		cnt += 1
#	New_Coeffs = SS_data
	return New_Coeffs

#function to scale spectra:
def Scale_Spect(Spectrum,scaling_factor):
	Scaled = []
	for Wave in Spectrum:
		Int_scaled = float(Wave[1]) * scaling_factor 
		W_scaled = [Wave[0],Int_scaled]
		Scaled.append(W_scaled)
	return Scaled

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

# function to calculate the deviation between two spectra:
def Compare_Spectra(Spect1,Spect2):  
	Vprint(1, "\nCalculating spectral Similarity:")
	ref_num = len(Spect1)
	matches = 0
	Perr = 0.0
	Pdev = 0.0
	Spect1_Amp = 0.0
	RESULTS = [[],[]]
	MISS = []
	for P1 in Spect1:
# 	find matching Point if it exists:
		X1,Y1 = P1
		found = 0
		P2 = [0.0,0.0,0]
		for Point in Spect2:
			if Point[0] == X1:
				P2 = [Point[0],Point[1],1]
				Vprint(4, "Ref.point: (%2.1f,%2.1f)" % tuple(P1), " Match: (%2.1f, %2.1f)" % (P2[0],P2[1]))
#	try intrapolation if it does not:
		if P2[2] == 0:
			P2 = Intrapolate(Spect2,X1)
			if P2[2] == 1:
				Vprint(4, " Ref. point :(%2.1f,%2.1f)" % tuple(P1), " Intrpolation:  (%2.1f, %2.f)" % (P2[0],P2[1]))
#	calculate deviation for matched points:
		if P2[2] == 1:	
			diff = P1[1]-P2[1]
			squares =  math.pow((P1[1]-P2[1]),2)
			Perr += math.fabs(diff)
			Pdev += squares
			Spect1_Amp += math.pow(P1[1],2)
			matches += 1
			Point_data = [X1,Y1,P2[1],squares,diff]
			RESULTS[0].append(Point_data)

		else:
			Vprint(2, " Ref. point skipped: (%2.1f,%2.1f)" % tuple (P1))
			MISS.append(P1)

#	Compute mean deviation and error:
	if matches != 0 and ref_num != 0:
		match_perc = float(matches)/ref_num*100
		RMSD = math.sqrt(float(Pdev)/matches)
		MUSE = float(Perr)/matches
		if Spect1_Amp != 0.0:
			NRMSD = math.sqrt(float(Pdev)/Spect1_Amp)
		else:
			NRMSD =0.0
	else:
		match_perc = 0.0
		RMSD = 0.0
		MUSE = 0.0
		NRMSD = 0.0
	RESULTS[1] = [matches,RMSD,MUSE,NRMSD,MISS]
		
	Vprint(2, "\nNumber matching entries found: %4d (%2.1f" % (matches,match_perc) + " %)")	
	Vprint(3, "RMSD: %1.3f\nMUSE: %1.3f\nNRMSD: %1.3f" % (RMSD,MUSE,NRMSD))

	return RESULTS

#function to sum or scale CD spectra:
def Modify_Spectra(Spect1,Spect2,sc_factor):
	Spect = []
	MISS = []
#	first Scale the target spectrum:
	Scaled_spectrum = Scale_Spect(Spect1,sc_factor)
#	check if Spect2 is not empty:
	if Spect2 != [] and Spect2 != "None":
#		add scaled target to the reference spectrum if both are provided:
		Vprint(2, "\nSpectra summed, scaling applied to target spectrum")
		matches = 0
		for P1 in Spect2:
# 		find matching Point if it exists:
			X1,Y1 = P1
			found = 0
			P2 = [0.0,0.0,0]
			for Point in Scaled_spectrum:
				if Point[0] == X1:
					P2 = [Point[0],Point[1],1]
					Vprint(4, "Ref.point: (%2.1f,%2.1f)" % tuple(P1), " Match: (%2.1f, %2.1f)" % (P2[0],P2[1]))
#		try intrapolation if it does not:
			if P2[2] == 0:
				P2 = Intrapolate(Scaled_spectrum,X1)
				if P2[2] == 1:
					Vprint(4, " Ref. point :(%2.1f,%2.1f)" % tuple(P1), " Intrpolation:  (%2.1f, %2.1f)" % (P2[0],P2[1]))
#		add up matched points
			if P2[2] == 1:
				Psum = [X1, (Y1+P2[1])]
				Spect.append(Psum)
			else:
				Vprint(2, " Ref. point skipped: (%2.1f,%2.1f)" % tuple (P1))
				MISS.append(P1)
			
		Mod_Spectrum =  [Spect, MISS]
	else:
#		return scaled spectrum if no reference was provided:
		Vprint(1, "\nTarget Spectrum scaled (%1.3f)" % sc_factor)
		Mod_Spectrum = [Scaled_spectrum, MISS]
	return Mod_Spectrum	

#function to estimate the SS error (based on the RMSD):
def Estimate_SSerror_linear(Rmsd,Error_Par):
	mf = Error_Par[0][1]
	Sm = Error_Par[0][2]
	Rf = Error_Par[1][1]
	Sr = Error_Par[1][2]
	SSmin = Error_Par[2][1]
	check = -1
 
#	calculate estimated SS error and typical uncertainty with weak error correlations
	SS_est  = float(Rmsd)/mf*100
	typ_err = float(SSmin)/math.sqrt(2)*100
	if SS_est < 0.0:
		check = -2
		SS_est = 0.0
		Vprint(2, "Warning the estimated model error is smaler than 0%, please double check the error parameters!")
	elif SS_est > 100.0:
		check = -3
		SS_est = 100.0
		Vprint(2, "Warning the estimated model error is larger than 100%, like due errors in the reference data!")

#	calculate upper and lower bounds of the SS_error assuming strong error correlations
	SS_up   = min(100.0, float(Rmsd+Rf+Sr)/(mf-Sm)*100)
	SS_low  = max(0.0, float(Rmsd-Rf-Sr)/(mf+Sm)*100)
	Data = [SS_est, typ_err, SS_up, SS_low, SSmin, check]
	return Data

#funtion to perform non-linear SS error estimation:
def Estimate_SSerror_nonlin(Rmsd, Error_Par, curve):	
#	define baisc variables:
	Data =[]
	R_min = Error_Par[0][1]
	R_max = Error_Par[0][2]	
	e_num = len(Error_Par)
	col_num = len(Error_Par[1])
	dSS_curve = []
	SD_curve = []

#	extract calibration curves:
	Vprint(5, "\nCalibration curve:")
	k1, k2 = 2*(curve-1)+1, 2*(curve-1)+2
	for j in range(1,e_num):
		Rmsdj, dSSj, SDj = Error_Par[j][0], Error_Par[j][k1], Error_Par[j][k2]
		Vprint(5,"%1.3f %1.3f  %1.3f" % (Rmsdj,dSSj, SDj))
		dSS_curve.append([Rmsdj, dSSj])
		SD_curve.append([Rmsdj, SDj])

#	estimate model error for a given Rmsd:
	dSS_fin = [Rmsd,0.0]
	SD_fin = [Rmsd,0.0]
	check = 0
#	Vprint(2, Rmsd, R_max, R_min, check)
	if Rmsd > R_min and Rmsd < R_max:
		dSS_fin = Intrapolate(dSS_curve, Rmsd)
		SD_fin = Intrapolate(SD_curve, Rmsd)
		check = 1
	elif Rmsd <= R_min:
		Vprint(2, "Warning, the spectral deviation (RMSD: %1.3f) is smaller than the calibration range! (Rmin: %1.3f)"%(Rmsd,R_min))
		for j in range(e_num-1):
			Rmsdj,dSSj,SDj = dSS_curve[j][0],dSS_curve[j][1], SD_curve[j][1]
			if R_min == Rmsdj:
				dSS_fin = [R_min, dSSj]
				SD_fin = [R_min, SDj]
				check = 2
	elif Rmsd >= R_max:
		Vprint(2, "Warning, the spectral deviation is (RMSD: %1.3f) is larger than the calibration range! (Rmax: %1.3f)"%(Rmsd,R_max))
		for j in range(e_num-1):
			Rmsdj,dSSj,SDj = dSS_curve[j][0],dSS_curve[j][1], SD_curve[j][1]
			if R_max == Rmsdj:
				dSS_fin = [R_max, dSSj]
				SD_fin = [R_max, SDj]
				check = 3


	Vprint(4, "Estimated model error: %1.3f +/- %1.3f"%(dSS_fin[1], SD_fin[1]))
	SS_est = float(dSS_fin[1])*100
	typ_err = float(SD_fin[1])*100
	SS_low = max(0.0, SS_est-2*typ_err)
	SS_up = min(100.0, SS_est+2*typ_err) 
	Data = [SS_est, typ_err, SS_up, SS_low, curve, check]

	return Data

	
#function to format output for writing
def Format_Output(Filenames,Main_data,Aux_data,mode):
	workdir, reference, target, library = Filenames
	matches, Rmsd, Muse, Nrmsd, W_extra, Error, scaling = Aux_data
	Output = ""
#	format header file info:	
	Header = "#SESCA CD calculation module:\n"
	Header += "#Workdir: %1s\n" % workdir
	if reference != "":
		Header += "#Reference file: %1s\n" % reference
	if target != "":
		Header += "#Target file: %1s\n" % target
	if library != "" and mode in [0,1,4]:
		Header += "#Basis set file: %1s\n" % library
	if mode != 2 and scaling != 1.0:
		Header += "#Scaling factor applied to target: %1.3f\n" % scaling
	Output += Header+"\n"

#	format Coefficients:
	if mode in [0,1,4]:
		W_info = "#Weighting factors for the calculated CD spectrum:\n"
		for entry in Main_data[0]:
			W_info += "#  %12s  :   %1.3f \n" % tuple(entry)
		Output += W_info+"\n"
	
#	format CD spectra:
	Spectra = ""
	if mode in [0,3,4]:
		Spectra += "#   wvlgth       Icalc\n"
	elif mode == 1:
		Spectra += "#   wvlgth       Iref        Icalc       dev         diff\n"
	elif mode == 2:
		Spectra += "#   wvlgth       Int1        Int2        dev         diff\n"
 
	for entry in Main_data[1]:
		col_num = len(entry)
		string = col_num*"   %6.3f   "+"\n"
#		print entry
		Spectra += string % tuple(entry)
	Output += Spectra

#	format deviation data:
	Devs = ""
	if Rmsd != "":
		Devs += "\n#matches"+9*" "+ "RMSD"+8*" " +"MUSE"+8*" "+"NRMSD\n#"
		Devs += 2*" " + "%5d"%matches + 8*" " + "%6.3f"%Rmsd + 6*" " + "%6.3f"%Muse+6*" "+"%6.3f\n"%Nrmsd
	Output += Devs

#	format error estimates:
	SS_est = ""
	if Error != []:
		print(Error)
		if Error[5] < 0:
			SS_est += "\n#Estimated SS error of the protein model (linear estimate):\n"
			SS_est += "#dSS-est (assuming weak error correlations):  %2.2f +/- %2.2f" % (Error[0],Error[1]) + " %\n"
			SS_est += "#dSS-up & dSS-low (assuming strong corr.):    %2.2f  -  %2.2f" % (Error[2],Error[3]) + " %\n"
		else:
			SS_est += "\n#Estimated SS error of the protein model (calibration curve %1d):\n"%Error[4]
			SS_est += "#dSS-est (and standard deviation):  %2.2f +/- %2.2f" % (Error[0],Error[1]) + " %\n"
			SS_est += "#dSS-up & dSS-low (conf. interval): %2.2f  -  %2.2f" % (Error[2],Error[3]) + " %\n"
		if Error[5] in [2,-2]:
			SS_est += "#Warning, spectrum deviation is smaller than the expected range!"
		elif Error[5] in [3,-3]:
			SS_est += "#Warning, spectrum deviation is larger than the expected range!\n"
			SS_est += "#Error estimate is likely unreliable!\n" 
			
	Output += SS_est

#	format missing wavelengths:
	Missed_wave = ""
	if W_extra != []:	
		W_num = len(W_extra)
		Missed_wave += "\n#Experimental entries with no match found: %1d\n" % W_num
		for entry in W_extra:
			Missed_wave += "# %2.3f     %2.3f\n" % tuple(entry)
	Output += Missed_wave	
			
	return Output


# Main function for script execution:
def CDpred_Main(Args):
#	set run parameters:
	Inp_files, Out_files, Param, failmark, verbosity = Args
	ss_file, lib_file, ref_file = Inp_files
	out_file = Out_files[0]
	mode, L_range, sc_factor, norm, nrmsd, error = Param
	Set_verb(Args[4])

	Vprint(1, "\nSESCA CD_calc module")
	Main_Data = []	
	Aux_Data = ["","","","",[],[],sc_factor]
	Ref1 = []
	Ref2 = []
	Lib1 = []
	Input_Data = [Ref1, Ref2, Lib1]

#	check if somebody forgot the @tar flag for modes 0 or 3:
	if mode in [0,3] and ref_file != "" and ss_file == "":
		ss_file = ref_file
		ref_file = ""

#	read in structure/secondary spectrum file:
	if mode in [0,1,4] and failmark == 0:
		Ref2 = Read_Struct_file(ss_file)
		if mode == 4:
			Ref2 = Divide_Coeffs(Ref2)
	elif failmark == 0:
		Ref2 = Read_Spectrum_file(ss_file,L_range)
#	stop if no target is found:
	if Ref2 == "None" and failmark == 0:
		failmark =4
	else:
		Input_Data[1] = Ref2
		
#	read in basis set parameters:
	if failmark == 0:
		Basis = Read_BS(lib_file,L_range)
#		stop if basis set is not available when needed
		if mode in [0,1,4]:
			if Basis == "None":
				failmark = 5
			Input_Data[2] = Basis
#			print Basis

#	check the requested error estimation:
	if failmark == 0 and mode ==  1 and error != 0:
		if Basis[3] == []:
			Vprint(2, "\nNo error parameters were found!")
			error = 0
		elif len(Basis[3][0]) == 3:
			Vprint(2, "\nUsing linear error estimation")
			error = -1
		elif len(Basis[3][0]) == 4:
			curve_num = Basis[3][0][3]
			if error > 0 and error <= curve_num:
				Vprint(2, "\nUsing calibration curve %1d for estimating model error"%error)
			else:
				Vprint(1, "\nError, invalid error calibration (%1d) request, (%1d available)"%(error, curve_num))
				failmark = 1
					
	
#	read in reference spectrum file:
	if mode != 0 and failmark == 0:
		Ref1 = Read_Spectrum_file(ref_file,L_range)
		if mode in [3,4] and Ref1 == "None":
			pass
		elif Ref1 == "None" and failmark == 0:
			failmark = 3	
		else:
			Input_Data[0] = Ref1

#       Parse Input information:
	if failmark != 0:
		print("\nError while reading input files, script stops!")
		print(Usage)
		sys.exit(failmark)

#	Calculate theoretical spectra:
	if mode in [0,1,4]:
		if norm == 1:
			Vprint(1, "\nNormalizing coefficients:")
			Coeffs = Normalize_Coeffs(Ref2)
		else:
			Coeffs = Ref2			
		Vprint(1, "\nComputing theoritical spectrum:")
		if mode in [0,1]:
			CD_calc = Generate_CD(Coeffs,Basis)
		elif mode in [4]:
			Vprint(1,Coeffs)
			CD_calc = Compute_CD(Coeffs,Basis[2])
		
		if sc_factor != 1.0:
			Vprint(1, "\nScaling calculated spectrum by %1.3f" % sc_factor)
			Spectrum_new = Scale_Spect(CD_calc[1],sc_factor)
			CD_calc[1] = Spectrum_new

#	Organize results according to mode
	if mode == 0:
#		just get the calculated spectrum:
		Main_Data = CD_calc

	elif mode == 1:
#		compare reference and calc. spectrum
		Comp_summary = Compare_Spectra(Ref1,CD_calc[1])
		Main_Data = [CD_calc[0], Comp_summary[0]]
		Aux_Data = Comp_summary[1]
		RMSD = Aux_Data[1]
		if Basis[3] != [] and RMSD != "" and error != 0:
#		compute estimated error
			if error == -1:
				Error_data = Estimate_SSerror_linear(RMSD,Basis[3])
			else:
				Error_data = Estimate_SSerror_nonlin(RMSD,Basis[3],error)
#			print "Error data:",Error_data
			Aux_Data.append(Error_data) 
		else:
			Aux_Data.append([])
		Aux_Data.append(sc_factor)

	elif mode == 2:
#		compare two CD spectra:
		Comp_summary = Compare_Spectra(Ref1,Ref2)
		Main_Data = [[],Comp_summary[0]]
		Aux_Data =  Comp_summary[1]
		Aux_Data.append([])
		Aux_Data.append(sc_factor)
	elif mode == 3:
#		sum and/or scale provided spectra:
		Summed_spectra = Modify_Spectra(Ref2,Ref1,sc_factor)
		Main_Data = [[],Summed_spectra[0]]
		Aux_Data[4] = Summed_spectra[1]
#		Vprint(4, "\nRef1:\n",Ref1,"\nRef2:\n",Ref2,"\nscaling:",sc_factor)		
	elif mode == 4:
#		recalculate CD spectrum from weights:
		Vprint(2,"\nCoefficients:",Coeffs)
		Main_Data = [Coeffs, CD_calc]	
			      	
	else:
		Main_Data = Input_Data

#	writing output:
	Filenames = [workdir,ref_file,ss_file,lib_file]
	if out_file != "":
		Vprint(4,Aux_Data)
		Output_data = Format_Output(Filenames, Main_Data, Aux_Data, mode)
		Output = Output_data.encode("ascii")
		o = open(out_file,"wb")
		o.write(Output)
		o.close()
			
	return Main_Data	

# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
#	set arguments:
	Inp_files, Out_files, Param, failmark, verbosity = Custom_Args
	out_file = Out_files[0]
	Set_verb(Custom_Args[4])
	Vprint(2, "\nRun parameters:\n", Custom_Args)

#       executing main code
	Data_main = CDpred_Main(Custom_Args)
#       print "\nMain data:\nRef1:\n",Data_main[0],"\nRef2:\n",Data_main[1]
#	print "\nLib1:\n",Data_main[2]

#	print run-time message
	ftime = time.time()
	runtime = ftime-stime
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",out_file)
else:
	Vprint(2, "SESCA CD prediction module (SESCA_pred.py)")	
