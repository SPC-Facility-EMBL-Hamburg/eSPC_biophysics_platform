#!/usr/bin/env python
#

import sys
import os
import math
import time

workdir = os.getcwd()
stime = time.time()


#print workdir+"\n"
# ["spect","struct","map","mode","range","iter","sigma","write","mat","verb"]
usage0  = "*******************************************************************\n"
usage   =  "SESCA basis set solver module for deriving new basis sets.\n"
usage   += "usage: SESCA_solver.py <spectrum_list> <SS list> <map_file> @flag <argument>\nPossible command flags are:\n"
usage   += "   @spect <file_list> specify file list for reference CD spectra\n"
usage   += "   @struct <file_list> specify file list for reference SS compositions\n"
usage   += "   @map  <file>  specify file containing combination matrix / assignment\n"
usage   += "   @mode <int> specify calculation method, (default: 1)\n      0 - Simplex mode, linear fit\n"
usage   += "      1 - Linear solver mode (linear fit, requires numpy)\n      2 - Simplex mode, Rmsd-based fit\n" 
usage   += "   @range <float,float> set minimum and maximum wavelenghts for the basis spectra (also limited by the spectra)\n"
usage   += "   @iter <int> set the maximum number of iterations (for simplex mode, default: 2000)\n"
usage   += "   @write <0,file> specify Basis set file name, 0 - do not write basis set (default: Basis_spectra.dat\n"
usage   += "   @mat <0,file> specify filename writing data matrices, 0 - do not write matrices (default: Matrices.dat)\n"
usage   += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"
Usage = usage0 + usage + usage0

#Specify SESCA files and directories here:
#################################################
SESCA_Dir =     "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097"
SESCA_scripts =  os.path.join(SESCA_Dir,"scripts")
##################################################
#the code from here should not be changed

#default parameters:
spect_file = ""
Spect_Mat = []
struct_file = ""
Struct_Mat = []
map_file = ""
Comb_Mat = []
out_file = "Basis_spectra.dat"
mat_file = "Matrices.dat"
mode = 1
Maxiter = 2000
L_range = ["",""]
B_param = [0.20, 5.0]
IMPORTS = []

failmark = 0
verbosity = 1

Input_files = [spect_file, struct_file, map_file, Spect_Mat, Struct_Mat, Comb_Mat]
Output_files = [out_file, mat_file]
Param = [mode, L_range, Maxiter, B_param]

Def_Args = [Input_files, Output_files, Param, failmark, verbosity]


#function definitions:
#function to pass on defaults:
def Pass_Defaults():
        return Def_Args

def Pass_imports():
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

#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
#	New_Args = [Input_files, Output_files, Param, failmark, verbosity]
	New_Args = Pass_Defaults()
	FLAGS = ["spect","struct","map","mode","range","iter","sigma","write","mat","verb"]
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
		elif flag == "spect":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "struct":
			New_Args[0][1] = arg
			flag = ""	
		elif flag == "map":
			New_Args[0][2] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "mat":
			if arg == "0":
                                New_Args[1][1] = ""
			else:
				New_Args[1][1] = arg
			flag = ""
		elif flag == "mode":
			New_Args[2][0] = int(arg)
			flag = ""
		elif flag == "range":
			try:
				parts = arg.split(",")
				new_range = []
				for number in parts:
					value = float(number)
					new_range.append(value)
				if new_range[0] <= new_range[1]:
					New_Args[2][1] = new_range
				else:
					New_Args[2][1] = [new_range[1],new_range[0]]
			except Exception:
				Vprint(1, "@range only takes two comma-separated float arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "iter":
			try:
				New_Args[2][2] = int(arg)
			except Exception:
				Vprint(1, "@iter only takes integer arguments")
				New_Args[3] = 1
			flag = ""
	
		elif flag == "sigma":
			try:
				parts = arg.split(",")
				new_range = []
				for number in parts:
					value = float(number)
					new_range.append(value)
					New_Args[2][3] = new_range
			except Exception:
				Vprint(1, "@sigma only takes two comma-separated float arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "verb":
			New_Args[4] = int(arg)
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0][0] == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[0][1] == "" and acnt == 1:
			New_Args[0][1] = arg
			flag = ""
		elif flag == "" and New_Args[0][2] == "" and acnt == 2:
			New_Args[0][2] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args


#function to load non-standard modules:
def Import_Custom():
	IMPORTS = []
	try:
		globals()["np"] = __import__("numpy")
		IMPORTS.append("numpy")
	except ImportError:
		Vprint(2,"\nWarning, Could not import numpy module \nThis module is dependent on that packages!")
	try:	
		globals()["spo"] = __import__("scipy.optimize", globals(),locals(),["minimize"])
		IMPORTS.append("minimize")
	except ImportError:
		Vprint(2,"\nminimizer from scipy.optimize not found, importing SESCA simplex module...")
		try:
			if not SESCA_scripts in sys.path:
				sys.path.append(SESCA_scripts)
			globals()["min2"] = __import__("SESCA_min")
			IMPORTS.append("minimize2")
		except ImportError:
			Vprint(2,"\nWarning, Could not import simplex module!")

	return IMPORTS



#function to read file lists:
def Read_List(infile):
	if infile != "" and os.path.isfile(infile) == False:
		return "None"
	FILES = []
	Vprint(2, "\nReading file list:")
	r = open(infile,"rb")
	for line in r:
		line2 = str(line.decode("ascii")).strip("\n")
#		determine source directory
		if line2.startswith("#Dir="):
			exp_path = line2.split()[1]
			if exp_path[-1] == "/":
				exp_path = exp_path[:-1]
			Vprint(3, "Path used:",exp_path) 
#		read in file names
		if not line2.startswith("#") or line2 == "":
			entry = line2.strip()
			if exp_path == "./":			
				full_entry = entry
			else:
				full_entry = os.path.join(exp_path,entry)
#			add full paths to file list
			Vprint(3,entry)
			FILES.append(full_entry)
	r.close()
	return FILES
			
#function to read in spectral information:
def Read_Spectrum_File(File,L_range):
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



#function to read structural information:
def Read_Struct_File(File):
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
	i.close()
	return Classes


#function to read combination matrices:
def Read_Comb_Mat(File):
	flag = 0
	Classes = []
	Codes = []
	MATRIX = []
	Vprint(1,"\nReading Combination File:",)
	if os.path.isfile(File) == False:
		Vprint(1, "\nFile not found")
		return "None"

	Vprint(1, File)
	f = open(File,"rb")
	for line in f:
		line2 = str(line.decode("ascii")).strip("\n")		
		if line2.startswith("#Basis/"):
			flag = 1
			parts = line2.split()
			Codes = parts[1:]
			Vprint(4,"Library combination map read: Codes,",Codes)
		elif line2.startswith("#") and flag == 1:
			try:
				first, second = line2.split(":")
				name = first.split()[1]
				values = second.split()
				for i  in range(len(values)):
					val = float(values[i])
					values[i] = val
				MATRIX.append(values)	
				Classes.append(name)
				Vprint(4,"Library combination map read: %1s,"%name, values)
			except Exception:
				Vprint(2,"Library, unable to process:",line2)

	f.close()
	cnum = len(MATRIX)
	knum = len(MATRIX[0])
	Vprint(3, "Combination Matrix read with  %1d X %1d dimensions"%(cnum,knum))
	COMBMAT = [MATRIX,Codes,Classes,cnum,knum]

	return COMBMAT

#Function to generate unity combination matrix:
def Generate_Unit_CM(Codes):
	knum = len(Codes)
	Classes = []
	MATRIX = []
	Vprint(2,"Assigning classes to:",Codes)
	for cnt in range(knum):
#	for each structutral element create a class
		Class_label = "BS%1s" % (cnt+1)
		Classes.append(Class_label)
		TEMP = []
#	to each class add this one structural element, set everything else to zero:
		for kcnt in range(knum):
			if cnt == kcnt:
				TEMP.append(1.0)
			else:
				TEMP.append(0.0)
		MATRIX.append(TEMP)
		Vprint(4, "Assignment for Class: %1s"%Class_label,TEMP)
	cnum = len(MATRIX)
	Vprint(3, "Combination matrix generated with  %1d X %1d dimensions"%(cnum,knum))
	COMBMAT = [MATRIX,Codes,Classes,cnum,knum]
	
	return COMBMAT
		 
	

#function to generate spectrum matrix:
def Generate_SpectMat(FILES, L_range):
	MATRIX = []
	Waves = []
	if FILES == []:
		return "None"

#	cycle through spectrum files:
	cnt = 0
	lnum = 0
	Waves0 = []
	Vprint(1, "\nGenerating spectrum matrix:")
	for File in FILES:
#		read spectrum
		Spectrum = Read_Spectrum_File(File,L_range)
		if not Spectrum in ["None",[]]:
			lcnt = 0
			for entry in Spectrum:
#				Vprint(4, entry)
				Lambda, Int = entry
#				collect intensities for compatible wavelenghts:
				if cnt == 0:
					Waves0.append([Lambda,[Int]])
					lcnt +=1
				else:
					found = 0
					for Wave in Waves0:
						if Lambda == Wave[0]:
							Wave[1].append(Int)
							found = 1	
							lcnt += 1
							break
					if found == 0:
						Vprint(2, "Warning! unexpected wavelength (%1.1f),  does not match previous spectra"%Lambda)
			cnt += 1
		else:
			Vprint(2,"Warning: Spectrum file was not read:\n%1s"%File)

#	process wavelengths:
	Vprint(2, "\nProcessing Spectra:")
	pnum = cnt
	for j in range(pnum):
		MATRIX.append([])
	lcnt = 0
	for Wave in Waves0:
		inum  = len(Wave[1])
		if inum == pnum:
			Vprint(3, "Wavelength confirmed: %2.1f nm"%Wave[0])
			Waves.append(Wave[0])
			for j in range(pnum):
				MATRIX[j].append(Wave[1][j])
			lcnt += 1
		else:
			Vprint(2, "Warning: Wavelength %2.1f nm rejected"%Wave[0])
			
	pnum = len(MATRIX)
	lnum = len(Waves)
	Vprint(3,"\nSpectrum matrix generated with %1d X %1d dimensions"%(pnum,lnum))	
			
	return [MATRIX, Waves, pnum, lnum]

#function to generate structure matrix:
def Generate_StructMat(FILES):
	MATRIX = []
	Labels = []
	knum = 0
	pnum = 0
	if FILES == []:
		return "None"


	cnt = 0
	for File in FILES:
#		read structural informaion
		Structure = Read_Struct_File(File)
		if not Structure in ["None",[]]:
			kcnt = 0
			TEMP = []
#			Vprint(3,File)
			for entry in Structure:
#				Vprint(4, entry)
				Label, Perc = entry
#				Make sure Labels match:
				if cnt == 0:
					Labels.append(Label)
					TEMP.append(Perc)
					kcnt +=1			
				elif Label == Labels[kcnt]:
					TEMP.append(Perc)
					kcnt +=1			
				else:
					Vprint(2, "Warning! unexpected SS element (%1s)  does not match previous label (%s)"%(Label,Labels[kcnt]))
#			add Structural info to the matrix:
			if cnt == 0:
				knum = len(Labels)
				MATRIX.append(TEMP)				
				Vprint(3, "Structual data with %1d fractions added to structure matrix"%knum)
			elif cnt != 0 and len(TEMP) == knum:
				MATRIX.append(TEMP)				
				Vprint(3, "Struct. data with %1d fractions added to structure matrix"%knum)
			else:
				Vprint(2,"Warning: Structure %1d with %1d fractions ignored!" % (cnt,len(TEMP)))

		else:
			Vprint(2,"Warning: Structure file was not read:\n%1s"%File)
		cnt += 1

	pnum = len(MATRIX)
	Vprint(3,"\nStructure matrix generated with %1d X %1d dimensions"%(pnum,knum))	
	Struct_Data = [MATRIX, Labels, pnum, knum]
	
	return Struct_Data 

#function to generate the coefficient matrix:
def Generate_Coeffs(Struct_data, Comb_data):
	Vprint(2, "\nCalculating Coefficients:")
	Coeff_Data = []
#	check that we have data
	if Struct_data == [] or Comb_data == []:
		Vprint(1,"Error, cannot calculate coefficients!")
		return "None"
#	check that structure elements match
	if Struct_data[3] != Comb_data[4]:
		Vprint(2, "Warning, Structural information is defined by %1d elements, and combination matrix is defined by %1d"%(Struct_data[3],Comb_data[4]))

#	map coding to Structural info:
	Codes = Comb_data[1]
	Labels = Struct_data[1]
	Mapping = []
	knum = len(Codes)
	Vprint(3, "Mapping structural elements:")
	for i in range(knum):
		match = -1
		Code = Codes[i]
		for cnt in range(len(Labels)):
			if  Code == Labels[cnt]:
				match = cnt
		if match == -1:
			Vprint(1,"Error, cannot match %1s from the combination matrix to any structural element"%Code)
			return "None"
		else:
			Mapping.append(match)
			Vprint(4,"%1s: %1d"%(Code,match))	
	
#	determine BS coefficients for each protein: 
	Vprint(2, "\nGenerating Coefficient matrix:")
	Classes = Comb_data[2]
	Comb_Mat = Comb_data[0] 
	cnum  = len(Classes)
	Struct_Mat = Struct_data[0]
	pnum = len(Struct_Mat)
	Coeff_Matrix = []
	for j in range(pnum):
		Coeffs = []
		for i in range(cnum):
			Cji = 0.0
			Wsum = 0.0
			for k in range(knum):
				m = Mapping[k]
				Wjk = Struct_Mat[j][m]
				Aik = Comb_Mat[i][k]
				Cji += (Wjk/100)*Aik
			Coeffs.append(Cji)
		Coeff_Matrix.append(Coeffs)
		Vprint(4, "Basis spectrum coeffs for entry %1d:"%(j+1),cnum*"%1.3f, " % tuple(Coeffs))
	Coeff_Data = [Coeff_Matrix, Classes, pnum, cnum] 
	Vprint(3,"\nCoefficient matrix generated with %1d X %1d dimensions"%(pnum,cnum))	
		
	return Coeff_Data

#function to solve the linear equation system with numy:
def Solve_Basis_numpy(Spect_data,Coeff_data):
	Basis_spect = []
	Spect_Mat = Spect_data[0]
	Coeff_Mat = Coeff_data[0]
	Waves = Spect_data[1]
	Classes = Coeff_data[1]
	Cji = np.array(Coeff_Mat)
	Ct = Cji.transpose()
	pnum = Spect_data[2]
	lnum = Spect_data[3]
	cnum = Coeff_data[3]
#	Basic equation at wavelength l Sjl = Sum(Cji*Bil)
#	if Sl is the vector for all Sjl ata  given l

#	we need as many basis spectra as classes:
	for i in range(cnum):
		Basis_spect.append([])

	Vprint(1, "\nSolving Basis specra with numpy:")
	Vprint(4,len(Spect_Mat),pnum)
	for l in range(lnum):
#		collect all Spectrum intensities at wavelength l
		Sl = []
		for j in range(pnum):
			Sjl = Spect_Mat[j][l]
			Sl.append(Sjl)
		Sl = np.array(Sl)

#		multiply both sides by Ct
		CtxCji = np.dot(Ct,Cji)
		CtxSl  = np.dot(Ct,Sl)
#		solve the equation system for Bl
		Bl = np.linalg.solve(CtxCji, CtxSl)
		Vprint(4,"l=%1d, %2.1f nm"%(l,Waves[l]),Bl)

#		Add Bil to the results:
		for i in range(cnum):
			Bil = Bl[i]
			Basis_spect[i].append(Bil)
		
	Basis_Mat = [Basis_spect,Classes,Waves,cnum,lnum]
	return Basis_Mat

#function to calculate  RMSD for the reference set
def Calc_Rmsd_set(Bints, Set_Param):
	Rmsd_set = 0.0
	Sints, Coeffs = Set_Param
	pnum = len(Coeffs)
	cnum = len(Coeffs[0])
	pnum2 = len(Sints)
#	make sure we have all the numbers:
	if pnum != pnum2:
		Vprint(1, "Error: dimension mismatch!\n %1d spectral intensities, but %1d set of coefficients" % (pnum2, pnum))
		return "None"
	if pnum == 0:
		Vprint(1, "Error: the number reference proteins is zero")
		return "None"

#	Calculate RMSD:
	Sum = 0.0
	for j in range(pnum):
		Sjcalc = 0.0
		Sjexp  = Sints[j]
		for i in range(cnum):
			Sjcalc += Bints[i] * Coeffs[j][i]
		Sum += ((Sjexp - Sjcalc)**2) / pnum
	Rmsd_set = math.sqrt(Sum)
	return Rmsd_set

def Calc_Diff_set(Bints, Set_Param):
	Diff_set = 0.0
	Sints, Coeffs = Set_Param
	pnum = len(Coeffs)
	cnum = len(Coeffs[0])
	pnum2 = len(Sints)
#	make sure we have all the numbers:
	if pnum != pnum2:
		Vprint(1, "Error: dimension mismatch!\n %1d spectral intensities, but %1d coefficients" % (pnum2, pnum))
		return "None"
	if pnum == 0:
		Vprint(1, "Error: the number reference proteins is zero")
		return "None"

#	Calculate absolute difference:
	for j in range(pnum):
		Sjcalc = 0.0
		Sjexp  = Sints[j]
		for i in range(cnum):
			Sjcalc += Bints[i] * Coeffs[j][i]
		Diff_set += math.fabs(Sjexp - Sjcalc) / pnum
	
	return Diff_set

def Calc_Bayes_set(Bints, Set_Param):
	Bexp = 0.0
	Sints, Coeffs, Bparam = Set_Param
	sm, s0 = Bparam
	pnum = len(Coeffs)
	cnum = len(Coeffs[0])
	pnum2 = len(Sints)
#	make sure we have all the numbers:
	if pnum != pnum2:
		Vprint(1, "Error: dimension mismatch!\n %1d spectral intensities, but %1d coefficients" % (pnum2, pnum))
		return "None"
	if pnum == 0:
		Vprint(1, "Error: the number reference proteins is zero")
		return "None"

#	Calculate the exponent of the basis set being correct:
	Sum = 0.0
	Norm = 0.0
	for j in range(pnum):
		Sjcalc = 0.0
		Sjexp  = Sints[j]
		Ujl = sm * math.fabs(Sjexp) +s0 
		for i in range(cnum):
			Sjcalc += Bints[i] * Coeffs[j][i]
		Sum += ((Sjexp - Sjcalc) / (2*Ujl))**2
		Norm += math.sqrt(2*3.14)*Ujl
	Bexp = math.sqrt(Sum/pnum)
	return Bexp
	





#	define minimization parameters:
def Eval_Function(Vars, Fparam):
	mode = Fparam[0]
	Sints = Fparam[1]
	Coeffs = Fparam[2]
	Value = 10000.0

	if mode == 0:
#		absolute diference based calculation
		Value = Calc_Diff_set(Vars, [Sints, Coeffs])	
	elif mode == 2:
#		RMSD-based calculation
		Value = Calc_Rmsd_set(Vars, [Sints, Coeffs])
	elif mode == 3:
#		probability based calculation:
		Bparam = Fparam[3]
		Value = Calc_Bayes_set(Vars, [Sints, Coeffs, Bparam])

	if Value == "None":
		Value = 10000.0
	
	return Value	

#function to solve basis spectra through RMSD fitting:
def Solve_Basis_simplex(Spect_data, Coeff_data, Param):
	Basis_Mat = []	
	Spect_Mat = Spect_data[0]
	Coeff_Mat = Coeff_data[0]
	Waves = Spect_data[1]
	Classes = Coeff_data[1]
	pnum = Spect_data[2]
	lnum = Spect_data[3]
	cnum = Coeff_data[3]
	mode = Param[0]
	IMPORTS = Param[1]
	Maxiter = Param[2]

#	Set up basis spectra
	Basis_spect = []
	Bl_init = []
	for i in range(cnum):
		Basis_spect.append([])
		Bl_init.append(10.0)

#	search basis spectra wavelength by wavelength
	for l in range(lnum):
		Wave = Waves[l]
#		determine reference intensities:
		Int_l = []
		for j in range(pnum):
			Sjl = Spect_Mat[j][l]
			Int_l.append(Sjl)
		Vprint(4,Wave, Int_l)	
#		the inital guess of the Basis set intensities:
		if l == 0:
			Bl_0 = Bl_init
		else:
			Bl_0 = Bl_prev
#		the evaluation function parameters:
		Fparam = [mode, Int_l, Coeff_Mat]
		if mode == 3:
			Fparam.append(Param[3])

#		check if the evalution function works:
		Eval0 = Eval_Function(Bl_0,Fparam)
		Vprint(3, "\nInitial function value of %1.3f at:"%Eval0, Bl_0)

#		use minimizer scpy simplex minimizer:
		if "minimize" in IMPORTS:
			Vprint(1, "Simplex minimization with Scipy.optimize module: %2.1f" % Wave)
			adaptive,tol,disp = [True, 1e-10, False]
			if verbosity >= 3:
				disp = True
			Min_results = spo.minimize(Eval_Function, Bl_0 , Fparam, method='nelder-mead',
			options={"ftol":tol,  "disp":disp, "maxiter":Maxiter, "adaptive":adaptive})
			Bl_fin = Min_results.x

#		or the built-in SESCA simplex minimizer:
		elif "minimize2" in IMPORTS:
			Vprint(1, "Simplex minimization with SESCA_min module: %2.1f" % Wave)
			adaptive,tol,disp = [1, 1e-10, 0]
			if verbosity < 2:
				min2.Set_verb(0)
			elif verbosity >= 3:
				disp = 1	
			Min_param = [adaptive,tol, Maxiter,disp]
			NM_param = min2.Pass_NMparam()
			Vprint(4,Min_param)
			Min_results = min2.Iterate_Min(Eval_Function, Bl_0, Fparam, Min_param, NM_param)
			if Min_results != None:
				Bl_fin = Min_results[2]
			else:
				Bl_fin = Bl_0
				Vprint(1, "Minimization failed, Bil set to:",Bl_0)
		else:
			Vprint(1, "Minimizers disabled, Bil set to:",Bl_0)
			Bl_fin = Bl_0


		Vprint(4,"l=%1d, %2.1f nm"%(l,Waves[l]),Bl_fin)
		Bl_prev = []
#		Add Bil to the results:
		for i in range(cnum):
			Bil = Bl_fin[i]
			Basis_spect[i].append(Bil)
			Bl_prev.append(Bil)
		
	
	Basis_Mat = [Basis_spect,Classes,Waves,cnum,lnum]
	return Basis_Mat 

#function to print basis set:
def Format_Output(Inputs, Comb_data, Basis_data):
	spect_file, struct_file, map_file, pnum = Inputs
	Basis_spectra, Classes, Waves, cnum, lnum = Basis_data
	Comb_Mat, Codes, knum = Comb_data[0],Comb_data[1],Comb_data[4]	
	Output = ""

#	Add header about input files
	header =  "#Derived Basis Spectrum set for calculating CD spectra\n"
	header += "#Spectrum files: %1s\n#Structure files: %1s\n" % (spect_file, struct_file)
	header += "#Number of decomposed  spectra: %1d\n\n" % pnum
	Output += header

#	Add combination matrix:
	Comb = "#Combination matrix for Structural Properties:\n"
#	prepare matrix header
	Chead = "#Basis/    "
	for Label in Codes:
		if len(Label) < 10:
			Chead += " "+Label+(10-len(Label))*" "
		else:
			Chead += " "+Label+" "
	Comb += Chead+"\n"
#	Add assignment coeffs for each class
	for i in range(cnum):
		String = "# %6s : " % Classes[i]
		for Aki in Comb_Mat[i]:
			String += "  %2.4f   " % Aki
		Comb += String+"\n"
	Output += Comb+"\n"  		


#	Add basis spectra:
	Basis = "\n#Detemrined Basis spectra:\n"
#	preapre header for the basis spectra
	Bhead1 = "#|  Wavelength      "
	Bhead2 = "#      nm       "
	for i in range(cnum):
		bnum = len(Classes[i])
		if bnum < 15:
			Bhead1 += Classes[i] + (15-bnum)*" "
		else:
			Bhead1 += Classes[i] + " "
		Bhead2 += "  CD(kMRE)  "+3*" "
	Basis += Bhead1+"\n"+Bhead2+"\n"
#	write basis spectra:
	for l in range(lnum):
		String = "     %3.1f     " % Waves[l]
		for i in range(cnum):
			String += "   %6.4f   " % Basis_spectra[i][l] + 2*" "
		Basis += String+"\n"
	Output += Basis

	return Output

def Format_Matrices(Spect_data, Struct_data, Coeff_data):
	Output = ""

#	write out spectrum matrix:
	Spect_Mat, Waves, pnum, lnum = Spect_data
#	first the header:
	Shead1 = "#Spectrum Matrix:\n#Wvl    "
	for j in range(pnum):
		Shead1 += 3*" "+"%1s"%(j+1)+(5-len(str(j)))*" "
	Spect = Shead1+"\n"
#	then the spectra:
	for l in range(lnum):
		String = " %2.1f  " % Waves[l]
		for j in range(pnum):
			String += " %6.3f " % Spect_Mat[j][l]
		Spect += String+"\n"
	Output += Spect+"#&\n"

#	write out Structural data:
	Struct = ""
	Struct_Mat, Classes, pnum, cnum = Struct_data
#	add header with SS elements:
	Shead2 = "#Structure Matrix:\n#   Wjk   "
	for k in range(cnum):
		length = len(Classes[k])
		if length < 13:
			Shead2 += "%1s"%Classes[k] + (13-length)*" "
		else:
			Shead2 += " %1s " % Classes[k]
	Struct += Shead2+"\n"
#	then the matrix:
	for j in range(pnum):
		String = " %6d "%(j+1)
		for k in range(cnum):
			String += " %6.1f " % Struct_Mat[j][k]+5*" "
		Struct += String+"\n"
	Output += Struct+"#&\n"

#	write out Coefficient data:
	Coeffs = ""
	Coeff_Mat, Classes, pnum, cnum = Coeff_data
	print(Coeff_data)
#	add header with SS classes:
	Chead = "#Coefficient Matrix:\n#   Cji    "
	for i in range(cnum):
		length = len(Classes[i])
		if length < 13:
			Chead += "%1s"%Classes[i] + (13-length)*" "
		else:
			Chead += " %1s " % Classes[i]
	Coeffs += Chead+"\n"
#	then the matrix:
	for j in range(pnum):
		String = " %6d "%(j+1)
		for i in range(cnum):
			String += " %6.3f " % Coeff_Mat[j][i]+5*" "
		Coeffs += String+"\n"
	Output += Coeffs+"#&\n"

	return Output		
 

# Main function for script execution:
def Script_Main(Args):
#	set run parameters:
	Input_files, Out_files, Param, failmark,verbosity = Args
	spect_file, struct_file, map_file, Spect_Data, Struct_Data, Comb_Data = Input_files
	out_file, mat_file = Out_files
	mode, L_range, Maxiter, B_param = Param
	Set_verb(Args[4])

	IMPORTS = Import_Custom()
	Vprint(2, "Imported modules:",IMPORTS)
	
#	execute main code:
	Main_Data = []

#	read in spectrum files:
	Spect_Mat = []
	Waves = []
	if Spect_Data == [] and spect_file != "" and failmark == 0:
		Exp_List = Read_List(spect_file)
		Spect_Data = Generate_SpectMat(Exp_List, L_range)

	if not Spect_Data in [[],"None"] and failmark == 0: 
		Spect_Mat = Spect_Data[0]
		Waves = Spect_Data[1]
	elif failmark == 0:
		Vprint(1, "\nError: Spectral information not available!")
		failmark = 2

#	read in structural data
	Struct_Mat = []
	Labels = []
	if Struct_Data == [] and struct_file != "" and failmark == 0:
		Struct_List = Read_List(struct_file)
		Struct_Data = Generate_StructMat(Struct_List)
	
	if not Struct_Data in [[],"None"] and failmark == 0: 
		Struct_Mat = Struct_Data[0]
		Labels = Struct_Data[1]
	elif failmark == 0:
		Vprint(1, "\nError: Structural information not available!")
		failmark = 2

#	read in combination matrix
	Comb_Mat = []
	Classes = []
	if Comb_Data == [] and map_file != "" and failmark == 0:
		Comb_Data = Read_Comb_Mat(map_file)
#	or generate a new if no data was given
	elif Comb_Data == [] and map_file == "" and Labels != [] and failmark == 0:
		Comb_Data = Generate_Unit_CM(Labels)

	if not Comb_Data in [[],"None"] and failmark == 0: 
		Comb_Mat = Comb_Data[0]
		Codes = Comb_Data[1]
		Classes = Comb_Data[2]
	elif failmark == 0:
		Vprint(1, "\nError: Could not generate combination matrix!")
		failmark = 2
						
#	calculate the Coefficient Matrix:
	Coeff_Data = []
	Coeff_Mat = []
	if Comb_Data != [] and Struct_Data != []:
		Coeff_Data = Generate_Coeffs(Struct_Data, Comb_Data)
	if failmark == 0 and not Coeff_Data in [[],"None"]:
		Coeff_Mat = Coeff_Data[0]
	elif failmark == 0:
		Vprint(1, "\nError: Cannot calculate coefficient matrix!")
		failmark = 3

#	write out matrices for inspection:
	if mat_file != "" and failmark == 0:
		Vprint(2,"Writing out generated matrices:")
		Output2 = Format_Matrices(Spect_Data, Struct_Data, Coeff_Data)
		o = open(mat_file,"wb")
		o.write(Output2.encode("ascii"))
		o.close()

#	check if required modules are imported:
	if mode == 1 and not "numpy" in IMPORTS and failmark == 0:
		Vprint(1, "\nError: numpy solver was called (@mode 1), but numpy could not be imported!") 
		failmark = 4
	if mode == 2 and not ("minimize" in IMPORTS or "minimize2" in IMPORTS)  and failmark == 0:	
		Vprint(1, "\nError: Simplex solver was called (@mode 2), but scipy could not be imported!") 
		failmark = 4

#	Stop the script if there were error while reading the input files:
	if failmark != 0:
		print("\nError while reading input files, script stops!")
		print(Usage)
		sys.exit(failmark)

#Solving basis spectum intensities:
	if mode == 1 and not Coeff_Data in [[],"None"]:
#	use numpy linalg to calculate basis spectra:
		Basis_Data = Solve_Basis_numpy(Spect_Data,Coeff_Data)

	if mode in [0, 2, 3] and not Coeff_Data in [[],"None"]:
#	use RMSD based simplex solver to calculate basis spectra:
		Param = [mode, IMPORTS, Maxiter]
		if mode == 3:
			Param.append(B_param)
		Basis_Data = Solve_Basis_simplex(Spect_Data, Coeff_Data, Param)

		
#	write out basis set:
	if out_file != "" and not Basis_Data in [[],"None" ]:
		Inputs = (spect_file, struct_file, map_file, Coeff_Data[2])
		Output = Format_Output(Inputs, Comb_Data, Basis_Data)
		o = open(out_file,"wb")
		o.write(Output.encode("ascii"))
		o.close()
	elif out_file != "":
		Vprint(1, "Basis spectra unsolved! Cannot print results.")
		sys.exit(5)


	
	Main_Data = [Comb_Data,Basis_Data,Coeff_Data]


	return Main_Data
	



# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Vprint(2, "\nRun parameters:\n", Custom_Args)
	Input_files, Output_files, Param, failmark, verbosity = Custom_Args
	out_file = Output_files[0]

#       executing main code
	Data_main = Script_Main(Custom_Args)
	Vprint(5, "\nMain data:", Data_main)

#       print run-time messages
	ftime = time.time()
	runtime = ftime-stime
	outfiles = ""
	for File in Output_files:
		if File != "":
			outfiles += " %1s," % File
	outfile = outfiles[:-1]	
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",outfiles)
else:

	IMPORTS = Import_Custom()
	if "numpy" in IMPORTS or "minimize" in IMPORTS or "minimize2" in IMPORTS:
		print("SESCA Basis set solver module (SESCA_solver.py")
	else:
		raise ImportError("Failed to load dependencies")
