#!/usr/bin/env python
#


import sys
import os
import math
import time


stime = time.time()
workdir = os.getcwd()

#Python script designed to scale the intensity of a target spectrum, to match that of reference spectrum
usage0  = "*******************************************************************\n"
usage  = "SESCA module for scaling CD spectrum intensities\nusage: Scale_spect.py <target_file> (<reference_file>) @flag <argument>\nPossible command flags are:\n"
usage += "   @ref  <ref_file> specify reference file (CD spectrum to match, optional)\n   @tar  <target_file> specify target file (CD spectrum to be scaled)\n"
usage += "   @write <0 / out_file> specify output file (0 - dont print output file, default: Spectrum_scaled.out\n"
usage += "   @SF0 <float> use starting scaling factor for the target CD-spectrum, if iter >0 is set, the scaling factor will be optimized (default: 1.0)\n"
usage += "   @range <float,float> limit minimization range to work in (default: from 0.2 to 5.0)\n"
usage += "   @iter <int> set the maximum number of iterations durun the SF-search (default: 5000)\n"
usage += "   @mult <int> set Lagrange multiplier for forcing the SF-range (default: 100)\n"
usage += "   @norm <1,0> set normalization for for the mean CD intensity (default: 0 - off )\n"


Usage = usage0 + usage + usage0 

#Specify SESCA files and directories here:
#################################################
SESCA_Dir =     "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097"
SESCA_scripts =  os.path.join(SESCA_Dir,"scripts")
##################################################
#the code from here should not be changed


#default parameters:
ref_file = ""
tar_file = ""
out_file = "Spectra_scaled.out"
SF_range = [0.2,5]
SF_start = 1.0
norm = 0
Lambda = 100
Maxiter = 5000
failmark = 0
verbosity = 1

Input_files = [ref_file, tar_file]
Output_files = [out_file]
Scale_param =  [SF_start, SF_range, norm, Lambda, Maxiter]

Def_Args = [Input_files, Output_files, Scale_param, failmark, verbosity]

#define necessary functions:
def Pass_defaults():
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

#read arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(3, "Reading in %1d arguments:\n" % argnum, Args)
#	getting default arguments:
	New_Args = [Input_files, Output_files, Scale_param, failmark, verbosity]
#	Input_files = [ref_file, tar_file]
#	Output_files = [out_file]
#	Scale_Param =  [SF_start, SF_range, norm, Lambda, Maxiter]
	

	FLAGS = ["ref","tar","write","SF0","range","mult","iter","verb","norm"]
	flag = ""
	acnt = 0
	for arg in Args:
		if arg.startswith("@"):
			flag = arg.strip("@")
			if not flag in FLAGS:
				Vprint(1, "Unknown Flag:",flag)
				New_Args[3] = 1
		elif flag == "ref":
			New_Args[0][0] = arg
			flag = ""
		elif flag == "tar":
			New_Args[0][1] = arg
			flag = ""
		elif flag == "write":
			if arg == "0":
				New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "SF0":
			try:
				New_Args[2][0] = float(arg)
			except Exception:
				Vprint(1, "@SF0 only takes float arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "range":
			try:
				parts = arg.split(",")
				new_range = []
				for number in parts:
					value = float(number)
					new_range.append(value)
				New_Args[2][1] = new_range
			except Exception:
				Vprint(1, "@range only takes two comma-separated float arguments")
				New_Args[3] = 1
			flag = ""	
		elif flag == "norm":
			if arg in ["1","0"]:
				New_Args[2][2] = int(arg)
			else:
				Vprint(1, "@norm only takes 1 and 0 as arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "mult":
			try:
				New_Args[2][3] = float(arg)
			except Exception:
				Vprint(1, "@mult only takes float arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "iter":
			try:
				New_Args[2][4] = int(arg)
			except Exception:
				Vprint(1, "@iter only takes integer arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "verb":
			try:
				New_Args[4] = int(arg)
			except Exception:
				Vprint(1, "@verb only takes integer arguments")
				New_Args[3] = 1
			flag = ""
#setting defaults if no command flag is used:
		elif flag == "" and New_Args[0][1] == "" and acnt == 0:
			New_Args[0][1] = arg
			flag = ""
		elif flag == "" and New_Args[0][0] == "" and acnt == 1:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1][0] in ["","Spectra_scaled.out"] and acnt == 2:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1, "Unknown argument:", arg)
			New_Args[3] = 1

		acnt += 1

	return New_Args

#function to load non-standard modules:
def Import_Custom(*args):
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
			Vprint(2,"\nWarning, Could not import spimplex module!")

	return IMPORTS



#function to read in spectra:
def Read_Spectrum(Infile):
        Vprint(1, "\nReading in Spectrum file:",)
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
#               getting spectrum information
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


# fuction definition for chi^2 function (1/N*Sum[(Iexp(i)-C*Icalc(i))/sigma_exp(i)]^2:
def Func_chi_sq(Args):
        if len(Args) == 1:
                global RESULTS
                Data = RESULTS
                Const = Args[0]
        elif len(Args) == 2:
                Data = Args[0]
                Const = Args[1]
        else:
                Vprint(1, "Argument error, this functions takes 1 or 2 variables")
                return "None"
        if Data != []:
                cnt = 0
                Sum = 0
                for entry in Data:
                        if entry[3] == 0.0:
                                sigma = 0.001
                        else:
                                sigma = entry[3]
                        Sum += ((entry[1]-Const*entry[2])/sigma)**2
                        cnt += 1
                Chi_sq = Sum/cnt
        else:
                Vprint(1, "Chi_sq: Error! Empty data array!")
                Chi_sq = 0.0
        return Chi_sq


#function to intrapolate spectrum intensity at wavelength X
def Intrapolate(Spectrum, X):
			lower = ["",0.0,0.0]
			upper = ["",0.0,0.0]
			entry = [X,0.0,0]
			for point in Spectrum:
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

#function to parse spectra:
def Parse_Spectra(Spect1,Spect2):  
	Vprint(1, "\nParsing Spectra:")
	ref_num = len(Spect1)
	matches = 0
	Perr = 0.0
	Pdev = 0.0
	RESULTS = [[],[],[],[]]
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
			matches += 1
			RESULTS[0].append(X1)
			RESULTS[1].append(Y1)
			RESULTS[2].append(P2[1])
			
		else:
			Vprint(2, " Ref. point skipped: (%2.1f,%2.1f)" % tuple (P1))
			MISS.append(P1)

#	Compute mean deviation and error:
	if matches != 0 and ref_num != 0:
		match_perc = float(matches)/ref_num*100
		RMSD = math.sqrt(float(Pdev)/matches)
		MUSE = float(Perr)/matches
	else:
		match_perc = 0.0
		RMSD = 0.0
		MUSE = 0.0
	RESULTS[3] = [matches,RMSD,MUSE,MISS]		
	Vprint(3, "\nNumber matching entries found: %4d (%2.1f" % (matches,match_perc) + " %)")	
	Vprint(3, "RMSD: %1.3f\nMUSE: %1.3f" % (RMSD,MUSE))
	return RESULTS

#function to calculate root mean squared deviation between the spectra:
def Calc_RMSD(Const,Ints):
	Int1, Int2 = Ints
	point_num = len(Int1)
	Sum = 0.0
	for cnt in range(point_num):
		Sum += (Int1[cnt]- Const*Int2[cnt])**2
	if point_num != 0:
		Rmsd = math.sqrt(float(Sum)/point_num)
	else:
		Rmsd = 0.0

	return Rmsd

#function to calculate the amplitude normalized RMSD between the spectra:
def Calc_NRMSD(Const,Ints):
	Int1, Int2 = Ints
	point_num = len(Int1)
	Sum = 0.0
	Norm = 0.0
	for cnt in range(point_num):
		Sum += (Int1[cnt]- Const*Int2[cnt])**2
		Norm += (Const*Int2[cnt])**2
#		Norm += math.fabs(Const*Int2[cnt])
	if point_num != 0:
		Rmsd = math.sqrt(float(Sum)/point_num)
#		NormF = math.sqrt(float(Norm)/point_num)
		NormF =  float(Norm)/point_num
#		NormF = Const*math.fabs(max(Int2)-min(Int2))
		if NormF != 0.0:
			NRMSD = Rmsd/NormF
		else:
			NRMSD = Rmsd
	else:
		NMRSD = 0.0

	return NRMSD

#function to scale spectra:
def Scale_Spect(Spectrum,scaling_factor):
        Scaled = []
        for Wave in Spectrum:
                Int_scaled = float(Wave[1]) * scaling_factor
                W_scaled = [Wave[0],Int_scaled]
                Scaled.append(W_scaled)
        return Scaled

#modified RMSD function to minimize scaling factors:
def Eval_Function(Consts,Args):
	Const = Consts[0]
	Int1,Int2,SF_range,Lambda, norm = Args
	Full_function = 1000.0
#	Adding restraints on the scaling factor
	LG_Mult = 0.0
#	print SF_range,Lambda,Const
	if SF_range[0] != "" and Const < SF_range[0]:
		LG_Mult += Lambda * (SF_range[0] - Const)
	if SF_range[1] != "" and Const > SF_range[1]:
		LG_Mult += Lambda * (Const - SF_range[1])
		 
	if norm == 0:
		Full_function = Calc_RMSD(Const,[Int1,Int2]) + LG_Mult
	else:
		Full_function = Calc_NRMSD(Const,[Int1,Int2]) + LG_Mult
#	print (Full_function-LG_Mult),LG_Mult

	return Full_function

#function to optimize scaling:
def Find_Scaling(SF0, Ints, SF_range, Lambda, Maxiter,norm,IMPORTS):	
	Eval_args = [Ints[1],Ints[2],SF_range,Lambda,norm]
#	use a simplex minimizer to find SF-s with the lowest RMSD
	Vprint(2,IMPORTS)
	if Maxiter != 0 and "minimize" in IMPORTS:
		Vprint(1, "Simplex minimization with Scipy.optimize module:")
		adaptive,tol,disp = [False, 1e-15, True]
		if verbosity < 3:
			disp = False
		Min_results = spo.minimize(Eval_Function, [SF0], Eval_args, method='nelder-mead',
		options={"ftol":tol,  "disp":disp, "maxiter":Maxiter, "adaptive":adaptive})
		SF_fin = Min_results.x[0]
	elif Maxiter != 0 and "minimize2" in IMPORTS:
		Vprint(1, "Simplex minimization with SESCA_min module:")
		adaptive,tol,disp = [0, 1e-15, 1]
		if verbosity < 3:
			disp = 0
			min2.Set_verb(0)
		Min_param = [adaptive,tol, Maxiter,disp]
		NM_param = min2.Pass_NMparam() 
		Vprint(4,Min_param)
		Min_results = min2.Iterate_Min(Eval_Function, [SF0], Eval_args, Min_param, NM_param)
		if Min_results != None:
			SF_fin = Min_results[2][0]
		else:
			SF_fin = SF0
			Vprint(1, "Reference Scaling factor set to: %1.1f"%SF0)

	else:
		SF_fin = SF0	

#	recalculate RMSD
	RMSD_fin = Calc_RMSD(SF_fin,[Ints[1],Ints[2]])
	Vprint(1, "optimum scaling factor: %1.3f with RMSD: %1.3f" %(SF_fin, RMSD_fin))
#	rescale target spectrum
	Tar_Spect = []
	for l in range(len(Ints[0])):
		Wave = [ Ints[0][l], Ints[2][l] ]
		Tar_Spect.append(Wave)
	Scaled_Target = Scale_Spect(Tar_Spect,SF_fin)

	Scaled_results = [SF_fin, RMSD_fin, Scaled_Target]
	return Scaled_results

#function to generated formatted output:	
def Format_Output(Filenames,Data):
		Output = ""
		Header  = "#SESCA spectrum Scaling module:\n"
		Header += "#Workdir: %1s\n#Reference file: %1s\n#Target file: %1s\n" % tuple(Filenames)
		Header += "#Scaling factor applied to target: %1.3f\n" % Data[0]
		Output += Header+"\n"

		Spectrum = "#   wvlght       Icalc\n"
		string = 2*"   %6.3f   "+"\n"      
		for entry in Data[2]:
			Spectrum += string % tuple(entry)

		Output += Spectrum
		return Output

#Main function for script execution:
def Scale_Main(Args):
	Input_files, Output_files, Scale_param, failmark, verbosity = Args
	ref_file, tar_file = Input_files
	out_file = Output_files[0]
	SF_start, SF_range, norm, Lambda, Maxiter = Scale_param
	Set_verb(Args[4])

	IMPORTS = Import_Custom()
#	read in spectra:
	Tar_Spect = []	
	if tar_file != "" and failmark == 0:
		Tar_Spect = Read_Spectrum(tar_file)
		if Tar_Spect == "None":
			failmark = 2
	elif failmark == 0:
		Vprint(1, "No target file was provided!")
		failmark = 2

	Ref_Spect = []
	if ref_file != "" and failmark == 0:
		Ref_Spect = Read_Spectrum(ref_file)
		if Ref_Spect == "None":
			failmark = 3

	if not "minimize" in IMPORTS and not "minimize2" in IMPORTS and Maxiter != 0:		
		print("\nError, missing dependencies! Scipy.optimize or SESCA_min are required to optimize the scaling factors, script stops!")
		failmark = 6

#	check if the spectra were read correctly:
	if failmark != 0:
		print("\nError while reading input files, script stops!")
		print(Usage)
		sys.exit(failmark)


#	set initial scaling factor
	if SF_start != "":
		SF0 = SF_start
	else:
		SF0 = 1.0

#	handle simple scaling problems:
	if Ref_Spect == [] and Tar_Spect != []:
		Scaled_Spect = Scale_Spect(Tar_Spect,SF0)
		Scale_Data = [SF0, 0.000, Scaled_Spect]
#	fit scaling to reference:
	else:
#		parse spectra:
		Parsed_data = Parse_Spectra(Ref_Spect,Tar_Spect)
		Ints = [Parsed_data[1],Parsed_data[2]]
		RMSD0 = Calc_RMSD(1.0,Ints)
#		print Parsed_data

#		Set minimzation parameters:
		Eval_args = [Ints[0],Ints[1],SF_range,Lambda,norm]
		test_Eval = Eval_Function([SF0],Eval_args)
		Vprint(2, RMSD0, test_Eval,SF0)

#		perform RMSD minimization to find the optimal scaling factor:
		Scale_Data = Find_Scaling(SF0,Parsed_data,SF_range,Lambda,Maxiter,norm,IMPORTS)
	
	#write out result if necessary:
	Filenames = [workdir,ref_file,tar_file]	
	if out_file != "":
		Output_data = Format_Output(Filenames,Scale_Data)
		Output = Output_data.encode("ascii")
		o = open(out_file,"wb")
		o.write(Output)
		o.close()

	return Scale_Data
	


#executing standalone script:
if __name__ == "__main__":
#	handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Input_files, Output_files, Scale_param, failmark, verbosity = Custom_Args
#	ref_file, tar_file = Input_files
#	out_file = Output_files
#	SF_start, SF_range, norm, Lambda, Maxiter = Scale_param
	
	Set_verb(Custom_Args[4])


	Vprint(1, "\nScaling module\nRun parameters:\n", Custom_Args)
#	executing main code:
	Data_main = Scale_Main(Custom_Args)
#	print "\nMain_data:\n",Data_main

#       print run-time message
	ftime = time.time()
	runtime = ftime-stime
	outfile = Custom_Args[1][0]
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",outfile)
else:
#	print(IMPORTS)
#	if module is called but dependencies are missing, raise an import error
	IMPORTS = Import_Custom()
	if  ("numpy" in IMPORTS and "minimize" in IMPORTS) or "minimize2" in IMPORTS:
		Vprint(1,"SESCA Spectrum Scaling module (SESCA_scale.py)")
	else:	
		raise ImportError("Failed to load dependencies")
