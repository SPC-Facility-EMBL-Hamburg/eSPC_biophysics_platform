#!/usr/bin/env python
#

import sys
import os
import math
import time
import random

workdir = os.getcwd()
stime = time.time()


usage0  = "*******************************************************************\n"
usage   = "SESCA Bayesian secondary structure (SS) estimation module\n"
usage   += "Usage: SESCA_bayes.py <spectrum_file>  @flag <argument>\n"
usage   += "Possible command line flags are:\n"
usage   += "   @spect <file> specify spectrum file (default: None)\n"
usage   += "   @lib   <file> specify basis set file (default: DS-dT)\n"
usage   += "   @par   <auto, file> specify Bayesian parameter file (default: auto)\n"
usage   += "   @corr  <file> specify baseline/sidechain correction file (default: None)\n"
usage   += "   @write <0, file> specify main output file name, 0 - no output (default: 'SS_est.out')\n"
usage   += "   @data  <0, file> specify sampling data file name, 0 - no output (default: None)\n"
usage   += "   @proj  <0, file> specify 2D projection file name, 0 - no output (default: None)\n"
usage   += "   @pdim  <int,int> specify 2D protjection along SS classes a and b, (default 1,2)\n" 
usage   += "   @iter  <int> set the maximum number of iterations during the SS-search (default: 500)\n"
usage   += "   @size  <int> set the number inital SS compositions for the SS search (default: 100)\n"
usage   += "   @discard <float> set the fraction of SS composition discarded as burn in (default 0.0)\n"
usage   += "   @scale <0,1> allow input spectrum re-scaling, 1 - on, 0 -off (default: 1)\n"
usage   += "   @verb  <int> set verbosity level 0 - 5 (default: 1)\n"
Usage = usage0 + usage + usage0

#Specify SESCA files and directories here:
################################################# 
SESCA_Dir =     "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097"
SESCA_scripts =  os.path.join(SESCA_Dir,"scripts")
SESCA_lib =      os.path.join(SESCA_Dir,"libs")

##################################################
#the code from here should not be changed

#default parameters:
spect_file = ""
basis_file = os.path.join(SESCA_lib,"Map_DSSP-T.dat")
param_file = "auto"
corr_file = ""
out_file = "SS_est.out"
dist_file = ""
data_file = ""
discard = 0.05
max_iter = 500
sample_size = 100
scaling = "auto"
project_dim = []
failmark = 0
verbosity = 2

Input_files =  [spect_file, basis_file, param_file, corr_file]
Output_files = [out_file, dist_file, data_file]
Dev_Par =      [["Gauss", 1.3, 0.9, 0.0, 12.0, 0.3], []]
Scale_Par =    [["Uniform", 0.2, 5.0, 0.2], []]
Change_Par =   [["Gauss",0.05, 0.1, 0.0, 0.2, 0.01],[]]
SS_Par0 =      [["Uni_ND",3, [[0.0, 1.0, 0.1]]],[],[]]
Error_Par =    [Dev_Par, Scale_Par, Change_Par]
Run_Par =      [max_iter, sample_size, scaling, discard, project_dim]
Param = [SESCA_Dir, Error_Par, SS_Par0, Run_Par]

Def_Args = [Input_files, Output_files, Param, failmark, verbosity]


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

#function to import custom modules:
def Import_Custom():
	IMPORTS = []
#	import predefined modules: 
	if not SESCA_scripts in sys.path:
		Vprint(2, "Using SESCA modules in: %1s" %SESCA_Dir)
		sys.path.append(SESCA_scripts)
	
	try:
#		SESCA modules
#		SESCA_pred:
		globals()["Spred"] = __import__("SESCA_pred")
		IMPORTS.append("Pred")
#		SESCA_scale
		globals()["Sscale"] = __import__("SESCA_scale")
		IMPORTS.append("Scale")		
	except Exception:
		Vprint(1,"Cannot load requested modules")
		Vprint(1, "Please double check specified SESCA directory")

	try:
#		SESCA_main:
		globals()["Smain"] = __import__("SESCA_main")
		IMPORTS.append("Main")
	except Exception:
		Vprint(1,"Cannot import SESCA_main module, default basis sets disabled")
	try:
#		numpy module to speed up rescaling:
		Vprint(2, "\nLoading Modules:")	
		globals()["np"] = __import__("numpy")
		IMPORTS.append("numpy")
	except Exception:
		Vprint(1,"Warning: Importing  numpy failed!")


	return IMPORTS


#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
	New_Args = [Input_files, Output_files, Param, failmark, verbosity]
	FLAGS = ["spect","lib","par","corr","write","proj","data","pdim","iter","size","scale","discard","verb"]
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
#		change default input files:
		elif flag == "spect":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "lib":
			New_Args[0][1] = arg
			flag = ""	
		elif flag == "par":
			New_Args[0][2] = arg
			flag = ""	
		elif flag == "corr":
			New_Args[0][3] = arg
			flag = ""	
#		change default output files:
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "proj":
			if arg == "0":
                                New_Args[1][1] = ""
			else:
				New_Args[1][1] =arg
			flag = ""
		elif flag == "data":
			if arg == "0":
                                New_Args[1][2] = ""
			else:
				New_Args[1][2] = arg
#		change default run parameters:
		elif flag == "pdim":
			try:
				parts = arg.split(",")
				X,Y = int(parts[0])-1, int(parts[1])-1
				New_Args[2][3][4] = [X, Y]
			except Exception:
				Vprint(1, "@pdim only takes two comma separated integers")
				New_Args[3] = 1
			flag = ""
		elif flag == "discard":
			try:
				fraction = float(arg)
				if fraction >= 0.0 and fraction < 1.0:
                                	New_Args[2][3][3] = fraction
				else:
					New_Args[3] = 1
					Vprint(1, "@discard only takes fractions between 0 and 1")
			except Exception:
				Vprint(1, "@discard only takes fractions between 0 and 1")
				New_Args[3] = 1
			flag = ""
		elif flag == "iter":
			New_Args[2][3][0] = int(arg)
			flag = ""
		elif flag == "size":
			New_Args[2][3][1] = int(arg)
			flag = ""
		elif flag == "scale":
			if arg == "0":
                                New_Args[2][3][2] = ""
			elif arg == "auto":
				New_Args[2][3][2] = "auto"
			else:
				New_Args[2][3][2] = float(arg)
			flag = ""
		elif flag == "verb":
			New_Args[4] = int(arg)
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0] == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1] == "" and acnt == 1:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args

#imported function to read basis spectra
def Read_BS_File(basis_file, BS_List):
	BS_info = []
	BS_param = []
	BS_Data = []
#	check if basis set in on the default list:
	if basis_file in BS_List:
		Vprint(2, "\nUsing default basis set: %1s"%basis_file) 
		BS_info = BS_List[basis_file]
		Vprint(2, BS_info)
		BB_file = os.path.join(SESCA_lib,BS_info[1])
	else:
		BB_file = basis_file
		if os.path.isfile(BB_file) == True:
			Vprint(2,"\nUsing basis set: %1s"%BB_file)
		else:
			Vprint(1,"\nBasis set file not found: %1s"%BB_file)
			return BS_Data

#	read in basis set paramters:
	BS_param = Spred.Read_BS(BB_file, ["",""])
	if BS_param != []:
		b_num = BS_param[4]
		BS_Data = [BB_file, b_num, BS_param, BS_info]
	else:
		Vprint(1, "Unable to read basis set file: %1s"%BB_file)

	return BS_Data

#function to read in error parameters for the SS estimation:
def Read_BayesPar(File, Par0, SS0):
	New_pars = [ [[],[]], [[],[]], [[],[]] ]
	New_SS = [[],[],[]]
	

#	check if file exists:
	Vprint(1, "\nReading in parameter file:",)
	if os.path.isfile(File) == False:
		Vprint(1, "File not found")
		return None

	Vprint(1, File)
	f =open(File,"rb")
	flag = ""
	for line in f:
		line2 = str(line.decode("ascii")).strip("\n")
#		read file line by line, mark lines based on keywords:
		if line2.startswith("#|Scale Dist"):
			flag = "scale"
		elif line2.startswith("#|NonSS Dist"):
			flag = "cont"
		elif line2.startswith("#|SS Dist"):
			flag = "struct"
		elif line2.startswith("#|SS change"):
			flag = "change"
		elif line2.startswith("#|SS project"):
			flag = "project"
		elif line2.startswith("&"):
			flag = ""
		elif flag != "" and line2 == "":
			flag = ""


		elif not line2.startswith("#") and line2 != "" and flag != "":
			parts = line2.split()
			pnum = len(parts)
#			read lines that request guassian distribution, with mean, width, Xmin, and Xmax 
			if parts[0] == "Gauss":
				try:
					Temp = [parts[0], float(parts[1]), float(parts[2]), float(parts[3]), float(parts[4]), float(parts[5])]
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []
#			read lines that request a uniform distribution within Xmin-Xmax
			elif parts[0] == "Uniform":
				try:
					Temp = [parts[0], float(parts[1]), float(parts[2]), float(parts[3])]
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []
#			read lines that request uniform N dimensional distribution within Xmin-Xmax and step Xstep in each dim:
			elif parts[0] in ["Uni_ND","Custom_ND"]:
				try:
					Temp = [parts[0], int(parts[1]), float(parts[2]), float(parts[3]), float(parts[4])]
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []	
#			read lines that request a custom 1D distribution within Xmin-Xmax and step Xstep
			elif parts[0] == "Custom_1D":
				try:
					Temp = [parts[0], float(parts[1]), float(parts[2]), float(parts[3])]
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []
#			read lines that request a custom 2D distribution within Xmin-Xmax and step Xstep in each dim.
			elif parts[0] == "Custom_2D":
				try:
					Temp = [parts[0], 2,[[float(parts[1]), float(parts[2]), float(parts[3])], [float(parts[4]), float(parts[5]), float(parts[6])]]]
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []
#			read lines that define class mapping for SS classes
			elif parts[0] == "SS_map":
				try:
					Temp = [parts[0], []]
					for j in range(1,3):
						mapped = int(parts[j])
						Temp[1].append(mapped)
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []
#			 read lines from discrete distribution / error parameters
			else:
				try:
					Temp = []
					for j in range(pnum):
						Temp.append(float(parts[j]))
				except Exception:
					Vprint(2, "Warning could not read line:",line2)
					Temp = []

#			assign data as appropriate based on the current keyword 
			if Temp != []:
				if flag == "cont" and pnum > 3:
					New_pars[0][0] = Temp
				elif flag == "cont" and pnum <= 3:
					New_pars[0][1].append(Temp)
				elif flag == "scale" and pnum > 2:
					New_pars[1][0] = Temp
				elif flag == "scale" and pnum <= 2:
					New_pars[1][1].append(Temp)
				elif flag == "change" and pnum > 2:
					New_pars[2][0] = Temp
				elif flag == "change" and pnum <= 2:
					New_pars[2][1].append(Temp)
				elif flag == "struct" and Temp[0] in ["Uni_ND", "Custom_ND"] and pnum > 3:
					Bound_1d = [Temp[2], Temp[3], Temp[4]]
					if New_SS[0] == []:
						New_SS[0] = [Temp[0], Temp[1], [Bound_1d] ]
					else:
						New_SS[0][2].append(Bound_1d)
				elif flag == "struct" and pnum >= 3:
					New_SS[1].append(Temp)
				elif flag == "project" and Temp[0] == "SS_map":
					New_SS[2] = Temp

	 
	f.close()
#	if the parameter file did not contain data, fill out with defaults:
	blocks = len(New_pars)
	for i in range(blocks):
#		check Error parameters separately:
		if New_pars[i] == [[],[]]:
			New_pars[i] = Par0[i]
#	check general SS parameter:
	if New_SS == [[],[],[]]:
		New_SS = SS0
	elif New_SS[0] == [] and New_SS[1] == []:
		New_SS[0] = SS0[0]

		

	Vprint(4,"\nError parameters:\n", New_pars)
	Vprint(4,"\nSS parameters:\n", New_SS)
		
	return [New_pars, New_SS]

#functior to generate a discretized 1D gaussian distribution:
def Gauss_1D(X,Param):
	Mean, SDev, Int = Param
	Fx = Int/math.sqrt(2*math.pi*SDev**2) * math.exp(-0.5*((X-Mean)/SDev)**2)
	return Fx

#function to calculate 2D Gaussian:
def Gauss_2D(Point, Param):
	Int, Mx, Sx, My, Sy = Param
	X,Y = Point[0], Point[1]
	Fxy = Int/(2*math.pi*Sx*Sy) * math.exp(-0.5*((X-Mx)/Sx)**2 -0.5*((Y-My)/Sy)**2)
	return Fxy

def Gauss_Dist(Mean, SDev, Bounds):
	Distribution = []
	Fmin, Fmax, Fstep, Fint  = Bounds
#	set up boundaries if they were not defined
	if Fmin == "":
		Fmin = Mean-5*SDev
	if Fmax == "":
		Fmax = Mean+5*SDev
	if Fstep == "":
		Fstep = float(Fmax-Fmin)/50
	if Fint == "":
		Fint = 1.0
	Pars = [Mean,SDev,Fint]
#	determine values at bin xi (left-stairs):
	xi = Fmin
	while xi < Fmax:
		Fxi = Gauss_1D(xi,Pars)
		Fxj = Gauss_1D((xi+Fstep),Pars)
		Fx_bin = 0.5*(Fxi+Fxj)
		x_bin = xi + 0.5*Fstep
		Distribution.append([x_bin, Fx_bin])	
		xi += Fstep
		 
	return Distribution



#project down an ND distribution to a lower dimension:
def Project_ND(ND_Data, Low_Dims):
	low_dim = len(Low_Dims)
	Vprint(2, "\nProjecting to coordinates:",Low_Dims)
#	determin if the ND distribution has extra info:
	if len(ND_Data) <= 3:	
		high_dim = len(ND_Data[1][0])-1
		if len(ND_Data) == 3:
			Param = ND_Data[0]
	else:
		high_dim = len(ND_Data[0])-1
		Param = ["None"]

	bin0 = len(ND_Data)
#	check that low dim is larger then high dim:
	if low_dim > high_dim:
		Vprint(1,"cannot project to higher dimension")
		return LD_Data

#	prepare LD array:
	Vprint(4, "Projecting %1dD data to %1d dimensions:"%(high_dim, low_dim))
	LD_Data = []
	Bin_codes = {}
	added = 0
	for hi_bin in ND_Data[1]:
		Vprint(5, hi_bin)
#		project down hi-D bin numbers
		low_bin = []
		bin_code = ""		
		for d in range(high_dim):
			if d in Low_Dims:
				bin_d = hi_bin[d]
				low_bin.append(bin_d)
				bin_code += "%1d,"%bin_d
		low_bin.append(0.0)
#		add new bin if this is first projection:
		if not low_bin in LD_Data:
			LD_Data.append(low_bin)
			Bin_codes[bin_code] = added
			added += 1

#	Fill LD array:
	for hi_bin in ND_Data[1]:
#		find LD bin number
		low_code = ""
		for d in range(high_dim):
			if d in Low_Dims:
				low_code += "%1d,"%hi_bin[d]
		low_bin_num = Bin_codes[low_code]
#		add probability to the LD bin:
		P_hi = hi_bin[high_dim]
		LD_Data[low_bin_num][low_dim] += P_hi

#	get highest and lowes probabilities:
	Vprint(5, "Projected data:")
	Pmin, Pmax = ["",""]
	for entry in LD_Data:
		Vprint(5, entry)
		P_low = entry[low_dim]
		if Pmin == "" or Pmin > P_low:
			Pmin = P_low
		if Pmax == "" or Pmax < P_low:
			Pmax = P_low
	low_info = [low_dim, added, Pmin,Pmax]
	Bin_codes["matpar"] = low_info
	Vprint(4, low_info)

#	copy array info if it was provided:
	LD_Param = ["Custom_ND", low_dim, []]
	if Param[0] in ["Uni_ND", "Data_ND"]:
		Min, Max, Step = Param[2],Param[3],Param[4]
		for d in Low_Dims:
			LD_Param[2].append([Min,Max,Step])
	if Param[0] in ["Custom_ND"]:
		for d in Low_Dims:
			LD_Param[2].append(Param[2][d])
			

	return [LD_Param, LD_Data, Bin_codes]


#create 1D histogram:
def Histogram_Data(Data,Dparam):
	Dmin, Dmax, Dstep,Dnorm= Dparam
	Histogram = []
	data_num = len(Data)
	if Dnorm == 1:
		Hadd = 1.0/(data_num*Dstep)
	else:
		Hadd = 1.0
	cnt = Dmin
	Vprint(2, "\nGenerating Histogram")
#	generate bins:
	while cnt <= Dmax:
		Histogram.append([cnt,0.0])	
		cnt += Dstep

#	fill bins:
	for Pi in Data:
		bnum = int(math.floor(float(Pi-Dmin)/Dstep))
		Histogram[bnum][1] += Hadd	
		Vprint(4, "%1.3f assigned to bin %1d (%1.3f)"%(Pi,bnum,Histogram[bnum][0]))

	return Histogram

#create 2D heatmap:
def Heatmap_Data(Data, Dparam):
	Vprint(2, "Generating heatmap")
	Output = ""
	Param = Data[0][2]

#	prepare file header:
	Files, Run_P, Labels = Dparam
	Header0 =  "#SS heatmap generated by SESCA_bayes:\n"
	Header0 += "#Workdir:  %10s\n" % Run_P[5]
	Header0 += "#Spectrum file:        %10s\n"%Files[0]
	Header0 += "#Basis set file:       %10s\n"%Files[1]
	if Files[2] != "":
		Header0 += "#Parameter file:       %10s\n" % Files[2] 
	if Files[3] != "":
		Header0 += "#Correction file:      %10s\n" % Files[3]		
	Header0 += "#Initial sample size: %1d\n#Monte Carlo steps:     %1d\n"% (Run_P[1],Run_P[0])
	if Run_P[3] != 0:
		Header0 += "#Fraction of discarded samples:   %1.3f\n"%Run_P[3]
	if Run_P[2] == 0:
		Header0 += "#Spectrum scaling was disabled\n"
	Header0 += "\n#Projection along SS classes:\n"
	Header0 += "#%10s fraction (vertical)\n"%Labels[0]
	Header0 += "#%10s fraction (horizontal)\n"%Labels[1]
	Header0 += "#Bin probabilities are expressed in percents\n"
	Output += Header0
	
#	determine bin labels:
	Labels = [[],[]]
	Vprint(4,Data[0])
	for d in range(2):
		Min, Max, Step = Param[d][0], Param[d][1], Param[d][2]
		if Min != "" and Max != "" and Step != "":
			bin_value = Min + 0.5*Step
			while bin_value <= Max:
				Labels[d].append(bin_value)
				bin_value += Step
			Vprint(3, Labels[d])

		else:
			Vprint(1,"Missing parameters for heatmap")
			return Output

#	generate heatmap
	Heatmap = []
	bin_num1 = len(Labels[0])
	bin_num2 = len(Labels[1])
	for j in range(bin_num1):
		Tmp = []
		for k in range(bin_num2):
			Tmp.append(-0.0)
		Heatmap.append(Tmp)

# 	fill viable coordinates:
	for entry in Data[1]:
		j, k, Pjk = entry
		Heatmap[j][k] = Pjk * 100

#	print heatmap:	
	Header = "#| Hist "
	Map = ""
	for k in range(bin_num2):
		bin_val = Labels[1][k]
		Header +=  " %6.3f " % bin_val

	for j in range(bin_num1):
		bin_val = Labels[0][j]
		New_line = " %6.3f " % bin_val
		for entry in Heatmap[j]:
			New_line += " %6.3f " % entry
		Map += New_line+"\n"
	Output += Header+"\n"+Map+"\n"
		 
	return Output
	
#function to complete error parameters based on the current data:
def Complete_ErrDist(Err):
#	generate distribution based on parameters:
	if Err[0] != [] and Err[1] == []:
		Vprint(3, "\nGenerating error distribution:")
#		create 1D Gaussian distribution with mean and SD within a given boundary:
		if Err[0][0] == "Gauss":
			Mean, SD, Bounds = Err[0][1], Err[0][2], [Err[0][3], Err[0][4], Err[0][5], ""]
			Dist = Gauss_Dist(Mean,SD,Bounds)
			Err[1] = Dist
#		create uniform distribution within a given boundary:
		if Err[0][0] == "Uniform":
			Min, Max, Step = Err[0][1], Err[0][2], Err[0][3]
			bins = float(Max-Min)/Step
			Value = float(1.0)/(bins+1)
			cnt = Min
			while cnt <= Max:
				Point = [cnt, Value]
				Err[1].append(Point)
				cnt += Step
	if Err[1] != [] and Err[0] == []:
#		check dimensionality:
		dim = len(Err[1][0])-1
		if dim == 1:
			Temp = ["Custom_1D"]
		elif dim == 2:
			Temp = ["Custom_2D",2, []]
		elif dim >= 3:
			Temp = ["Custom_ND",dim, []]	
#		get basic step parameters on custom probability dist.
		for d in range(dim):
			Min = ""
			Max = ""
			Step = ""
			P_Set = set()
			for P in Err[1]:
				P_Set.add(P[d]) 
				if Min == "" or Min > P[d]:
					Min = P[d]
				if Max == "" or Max < P[d]:
					Max = P[d]
			Vprint(4,P_Set)
			P_vals = list(P_Set)
			P_vals.sort()
			Vprint(5, P_vals)
			Step = P_vals[1] - P_vals[0]
			Dmin = Min - 0.5*Step
			Dmax = Max + 0.5*Step
			if Temp[0] in "Custom_1D":
				Temp.append([Dmin, Dmax, Step])
			elif Temp[0] in ["Custom_2D", "Custom_ND"]:
				Temp[2].append([Dmin, Dmax, Step])


		Err[0] = Temp
		Vprint(4, "Custom distribution parameters:",Err[0])
			

#	print distrtibution histogram:
	dim = len(Err[1][0])-1
	format_string = dim*" %1.3f "+ " %1.3e "
	for Point in Err[1]:
		Vprint(5, format_string%tuple(Point))
	
	return Err

#function to generate an N-dimensional container with N bins along each dimension:
def Generate_NDmat(dim, bin_num):
	Temp_prev = [[0.0]]
	for N in range(dim):
#		expand array dimension by one:
		Temp = []
		for j in range(bin_num):
			for entry in Temp_prev:
#				add new coordinate:
				new_entry = [j]
#				add old coordinates:
				for coord in entry:
				 	new_entry.append(coord)
				Temp.append(new_entry)
		Temp_prev = Temp
	ND_mat = Temp_prev

	return ND_mat


#function for generating uniform distribution:
def SS_Index(SS_par):
	Nd_hist = []
	Codes = {}
	SS_Mat = []
	SS_bounds = []
	added = 0
	dim,Min,Max,Step = [0,"","",""]
#	check if boundaries were supplied:
	if SS_par[0] == []:
		Vprint(1, "No boundaries defined for prior SS distribution!")
		return Nd_hist

	if SS_par[0][0] in ["Uni_ND", "Custom_ND"]:
#	check dimension and boundaries:
		Type = SS_par[0][0]
		dim = SS_par[0][1]
		b_num = len(SS_par[0][2])
		SS_bounds = [Type, dim, []]
		if b_num == dim:
#		copy boundaries in each dimension if provided:
			for d in range(dim):
				Min,Max,Step = SS_par[0][2][d]
#				make a basic sanity check on the boundaries defined:
				if Min < Max and Min >= 0.0 and Max <= 1.0 and Step > 0.0:
					SS_bounds[2].append([Min,Max,Step])
				else:
					Vprint(1, "Wrong boundary parameters:",SS_par[0][2][d])
					return Nd_hist
		elif b_num == 1 and dim > b_num:
#		add boundaries in every dimension if only one was supplied:
			Min,Max,Step = SS_par[0][2][0]
			if Min < Max and Min >= 0.0 and Max <= 1.0 and Step > 0.0:
				for d in range(dim):
					SS_bounds[2].append([Min,Max,Step])
			else:
				Vprint(1, "Wrong boundary parameters:",SS_par[0][2][d])
		else:
#		quit if the dimensions are off:
			Vprint(1,"The number of boundary entries (%1d) does not match the number of dimensions (%1d)!"%(b_num,dim))
			return Nd_hist
	
	else:
#		quit if we dont have the right type:
		Vprint(1, "%1s is not an approriate type for the prior SS distribution, please check parameter file!"%SS_par[0][0])
		return Nd_hist			

	if SS_par[0][0] in ["Uni_ND"]:
#		determine bin number in 1 dimension:
		bin_num = int(math.ceil((Max-Min)/Step))
#		define possible sum of bin numbers, if the division is not precise:
		bin_max = 1.0
		bin_min = 1.0 - dim*Step


		Vprint(2, "Mesh size for %1dD SS array: %1d" % (dim,bin_num))
		ND_matrix = Generate_NDmat(dim,bin_num)
		Vprint(3,len(ND_matrix))

#		collect matix elements that refer to possible SS compostiions:
		for entry in ND_matrix:
			dist_sum = 0
			code = ""
			for d in range(dim):
				dist_sum += entry[d]*Step	
				code += "%1d,"%entry[d]
			if dist_sum <= (bin_max) and dist_sum >= (bin_min):
#				if the sum of SS fractions add up to 1, add them:
				SS_Mat.append(entry)
#				also add code for the helper function
				Codes[code] = added
				added += 1

#	normalize PDF and return matrix:
	if SS_par[0][0] == "Uni_ND" and added != 0:
		Vprint(2, "\nGenerating uniform SS distribution for %1d classes"%dim) 
		for entry in SS_Mat:
			entry[-1] = float(1.0)/added
			Vprint(4, entry)
		Mat_par = [dim, added, float(1.0)/added, float(1.0)/added]
		Codes["matpar"] = Mat_par
		Nd_hist = [SS_bounds, SS_Mat, Codes]

	elif SS_par[0][0] in ["Custom_ND"]:
		Vprint(2, "\nProcessing custom SS distribution for %1d classes"%dim) 
#		reduce bin coordinates to integer codes:
		for entry in SS_par[1]:
			new_entry = []
			for d in range(dim):
				Min,Max,Step = SS_bounds[2][d]		
				code = int(math.floor((entry[d]-Min)/Step))
				new_entry.append(code)
			new_entry.append(entry[dim])
			Vprint(4, "reduced bin entry:",new_entry)
			SS_Mat.append(new_entry)

#		generate indexing 
		Codes = Dist_List([SS_bounds, SS_par[1]])
		Mat_par = Codes["distpar"]
		Codes["matpar"] = Mat_par
		Nd_hist = [SS_bounds, SS_Mat, Codes]

#	print number of bins, minimum and maximum SS probability:
	Vprint(3, "\n%1dD matrix, bins: %1d, Pmin: %1.4f  Pmax: %1.4f" % tuple(Codes["matpar"]))
	Vprint(3, "\nBoundaries in each dimension are set:")
	for Bound in SS_bounds[2]:
		Vprint(3, "Min: %1.3f Max: %1.3f Step: %1.3f"%tuple(Bound))	
		
	return Nd_hist

#function to generate random SS compositions:
def Generate_randomSS(c_num):
	SS_comp = []
	Sum = 0.0
	for i in range(c_num):
		value = random.uniform(0.0, 100.0)
		SS_comp.append(value)
		Sum += value
	
	for i in range(c_num):
		SS_comp[i] = SS_comp[i]/Sum
	Vprint(5, SS_comp)
	return SS_comp

# helper function to determine the bin number of an SS composition:
def Bin_SScomp(SS, SS_dist):
#	get binning parameters:	
	label,dim = SS_dist[0][0], SS_dist[0][1]
	Bounds = []
	if label == "Uni_ND" and len(SS_dist[0][2]) == 1:
		Min,Max,Step = SS_dist[0][2][0]
		Bounds = []
		for d in range(dim):
			Bounds.append([Min,Max,Step])
	else:
		Bounds = SS_dist[0][2]

	Coords = []
	Code = ""
	for d in range(dim):
		Coeff = SS[d]
		Min,Max,Step = Bounds[d]
#		determine bin coordinates:
		Cnum_i = int(math.floor(float(Coeff-Min)/Step))
		if Coeff >= 1.0:	
			Cnum_i = Cnum_i-1
		Coords.append(Cnum_i)
		Code += "%1d,"%Cnum_i
#	acquire bin number based on the coordinates:
	try:	
		bin_num = SS_dist[2][Code]
		Vprint(5,"bin coordinates:",Coords,"bin number: %1d"%bin_num)
	except Exception:
		bin_num = ""
		Vprint(2, "Warning: bin coordinates cannot be detemined for:",Coords)
		
	return bin_num

#function to generate SS compositions to sample a distribution:
def Sample_SSdist(SS_dist, size):
	SS_samples = []
	dim, bin_num, Pmax, Pmin = SS_dist[2]["matpar"]
	Vprint(4, "Generating %1d SS compositions..."%size)
	cnt = 0
	while cnt < size:
#		generate random SS
		SS_j = Generate_randomSS(dim)
#		check which bin it would belong to:
		bin_num = Bin_SScomp(SS_j, SS_dist)
#		determine SS probability in the current distribution:
		if bin_num != "":
			Prob_SSj = SS_dist[1][bin_num][dim]
		else:
			bin_num = -1
			Prob_SSj = 0.0
#		accept or reject SS composition based on its probability:
		Random_val = random.uniform(0.0, Pmax)
		if Random_val <= Prob_SSj:
			Vprint(5, "SS comp accepted with R:",Random_val,"(Pssj:",Prob_SSj,")")
			SS_samples.append([bin_num, SS_j, Prob_SSj])
			cnt += 1
		else:
			Vprint(6, "SS comp rejected with R:",Random_val,"(Pssj:",Prob_SSj,")")


	return SS_samples

#function to calculate the most likely SS composition based on the sample:
def Average_SS(SS_sample):
	size = len(SS_sample)
	SS_0 = SS_sample[0][1]
	dim = len(SS_0)

#	Set up Avg SS array:
	Avg_SS = []
	Vprint(3, "\nAveraging SS samples:")
	for i in range(dim):
		Avg_SS.append([SS_0[i][0], 0.0, 0.0])

#	compute the averages:
	P_sum = 0
	for j in range(size):
		SS_j = SS_sample[j][1]
		P_sum += 1
		for i in range(dim):
			Avg_SS[i][1] += SS_j[i][1]
		

	for i in range(dim): 
		if P_sum != 0:
			Avg_SS[i][1] = Avg_SS[i][1] / P_sum
		else:	
			Avg_SS[i][1] = 0.0

#	compute standard deviations:	
	for j in range(size):
		SS_j = SS_sample[j][1]
		for i in range(dim):
			Avg_SS[i][2] += math.pow((Avg_SS[i][1] - SS_j[i][1]),2)
		
	for i in range(dim): 
		if P_sum != 0:
			Avg_SS[i][2] = math.sqrt(Avg_SS[i][2] / P_sum)
		else:	
			Avg_SS[i][2] = 0.0

#	print and return Avg SS
	Vprint(3,"Estimated SS composition:")
	for i in range(dim):
		Vprint(3,"%12s: %1.3f +/- %1.3f"%tuple(Avg_SS[i]))

	return Avg_SS


#function to update the SS dustribution:
def Update_SSdist(Points, SS_dist):
	dim, bin_num, Pmax, Pmin = SS_dist[2]["matpar"]
	Bounds = []
	if SS_dist[0][0] in ["Uni_ND","Custom_ND"] and len(SS_dist[0][2]):
		Min,Max,Step = SS_dist[0][2][0]
		Bounds = []
		for d in range(dim):
			Bounds.append([Min,Max,Step])
	else:
		Bounds = SS_dist[0][2]	

	Vprint(3, "\nUpdating SS probability distribution:")

#	Create temporary array for probability averaging:
	Temp_SS = []
	for n in range(bin_num):
		Temp_SS.append(0.0)

#	Increase posterior probabilities for each visited bin:
	for Point_j in Points:
		ID_j = Point_j[0]
		SS_j = Point_j[1]
		Post_j= Point_j[5]
		Temp_SS[ID_j] += 1.0

#	Update the probability of visited bins:
	P_sum = 0.0
	new_SSdist = []
	New_Codes = {}
	for n in range(bin_num):
		b_curr = SS_dist[1][n]
		b_new = []
		P_curr = b_curr[dim]
		new_code = ""
#		add coordinates:	
		for i in range(dim):
			b_new.append(b_curr[i])
			new_code += "%1d," % b_curr[i]
#		add prbability
		P_new = Temp_SS[n]
		Vprint(5, "bin %1d probability updated:"%(n+1),P_new," (",P_curr,")") 
#		or set to zero if the bin was not visited:
		b_new.append(P_new)
		new_SSdist.append(b_new)
		New_Codes[new_code] = n
		P_sum += P_new

#normalize probabilities:
	P_Bounds = ["", ""]
	if P_sum != 0.0:
		for n in range(bin_num):
			b_new = new_SSdist[n]
			b_new[dim] = b_new[dim]	/ P_sum
			Vprint(5,b_new)
#			get minimum and maximum probability:
			P_new = b_new[dim]
			if P_Bounds[0] == "" or P_Bounds[0] > P_new:
				P_Bounds[0] = P_new
			if P_Bounds[1] == "" or P_Bounds[1] < P_new:
				P_Bounds[1] = P_new
	Vprint(4, "new probability boundaries:",P_Bounds)

#update Coding block:
	new_par = [dim, bin_num, P_Bounds[1], P_Bounds[0]]
	New_Codes["matpar"] = new_par 
	New_Type = ["Custom_ND", dim, Bounds]

	New_Dist = [New_Type, new_SSdist, New_Codes]

	return New_Dist

#function to update scalinf factor distributions:
def Update_SFdist(Points, old_Dist):
	New_Dist = []
	Vprint(2, "Updating Scaling factor distribution:")

#	extracting binning parameters:
	Vprint(5, old_Dist)	
	Type = ["Custom_1D","","",""]
	if old_Dist[0][0] == "Gauss":
		Min, Max, Step = old_Dist[0][3], old_Dist[0][4], old_Dist[0][5]
	elif old_Dist[0][0] in ["Custom_1D","Uniform"]:
		Min, Max, Step = old_Dist[0][1], old_Dist[0][2], old_Dist[0][3]
	Type = ["Custom_1D",Min,Max,Step]
	Vprint(4, Type)
	New_Dist = [Type, []]

#	start new distribution:	
	bin_num = len(old_Dist[1])
	for entry in old_Dist[1]:
		X_val = entry[0]
		New_Dist[1].append([X_val,0.0])	
	
#	fill up distibution:
	p_num = len(Points)
	Stats = [0.0, 0.0]
	if p_num != 0 and bin_num != 0:
		for Pj in Points:
#			determine which bin:
			SF_j = Pj[2]
			bin_ID = int(math.floor((SF_j-Min)/Step))
			Vprint(4,"\nSF_j: %1.3f bin: %1d,"%(SF_j,bin_ID))

#			increase bin probability:
			if bin_ID >= 0 and bin_ID < bin_num:
#				New_Dist[1][bin_ID][1] += 1.0/(Step*p_num)
				New_Dist[1][bin_ID][1] += 1.0/(p_num)
			else:
				Vprint(2,"Warning, Scaling factor skipped, out of range")

#			determine average:
			Stats[0] += SF_j/p_num

		Dev = 0.0
		for Pj in Points:
			SF_j = Pj[2]
			Dev += (SF_j - Stats[0])**2
		Stats[1] = math.sqrt(Dev / p_num)
	New_Dist.append(Stats)

	Vprint(4, New_Dist)

	return New_Dist

#function for rejection sampling from 1D distribution function:
def Rejection_1D(Dist):
#check	boundaries:
	if Dist[0][0] == "Gauss":
		Min, Max, Step, Curve = Dist[0][3], Dist[0][4], Dist[0][5],  Dist[1]
	elif Dist[0][0] in ["Custom_1D", "Uniform"]:
		Min, Max, Step, Curve = Dist[0][1], Dist[0][2], Dist[0][3],  Dist[1]
	else:
		Vprint(1, "Warning, Unknown Function type:", Dist[0]) 
	Move_Accept = ""
#determine maximum probability:
	Pmax = 0.0
	for P in Curve:
		if P[1] > Pmax:
			Pmax = P[1]
	cnt = 0
#randomly select bin, and test probability:
	while cnt < 1:
		Move, R1 = random.uniform(Min,Max), random.uniform(0.0,Pmax)
		Move_bin = int(math.floor((Move-Min)/Step))
		Move_prob = Curve[Move_bin][1]
		if R1 < Move_prob:
			Move_Accept = Move 
			cnt += 1
	return Move_Accept


def Update_SS(SS0, Step_dist):
#	first Copy coefficients:
	SS_new = []
	dim = len(SS0)
	for i in range(dim):
		Label, Coeff = SS0[i]
		SS_new.append([Label,Coeff])

#	then select two random classes:
	Shift = []
	cnt = 0
	while cnt < 2:
		Pick = random.randrange(dim)
		Candidate = [Pick,SS_new[Pick][1]]
		if not Candidate in Shift:
			Shift.append(Candidate)
			cnt += 1
#	switch roles if the first condidate has zero fraction and the second doesnt:
	if Shift[0][1] == 0.0 and Shift[1][1] != 0:
		Shift = [Shift[1], Shift[0]]	
#	determine the fraction of SS to be moved by rejection sampling:
	Move_Accept = Rejection_1D(Step_dist)
#	make sure we dont move a larger fraction than there is in the class:	
	Move_final, i,j = min(Shift[0][1], Move_Accept), Shift[0][0], Shift[1][0]
#	change class coefficents:
	SS_new[i][1] = SS_new[i][1] - Move_final
	SS_new[j][1] = SS_new[j][1] + Move_final
	Vprint(4,"SS composition updated: %1.3f moved %12s -> %12s"%(Move_final, SS_new[i][0], SS_new[j][0]))
	Vprint(5, SS0,"\n",SS_new)	
	return SS_new

#function to determine a new scaling factor:
def Update_SF(SF0, Step_dist, Bounds):
	SF_new = SF0
#	change the current value by a random number:
	Change = Rejection_1D(Step_dist)
	sign = random.choice([1.0, -1.0])
	SF_new = SF0 + sign*Change
	if SF_new < Bounds[0] or SF_new > Bounds[1]:
		SF_new = SF0

	return SF_new
	
#function to determine prior probability for scaling factors:
def Get_SFprior(SF, SF_dist):
#	Calculate Scaling factor probability:
	if SF_dist[0][0] == "Gauss":
		Min, Max, Step, Curve = SF_dist[0][3], SF_dist[0][4], SF_dist[0][5], SF_dist[1]
	elif SF_dist[0][0] in ["Custom_1D","Uniform"]:
		Min, Max, Step, Curve = SF_dist[0][1], SF_dist[0][2], SF_dist[0][3], SF_dist[1]
#	make sure we are not out of range:
	if SF < Min or SF > Max:
		SF_prob = 1.0e-15
		Vprint(5, "Scaling factor %1.3f, Out of range: (%1.3f, %1.3f)"% (SF,Min,Max),"\nSetting P(SF)= e-15")
	else:
		SF_bin = int(math.floor((SF - Min)/Step))
		SF_prob = Curve[SF_bin][1]
		Vprint(5, "Scaling factor %1.3f, Probability:"% SF, Curve[SF_bin])
	
	return  SF_prob


# helper function to sort intensities for spectrum rescaling:
def Prep_Spectra(Target, Calculated):
#	match spectra using the SESCA_pred
	Matched_Spectra = Spred.Compare_Spectra(Target, Calculated)
	Ints = [[],[],[]]
#	take wavelengths and intensities
	for Point in Matched_Spectra[0]:
		Ints[0].append(Point[0])
		Ints[1].append(Point[2])
		Ints[2].append(Point[1])

	return Ints

#helper function to automaticaly re-scale CD spectra 
def Rescale_spectra(Exp_Spect, Calc_Spect, scaling):	
#	define paramteres for automated intensity scaling:
	if scaling == "auto": 
		Scale_Ints = Prep_Spectra(Exp_Spect, Calc_Spect)
		SF0,SF_range,Lambda,Scale_iter,norm = [1.0,[0.2,5], 100, 1000, 0]
		Sscale.Set_verb(0)
		Scale_imports = Sscale.Import_Custom()
		Scaled_data = Sscale.Find_Scaling(SF0, Scale_Ints, SF_range, Lambda, Scale_iter,norm,Scale_imports)
		Exp_SF = Scaled_data[0]
		Scaled_Spect = Scaled_data[2]
#	if provided, apply requested scaling factor:
	elif not scaling in [0,"",1.0]:
		Exp_SF = scaling
		Scaled_Spect = Spred.Scale_Spect(Exp_Spect, Exp_SF)
#	otherwise, do not scale
	else:
		Exp_SF = 1.0
		Scaled_Spect = Exp_Spect

	return [Exp_SF, Scaled_Spect]
	
#function to calculate Spectrum likelihoods:		
def Spect_Likelihood(Exp_Spect, Corr_Spect, scaling, SS_Coeffs, Basis_Spect, Error_dist):
#		calculate theoretical CD spectrum
		Spred.Set_verb(0)
		Pred_Spect = Spred.Compute_CD(SS_Coeffs,Basis_Spect)

#		Add (Side chain) correction to theoretical spectrum:
		if Corr_Spect != [] :
			Calc_Spect = Spred.Modify_Spectra(Corr_Spect, Pred_Spect, 1.0)[0]
		else:
			Calc_Spect = Pred_Spect		 
		
#		Re-scale the experimental spectrum:
		SF_j, Scaled_Spect = Rescale_spectra(Exp_Spect, Calc_Spect, scaling)									

#		calculate deviation from measured spectrum:
		Dev_j = Spred.Compare_Spectra(Scaled_Spect, Calc_Spect)
		Rmsd_j = Dev_j[1][1]	

#		Determine the bounds and dimensions for the joint probability function:
		R_bounds = []
		if Error_dist[0][0] == "Gauss":
			R_bounds  = [[Error_dist[0][3], Error_dist[0][4], Error_dist[0][5] ]]
		elif Error_dist[0][0] in ["Custom_1D", "Uniform"]:
			R_bounds = [[Error_dist[0][1], Error_dist[0][2], Error_dist[0][3] ]]
		elif Error_dist[0][0] in ["Gauss_2D", "Custom_2D", "Custom_ND"]:
			R_bounds = Error_dist[0][2]
		dimension = len(R_bounds)
		Vprint(4, "Error bounds:",R_bounds)

#		determine coordinates:
		R_bin, S_bin, penalty =  "","",0
#		look up nonSS coordinates:
		R_min,R_max, R_step = R_bounds[0]
		if Rmsd_j <= R_min:
			R_bin = 0
			penalty += int(math.ceil(float(R_min-Rmsd_j)/R_step))	
		elif Rmsd_j < R_max:
			R_bin = int(math.floor(float(Rmsd_j - R_min)/R_step))
		else:
#			if Rmsdj is too large, give an icreasing penalty for each increment
			R_bin = int(math.ceil(float(R_max - R_min)/R_step))-1
			penalty += int(math.ceil(float(Rmsd_j-R_max)/R_step))
		if dimension == 1:
#			look up likelihood, based on the nonSS distribution:
			Rj, L_j = Error_dist[1][R_bin]
			Vprint(5, "Rmsd: %1.3f, Likelihood:"%Rmsd_j,Error_dist[1][R_bin])

		if dimension == 2:
#		for 2D probability distributions, determine scaling coordinate:
			S_min, S_max, S_step = R_bounds[1]
			if SF_j <= S_min:
				S_bin = 0
				penalty += int(math.ceil(float(S_min-SF_j)/S_step))	
			elif SF_j < S_max:
				S_bin = int(math.floor(float(SF_j - S_min)/S_step))
			else:
#			if SFj is too large, give an icreasing penalty for each increment
				S_bin = int(math.ceil(float(S_max - S_min)/S_step))-1
				penalty += int(math.ceil(float(SF_j-S_max)/S_step))

			code = "%1d,%1d," % (R_bin,S_bin)		
			Vprint(5, "Rmsdj: %1.3f, SFj: %1.3f, bins: %1s"%(Rmsd_j,SF_j, code[:-1]))
			try:
				P_bin = Error_dist[2][code]
				Rj, Sj, L_j = Error_dist[1][P_bin]				
				Vprint(5, "Likelihood:",Error_dist[1][P_bin])
			except Exception:
				Vprint(2, "Warning: Could not determine likelihood for:")
				Vprint(2, "Rmsdj: %1.3f, SFj: %1.3f, bins: %1s"%(Rmsd_j,SF_j, code[:-1]))
				L_j = Error_dist[2]["distpar"][2]
				Vprint(2, "Assigned probabilty: %1.3e"%L_j)
				

#		apply cumulative penalties for being out of range:
		if penalty != 0:
			L_mult = math.pow(10,(-1*penalty))
			Vprint(5, "%1.1e Likelihood penalty applied for error parameters out of range"%L_mult)
			L_j = L_j*L_mult

		return [Rmsd_j, SF_j, L_j]

#Helper function for Listing points in a distribution:
def Dist_List(Dist):
	Codes = {}
	
	Bounds = []
	if Dist[0][0] in ["Custom_2D","Custom_ND", "Gauss_2D"]:
		Bounds = Dist[0][2]
	elif Dist[0][0] in ["Uniform","Custom_1D"]:
		Bounds = [[Dist[0][1], Dist[0][2], Dist[0][3] ]]
	elif Dist[0][0] == "Gauss":
		Bounds = [[Dist[0][3], Dist[0][4], Dist[0][5] ]]

	d_num = len(Bounds)
	p_num = len(Dist[1])
	Vprint(4,"Dist parameters:\n",Dist[0])
	Vprint(4,"Dimensions: %1d, entry number: %1d"%(d_num, p_num))
	

#	go over each point:
	Vprint(3,"Sorting Points in distribution:")
	Pmin = ""
	Pmax = ""
	for j in range(p_num):
#		for each point add code based on coordinates in the distribution:
		code = ""
		Pj = Dist[1][j]
		for d in range(d_num):
			Min,Max,Step = Bounds[d]
			bin_d = int(math.floor((Pj[d]-Min)/Step))
			code += "%1d,"%bin_d
		Codes[code] = j

#		determine minimum and maximum for the distribution values:		
		Pval = Pj[d_num]
		if Pmin == "" or Pval < Pmin:
			Pmin = Pval
		if Pmax == "" or Pval > Pmax:
			Pmax = Pval
#	add distribution parameters to codes:	
	Dist_param = [d_num, p_num, Pmin, Pmax]
	Codes["distpar"] = Dist_param
	Vprint(4, "Distribution parameters:",Dist_param)	

	return Codes
	
#function to discard a franction of the inital points:
def Discard_init(Points, fraction):
	Filtered = []
	point_num = len(Points)
	if fraction >= 0 and fraction < 1:
		dis_num = int(math.ceil(point_num*fraction))
		Vprint(3, "Discarding first %1d points (%1.3f)"%(dis_num,fraction))
		for j in range(point_num):
			Pj = Points[j]
			if j >= dis_num:
				Filtered.append(Pj)
	else:
		Vprint(1, "invalid fraction to discard points.")
		return Points
	
	return Filtered	

# function to do MC-Hastings refinement of the SS distribution:
def Iterate_SS_dist(Exp_Spect, Corr_Spect, SS0, BS_data, Param):
	Max_iter, Samp_size, scaling, discard = Param[0], Param[1], Param[2], Param[3]
	Error_dist = Param[4][0]
	Scale_dist = Param[4][1]
	Change_Par = Param[4][2]
#	Vprint(2, SS0)
	dim, bin_num, Pmin, Pmax = SS0[2]["matpar"]
	Basis_Spect = BS_data[2][2]

#	add list to the Error distribution for quick lookup:
	Point_codes = Dist_List(Error_dist)
	P_dim, P_num, P_min, P_max = Point_codes["distpar"]
	Error_dist.append(Point_codes) 

#	start-up phase: draw SS comositions from a distribution:
	Vprint(2,"\nRefining SS probability distribution:")
	Vprint(2,"MC steps: %1d\nInitial sample size: %1d\nSS dimensions: %1d\nPmax: %1.3e"%(Max_iter, Samp_size, dim, Pmax))
	Points_all = []
	SS_dist = SS0
	Sample_prev = Sample_SSdist(SS_dist, Samp_size)
	SS_est_prev = []
	Points_prev = []
	Accept_all =  [0.0,[]]
#	prepare initial points:
	for j in range(Samp_size):
		N_j, SS_j, Prior_SSj, = Sample_prev[j]
		SF_j, Rmsd_j, Post_SSj = 1.0, 0.0, 0.0
		Point_j = [N_j, [], SF_j, Rmsd_j, Prior_SSj, Post_SSj]		
#	format coefficients:
		for i in range(dim):
			Label = Basis_Spect[(i+1)][0] 
			Point_j[1].append([Label,SS_j[i]])
		Points_prev.append(Point_j)
	
#	Sampling phase: iteratively improve SS distribution:
	for m in range(Max_iter):
		Vprint(3, "\n\n\nMonte Carlo step: %1d"%(m+1))

#		predict CD spectra for each sample, determine likelihood:
		Vprint(4, "\nUpdating SS samples:")
		Points = []
		accept_m = 0.0
		for j in range(Samp_size):
			Vprint(4,"\nSS comp. %1d, Step %1s"%(j+1, m+1))

#			Change coefficients:
			Coeff_prev = Points_prev[j][1]
			Coeff_j = Update_SS(Coeff_prev, Change_Par)
			Vprint(4,"Coeff_prev:",Coeff_prev)
	
#			update bin number:
			SS_coords = []
			for i in range(dim):
				SS_coords.append(Coeff_j[i][1])
			N_j = Bin_SScomp(SS_coords, SS_dist)
			if N_j != "":
#				determine SS prior if bin is found
				Prior_SSj = SS_dist[1][N_j][dim]
				Vprint(4, "Coeff_curr:", Coeff_j, "\nbin_num:",N_j, "bin:", SS_dist[1][N_j])	
			else:
#				give low probabilty if not found
				N_j = -1
				Prior_SS_j = 0.1*float(P_min)
				Vprint(2, "Warning, no bin found for Coefficients:",Coeff_j) 
			SF_j0 = scaling

#			compute the likelihood of the observed CD spectrum based on Coeff_j:
			Rmsd_j, SF_j, L_j = Spect_Likelihood(Exp_Spect, Corr_Spect, SF_j0, Coeff_j, Basis_Spect, Error_dist)
			Vprint(4, "Rmsd_j= %1.3f, SF_j= %1.3f, Likelihood= %1.3e"%(Rmsd_j, SF_j, L_j))			
			if m == 0:
#				in the first cycle, also calculate the likelihood of initial SS compositions:
				Coeff_prev = Points_prev[j][1]
				Prior_prev = Points_prev[j][4]
				Rmsd_prev, SF_prev, L_prev = Spect_Likelihood(Exp_Spect, Corr_Spect, scaling, Coeff_prev, Basis_Spect, Error_dist)
				Post_prev = L_prev * Prior_prev
				Points_prev[j][2] = SF_prev
				Points_prev[j][3] = Rmsd_prev
				Points_prev[j][5] = Post_prev
				
			#compute the posterior probabilty for this SS:
			Post_SSj = L_j * Prior_SSj
			Post_prev = Points_prev[j][5]
			Vprint(4,"New posterior: %1.3e, Old posterior: %1.3e"%(Post_SSj,Post_prev))

#			apply Metropolis criterion to the SS move:
			accept = 0
			if Post_prev != 0.0:
				P_ratio = min(1.0, Post_SSj / Post_prev)
			else:
				P_ratio = 1.0
			R_test = random.uniform(0.0,1.0)
			if R_test <= P_ratio:
				accept = 1

#			add either the new or the old SS composition to the Markov chain:
			if accept == 1:
				New = [N_j, Coeff_j, SF_j, Rmsd_j, Prior_SSj, Post_SSj]
				Points.append(New)
				Vprint(4, "Chain %1d,  new point accepted Pacc(%1.3f)"% (j, P_ratio))
				Vprint(4, "bin: %1d, SF: %1.3f, Rmsd: %1.3f, Prior: %1.3e, Posterior: %1.7e\n" % (New[0], New[2],New[3], New[4], New[5]))
				accept_m += float(1.0)/Samp_size 
			else:
				Old = []
				for entry in Points_prev[j]:
					Old.append(entry)
				Points.append(Old)
				Vprint(4, "Chain %1d, new point rejected Pacc(%1.3f), falling back to"%(j,P_ratio))
				Vprint(4, "bin: %1d, SF: %1.3f, Rmsd: %1.3f, Prior: %1.7e, Posterior: %1.7e\n" % (Old[0], Old[2],Old[3], Old[4], Old[5]))
		

#		Add points to the Sample:
		for Point in Points:
			Points_all.append(Point)
		Points_prev = Points

#		Calculate the average SS from the sample:
		SS_est_m = Average_SS(Points_all)
		SS_est_prev = SS_est_m

#		Calculate acceptance ratios:
		Accept_all[1].append(accept_m)
		Accept_all[0] += accept_m/Max_iter
		Vprint(3,"\nAcceptance ratio in step: %1.3f"%accept_m)



#	Evaluation pahse: discard the initial fraction of points:
	Points_final = Discard_init(Points_all, discard)

#	calculate final SS estimate:
	SS_est_final = Average_SS(Points_all)
#	print acceptance ratio:
	Vprint(3,"\nOverall acceptance ratio: %1.3f"%Accept_all[0])
#	update the SS probabilty distribution:
	SS_dist_final = Update_SSdist(Points_final, SS_dist)
#	calculat scaling factor distribution:
	Scale_dist_final = Update_SFdist(Points_final, Scale_dist)
	
	Final_SS = [SS_est_final, SS_dist_final, Scale_dist_final, Points_final, Accept_all]

	return Final_SS

#funtion to format output data:
def Format_SSdata(SS_output, Aux_data):
	SS_est = SS_output[0]
	SS_dist = SS_output[1]
	SF_dist = SS_output[2]
	dim =len(SS_est)
	Binpar = SS_dist[0][2]
	Files, Param = Aux_data
	max_iter, sample_size, scaling, discard, project,workdir = Param 
	Output = ""

#	header for output file:
	Header = "#Bayesian SS estimation module:\n"
	Header += "#Workdir:  %10s\n" % Param[5]
	Header += "#Spectrum file:       %10s\n" % Files[0]
	Header += "#Basis set file:      %10s\n" % Files[1]
	if Files[2] != "":
		Header += "#Parameter file:      %10s\n" % Files[2] 
	if Files[3] != "":
		Header += "#Correction file:     %10s\n" % Files[3]		
	Header += "#Initial sample size: %3d\n#Monte Carlo steps:   %3d\n"% (Param[1],Param[0])
	if Param[3] != 0:
		Header += "#Fraction of discarded samples:     %1.3f\n"%Param[3]
	if Param[2] == 0:
		Header += "#Spectrum scaling was disabled\n"
	
	Output += Header+"\n"


#	header for SS distribution:
	SS_map = "#Discretized SS dsitribution map:\n#  i      SS Class  \n"
	Hstring2 = "#bin: "
	for i in range(dim):
		Hstring = "#  %1d  : %12s \n" % (i, SS_est[i][0])
		SS_map += Hstring
		Hstring2 += "   %1d   " % i
	SS_map += Hstring2+"  P(bin)\n"

#	print SS probability map:
	Bstring = 3*"  "+ dim*" %1.3f "+" %1.3e\n"
#	sort bins based on probability:
	Bins_sorted = sorted(SS_dist[1],reverse= True, key= lambda Bin : Bin[dim])
#	Bins_sorted = SS_dist[1]
	for Bin in Bins_sorted:
		data = Convert_index(Bin, Binpar)
		SS_map += Bstring % tuple(data)
	Output += SS_map+"\n"

#	print SF distribution:
	SF_map = "\n#Scaling factor distribution:\n#  SFj    P(SFj)\n"
	for entry in SF_dist[1]:
		SF_map += "%6.3f   %2.3f\n"%tuple(entry)

	SF_map += "\n#Avg. Scaling factor : %1.4f +/- %1.4f\n" % tuple(SF_dist[2])
	Output += SF_map

#	format SS composition
	SS_comp = "#Weighting factors for the most likely SS composition (Avg +/- SD):\n"	
	for entry in SS_est:
		W_string = "# %12s  :   %1.4f   +/- %1.4f\n" % tuple(entry)
		SS_comp += W_string
	Output+= SS_comp+"\n"

	return Output	

#function to print the sorted set of points sampled:
def Format_Points(Points):
	Output = ""
	dim = len(Points[0][1])
#	sort points by RMSD:
	Points.sort(key=lambda Bin: Bin[5], reverse= True)
	Header = "#    Nj   SFj   Rmsd_j    Pri_j    Post_j      Coeff_j\n"
	Output += Header
	line_prev = ""
	for Pj in Points:
		pstring = " %6d   %1.3f  %1.3f  %1.3e  %1.3e  [ "%(Pj[0], Pj[2],  Pj[3], Pj[4], Pj[5])
		Vprint(5, Pj[1])
		for Coeff in Pj[1]:
			pstring += " %1.3f " % Coeff[1]	
		line = pstring + " ]\n"
#		add non-redundant lines:
		if line != line_prev:
			Output += line
		line_prev = line      	

	return Output

def Format_Dist1D(Points):
	Output = "# X     P(x)\n"
	for Pj in Points:
		Output += " %6.3f  %1.3e\n" % tuple(Pj)

	return Output

#function to transform SS matrix indeces into SS coefficients:
def Convert_index(SS_bin, Binpars):
	dim = len(SS_bin)-1
	new_bin = []
	for i in range(dim):
		Min,Max,Step = Binpars[i]
		Coeff = Min + (SS_bin[i]+0.5)*Step
		new_bin.append(Coeff)
	new_bin.append(SS_bin[dim])

	return new_bin
	

#Main function for script execution:
def SSbayes_Main(Args):
#	set run parameters:
	Input_files, Out_files, Param, failmark,verbosity = Args
	Set_verb(Args[4])
	spect_file, basis_file, par_file, corr_file = Input_files
	out_file, dist_file, data_file = Out_files
	Sesca_path, Error_Par0, SS_Par0, Run_Par = Param
	max_iter, sample_size, scaling, discard, Project_dim = Run_Par
	
#	import required modules:
	Vprint(2, "\nLoading required modules:")
	IMPORTS = Import_Custom()
	if "Pred" in IMPORTS:
		Spred.Set_verb(Args[4])
	else:
		Vprint(1,"Essential modules missing: SESCA_pred.py")
		failmark = 1
	if "Scale" in IMPORTS:
		Sscale.Set_verb(Args[4])
	else:
		Vprint(1,"Essential modules missing: SESCA_scale.py")

#	get default basis set list:
	BS_List = []
	if "Main" in IMPORTS:
		Smain.Set_verb(Args[4])
		Vprint(2, "\nDefault Basis sets:")
		BS_List = Smain.Pass_BSdefs()
	else:
		Vprint(2, "\nNo default Basis sets were found!")

#	execute main code:
	Main_Data = []

#	read basis set parameters:
	BS_Data = [basis_file, 0, [], []] 
	if basis_file != "":
		BS_Data = Read_BS_File(basis_file, BS_List)
		if BS_Data == []:
			failmark = 2
			bs_num = 0
			Vprint(1, "%1s\nNo Basis set parameters were read!"%basis_file)
		else:
			bs_num = BS_Data[1]
			Vprint(2, "%1d Basis spectra were read from basis set"%bs_num) 
	else:
		failmark = 2
		bs_num = 0
		Vprint(1, "%1s\nBasis set was not found!"%basis_file)

#	read CD spectrum:
	CD_spect = []
	if spect_file != "" and failmark == 0:
		L_range = ["",""]
		CD_spect = Spred.Read_Spectrum_file(spect_file,L_range)
		if CD_spect in ["None",[]]:
			Vprint(1,"Failed to read spectrum:",spect_file)
			failmark = 3
	elif failmark == 0:
		Vprint(1,"Error: no CD spectrum was found!")
		failmark = 3

#	read Correction spectrum:
	Corr_Spect = []
	if corr_file != "" and failmark == 0:
		L_range = ["",""]
		Corr_Spect = Spred.Read_Spectrum_file(corr_file,L_range)
		if Corr_Spect in ["None",[]]:
			Vprint(1,"Failed to read spectrum:",corr_file)
			failmark = 4
		
#	read in Bayesian parameters:
	if par_file == "auto" and failmark == 0:
#	use default files without and with side chain correctrions:
		if corr_file == "":
			par_file = os.path.join(SESCA_lib,"Bayes_2D-noSC.dat")
		else: 
			par_file = os.path.join(SESCA_lib,"Bayes_2D-SC.dat")
	if par_file != "" and failmark == 0:
#	if the user specified a parameter file, use that instead:
		Error_Par1, SS_Par1 = Read_BayesPar(par_file, Error_Par0, SS_Par0)
	else:
#	if no parameter file was specified, fall back to hard-coded defaults:
		Error_Par1, SS_Par1 = Error_Par0, SS_Par0

#	make sure the SS prior dimensions match the basis set size:
	ss_dim = 0
	if SS_Par1[0] != [] and failmark == 0:
		ss_dim = SS_Par1[0][1]	
	elif SS_Par1[0] == [] and SS_Par1[1] != [] and failmark == 0:
		Vprint(2, "\nDetermining boundaries for SS Prior")
		SS_Par1   = Complete_ErrDist(SS_Par1)
		ss_dim = SS_Par1[0][1]
	elif failmark == 0:
		Vprint(1, "No SS prior read from Bayesian SS parameter file")
		failmark = 5		

	if bs_num != ss_dim:
		Vprint(2, "Warning: the dimension of the SS prior (%1d) does not match the number SS classes in the basis set (%1d)"%(ss_dim,bs_num))
		if SS_Par1[0][0] in ["Uni_ND"]:
#			reset dimensions if the disrtibution is uniform		
			Vprint(2, "Setting SS dimensions to %1d"%bs_num)
			SS_Par1[0][1] = bs_num
		else:
#			stop otherwise
			Vprint(1,"Error, cannot reset SS dimensions, please set an approriate SS prior in the parameter file!")
			failmark = 6

#	generate prior SS probabilities:
	Vprint(3, "\nProcessing Prior SS distribution:")
	SS_dist0 = SS_Index(SS_Par1)
	if SS_dist0 == []:
		Vprint(1,"Error: SS prior cannot be processed!")
		failmark = 6

#	stop if input parameters are not OK:
	if failmark != 0:
		Vprint(1,Usage)
		sys.exit(failmark)

#	complete error parameters:
	Dev_Par, Scale_Par, Change_Par = Error_Par1
	Vprint(3, "\nProcessing Non-SS deviation distribution:")	
	Dev_Par   = Complete_ErrDist(Dev_Par)
	Vprint(3, "\nProcessing Scaling Factor distribution:")	
	Scale_Par = Complete_ErrDist(Scale_Par)
	Vprint(3, "\nProcessing SS step distribution:")	
	Change_Par = Complete_ErrDist(Change_Par)
	Error_Par = [Dev_Par, Scale_Par, Change_Par]


#	get projection parameters for heat maps:
	if SS_Par1[2] != [] and Project_dim == []:
		print(SS_Par1[2])
		Project_dim = [(SS_Par1[2][1][0])-1, (SS_Par1[2][1][1])-1]
	elif Project_dim == []:
		Project_dim = [0,1]
		
	Vprint(3, "Projecting coordinates to classes: %1d, %1d"%(Project_dim[0]+1, Project_dim[1]+1))
		
#	Iteratively refine the SS distribution to match the CD spectrum:
	Run_param = [max_iter, sample_size, scaling, discard, Error_Par]
	Final_SS = Iterate_SS_dist(CD_spect, Corr_Spect, SS_dist0, BS_Data, Run_param)

#	write out results:
	Input_final = [spect_file, BS_Data[0], par_file, corr_file]
	Run_final = [max_iter, sample_size, scaling, discard, Project_dim, workdir]
	Aux_data = [Input_final, Run_final]
	Main_Data = Final_SS
	if out_file != "":
		Output = Format_SSdata(Final_SS, Aux_data)
		o = open(out_file, "wb")
		o.write(Output.encode("ascii"))
		o.close()

#	generating reduced coordinate heatmap:
	if dist_file != "" and Project_dim != []:
#		project probabilty density onto selected dimensions:		
		Projected_2D = Project_ND(Final_SS[1],Project_dim)

#		determine projection labels:
		Labels = []
		for dim in Project_dim:
			label = Final_SS[0][dim][0]
			Labels.append(label)
		Aux_data.append(Labels)

#		generate and write heatmap data:	
		Heatmap = Heatmap_Data(Projected_2D, Aux_data)
		h = open(dist_file,"wb")
		h.write(Heatmap.encode("ascii"))
		h.close() 

	if data_file != "":
		Point_output = Format_Points(Final_SS[3])
		d = open(data_file,"wb")
		d.write(Point_output.encode("ascii"))
		d.close()			

	return Main_Data
	
# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Set_verb(Custom_Args[4])
	Vprint(4, "\nRun parameters:\n", Custom_Args)
	Input_files, Output_files, Param, failmark, verbosity = Custom_Args
	out_file = Output_files[0]

#       executing main code
	Data_main = SSbayes_Main(Custom_Args)
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
	Vprint(2, "SESCA Bayesian SS estimation module (SESCA_bayes.py)")

