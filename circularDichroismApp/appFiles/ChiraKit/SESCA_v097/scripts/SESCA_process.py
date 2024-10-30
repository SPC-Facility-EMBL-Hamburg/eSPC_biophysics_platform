#!/usr/bin/env python
#

import sys
import os
import math
import time

workdir = os.getcwd()
stime = time.time()


#print workdir+"\n"

FLAGS = ["list","Clist","spect","corr","Uinp","Uout","col","wl","conc","Mconc","res","mr","path","write","data","verb"]
usage0  = "*******************************************************************\n"
usage   = "SESCA_module for processing and reformatting CD spectra\n"
usage  += "Usage:\n\nSESCA_process.py @flag <argument> (@flag <argnument> ...)\n\n"
usage  += "Possible command flags are:\n"
usage  += "   @spect <file> (<file> ...)  provide input spectra to process (default: None)\n"
usage  += "   @corr  <file> (<file> ...)  provide baseline spectra to process (default: None)\n"
usage  += "   @list  <file> provide file list for spectrum files (default: None)  \n"
usage  += "   @Clist <file> provide file list for baseline files (default: None)\n"
usage  += "   @write <0,file> set output file name (default: CD_processed.dat)\n"
usage  += "   @data  <0,file> set filename for individual spectra, (default: None)\n"
usage  += "   @Uinp  <str>  set input  unit (mdeg, dAbs, MRE, kMRE, dEps, default: mdeg)\n"
usage  += "   @Uout  <str>  set output unit (mdeg, dAbs, MRE, kMRE, dEps, default: mdeg)\n"
usage  += "   @col   <int,int> set wavelength and intensity columns in input (default: 1,2)\n"
usage  += "   @wl    <float,float,float> set wavelength range (start,end,step in nm, default: auto)\n"
usage  += "   @conc  <float> set protein concentration for unit conversion (in mol/L, default: 0.0) \n"
usage  += "   @path  <float> set cuvette pathlength for unit conversion (in cm, default: 1.0)\n"
usage  += "   @res   <int>   set the number of residues for unit conversion (default: 1)\n"
usage  += "   @Mconc <float> set protein mass concentraion (in g/L, requires Mr, default: 0.0)\n"
usage  += "   @mr    <float> set protein molecular mass (in g/mol default: 1.0) \n"
usage  += "   @verb  <int>   set verbosity level(0 - silent, 5 - debug mode, default: 3)\n"
Usage = usage0 + usage + usage0

#default parameters:
#input and output files:
list_file = ""
corr_list = ""
SPECT_LIST = []
CORR_LIST = []
out_file = "CD_processed.dat"
dat_file = ""
Separator = ""
Decimal = ","

#spectrum parameters:
in_format = "mdeg"
out_format = "mdeg"
CD_Units = ["mdeg","dAbs","MRE","kMRE","dEps"]
Cols = [1,2]
Waves = ["","",""]
CD_Par = [in_format, out_format, Cols, Waves]

#Sample parameters:
Conc0  = 0.0
ResNum = 1
path = 1.0
ConcM = 0.0
MolW = 1.0
Ext = []
Sample_Par = [Conc0, path, ResNum, ConcM, MolW] 



failmark = 0
verbosity = 2

Input_files = [list_file, corr_list]
Output_files = [out_file, dat_file]
Param = [CD_Par, Sample_Par]

Def_Args = [Input_files, Output_files, Param, failmark, verbosity]


#function definitions:
#function to pass on defaults:
def Pass_Defaults():
        return Def_Args

#function to control verbosity:
def Vprint(level, *Messages):
	global verbosity
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
	try:
		Vprint(2, "Loading Modules:")	
		globals()["np"] = __import__("numpy")
		IMPORTS.append("numpy")

	except Exception:
		Vprint(1,"Cannot load requested modules")
	
	return IMPORTS


#function to read in arguments:
def Read_Args(Args, Defs= Def_Args, Flags= FLAGS, Units= CD_Units):
	argnum = len(Args)
	global SPECT_LIST
	global CORR_LIST
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
	New_Args = Defs
	flag = ""
	acnt = 0
#	processing new passed arguments: 
#FLAGS = ["list","Clist","inp","corr","Uinp","Uout","col","wl","conc","Mconc","res","mr","path","write","verb"]
	Vprint(2, "Recognized flags:")
	for arg in Args:
		if arg.startswith("@"):
			flag = arg.strip("@")
			if flag in Flags:
				Vprint(2, flag)
			else:
				Vprint(1,"Unknown flag:",flag)
				New_Args[3] = 1
#		handle I/O parameters:
		elif flag == "list":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "Clist":
			New_Args[0][1] = arg
			flag = ""	
		elif flag == "spect":
			SPECT_LIST.append(arg)
		elif flag == "corr":
			CORR_LIST.append(arg)
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "data":
			if arg == "0":
                                New_Args[1][1] = ""
			else:
				New_Args[1][1] = arg
			flag = ""
#		handle spectrum parameters:
		elif flag == "Uinp":
			if arg in Units:
				New_Args[2][0][0] = arg
			else:
				Vprint(1,"Unrecognized input format!\nAcceptable formats are: %1s"%str(Units))
				New_Args[3] = 1	
			flag = ""
		elif flag == "Uout":
			if arg in Units:
				New_Args[2][0][1] = arg
			else:
				Vprint(1,"Unrecognized input format!\nAcceptable formats are: %1s"%str(Units))
				New_Args[3] = 1	
			flag = ""
		elif flag == "col":
			try:
				parts = arg.split(",")
				X,Y = int(parts[0]), int(parts[1])
				New_Args[2][0][2] = [X,Y]
				if len(parts) > 2:
					E = int(parts[2])
					New_Args[2][0][2] = [X,Y,E]
			except Exception:
				Vprint(1,"@cols takes comma separated integers (for wavelengths and intensities)")
				New_Args[3] = 1
			flag = ""
		elif flag == "wl":
			try:
				parts = arg.split(",")
				Min,Max,Step = float(parts[0]), float(parts[1]), float(parts[2])
				New_Args[2][0][3] = [Min, Max, Step]
			except Exception:
				Vprint(1,"@wl takes comma separated foats (Wli,Wlf,dWl)")
				New_Args[3] = 1
			flag = ""
#		handle sample parameters:
		elif flag == "conc":
			try:
				print(arg)
				New_Args[2][1][0] = float(arg)
			except ValueError:
				Vprint(1, "@conc takes only floats as argument")
				New_Args[3] = 1
			flag = ""
		elif flag == "path":
			try:
				New_Args[2][1][1] = float(arg)
			except ValueError:
				Vprint(1, "@path takes only floats as argument")
				New_Args[3] = 1
			flag = ""
		elif flag == "res":
			try:
				New_Args[2][1][2] = int(arg)
			except ValueError:
				Vprint(1, "@res takes only floats as argument")
				New_Args[3] = 1
			flag = ""
		elif flag == "Mconc":
			try:
				New_Args[2][1][3] = float(arg)
			except ValueError:
				Vprint(1, "@Mconc takes only floats as argument")
				New_Args[3] = 1
			flag = ""
		elif flag == "mr":
			try:
				New_Args[2][1][4] = float(arg)
			except ValueError:
				Vprint(1, "@mr takes only floats as argument")
				New_Args[3] = 1
			flag = ""
		elif flag == "verb":
			New_Args[4] = int(arg)
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0][0] == "" and acnt == 0:
			SPECT_LIST.append(arg)
			flag = ""
		elif flag == "" and New_Args[1][0] == "" and acnt == 1:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args

#function read to file lists:
def Read_List(File):
	In_Files = []
	Dir = "./"
	
#	make sure file is there:
	if os.path.isfile(File) == False:
		Vprint(1, "Cannot file input file: %1s"%File)
		return [In_files, Frames]

	Vprint(2, "Reading list file:")
	f= open(File, "rb")
	flag = "traj"
	for line0 in f:
		line = str(line0.decode("ascii")).strip("\n")
#		read in source directory:
		if line.startswith("#Dir="):
			Dir = line.split()[1]
			flag = "traj"
			Vprint(2, "Dir: %1s"%Dir)
#		ignore comments
		elif line.startswith("#"):
			pass
#		collect files, join and normalize paths 
		elif flag == "traj":
			parts = line.split()
			for ex in parts:
				path = os.path.join(Dir, ex)
				npath = os.path.normpath(path)	
				In_Files.append(npath)
				Vprint(3, npath)
		elif flag != "" and line == "":
			flag = ""
	
	return  In_Files

#function to read input spectra:
def Read_Spectra(File= "", Cols= Cols, decimal= Decimal, Sep= Separator):
	Data = []

#	check file
	if File == "" or os.path.isfile(File) == False:
		Vprint(2, "Cannot locate file: %1s"%File)
		return Data

#	set data dimension:
	d_num = len(Cols)
	for i in Cols:
		Data.append([])

#	open file:
	Vprint(2, "\nReading File: %1s"%File)
	f = open(File,"rb")
	for line0 in f:
		line = str(line0.decode("ascii")).strip('\n')
		if line != "" and not line[0] in ["#","&",";"]:
			Vprint(4, "",str(line))
#			separate data by the given separator (white space by default)
			if Sep == "":
				Pj = line.split()
			else:
				Pj = line.split(Sep)

#			try to read out data in the appropriate columns
			try:
#			if 1 == 1:
				Temp = []
				for i in range(d_num):
					ndx = Cols[i]-1
#					if the decimal is a dot, just convert number
					if decimal == ".":
						value = float(Pj[ndx])
#					if the decimal is not a dot, replace decimal with dot
					else:
						X_string = ""
						Vprint(5,Pj, ndx, Cols)
						X_val = Pj[ndx].split(decimal)
						for X in X_val:
							X_string += X+"."
						value = float(X_string[:-1])
					Temp.append(value)

#				if we collected all data, add them to data array:
				Vprint(3, Temp)
				for i in range(d_num):
					Data[i].append(Temp[i])

#			skip line if we ran into trouble:
			except Exception:
				Vprint(3, "","Skipped line: %1s"%line)
		
		
	return Data	

#helper function to calculate statistics:
def Get_Stat(Data):
	Mean, St_dev, Pmin, Pmax = [0.0, 0.0, "", ""]
	d_num = len(Data)

#	go over all points:
	if d_num != 0:
		for Py in Data:
#			calculate mean
			Mean += float(Py)/d_num
#			track minimum and maximum
			if Pmin	== "" or Py < Pmin:
				Pmin = Py
			if Pmax == "" or Py > Pmax:
				Pmax = Py

#		compute standard deviation:
		for Py in Data:
			St_dev += (Py-Mean)**2 / d_num
		St_dev = math.sqrt(St_dev)
			
	return [Mean, St_dev, Pmin, Pmax]	

def Get_WL_Range(data= [], Wl_Par= Waves):
	Wave_final = []
	d_num = len(data)
	Vprint(2,"\nMerging wavelength information:")

	Wave0 = []
	Lmin = 0.0
	Lmax = 0.0
#	determine spectral range :
	if Wl_Par != ["","",""]:
		Wave0 = [Wl_Par[0], Wl_Par[1], Wl_Par[2]]	
		if Wave0[0] <= Wave0[1] and Wave0[2] != 0.0:
			Wave0[2] = math.fabs(Wl_Par[2])	
			Lmin, Lmax = Wave0[0], Wave0[1]
		elif Wave0[0] > Wave0[1] and Wave0[2] != 0.0:
			Wave0[2] = -1.0*math.fabs(Wl_Par[2])
			Lmin, Lmax = Wave0[1], Wave0[0]
		else:
			Vprint(2,"Warning, improper wavelength parameters: %1s"%Wl_Par)
			return Wave_final
		Vprint(2,"Requested wavelength range: %2.1f - %2.1f (step %1.2f)"%tuple(Wave0))

#		check if requested wavelength match the data:
		for i in range(d_num):
#			get wavelength statistics:
			Stat_i = Get_Stat(data[i][0])
			Min_i, Max_i = Stat_i[2], Stat_i[3]
			Step_i = data[i][0][1] - data[i][0][0]
			w_num = len(data[i][0])
			Sum_i = [i, Min_i, Max_i,Step_i, w_num]
			if Min_i > Lmin:
				Lmin = Min_i
			if Max_i < Lmax:
				Lmax = Max_i
		if Wave0[2] > 0:
			Wave0[0], Wave0[1] = Lmin, Lmax
		if Wave0[2] < 0:
			Wave0[0], Wave0[1] = Lmax, Lmin
		Vprint(2,"Possible wavelength range: %2.1f - %2.1f (step %1.2f)"%tuple(Wave0))
			

	else:
		Wave0 = ["","",""]
		for i in range(d_num):
#			get wavelength statistics:
			Stat_i = Get_Stat(data[i][0])
			Min_i, Max_i = Stat_i[2], Stat_i[3]
			Step_i = data[i][0][1] - data[i][0][0]
			w_num = len(data[i][0])
			Sum_i = [i, Min_i, Max_i,Step_i, w_num]
			Vprint(4, "Set statistics: %1s"%Sum_i)
			if Wave0[0] == "" or Wave0[0] < Min_i:
				Wave0[0] = Min_i
			if Wave0[1] == "" or Wave0[1] > Max_i:
				Wave0[1] = Max_i
			if Wave0[2] == "" or Wave0[2] > Step_i:
				Wave0[2] = Step_i
#			exchange Min/Max if step is negative:
		if Wave0[2] < 0.0:
			Wave0 = [Wave0[1], Wave0[0], Wave0[2]]
		Vprint(2,"Determined wavelength range %2.1f - %2.1f (step %1.2f)"%tuple(Wave0))	
#	generate wavelengths:
	Wl = Wave0[0]
	if Wave0[2] > 0.0:
		while Wl <= Wave0[1]:
			Wl = round(Wl,3)
			Wave_final.append(Wl)
			Wl += Wave0[2]	
	else:
		while Wl >= Wave0[1]:
			Wl = round(Wl,3)
			Wave_final.append(Wl)
			Wl += Wave0[2]
	Vprint(4, Wave_final)	
	
	return Wave_final

#helper function to interpolate points:
def Interpolate_Y(Xi= [], Yi= [], X0= 0.0):
	Y0 = 0.0
	lower = ["",0.0,0.0]
	upper = ["",0.0,0.0]
	P0 = [X0, Y0, 0]
	p_num = len(Xi)
#	cycle over the Points:
	for j in range(p_num):
		X_dist = math.fabs(Xi[j]-X0)
#		return value if we have an exact match:
		if X0 == Xi[j]:
			Y0 = Yi[j]
			P0 = [X0, Y0, 1]
			Vprint(5, "match found! Xi= %2.1f, Yi= %1.3f"%(Xi[j],Yi[j]))
			break 
#		find lower neightbor:
		if X0 > Xi[j] and (lower[0] == "" or lower[0] > X_dist):
			lower = [X_dist, Xi[j], Yi[j]]
#		find upper neightbor:
		if X0 < Xi[j] and (upper[0] == "" or upper[0] > X_dist):
			upper = [X_dist, Xi[j], Yi[j]]

#	interpolate if no match was found:
	if P0[2] == 0 and lower[0] != "" and upper[0] != "":
		div = math.fabs(upper[1]-lower[1])
		Wup, Wlo = lower[0]/div, upper[0]/div
		Y0 = Wup*upper[2] + Wlo*lower[2]
		P0 = [X0, Y0, 2]
		Vprint(5, "(%2.1f, %1.3f)value interpolated from:"%(X0,Y0))
		Vprint(5, "   Plo(%2.1f, %1.3f)"%(lower[1],lower[2]))
		Vprint(5, "   Pup(%2.1f, %1.3f)"%(upper[1],upper[2]))
	
	if P0[2] == 0:
		Vprint(4, "Could not interpolate value for: %2.1f"%X0)	
	

	return P0

def Merge_Spectra(data= [], wave= []):
	Spectra_Reg= []
	d_num = len(data)
	w_num = len(wave)
	Vprint(2, "\nProcessing Spectra:")
	
#	compute CD intensities for each wavelength:
	Int_mean = []
	Int_sdev = []
	Int_all = []
	for l in range(w_num):
		Wave_l = wave[l]
		Temp = [Wave_l, 0.0, 0.0, [], []]
#		interpolate spectra:
		for j in range(d_num):
			X_j,Y_j = data[j][0], data[j][1]
			Int_jl = Interpolate_Y(Xi= X_j, Yi= Y_j, X0= Wave_l)
			Vprint(4,"   ",Int_jl)
			Temp[3].append(Int_jl[1])
			if len(data[j]) > 2:
				E_j = data[j][2]
				Err_jl = Interpolate_Y(Xi= X_j, Yi= E_j, X0= Wave_l)
				Temp[4].append(Err_jl[1])
#		compute intensity statistics:
		Stats_l = Get_Stat(Temp[3])		
		Temp[1] = round(Stats_l[0],3)
		Temp[2] = round(Stats_l[1],3)
#		add uncertainty from spectra:
		if Temp[4] != []:
			Vprint(5,"averaging spectrum uncertainty:")
			Vprint(5, Temp[4])
			Estats_l = Get_Stat(Temp[4])
			Err_av = Estats_l[0]
			Err_fin = math.sqrt((Temp[2]+Err_av)**2)
			Vprint(5,Temp[2],Err_av,Err_fin)
			Temp[2] = Err_fin
#		save results:
		Vprint(3,Temp)
		Int_mean.append(Temp[1])	
		Int_sdev.append(Temp[2])
		Int_all.append(Temp)	
				
	Spectra_Reg = [Int_mean, Int_sdev, Int_all]

	return Spectra_Reg

#function to manipulate spectrum intensities:
def Modify_Spectra(tar= [], ref= [], op= "subtract"):
	Modified = []
	t_num = len(tar)
	r_num = len(ref)
	if op in ["add","subtract","rms"]:
		if t_num != r_num:
			Vprint(2,"Mismatching data length, cannot perform operation %1s"%op)
			return Modified
		else:
			Vprint(2,"Spectrum operation: %1s"%op)
			for i in range(t_num):
				if op == "add":
					value_i = round(tar[i]+ref[i],3)
				elif op == "subtract":
					value_i = round(tar[i]-ref[i],3)
				elif op == "rms":
					value_i = round(math.sqrt((tar[i]+ref[i])**2),3)
				Modified.append(value_i)
	if op in ["multiply","divide"]:	
		if r_num != 1:
			Vpirnt("Unexpected parameters for operation %1s"%op)
		else:
			Vprint(2, "Spectrum operation: %1s"%op)
			for i in range(t_num):
				if op == "multiply":
					value_i = ref[0]*tar[i]
				if op == "divide" and ref[0] != 0.0:
					value_i = float(tar[i])/ref[0]
				Modified.append(value_i)
	Vprint(4, Modified)
	return Modified

def Extract_Processed(data= [], base= []):
	Ints = []
	Labels = ["Individual input spectra:\n#","wvlgth"]
	w_num = len(data)
	d_num = len(data[0][3])
	b_num = 0
	if base != []:
		b_num = len(base[0][3])

#	prepare data arrays and labels:
	cnt = 0
	for i in range(d_num):
		cnt += 1
		Labels.append("Int%1d "%cnt)
	for i in range(b_num):
		cnt += 1
		Labels.append("Base%1d"%cnt)

#	collect wavelengths:
	Wl = []
	for j in range(w_num):
		Wl.append(data[j][0])
	Ints.append(Wl)
		
#	collect intensities:
	for i in range(d_num):
		Temp_i = []
		for j in range(w_num):
			value_ji = data[j][3][i]
			Temp_i.append(value_ji)	
		Ints.append(Temp_i)
#	collect baseline intensities:
	if base != []:
		for i in range(b_num):
			Temp_i = []
			for j in range(w_num):
				value_ji = base[j][3][i]
				Temp_i.append(value_ji)
			Ints.append(Temp_i)

	return [Ints, Labels]

#function to calculate CD unit conversion matrix based on sample parameters:
def Compute_ConvMatrix(Cmol= Conc0, plength= path, Rnum= ResNum, Cmass= ConcM, Mw= MolW):
	#we will convert units exactly in this order: 
	CD_Units = ["mdeg","dAbs","MRE","kMRE","dEps"]
	ConvMat = []
	Index = {}
#	conversion between dESP and theta:
	T2E = 3298.2

#	determine concentraion:
	Cprot = 0.0
	if Cmol != 0.0:
		Cprot = Cmol
	elif Cmass != 0.0 and Mw != 0.0:
		Cprot = Cmass/Mw
	else:
		Cprot = 1.0

#	add generate index dictionary:
	u_num = len(CD_Units)
	for i in range(u_num):
		unit_i = CD_Units[i]
		Index[unit_i] = i

#	Generate Matrix:
	for i in range(u_num):
		Row = []
		for i in range(u_num):
			Row.append(1.0)
		ConvMat.append(Row)

#	Fill up conversion matrix:
#	mdeg to dAbs:
	ConvMat[0][1] = 1.0/(T2E * 10.0)
	ConvMat[1][0] = (T2E * 10.0)
#	mdeg to MRE:
	ConvMat[0][2] = 0.1/(Cprot*plength*Rnum)
	ConvMat[2][0] = 1.0/ConvMat[0][2]
#	mdeg to kMRE:
	ConvMat[0][3] = ConvMat[0][2]/1000.0		
	ConvMat[3][0] = 1.0/ConvMat[0][3]
#	mdeg to dEPS:
	ConvMat[0][4] = ConvMat[0][2]/T2E
	ConvMat[4][0] = 1.0/ConvMat[0][4]
#	dAbs to dEPS:
	ConvMat[1][4] = 1.0/(Cprot*plength*Rnum)
	ConvMat[4][1] = 1.0/ConvMat[1][4]
#	dAbs to MRE:
	ConvMat[1][2] = ConvMat[1][4]*T2E
	ConvMat[2][1] = 1.0/ConvMat[1][2]
#	dAbs to kMRE:
	ConvMat[1][3] = ConvMat[1][2]/1000.0
	ConvMat[3][1] = 1.0/ConvMat[1][3]
#	MRE to kMRE:
	ConvMat[2][3] = 0.001
	ConvMat[3][2] = 1.0/ConvMat[2][3]
#	MRE to dEPS:
	ConvMat[2][4] = 1.0/T2E
	ConvMat[4][2] = 1.0/ConvMat[2][4]
#	kMRE to dEPS:
	ConvMat[3][4] = 1000.0/T2E	
	ConvMat[4][3] = 1.0/ConvMat[3][4]

#	print matrix:
	Vprint(5,"CD unit conversion Matrix:")
	Vprint(5,CD_Units)
	string = u_num*" %1.3e "
	for Row in ConvMat:
		Vprint(5,string%tuple(Row))

	return [ConvMat, Index]

#output formatting:
def Format_Output(data = [],labels= [], aux= []):
	Output = ""
	a_num = len(aux)

#	add file header
	Header = "#Preprocessed CD spectrum:\n"
	Header += "#source: %1s\n"%os.getcwd()
	if a_num >= 3:
		Uin, Uout = aux[2][0], aux[2][1]
		Cols = aux[2][2]
		Header += "#data columns:  %1s\n"%Cols
		Header += "#input format:  %1s\n"%Uin
		if Uin != Uout:
			Header += "#output format: %1s\n"%Uout
			Cmol, pathl, resN, Cmass, Mr = aux[3]
			if Cmol != 0.0:
				Cprot = Cmol
			elif Cmass != 0.0 and Mr != 0.0:
				Cprot = Cmass/Mr
			else:
				Cprot = 0.0
			Header += "#protein conc(M): %1.2e\n"%Cprot
			Header += "#residue number:  %1d\n"%resN
			Header += "#pathlength(cm):  %1.3f\n"%pathl

	if a_num >= 1:
		Header += "#spectrum files:\n"
		for File in aux[0]:
			Header += "#%1s\n"%File
	if a_num >= 2:
		Header += "#basline files:\n"
		for File in aux[1]:
			Header += "#%1s\n"%File	
	Output += Header+"\n"

#	add data block
	Data_Block = ""
#	add column labels:
	if labels != []:
		Data_Block += "# "
		for Label in labels:
			Data_Block += "%8s "%Label
		Data_Block += "\n"
#	add data row by row:
	d_num = len(data)
	w_num = len(data[0])
	Vprint(4, "dimensions: %1d,%1d"%(d_num,w_num))
	digit = 3
	if digit < 5:
		sep = 5-digit
	else:
		sep = 0
	var_string = "  %6."+"%1df"%digit+sep*" "
	for i in range(w_num):
		Row_i = "   %6.1f " % data[0][i]
		for j in range(1,d_num):
			Vprint(5, j,i)
#			Row_i += "  %6.3f  "%data[j][i]
			Row_i += var_string%data[j][i]
		Data_Block += Row_i+"\n"
	Output += Data_Block	
		
	return Output

# Main function for script execution:
def Script_Main(Args):
#	set run parameters:
	Input_files, Out_files, Param, failmark,verbosity = Args
	Set_verb(Args[4])
	list_file, list_corr = Input_files
	out_file, data_file = Out_files
	CD_Par, Sample_Par = Param

#	read input files:
	global SPECT_LIST
	global CORR_LIST

#	read in spectrum files:
	if list_file != "" and failmark == 0:
		Spectra = Read_List(list_file)
		for S in Spectra:
			SPECT_LIST.append(S)
#	read in basline files:
	if list_corr != "" and failmark == 0:
		Baselines = Read_List(list_corr)
		for B in Baselines:
			CORR_LIST.append(B)

#	check if we have input to work with:	
	if failmark == 0 and SPECT_LIST == []:
		Vprint(1, "\nNo input files!")
		failmark = 2

#	check if concentrations are provided for unit conversion:
	conv_flag = False
	if CD_Par[0] in ["mdeg","dAbs"] and CD_Par[1] in ["MRE","kMRE","dEps"]:
			conv_flag = True
	if CD_Par[1] in ["mdeg","dAbs"] and CD_Par[0] in ["MRE","kMRE","dEps"]:
			conv_flag = True
	if conv_flag == True and Sample_Par[0] == 0.0 and Sample_Par[3] == 0.0:
		Vprint(1, "\nNo concentraion was provided for unit conversion!")
		Vprint(1, "Please set sample parameters for unit conversion.")
		failmark = 3	

#	check input data:
	if failmark != 0:
		Vprint(1, Usage)
		sys.exit(failmark)

#	execute main code:	
	Main_Data = []
	spect_num = len(SPECT_LIST)
	base_num = len(CORR_LIST)

#	extract spectra:
	Spectra_Raw = []
	Vprint(2, "\nReading %1d spectrum files:"%spect_num) 
	for File in SPECT_LIST:
		Spectra_i = Read_Spectra(File= File, Cols= CD_Par[2])
		if Spectra_i != []:
			Spectra_Raw.append(Spectra_i)

	Vprint(2, "\nReading %1d baseline files:"%base_num)
	Baselines_Raw= []
	for File in CORR_LIST:
		Baseline_i = Read_Spectra(File, Cols= CD_Par[2])
		Baselines_Raw.append(Baseline_i)

#	process and average spectra:
	Wavelengths = []
	Spectra_full = []
	Sdev_full = []
	Processed_all = []
	if Spectra_Raw != []:
		Vprint(2,"\nAnalysing CD spectra:")
		Wavelengths = Get_WL_Range(Spectra_Raw, CD_Par[3])
		Spectra_processed = Merge_Spectra(Spectra_Raw, Wavelengths)
	else:
		Vprint(1, "\nNo wavelengths found in the provided spectra!")
		failmark = 3
		Vprint(1, Usage)
		sys.exit(failmark)

	if Baselines_Raw != []:
		Vprint(2,"\nAnalysing Basline spectra:")
		Base_processed = Merge_Spectra(Baselines_Raw, Wavelengths)
		Vprint(2,"\nPerforming Baseline corrections:")
		Spectra_full = Modify_Spectra(tar= Spectra_processed[0], ref= Base_processed[0], op= "subtract")
		Sdev_full = Modify_Spectra(tar= Spectra_processed[1], ref= Base_processed[1], op= "rms")
		Processed_all = Extract_Processed(Spectra_processed[2], Base_processed[2])
	elif Baselines_Raw == [] and Spectra_Raw != []:
		Spectra_full = Spectra_processed[0]
		Sdev_full = Spectra_processed[1]
		Processed_all = Extract_Processed(Spectra_processed[2])

#	convert spectra to a different unit if necessary:
	if CD_Par[0] != CD_Par[1]:
		Vprint(2, "\nConverting processed spectra from %1s to %1s units"%(CD_Par[0], CD_Par[1]))
		ConvMat, Uindex = Compute_ConvMatrix(*Sample_Par)
		ndx_i = Uindex[CD_Par[0]]
		ndx_o = Uindex[CD_Par[1]]
		Conv_factor = ConvMat[ndx_i][ndx_o]
		Vprint(3, "Conversion factor: %1.2e"%Conv_factor)
		Spectra_full = Modify_Spectra(tar=Spectra_full, ref= [Conv_factor], op= "multiply")
		Sdev_full = Modify_Spectra(tar=Sdev_full, ref= [Conv_factor], op= "multiply")
 

	Main_Data = [Wavelengths, Spectra_full, Sdev_full]

#	write output data:
	if out_file != "":
		Labels_main = ["wvlgth", "Int_avg","Int_sd"]
		Aux_main = [SPECT_LIST, CORR_LIST, CD_Par, Sample_Par]
		Output = Format_Output(data= Main_Data, labels= Labels_main, aux= Aux_main)
		o=open(out_file,"wb")
		o.write(Output.encode("ascii"))
		o.close()
#	write data file if requested
	if data_file != "" and Processed_all != []:
		Data_Out = Format_Output(data= Processed_all[0], labels= Processed_all[1], aux= Aux_main)
		o = open(data_file,"wb")
		o.write(Data_Out.encode("ascii"))
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
	Set_verb(Custom_Args[4])
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
	outfiles = outfiles[:-1]	
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",outfiles)
else:
	print("SESCA spectrum processing module")

