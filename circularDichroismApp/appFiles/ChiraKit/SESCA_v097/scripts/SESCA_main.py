#!/usr/bin/env python
#

import sys
import os
import math
import time

stime = time.time()
workdir = os.getcwd()
#print workdir

#SESCA main module to predict analyse the CD spectra and secondary structure of proteins:
usage0  = "*******************************************************************\n"
usage   = "Main SESCA module for automatic CD spectrum prediction and model validation\n"
usage  += "usage: SESCA_main.py <target_file> (<reference_file>)  <basis_set> @flag <argument>\nPossible command flags are:\n"
usage  += "   @ref  <ref_file> specify reference file (CD spectrum)\n   @pdb <target_file> specify structure file (single structure or trajectory in PDB format)\n"
usage  += "   @BB_file  <SS_file> specify secondary structure file (and skip SS preprocessing of the pdb file)\n"
usage  += "   @SC_file  <SC_file> specify sequence summary file (and skip sequence preprocessing by the SESCA_seq module, necessary for mixed basis set only)\n"
usage  += "   @lib  <BS_file> specify basis spectrum library (default is the DS-dT basis set, see @lib help for custom basis set options)\n"
usage  += "   @write <output_file> specify output file name (default: CD_comp.out)\n"
usage  += "   @range <float,float> limit wavelength range to work in (default: none)\n"
usage  += "   @scale <0/float>  use scaling factor for the calculated CD spectrum, 0 - no scaling (default: 1.0)\n"
usage  += "   @refscale <0/float/auto> use scaling factor for the reference CD spectrum, 0 - no scaling, auto - fit ref. intensity to match calculated spectrum (default: 0)\n"
usage  += "   @norm <0,1 / float> normalize structure composition to 100%, if a float other than 1.0 is provided, the calculated CD spectrum is scaled by that amount (default: 0 - off)\n"
usage  += '   @err <int> select calibration curve for model error estimation (found in the basis set), 0 - no error estimation (default: "auto")'
usage  += '   @prep "<string>" provide custom arguments for the structure preprocessor (depends on the basis set)\n'
usage  += '   @main "<string>" provide custom arguments directly to the main CD module (for development only)\n'
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"

Usage = usage0 + usage + usage0



#Specify SESCA files and directories here:
#################################################
#1- use windows mode, 0- use unix mode
win = 0
SESCA_Dir =     "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097"
#SESCA_Dir =      "D:\\tan\Gottingen\Programs\SESCA_dev"
#SESCA_Dir =      "/home/gnagy/Programs/SESCA_dev"
SESCA_scripts =  os.path.join(SESCA_Dir,"scripts")
SESCA_lib =	 os.path.join(SESCA_Dir,"libs")
SESCA_seq  =	 os.path.join(SESCA_scripts,"SESCA_seq.py")
SESCA_defBS =	 os.path.join(SESCA_lib,"Basis_sets.nfo")
DISICL_Dir =     os.path.join(SESCA_Dir,"DISICL")
DISICL_main =	 os.path.join(DISICL_Dir,"DISICL_main.py")
DSSP_ref  = 	 os.path.join(SESCA_scripts,"SESCA_dssp.py")
DSSP_Dir =       os.path.join(SESCA_Dir,"DSSP")
HBSS_Dir =	 os.path.join(SESCA_Dir,"HBSS")
HBSS_main =      os.path.join(HBSS_Dir,"HBSS_main.py")
if win == 0:
	DSSP_main = 	 os.path.join(DSSP_Dir,"dssp-2.0.4-linux-amd64")
elif win == 1:
	DSSP_main = 	 os.path.join(DSSP_Dir,"dssp-2.0.4-win32.exe")


##################################################
#the code from here should not be changed

Libhelp =  "\nDefault basis set library in: %1s\n" % SESCA_lib
Libhelp	+= "Defaults read from: %1s\n" % SESCA_defBS
Libhelp += "to specify custom basis sets, edit the default file or use the flags:\n"
Libhelp += '@lib "custom", in conjuction with\n@method <SS_method> to specify structure preprocessor,'
Libhelp += " valid options: (DS_det, DS_sim ,Dssp, Hbss, Seq)\n" 
Libhelp += '@BB_lib <BS_file> to specify backbone basis set file\n@SC_lib <BS_file/"None"> to specify side chain library (based on Seq by default)\n'


#defining default parameters:
pdb_file = ""
BB_file = ""
SC_file = ""
ref_file = ""
Inp_files = [pdb_file, BB_file, SC_file, ref_file]
lib_name = "DS-dT"
SS_method = ""
BB_lib = ""
SC_lib = ""
BS_data = [lib_name, SS_method, BB_lib, SC_lib]
out_file = ""
dec_file = ""
Out_files = [out_file, dec_file]
L_range = ["",""]
scale = 1.0
exp_scale = ""
err_curve = "auto"
main_args = []
prep_args = []
Param = [L_range, scale, exp_scale, main_args, prep_args, err_curve]
failmark = 0
verbosity = 1

#Def_Main = [Inp_files, BS_data, out_file, deconv, main_args, prep_args, L_range, scale, exp_scale, failmark, verbosity]
Def_Main = [Inp_files, BS_data, Out_files, Param, failmark, verbosity]


#necessary functions:
#function to pass defaults:
def Pass_Defaults():
	return Def_Main

def Pass_BSdefs():
	BS_param = BS_defaults(SESCA_defBS)
	return BS_param

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

#function to return current verbosity level:
def Print_verb():
	print("Current verbosity: %1d"%verbosity)

#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4,"Reading in %1d arguments:\n" % argnum, Args)
#	setting default arguments:
	Inp_files, BS_data, Out_files, Param, failmark, verbosity = Pass_Defaults()
	pdb_file, BB_file, SC_file, ref_file = Inp_files
	lib_name, SS_method, BB_lib, SC_lib = BS_data
	out_file, dec_file = Out_files
	L_range, scale, exp_scale, main_args, prep_args, error_curve = Param

#	modify defaults:
	New_Args = [Inp_files, BS_data, Out_files, Param, failmark, verbosity]
	FLAGS = ["pdb","BB_file","SC_file","ref","write","deconv","lib","method","BB_lib","SC_lib","main","prep","range","scale","refscale","norm","mode", "err", "verb"]
	METHODS = ["Dssp","DS_det","DS_sim","Hbss","Hbss_ext","Seq"]
	flag = ""
	acnt = 0
#	processing new passed arguments: 
	Vprint(2, "Recognized flags:")
	for arg in Args:
		if arg.startswith("@"):
			if not flag in ["prep","main"]:
				flag = arg.strip("@")
				if flag in FLAGS:
					Vprint(2, flag)
				else:
					Vprint(1, "Unknown flag:",flag)
					New_Args[4] = 1
			else:
				if arg.strip("@") in FLAGS:
					flag = arg.strip("@")
				elif flag == "main" and not arg.strip("@") in FLAGS:
					New_Args[3][3].append(arg)
				elif flag == "prep" and not arg.strip("@") in FLAGS:
					New_Args[3][4].append(arg)
#		basic I/O flags					
		elif flag == "pdb":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "BB_file":
			New_Args[0][1] = arg
			flag = ""	
		elif flag == "SC_file":
			New_Args[0][2] = arg
			flag = ""	
		elif flag == "ref":
			New_Args[0][3] = arg
			flag = ""
		elif flag == "write":
			if arg == "0":
				New_Args[2][0] = ""
			else:
				New_Args[2][0] = arg
			flag = ""
		elif flag == "deconv":
			if arg == "0":
				New_Args[2][1] = ""
			else:
				New_Args[2][1] = arg
			flag = ""
		elif flag == "verb":
			try:
				New_Args[5] = int(arg)
			except Exception:
				Vprint(1, "@verb only takes integer arguments")
				New_Args[4] = 1
			flag = ""
#		Basis set flags
		elif flag == "lib":
			New_Args[1][0] = arg
			flag = ""
		elif flag == "method":
			if arg in METHODS:
				New_Args[1][1] = arg
				New_Args[1][0] = "custom"
			else:
				Vprint(1, "Unreconized method: %1s" % arg)
				Vprint(1, "Method options are:",METHODS)
			flag = ""
		elif flag == "BB_lib":
			New_Args[1][2] = arg
			New_Args[1][0] = "custom"
			flag = ""
		elif flag == "SC_lib":
			New_Args[1][3] = arg
			New_Args[1][0] = "custom"
			flag = ""
#		Spectrum parameters:
		elif flag == "range":
			try:
				parts = arg.split(",")
				lower = float(parts[0])
				upper = float(parts[1])
				if upper >= lower:
					New_Args[3][0] = [lower,upper]
				else:
					New_Args[3][0] = [upper,lower]
				for entry in [ "@range", arg ]:
                        	        New_Args[3][3].append(entry)
			except Exception:
				Vprint(1, "@range only takes two comma-separated float arguments")	
				New_Args[4] = 1
			flag = ""
		elif flag == "scale":
			if arg == "0":
				New_Args[3][1] = 1.0
			else:
				try:
                                	New_Args[3][1] = float(arg)
				except Exception:
					Vprint(1, "@scale only takes float arguments or '0'")
					New_Args[4] = 1
				for entry in ["@scale",arg]:
                        	        New_Args[4].append(entry)
			flag = ""
		elif flag == "refscale":
			if arg == "0":
				New_Args[3][2] = ""
			elif arg == "auto":
				New_Args[3][2] = arg
			else:
				try:
                                	New_Args[3][2] = float(arg)
				except Exception:
					Vprint(1, "@refscale takes only 'auto', '0', or floats as arguments")
					New_Args[4] = 1	
			flag = ""
		elif flag == "norm":
			if arg in ["0","1"]:
				for entry in ["@norm",arg]:
                        	        New_Args[3][3].append(entry)
			else:
				Vprint(1, "@refscale takes only 'auto', '0', or floats as arguments")
				New_Args[4] = 1		
			flag = ""
		elif flag == "err":
			try:
				if arg == "auto":
					New_Args[3][5] = arg
				else:
                               		New_Args[3][5] = int(arg)
			except Exception:
				Vprint(1, "@err takes only 'auto' or integer arguments")
				New_Args[4] = 1	


#		Other control fags:
		elif flag == "mode":
			for entry in ["@norm",arg]:
                                New_Args[3][3].append(entry)
			flag = ""
			
#		allow custom arguments for SESCA modules and preprocessors:rgs[4] = int(ar
		elif flag == "main":
			parts = arg.split()
			for part in parts:
				New_Args[3][3].append(part)

		elif flag == "prep":
			parts = arg.split()
			for part in parts:
				New_Args[3][4].append(part)
		
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0][0] == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1] == "" and acnt == 1:
			New_Args[0][3] = arg
			flag = ""
		elif flag == "" and acnt == 2:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1, "Unknown argument:",arg)
			New_Args[4] = 1

		acnt += 1
	return New_Args


#function to import auxilary modules:
def Import_Custom():
	IMPORTS = []
#	importing SESCA core modules:
	if not SESCA_scripts in sys.path:
		Vprint(2,"Using SESCA modules in: %1s" %SESCA_scripts)
		sys.path.append(SESCA_Dir)

	try:
		Vprint(2,"\nLoading Modules:")
		globals()["Calc"] = __import__("SESCA_pred")
		IMPORTS.append("Pred")
		globals()["Seq"] = __import__("SESCA_seq")
		IMPORTS.append("Seq")
		globals()["Dssp"] = __import__("SESCA_dssp")
		IMPORTS.append("Dssp")	
	except Exception:
		Vprint(1,"\nCannot load basic modules!")
		Vprint(1,"Please double check your specified SESCA directories!")

#	importing HBSS modules:
	if not HBSS_Dir in sys.path:
		Vprint(2,"Using HBSS modules in: %1s" % HBSS_Dir)
		sys.path.append(HBSS_Dir)
	try:
		globals()["HBSS"] = __import__("HBSS_main")
		IMPORTS.append("Hbss")	
	except ImportError:
		Vprint(1,"HBSS classification module disabled!")

#	import auxilary modules
	try:
		globals()["Scale"] = __import__("SESCA_scale")
		IMPORTS.append("scale")	
	except ImportError:
		Vprint(2,"SESCA Scaling module disabled!")
	try:
		globals()["Deconv"] = __import__("SESCA_deconv")
		IMPORTS.append("deconv")	
	except ImportError:
		Vprint(2,"SESCA deconvolution module disabled!")

	return IMPORTS



#function to read in Basis set defaults:
def BS_defaults(libinfo):
#	check if file exists:
	if os.path.isfile(libinfo) == False:
		Vprint(1, "\nDefault library file %1s not found!" % libinfo)
		Vprint(1, "Please check your SESCA paths and the SESCA_defBS variable")
		return "None"

	Vprint(4, "%1s:" % libinfo)
	Basis_Sets = {}
	i = open(libinfo,"rb")
#	read in default options
	for line in i:
		line2 = str(line.decode("ascii")).strip("\n")
		if not line2.startswith("#") and line2 != "":
			parts = line2.split()
			part_num = len(parts)
			note = ""
			name,method,BB_lib,SC_lib,calib = parts[0],parts[1],parts[2],parts[3],parts[4]
			for j in range(5,part_num):
				note += str(parts[j])+" "
			Basis_Sets[name] = [method, BB_lib, SC_lib, calib, note]

	i.close()
#	return dictionary
	for entry in Basis_Sets:
		Vprint(2, entry)

	return Basis_Sets
	
#function handle structure processing:
def Handle_Struct(File,SS_method, prep_args, verbosity):
#	check if file exists:
	if os.path.isfile(File) == False:
		Vprint(1, "\nStructure file %1s not found!" % File)
		return "None"

#	getting base for input/output files
	parts = File.split(".")
	ext_num = len(parts[-1])
	Source =   os.path.split(File)[0]
	Filename = os.path.split(File)[-1] 
	namebase = Filename[:-1*(ext_num+1)]
	ext = Filename[-1*(ext_num+1):]
#	print namebase,ext

#	main command line formats for pre-processing:
	Method_Options = {
		"DS_det" : ["%1s  @traj %1s @protlib 1 > prep.log",[DISICL_main, Filename], "DISICL_pdet_%1s.out" ],	
		"DS_sim" : ["%1s  @traj %1s @protlib 2 > prep.log",[DISICL_main, Filename], "DISICL_psim_%1s.out" ],
#		"HbSS_e" : ["%1s  @traj %1s @protlib 2 > prep.log",[DISICL_main, File], "DISICL_psim_%1s.out" ],
#		"Dssp"   : ["%1s %1s %1s %1s > prep.log",[DSSP_ref, File, "Dssp_%1s.out"%namebase, DSSP_main], "Dssp_%1s.out"],
#		"Seq"    : ["%1s %1s > prep.log",[Prot_seq,File], "pseq.out"],
#		"Seq_2"  : ["mv %1s %1s >> prep.log",["pseq.out","Seq_%1s.out" % namebase], "Seq_%1s.out"]
			 }
#	getting extra arguments:

	Output = ""
#	Handle DSSP analysis:
	if SS_method == "Dssp":
		Vprint(1, "Processing File:",File)
		Output = "Dssp_%1s.out"%namebase
		Dssp_args = ["@dssp",DSSP_main,"@pdb",File,"@write",Output,"@verb", str(verbosity)]
		for arg in prep_args:
			Dssp_args.append(arg)
		Dssp.Set_win(win)	
		Dssp_Par = Dssp.Read_Args(Dssp_args)
		Dssp_Result = Dssp.Reformat_Main(Dssp_Par)

#	Handle HBSS analysis:	
	if SS_method == "Hbss" or SS_method == "Hbss_ext":
		Vprint(1, "Processing File:",File)
		Output = "Hbss_%1s.out"%namebase
		Hbss_args = ["@pdb",File,"@write",Output,"@verb", str(verbosity)]
		if SS_method == "Hbss":
			Hbss_args.append("@ext")
			Hbss_args.append("0")	
		for arg in prep_args:
			Hbss_args.append(arg)
		Hbss_Par = HBSS.Read_Args(Hbss_args)
		Hbss_Result = HBSS.HBSS_Main(Hbss_Par)
		global HB_Data
		HB_Data = []
		global SS_Data
		SS_Data = []

#	Handle Sequence analysis:
	elif SS_method == "Seq":
		Output = "Seq_%1s.out"%namebase	
		Seq_args = ["@pdb",File,"@write",Output,"@verb",str(verbosity)]
		Seq_Par = Seq.Read_Args(Seq_args)
		Seq_Result =Seq.Seq_Main(Seq_Par)
			
#	Handle DISICL analysis:
	elif SS_method in ["DS_sim","DS_det"]:
		Vprint(1, "Processing File:",File)
		prep_param = Method_Options[SS_method]
		Vprint(4, prep_param)
		Extra_args = ""
		if Source != "":
			Extra_args += " @path %1s " % Source
		for arg in prep_args:
			Extra_args += arg+" "
		command1 = prep_param[0] % tuple(prep_param[1]) + Extra_args
		Vprint(3, command1)
		os.system(command1)
		Output = prep_param[2] % namebase	
	
	if os.path.isfile(Output) == True:
		Vprint(2, "Structure processed in:",Output)
	else:
		Output = "None"

	return Output

#function to control CD predictions:
def Handle_CDpred(Files,Args):
	SS_file, SS_lib, SC_file, SC_lib,workdir,out_file = Files

#	first add basic backbone calculation parameters:
	Vprint(1, "Calculating backbone contributions:")
	Main_args = []
	BB_output = []
	if not SS_file in ["","None"]:
		for entry in  ["@tar",SS_file,"@lib",SS_lib,"@mode","0","@verb","1"]:
			Main_args.append(entry)
	if not SC_file in["","None"] and out_file != "":
		for entry in ["@write","BB_calc.out"]:
			Main_args.append(entry)	
	else:
		for entry in ["@write","0"]:
			Main_args.append(entry)
	for arg in Args:
		parts = arg.split()
		for entry in parts:
			Main_args.append(entry)
	Vprint(3, "Backbone calculation parameters:\n", Main_args)
#	format and pass on arguments to the CD_calc module:
	BB_args = Calc.Read_Args(Main_args)
	BB_output = Calc.CDpred_Main(BB_args)

#	then calclate side_chain contributions if needed:
	Side_args = []
	SC_output = []
	if not SC_file in ["","None"]:
		Vprint(1, "\nCalculating side-chain contributions:")
		Side_args = ["@tar",SC_file,"@lib",SC_lib,"@mode","0","@verb","1"]
		if out_file != "":
			for entry in ["@write","SC_calc.out"]:
				Side_args.append(entry)	
		for arg in Args:
			parts = arg.split()
			for entry in parts:
				Side_args.append(entry)
		Vprint(3, "Sidechain calculation parameters:\n",Side_args)
#		format and pass on arguments to the CD_calc module:
		SC_args = Calc.Read_Args(Side_args)
		SC_output = Calc.CDpred_Main(SC_args)

#	If both BB and SC calculations succeeded, add up spectra:
	Comb_args = []
	Sum_output = []
	Main_data = []
	if BB_output != [] and SC_output != []:
		Vprint(1, "\nCombining backbone and side chain spectra...")
		Sum_spectra = Calc.Modify_Spectra(BB_output[1],SC_output[1],1.0)[0]
		Sum_output = [BB_output[0], Sum_spectra]
		Vprint(5,"\nSum:",Sum_output)

	Data = ["Backbone_data:",BB_output,"Sidechain_data:",SC_output,"Combined_data:",Sum_output]
	return Data



def Compare_Spectra(Ref_data,Pred_data,Args):
	SS_file, SS_lib, SC_file, SC_lib,workdir,out_file, Ref_file, L_range, scale, exp_scale, err_curve = Args
#	Compare combined CD spectrum with the experimental one:
	Compared_prescale = Calc.Compare_Spectra(Ref_data,Pred_data[1])

#	scale reference spectra if required:
	if exp_scale == "auto":
#		find best scaling factor if auto was called:
		Ints = Extract_Ints(Compared_prescale)
		Ints_reversed = [Ints[0],Ints[2],Ints[1]]
		SF0,SF_range,Lambda,Maxiter,norm = [1.0,[0.2,5], 100, 1000, 0]
		Scale_imports = Scale.Import_Custom()
		Scaled_data = Scale.Find_Scaling(SF0, Ints_reversed, SF_range, Lambda, Maxiter,norm,Scale_imports)
#		compare scaled reference spectrum to the predicted spectrum
		exp_SF = Scaled_data[0]
		Scaled_Ref = Scaled_data[2]
		Compared_data = Calc.Compare_Spectra(Scaled_Ref,Pred_data[1])
		Vprint(2, "optimized ref. scaling factor: %1.3f" % exp_SF)

	elif not exp_scale in [1.0,""]:
#		scale and compare reference spectrum if a scaling factor was given
		exp_SF = exp_scale
		Scaled_Ref = Calc.Scale_Spect(Ref_data,exp_SF)
		Compared_data = Calc.Compare_Spectra(Scaled_Ref,Pred_data[1])
		Vprint(2, "used ref. scaling factor: %1.3f" % exp_SF)

	else:
#		pass data if no scaling is needed
		exp_SF = 1.0
		Compared_data = Compared_prescale
		
#	estimate model error if parameters are available:
	Error_pars = Calc.Read_BS(SS_lib,L_range)[3]
	Rmsd = Compared_data[1][1]
	if Error_pars != []:
#		determine the number error calibration curves:
		e_num = -1
		if len(Error_pars[0]) > 3:
			e_num = Error_pars[0][3]
	
		if e_num == -1:
#		use the linear estimate if no curves are available: 
			Error_est = Calc.Estimate_SSerror_linear(Rmsd,Error_pars)
		elif e_num > 0:
#		if curves are available, determine which one to use:	
			if err_curve == "auto" and (exp_SF == 1.0 or e_num < 2):
#			use the first (unscaled) curve if no scaling was done
				calib = 1
				Vprint(3, "Using default calibration curve (1)")
			elif err_curve == "auto" and exp_SF != 1.0 and e_num >= 2:
#			use the second (rescaled) curve if scaling was preformed
				calib = 2
				Vprint(3, "Using rescaled calibration curve (2)")
			elif err_curve < 0 or err_curve > e_num:
#			use the default (unscaled) curve if the requested curve is not available
				calib = 1
				Vprint(2, "Warning, requested calibration curve not available for the basis set!\nFalling back to default...")
			else:
#			use the requested curve otherwise:
				calib = err_curve
				Vprint(3, "Using calibration curve %1d"%calib)

#		use the chosen calibration curve:
			Error_est = Calc.Estimate_SSerror_nonlin(Rmsd,Error_pars,calib)

		Compared_data[1].append(Error_est)
	else:
		Compared_data[1].append([])
	Compared_data[1].append(scale)

	Main_data = [Pred_data[0],Compared_data[0],Compared_data[1],exp_SF]
	return Main_data

#function to extract the wavelength and intensities of compared spectra:		
def Extract_Ints(Compared_data):
	Ints = [[],[],[]]
	for Point in Compared_data[0]:
		Ints[0].append(Point[0]) 
		Ints[1].append(Point[1]) 
		Ints[2].append(Point[2])
	return Ints 
		


#Main function for script execution:
def SESCA_Main(Args):
#	handle arguments:
	Inp_files, BS_data, Out_files, Param, failmark, verbosity = Args
	pdb_file, BB_file, SC_file, ref_file = Inp_files
	lib_name, SS_method, BB_lib, SC_lib = BS_data
	out_file, dec_file = Out_files
	L_range, scale, exp_scale, main_args, prep_args, err_curve = Param
	Ref_spectrum = []
	BS_param = []

#	import required modules:
	Imports = Import_Custom()


#	check basis set defaults:
	Vprint(1, "\nLoading Default basis sets:")
	BS_param = BS_defaults(SESCA_defBS)	

#	check if user defined a basis set instead of a reference spectrum:
	if ref_file != "" and ref_file in BS_param and lib_name == "DS-dT":
		lib_name = ref_file
		ref_file = ""		
#	stop if user requested help on the basis sets:
	if lib_name == "help":
		print(Libhelp)
		sys.exit()

#	Set up basis set parameters:
	elif not lib_name in ["custom","help"] and lib_name in BS_param:
		Vprint(2, "\nBasis set code recognized:",lib_name)
		lib_param = BS_param[lib_name]
		Vprint(3, lib_param)
#		find SS classification method::
		if SS_method in ["","None"] and not lib_param[0] in ["","None"]:
			SS_method = lib_param[0] 
		else:
			Vprint(2, "Default SS method overwritten by:",SS_method)
#		determine backbone basis set:
		if BB_lib in ["","None"]:
			BB_lib = os.path.join(SESCA_lib,lib_param[1]) 
		else:
			Vprint (2,"Default bakcbone library overwritten by:",BB_lib)
#		determin side chain basis set:
		if SC_lib in ["","None"] and not lib_param[2] in ["","None"]:
			SC_lib = os.path.join(SESCA_lib,lib_param[2])
		elif SC_lib in ["","None"]:
			pass 
		else:
			Vprint(2, "Default side chain library overwritten by:",SC_lib)
#	if one of parameters were custom-defined, do not look up defaults:
	elif lib_name == "custom":
		Vprint(2, "custom library parameters requested, no defaults are used.")
	else:
		Vprint(1, "Basis set code not recognized!")

#	make sure we have a basis set:
	if BB_lib in ["","None"] and failmark == 0:
		failmark = 4
		Vprint(1, "Error while reading main basis set file. Script stops!")

#	read in reference spectrum (if given):
	if ref_file != "" and failmark == 0:
		Ref_spectrum = Calc.Read_Spectrum_file(ref_file,L_range)
		if Ref_spectrum in [ [],"None"]:
			failmark = 5
			Vprint(1, "Error while reading reference file. Script stops!")
			Ref_spectrum = []
#		stop if spetrum scaling is requested, but scling module is disabled:
		if  exp_scale == "auto" and not "scale" in Imports:
			failmark = 6
			Vprint(1, "Error: automatic reference spectrum scaling was requested, but scaling module is disabled!")

#	make sure no errors occured yet:
	if failmark != 0:
		print(Usage)
		sys.exit(failmark)

#	Acquire Structural data:
	Vprint(1, "\nCollecting Structural data:")
	Backbone_File = ""
	Sidechain_File = ""
#	collect SS data for CD perdictions:
	if BB_file != "":
#	take processed data if provided:
		SS_check = Calc.Read_Struct_file(BB_file)
		if not SS_check in ["None",[]]:
			Backbone_File = BB_file 
			Vprint(2, "Processed data file read:",BB_file)
	elif pdb_file != "" and not SS_method in ["","None"]:
#	otherwise process PDB structure:
		Vprint(2, "Structure preprocessor:", SS_method)
		Backbone_File = Handle_Struct(pdb_file,SS_method,prep_args, verbosity)
	else:
		Backbone_File = ""

#	Collect sequence data if side-chain correction are needed:
	if SC_file != "" and not SC_lib in ["","None"]:
#	take processed sequence data if provided:
		SC_check = Calc.Read_Struct_file(SC_file)
		if not SC_check in ["None",[]]:
			Sidechain_File = SC_file 
			Vprint(2, "Processed sequence data in:",SC_file)
	elif pdb_file != "" and not SC_lib in ["","None"]:
#	otherwise process PDB sequence:
		Vprint(2, "\nAdding sequence information:")	
		Sidechain_File = Handle_Struct(pdb_file,"Seq",[],verbosity)
	else:
		Sidechain_File = ""

#	check if processing went OK:
	if Backbone_File in ["","None"] and failmark == 0:
		Vprint(1, "Error while reading structure file. Script stops!")
		failmark = 3
	if Sidechain_File in ["","None"] and not SC_lib in ["","None"] and failmark == 0:
		Vprint(1, "Error during sequence proccesing. Script stops!")
		failmark = 4 
	if failmark != 0:
		print(Usage)
		return ["Failed",failmark]
		
#	predict the CD spectrum based on collected data:
	Vprint(1, "\nPredicting CD spectrum:")
	Files = [Backbone_File, BB_lib, Sidechain_File, SC_lib, workdir, out_file]
	CD_calc = Handle_CDpred(Files,main_args)
	if CD_calc[5] != []:
		Pred_data = CD_calc[5]
	else:
		Pred_data = CD_calc[1]
	
#	Compare to reference CD spectrum if one was given:
	if Ref_spectrum != []:
		Vprint(1, "\nComparing to reference:")
		Comp_args = [Backbone_File, BB_lib, Sidechain_File, SC_lib, workdir, out_file,ref_file,L_range, scale, exp_scale, err_curve]
		CD_comp = Compare_Spectra(Ref_spectrum,Pred_data,Comp_args) 	
		Main_Data = [CD_comp[0], CD_comp[1]]
		Aux_Data = CD_comp[2]
		exp_SF = CD_comp[3]
		print_mode = 1
#	otherwise just accept the predicted spectrum:
	else:
		Main_Data = Pred_data
		Aux_Data = ["","","","",[],[],scale]
		exp_SF = 1.0
		print_mode = 0

	DATA = [Main_Data,Aux_Data]

#	print output file if requested:
	if out_file != "":
		Proc_files,Proc_libs = Backbone_File, "\n#    %1s"%BB_lib
		if Sidechain_File != "":
			Proc_files += " + %1s" % Sidechain_File  
			Proc_libs  += "\n#    %1s" % SC_lib
		if exp_SF != 1.0:
			Proc_libs += "\n#Scaling factor applied to reference: %1.3f" % exp_SF  
		Filenames = [workdir,ref_file,Proc_files, Proc_libs]

		Output_data = Calc.Format_Output(Filenames,DATA[0],DATA[1],print_mode)
		Output = Output_data.encode("ascii")
		o = open(out_file,"wb")
		o.write(Output)
		o.close()

#	set scaling factor to correct value:
	Aux_Data[6] = exp_SF

	return 	DATA



#executing main script:
if __name__ == '__main__':
#	read in command line arguments:
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
                        Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Set_verb(Custom_Args[5])
	Out_files = Custom_Args[2]

#	execute main module:
	Vprint(1, "\nExecuting main SESCA module:")	
	Main_Data = SESCA_Main(Custom_Args)
	Vprint(5, "\n",Main_Data)	


#       print run-time message
	ftime = time.time()
	runtime = ftime-stime
	if Main_Data[0] == "Failed":
		Vprint(1, "Script failed,  runtime was %2.2f seconds"%runtime)
	else:
		outfile = " "
		for File in Out_files:
			if File != "":
				outfile += " %1s," %File
		outfile = outfile[:-1]
		Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
		Vprint(1, "SESCA finished sucessfully! Output written to:",outfile)

else:	
	Vprint(1, "SESCA main  module (SESCA_main.py)")

