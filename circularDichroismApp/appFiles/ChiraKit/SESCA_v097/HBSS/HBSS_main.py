#!/usr/bin/env python
#

import sys
import os
import math
import time

workdir = os.getcwd()
stime = time.time()


#print workdir+"\n"

usage0  = "*******************************************************************\n"
usage   = "\n"
usage   = "HBSS command module for secondary structure classification\n"
usage  += "usage: HBSS_main.py <pdb_file> <output_file> @flag <argument>\n"
usage  += "Possible command flags are:\n   @pdb specify protein trajectory in PDB format\n"
usage  += "          Pleas note that trajectory frames should be placed between lines starting with MODEL and ENDMDL respectively\n"
usage  += "   @hbond <0, file> specify hydrogen bond trajectory file (produced by HBSS_prep.py, default: 0 - determine hbonds)\n"
usage  += "   @SS_file <0, file> specify time series file for basic HBSS classes (produced by HBSS_basic.py, default: 0 - determine basic SS)\n"
usage  += "   @write <0, file> specify extended classification output file name ('0'- dont write, default: HBSS_main.out)\n"
usage  += "   @det <0,1> Request the avg number of residues and segments in the SS class printed, '0' - no, SS fractions ony (default: '1')\n"
usage  += "   @sum   <0, file> specify residue statistics summary file name ('0'- dont write,default: 0)\n"
usage  += "   @twist <0, file> specify twist angle output file name ('0'- dont write, default: 0)\n"
usage  += "   @ext <0,1> request extended classification based on beta-sheet twists ('0' - dont extend, '1' - extend, default: '1')\n"
usage  += "   @hb_mode <0,1> set H-bond determination mode ('0' - heavy atoms only, '1' - with explicit hydrogens, default: '0')\n"
usage  += '   @hb_args <" str"> pass arguments to the H-bond module (HBSS_prep.py), arguments should provided by a string between quotes,\n'
usage  += '            starting with a white space (eg. @hb_args " @verb 4", see HBSS_prep.py for possible arguments, default: None)\n'
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"

Usage = usage0 + usage + usage0



#Specify HBSS files and directories here:
##################################################################################
HBSS_dir =  "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097/HBSS"
HBSS_prep = os.path.join(HBSS_dir,"HBSS_prep.py")
HBSS_basic = os.path.join(HBSS_dir,"HBSS_basic.py") 
HBSS_extend =  os.path.join(HBSS_dir,"HBSS_extend.py")

##################################################################################
#The code from here should not be changed

#importing modules:
if not HBSS_dir in sys.path:
	print("Using HBSS_modules in:",HBSS_dir)
	sys.path.append(HBSS_dir)

try:
	print("Loading HBSS modules:")
	import HBSS_prep as Hbond
	import HBSS_basic as Basic
	import HBSS_extend as Extend
except Exception:
	print("Cannot load basic modules!")
	print("Please check your specified HBSS directory")

#default parameters:
pdb_file = ""
hbond_file = ""
SS_file = ""
out_file = "HBSS_main.out"
sum_file = ""
hbond_sum = ""
twist_sum = ""

ext_mode = 1
SS_det = 0
hbond_mode = 0 
hbond_args = ""
failmark = 0
verbosity = 1

Input_files = [pdb_file, hbond_file,  SS_file]
Output_files = [out_file, sum_file, hbond_sum, twist_sum]
Param = [ext_mode, hbond_mode, SS_det, hbond_args]

Def_Args = [Input_files, Output_files, Param, failmark, verbosity]
HB_Data = []
SS_Data = []

#function definitions:
#function to pass on defaults:
def Pass_Defaults():
        return Def_Args

#function to control verbosity:
def Vprint(level,*Messages):
        if level <= verbosity:
                for message in Messages:
                        print(message),
                print("")
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
	New_Args = [Input_files, Output_files, Param, failmark, verbosity]
	FLAGS = ["pdb","hbond","SS_file","write","sum","twist","ext","hb_sum","hb_mode","hb_args","det","verb"]
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
		elif flag == "pdb":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "hbond":
			New_Args[0][1] = arg
			flag = ""	
		elif flag == "SS_file":
			New_Args[0][2] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "sum":
			if arg == "0":
                                New_Args[1][1] = ""
			else:
				New_Args[1][1] = arg
			flag = ""
		elif flag == "hb_sum":
			if arg == "0":
                                New_Args[1][2] = ""
			else:
				New_Args[1][2] = arg
			flag = ""
		elif flag == "twist":
			if arg == "0":
                                New_Args[1][3] = ""
			else:
				New_Args[1][3] = arg
			flag = ""
		elif flag == "ext":
			if arg in ["0","1"]:
				New_Args[2][0] = int(arg)
			else:
				Vprint(1, "@ext takes only '0' and '1' as arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "hb_mode":
			if arg in ["0","1"]:
				New_Args[2][1] = int(arg)
			else:
				Vprint(1, "@ext takes only '0' and '1' as arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "det":
			if arg in ["0","1"]:
				New_Args[2][2] = int(arg)
			else:
				Vprint(1, "@det takes only '0' and '1' as arguments")
				New_Args[3] = 1
			flag = ""
		elif flag == "hb_args":
			New_Args[2][3] = arg.split()
			flag = ""
		elif flag == "verb":
			New_Args[4] = int(arg)
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0][0] == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1][0] in ["","HBSS_main.out"] and acnt == 1:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args


#function to execute the Hbond_module and process its data:
def Handle_Hbond(pdb_file, sumfile, mode, args, verb):
	Processed_Data = []
	Hbond_param = []
	Raw_data = []
	
	Vprint(2, "\nExecuting H-bond analysis:")
#	run the H-bond module and obtain data
	Hbond.Set_verb(verb)
	files = [pdb_file,"0","0"]
	if sumfile != "":
		files[1] = sumfile
	Hbond_param = ["@pdb",files[0],"@write",files[1],"@sum",files[2],"@mode",str(mode)]
	for arg in args:
		Hbond_param.append(arg)

	Run_param = Hbond.Read_Args(Hbond_param)
	Raw_data = Hbond.Prep_Main(Run_param)

	if Raw_data[0] in [[],"None"]:
		return "None"
	
#	identify the H-bonded pairs:
	for entry in Raw_data[0]:
		frame,resD,resA = entry[0],entry[1],entry[2]
		Processed_Data.append([frame,resD,resA])
	Chain_info = Obtain_ChainInfo(Raw_data[1])

	return [Processed_Data,Chain_info]

#recalculate Chain_code information:
def Obtain_ChainInfo(Chain_data):
	entry_num = len(Chain_data)-2
	res_total = Chain_data[-2]
	Chain_info = []
	for i in range(entry_num):
		Chain = Chain_data[i]
		chain_code = Chain[0]
		res_mod = Chain[2]
		new_codes = ""
		cmin = Chain[1][0][0]
		cmax = Chain[1][-1][1]
		old_codes = "%1d-%1d " %(cmin,cmax)
		for entry in Chain[1]:
			rmin,rmax = entry
			rmin_new = rmin + res_mod
			rmax_new = rmax + res_mod
			new_codes += "%1d-%1d," % (rmin_new,rmax_new)
		Chain_info.append([chain_code,old_codes,i,new_codes[:-1]])

	return Chain_info


#function to execute basic HBSS classification based on Hbond data:
def Handle_SSbasic(Inputs, Outputs, HB_data, Param, verb):
	pdb_file,hbond_file,SS_file = Inputs
	out_file, sum_file, hbond_sum, twist_sum = Outputs
	ext_mode, hbond_mode, SS_det, hbond_args = Param
	Processed_SS = []
	Raw_SS = []

	Vprint(2, "\nPerforming basic HBSS classification:")
#	pass on H-bond data:
	Basic.Set_verb(verb)
	Basic.Set_HBdata(HB_data)
#	set run parameters:
	if ext_mode == 1:
		files = ["0","0"]
	else:
		files = [out_file,sum_file]
	
	if hbond_file != "":
		files.append(hbond_file)
	else:
		files.append(pdb_file)
	Basic_args = ["@write",files[0],"@sum",files[1],"@inp",files[2]]
	if SS_det == 0:
		Basic_args.append("@SSdet")
		Basic_args.append("0")
	 
#	run the classification module:
	Basic_par = Basic.Read_Args(Basic_args)
	Raw_SS = Basic.Classify_Main(Basic_par)

#	extract relevant info:
	if not Raw_SS in [[],"None"]:
		frame_num, res_total, bond_num, Chain_codes = Raw_SS[2]
		Basic_info = [frame_num, res_total, Chain_codes]
		Processed_SS = [[],Basic_info]
		for entry in Raw_SS[0]:
			Processed_SS[0].append(entry)
	else:
		Processed_SS = "None"

#	purge the data
	Basic.Unset_HBdata()

	return Processed_SS

#function to calculate twist angles and execute the extended HBSS classification:
def Handle_SSext(Inputs, Outputs, HB_Data, SS_Data, SS_det, verb):
	pdb_file,hbond_file,SS_file = Inputs
	out_file, sum_file, hbond_sum, twist_sum = Outputs
	Extended_SS = []
	Raw_SSext = []

	Vprint(2, "\nExtending HBSS classification:")
#	pass on current data:
	Extend.Set_verb(verb)
	Extend.Set_HBdata(HB_Data)
	Extend.Set_SSdata(SS_Data)
	Extend_args = ["@pdb",pdb_file ]
	if hbond_file != "":
		Extend_args.append("@hbond")
		Extend_args.append(hbond_file)
	if SS_file != "":
		Extend_args.append("@SS_file")
		Extend_args.append(SS_file)	
	if out_file != "":
		Extend_args.append("@write")
		Extend_args.append(out_file)	
	if twist_sum != "":
		Extend_args.append("@twist")
		Extend_args.append(twist_sum)
	if SS_det == 0:
		Extend_args.append("@SSdet")
		Extend_args.append("0")

#	execute HBSS_extend module:
	Ext_par = Extend.Read_Args(Extend_args)
	Raw_SSext = Extend.Extend_Main(Ext_par)	
	if not Raw_SSext in ["None",[]]:
		for entry in Raw_SSext:
			Extended_SS.append(entry)
	else:
		Extended_SS = "None"

#	add residue statistics if requested:
	if sum_file != "" and Extended_SS != "None":
		if hbond_file != "":
			Inp = hbond_file
		else:
			Inp = pdb_file
		Stat_data = Basic.Residue_Stats(Extended_SS)
		Stat_Out = Basic.Format_Sum(Inp,Stat_data)
		o = open(sum_file,"wb")
		o.write(Stat_Out.encode("ascii"))
		o.close()

#	purge the data we used:
	Extend.Unset_HBdata()
	Extend.Unset_SSdata()	

	return Extended_SS	 

#Main function for script execution:
def HBSS_Main(Args):
#	set run parameters:
	Input_files, Output_files, Param, failmark,verbosity = Args
	pdb_file, hbond_file,  SS_file = Input_files
	out_file, sum_file, hbond_sum, twist_sum = Output_files
	ext_mode, hbond_mode, SS_det, hbond_args = Param
	Set_verb(Args[4])

	Main_Data = []
	global HB_Data
	global SS_Data

#	execute main code:
	
#	#check if input parameters are valid:
	if failmark == 0:
		if ext_mode == 1 and pdb_file == "":
			Vprint(1, "Error, a PDB structure is required for extended classifications")
			failmark = 1

		if ext_mode == 0 and twist_sum != "":
			Vprint(2, "Warning, twist summary is disabled for basic HBSS classification.")
		if hbond_file != "" and hbond_sum != "":
			Vprint(2, "Warning, H-bonds are read from %1s, no new summary is written."%hbond_file)
	
#	determining H-bond network:
	HB_Data = []
	if pdb_file != "" and hbond_file == "" and failmark == 0:
		HB_Data = Handle_Hbond(pdb_file,hbond_sum, hbond_mode, hbond_args, verbosity)

	elif hbond_file != "" and failmark == 0:
		HB_Data = Basic.Read_Hbonds(hbond_file)

	if HB_Data in ["None",[]] and failmark == 0:
		HB_data = []
		Vprint(1,"Error: H-bond analysis failed!")
		failmark = 2

#	Then, classify residues based on Hbonds:
	SS_Data = []
	if HB_Data != [] and SS_file == "" and failmark == 0:
		SS_Data = Handle_SSbasic(Input_files,Output_files, HB_Data, Param, verbosity)
		Vprint(3,"Calculated basic SS data:")
		Vprint(3, SS_Data)
	elif SS_file != "" and failmark == 0:
		Class_table = Extend.Pass_Table_basic()
		SS_Data = Extend.Read_SSdata(SS_file,Class_table)
		Vprint(3, "basic SS data read from: %1s"%SS_file)
	if SS_Data in ["None",[]] and failmark == 0:
		Vprint(1, "Error: Basic classification failed!")
		failmark = 3

#	Finally, if requested extend the classification based on twist angles:
	SS_ext = []
	if ext_mode == 1 and failmark == 0:
		SS_ext = Handle_SSext(Input_files, Output_files, HB_Data, SS_Data, SS_det, verbosity)

		if SS_ext in ["None",[]] and failmark == 0:
			Vprint(1, "Error: Extended classification failed!")
			failmark = 4

#	check if everything went OK, and pass on data if it was:	
	Main_Data = []
	if failmark != 0:
		Vprint(1, Usage)
		sys.exit(failmark)
	elif ext_mode == 1:
		Main_Data = SS_ext
	elif ext_mode == 0:
		Main_Data = SS_Data
	

	return Main_Data
	

# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
#	set new I/O parameters:
	Set_verb(Custom_Args[4])
	Vprint(2, "\nRun parameters:\n", Custom_Args)
	Input_files, Output_files, Param, failmark,verbosity = Custom_Args
	out_file, sum_file, hbond_sum, twist_sum = Output_files

#       executing main code
	Data_main = HBSS_Main(Custom_Args)
	Vprint(3, "\nMain data:")
	for entry in Data_main[1]:
		Vprint(3, entry)

#       print run-time message
	ftime = time.time()
	runtime = ftime-stime
	Outfiles = ""
	for File in Output_files:
		if File != "":
			Outfiles += " %1s," % File
	if Outfiles != "":
		Outfiles = Outfiles[:-1]
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",Outfiles)
else:
	Vprint(2,"HBSS classification command module(HBSS_main.py)")
#	global HB_Data
#	HB_Data = []
#	global SS_Data
#	SS_Data = []
