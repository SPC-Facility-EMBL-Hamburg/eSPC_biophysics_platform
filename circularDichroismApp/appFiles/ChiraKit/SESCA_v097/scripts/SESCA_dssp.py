#!/usr/bin/env python
#


import sys
import os
import math
import time


stime = time.time()
workdir = os.getcwd()
#print workdir


#script to analyze trajectories  and reformat secondary structure output using the dssp algorithm:
usage0  = "*******************************************************************\n"
usage   = "DSSP interface for trajectory analysis\nusage: dssp_reform.py <pdb_file> <output_file> <dssp_path> @flag <argument>\n"
usage  += "Possible command flags are:\n   @pdb <file> specify input structure  file (in pdb format)\n"
usage  += "   Pleas note that trajectory frames should be placed between lines starting with MODEL and ENDMDL respectively\n"
usage  += "   @write <file> specify output  file name (default: Dssp_traj.out)\n"
usage  += "   @dssp  <file> specify path to DSSP program (default: /usr/local/bin/dssp)\n"
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"

Usage = usage0 + usage + usage0


#default parameters:
#############################################
#these parameters should be changed on installation:
win = 0
dssp_main = "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097/DSSP/dssp-2.0.4-linux-amd64"
#dssp_main = "/home/gnagy/Programs/SESCA_dev/DSSP/dssp-2.0.4-linux-amd64"
#dssp_main = "D:\\tan/Gottingen/Programs/SESCA_dev/DSSP/dssp-2.0.4-win32.exe"
#dssp_main = "/usr/local/bin/dssp"
#############################################

#the code from here should not be changed
pdb_file = ""
out_file = "Dssp_traj.out"
failmark = 0
verbosity = 2
copy = ""
remove = ""


Def_Args = [dssp_main, pdb_file, out_file, failmark, verbosity]


DSSP_elements =	{
		"G" : ["3-Helix",      0 ],
		"H" : ["4-Helix",      1 ],
		"I" : ["5-Helix",      2 ],
		"S" : ["Bend",         3 ],
		"B" : ["Beta-Bridge",  4 ],
		"E" : ["Beta-strand",  5 ],
		"T" : ["Turn",         6 ],
		"X" : ["Unclassified", 7 ]
		}

DSSP_short = {
		"3-Helix"      : ["3H", "G", 0],
		"4-Helix"      : ["4H", "H", 1],
		"5-Helix"      : ["5H", "I", 2],
		"Bend"         : ["BE", "S", 3],
		"Beta-Bridge"  : ["BB", "B", 4],
		"Beta-strand"  : ["BS", "E", 5],
		"Turn"         : ["TU", "T", 6],
		"Unclassified" : ["UC", "X", 7],
		}

#function definitions:
#function to pass on defaults:
def Pass_Defaults():
	return Def_Args

def Pass_Dssp_dict():
	return DSSP_elements

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
	
#function to set Windows mode:
def Set_win(int):
	global win
	win = int


#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
#	 Set default arguments:
	dssp_main, pdb_file, out_file, failmark, verbosity = Pass_Defaults()

#	modify parameters:
	New_Args = [dssp_main, pdb_file, out_file, failmark, verbosity]
	FLAGS = ["dssp","pdb","write","verb"]
	flag = ""
	acnt = 0
#	processing new passed arguments: 
	Vprint(3, "Recognized flags:")
	for arg in Args:
		if arg.startswith("@"):
			flag = arg.strip("@")
			if flag in FLAGS:
				Vprint(3, flag)
			else:
				Vprint(1, "Unknown flag:",flag)
				New_Args[3] = 1
		elif flag == "dssp":
			New_Args[0] = arg
			flag = ""	
		elif flag == "pdb":
			New_Args[1] = arg
			flag = ""
		elif flag == "write":
			if arg == "0":
				New_Args[2] = ""
			else:
				New_Args[2] = arg
			flag = ""
		elif flag == "verb":
			try:
				New_Args[4] = int(arg)
			except Exception:
				Vprint(1, "@verb only takes integer arguments")
				New_Args[3] = 1
			flag = ""	
#		setting default files if no flags are provided:
		elif flag == "" and acnt == 0:
			New_Args[1] = arg
			flag = ""
		elif flag == "" and acnt == 1:
			New_Args[2] = arg
			flag = ""
		elif flag == "" and acnt == 2:
			New_Args[0] = arg
			flag = ""
		else:
			Vprint(1, "unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	
#	dssp_main, pdb_file, out_file, failmark, verbosity = New_Args
	return New_Args

#function to execute DSSP:
def Run_Dssp(path, infile, outfile, dssp_args,logfile):
#	execute DSSP on the structure file
	dssp_command = "%1s %1s %1s %1s" % (path,infile,outfile,dssp_args)
	if logfile != "":
		dssp_command += " > %1s" % logfile
	Vprint(2, dssp_command)
	os.system(dssp_command)

#	check if the output was printed	
	Dssp_raw = []
	if os.path.isfile(outfile) == True:
#	extract DSSP output from the file:
		d = open(outfile,"rb")
		flag = 0
		for line0 in d:
			line = str(line0.decode("ascii")).strip("\n")
			parts = line.split()
			if parts[0] == "#":
				flag = 1
			elif flag == 1 and not line.startswith("#"):
				OK = 0
				TMP = []
				try:
					num = int(parts[0])
					resnum = int(parts[1])
					chain = parts[2]
					rescode = parts[3]
					struct = parts[4]
					TMP = [num,resnum,chain,rescode,struct]
					OK = 1
				except ValueError:
					Vprint(3, "\nValueError, Trying alternative:")
					try:
						num = int(parts[0])
						resnum = int(parts[1][:-1])
						chain = parts[1][-1]
						rescode = parts[2]
						struct = parts[3]
						TMP = [num,resnum,chain,rescode,struct]
						OK = 1
					except ValueError:
						Vprint(2, 'Line could not be read in:\n',line)
				if TMP != []:
					Dssp_raw.append(TMP)
			elif parts == []:
				flag = 0
		d.close()	

	return Dssp_raw

#function to process DSSP data:
def Process_Dssp(Data,SS_dict,timecode):
	Dssp_Data = [[],[]]
	rcnt = 0
	sscnt = 0

#	sort and arrange SS_elements according to their number
	TMP = []
	for Element in SS_dict.values():
		TMP.append(Element)
		Vprint(3, Element)
	def sort_key(element):
		return element[1]

#	Set up Data matrix:
	Labels = sorted(TMP, key=sort_key)
	for entry in Labels:
		TMP = [entry[0],[],0]
		Dssp_Data[0].append(TMP)
		sscnt += 1
	
#	Fill Data matrix:		
	Res_prev = ["","",""]
	for entry in Data:
		num, resnum, chain, rescode, struct = entry
#		Identify residue classification:
		try:
			ID = SS_dict[struct]
		except KeyError:
			ID = SS_dict["X"]
		Vprint(4, entry, ID)

#		Add datapoint to matrix:
		Datapoint = (timecode,num)
		SSID = ID[1]
		Dssp_Data[0][SSID][1].append(Datapoint)
		Dssp_Data[0][SSID][2] += 1
		rcnt += 1

#		check residue coding:
		Res_current = [num,resnum,chain]
		if Res_prev == ["","",""]:
			Dssp_Data[1].append(Res_current)
		elif Res_prev[2] != Res_current[2]:
			Dssp_Data[1].append(Res_prev)
			Dssp_Data[1].append(Res_current)
		Res_prev = Res_current
#	add SS class number and residue count for statistics checks:
	Dssp_Data[1].append(Res_prev)
	Dssp_Data.append([1,sscnt,rcnt])
		
	return Dssp_Data	

#function to execute DSSP analysis for trajectories:
def Do_Traj(DSSP_param, pdb_file, raw_file, logfile):
	Dssp_path, Dssp_dict, Dssp_args = DSSP_param
	namebase  = "traj_%1d"
	tmp_pdb = namebase+".pdb"
	tmp_dssp = namebase+'.dssp'
	tmpdir = "tmp"
	Data_traj = []

#	delete temp. files, except in debug mode!
	if verbosity >= 5:
		clean = 0
	else:
		clean = 1

#	create temporary dir to perform DSSP:	
	start_file = os.path.join(workdir,pdb_file)
	if os.path.isfile(start_file) == False:
		Vprint(1, "Cannot find input file: %1s" % start_file)
		return "None"
	if os.path.isdir(tmpdir) == False:
		os.system("mkdir "+tmpdir)
	os.chdir(tmpdir)

#	separate trajectory into single-conformation files:
	cnt = 0
	i = open(start_file,"rb")
	for line0 in i:
		line = line0.decode("ascii").strip("\n")
		if line.startswith("MODEL"):
			if cnt != 0:
				f.close()
			cnt += 1
			Temp_file = tmp_pdb % cnt
			f = open(Temp_file,"wb")
			new_line = "MODEL     "+str(cnt)+"\n"
			f.write(new_line.encode("ascii"))
		elif cnt != 0 and not line.startswith("#"):
			f.write(line0)
	if cnt != 0:
		f.close()
		i.close()

	else:
# 		if no frames were found try reading in a single-structure file without "MODEL" labels:
		i.close()
		Temp_file = tmp_pdb % 1
		i = open(start_file,"rb")
		f = open(Temp_file,"wb")
		found = 0
		for line0 in i:
			line = line0.decode("ascii").strip("\n")
			if found == 0 and (line.startswith("ATOM") or line.startswith("HETATM")):
				found = 1
				new_line = "MODEL    1\n"
				f.write(new_line.encode("ascii")) 
				f.write(line0)
			elif not line.startswith("#"):
				f.write(line0)
		
		if found == 1:
			f.write("ENDMDL\n".encode("ascii"))
			cnt += 1
			f.close()

#		stop if the script couldn't find a single atom in the input file. 
		else:
			Vprint(1, "Error while reading pdb file: no lines with ATOM or HETATM were found.")
			Vprint(1, "Please double-check your input file!")
			f.close()
			i.close()
			return "None"
	
	frame_num = cnt

#	execute DSSP analysis  on each frame:
	for cnt in range(1,frame_num+1):	
		Vprint(1, "Processing frame %1d" % cnt)
#		run DSSP on the structure:
		tmp_in, tmp_out = [tmp_pdb%cnt, tmp_dssp%cnt]
		Dssp_raw = Run_Dssp(Dssp_path, tmp_in, tmp_out, Dssp_args, logfile)
#		save raw_data for frame 1 (for comaprison purposes)
		if cnt == 1:
			raw_up = os.path.join("..",raw_file)
			save_frame = "%1s %1s %1s" % (copy, tmp_out,raw_up)
			os.system(save_frame)

#		process raw DSSP data into time series:
		Dssp_proc = Process_Dssp(Dssp_raw, Dssp_dict,cnt)

#		Collect processed data over several frames:
		if cnt == 1:
			Data_traj.append(Dssp_proc[0])		
			Data_traj.append(Dssp_proc[1])
			Data_traj.append(Dssp_proc[2])
			
		else:
#			Check if residues codes are the same:
			if Data_traj[1] != Dssp_proc[1]:
				Vprint(2, "Warning, Residue data mismatch:")
				Vprint(2, "frame 1: %1s\nframe %1d: %1s\n" % (Data_traj[1],cnt, Dssp_proc[1]))
			Data_traj[2][-1] += Dssp_proc[2][-1]
			Data_traj[2][0]  += 1 	
			cnt = 0
			for Class in Dssp_proc[0]:
				Data_traj[0][cnt][2] += Class[2]
				for entry in Class[1]:
					Data_traj[0][cnt][1].append(entry)
				cnt += 1
	
#	Calculate SS statistics:
	Vprint(4,"\nSS-class statistics:")
	SS_stats = []
	Res_total = Data_traj[2][2]
	for Class in Data_traj[0]:
		SS_ID = Class[0]
#		calculate the percentage residues in given class
		if Res_total != 0:
			SS_perc = float(Class[2])/Res_total*100
		else:
			SS_perc = 0.0
#		add class data to the results
		SS_data = (SS_ID,SS_perc)
		Vprint(4,"%1s  %2.1f" % SS_data)
		SS_stats.append(SS_data)
	Data_traj.append(SS_stats)

#	Clean up temporary files if requested:
	os.chdir("..")	
	if clean == 1:
		tmp_log = os.path.join(tmpdir,logfile)
		save_log = " %1s %1s ." % (copy,tmp_log)
		os.system(save_log)
		clean_command = "%1s %1s "%(remove,tmpdir)
		os.system(clean_command)
		if win == 1:
			os.system("rmdir %1s"%tmpdir)

	return	Data_traj
	



#function to format results:	
def Format_Output(Files, Data):
	dssp_file, input_file, raw_file, output_file = Files
	Class_Data, Res_Coding, Res_Data, Stats_Data = Data
	Frame_num,SS_num,Res_tot = Res_Data
	Chain_num = int((len(Res_Coding))/2)
	global DSSP_short
	if Frame_num != 0:
		Res_num = Res_tot/Frame_num
	else:
		Res_num = 0

#	output header
	Output = ""
	Header =  "#Reformatted DSSP classification:\n"
	Header += "#DSSP program: %1s\n#Structure File: %1s\n#Raw DSSP output: %1s\n" % (dssp_file,input_file,raw_file)
	Header += "\n#Residue encoding (%1d total residues):\n#      numbers        Chain       residues\n" % Res_num
	for chain_cnt in range(Chain_num):
		num_start = Res_Coding[2*chain_cnt][0]
		res_start = Res_Coding[2*chain_cnt][1]
		chain = Res_Coding[2*chain_cnt][2]
		num_fin = Res_Coding[2*chain_cnt+1][0]
		res_fin = Res_Coding[2*chain_cnt+1][1]
		chain_string = "# %6d   - %6d :  %1s   %6d   - %6d\n" % (num_start,num_fin,chain,res_start,res_fin)
		Header += chain_string	
	Output += Header+"\n"

#	SS-class Time Series:
	for Class in Class_Data:
		Cl_codes = DSSP_short[Class[0]]
		Series = "#%1s (%1s) residues\n" % (Class[0], Cl_codes[0])	
		for entry in Class[1]:
			string = " %6d   %6d\n" % entry
			Series += string
		Output += Series+"&\n"

#	SS-class statistics:
	Stat_block = "#Secondary structure distribution by DSSP (%):\n"
	for Class in Stats_Data:
		string = "# %1s: %2.1f\n" % (Class[0],Class[1])
		Stat_block += string
	Output += Stat_block

	return Output

	
#Main function for script execution:
def Reformat_Main(Args):
	dssp_main, pdb_file, out_file, failmark, verbosity = Args
	Set_verb(Args[4])
	Dssp_stat = []

#	check if input files exist:
	if os.path.isfile(dssp_main) == False and failmark == 0:
		Vprint(1, "Error, DSSP executable not found:",dssp_main)
		failmark = 2
	 
	if (pdb_file == "" or os.path.isfile(pdb_file) == False) and failmark == 0:
		Vprint(1, "Error, input PDB file not found:",pdb_file)
		failmark = 3

#	select os command words based on the win flag:
	global copy
	global remove
	if win == 0:
		copy = "cp"
		remove = "rm -fr"
	elif win == 1:
		copy = "copy"
		remove = "del"
	elif failmark == 0:
		failmark = 1
		

	if failmark != 0:
		print(Usage)
		sys.exit(failmark)

#	test mode without trajectory handling, turn traj = 0 on only for debugging:	
	traj = 1
	if traj == 0:
#		run DSSP on the structure:
		tmp_in, tmp_out = [pdb_file, "tmp.dssp"]
		Dssp_raw = Run_Dssp(dssp_main, tmp_in, tmp_out, "","logfile" )

#		process raw DSSP data into time series:
		Dssp_proc = Process_Dssp(Dssp_raw, DSSP_elements,1)
		Dssp_proc.append([])
		Dssp_stat = Dssp_proc

	elif traj == 1:
#		execute trajectory analysis:
		tmp_out = "frame_1.dssp"
		Dssp_param = [dssp_main,DSSP_elements,""]
		Dssp_stat = Do_Traj(Dssp_param, pdb_file, tmp_out,"traj.log")


#	write results if necessary:
	if out_file != "" and not Dssp_stat in ["None",[]] :
		Files = [dssp_main, pdb_file, tmp_out,  out_file]
		Output_data = Format_Output(Files,Dssp_stat)
		Output = Output_data.encode("ascii") 
		o = open(out_file,"wb")
		o.write(Output)
		o.close()		

	return Dssp_stat


# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	dssp_main, pdb_file, out_file, failmark, verbosity = Custom_Args
	Set_verb(Custom_Args[4])
#       executing main code
	Vprint(1, "DSSP interface module\nRun parameters:\n", Custom_Args)
	Main_Data = Reformat_Main(Custom_Args)
#	print "\nMain data:\n",Main_Data	

	
#       print run-time message
	ftime = time.time()
	runtime = ftime-stime
	outfile = Custom_Args[2]
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",outfile)
else:
	Vprint(2,"SESCA-DSSP inteface module (SESCA_dssp.py)")
