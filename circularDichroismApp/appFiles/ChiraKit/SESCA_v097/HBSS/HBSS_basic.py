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
usage   = "HBSS module for basic secondary structure classification\nusage: HBSS_basic.py <input_file> <output_file> @flag <argument>\n"
usage  += "Possible command flags are:\n   @inp <file> specify hydrogen bond trajectory file (produced by HBSS_prep.py)\n"
usage  += "   @write <0, file> specify output  file name ('0'- dont write, default: HBSS_basic.out)\n"
usage  += "   @sum   <0, file> specify summary file name ('0'- dont write,default: 0)\n"
usage  += "   @SSdet <0,1> specifies if the avg number of residues and SS segments are printed, '0' - no, SS fractions only (default: 1)\n"
usage  += "   @hier <0,1> control hierarchical classification (priority is derermined by SS order, 0 - no, 1 - yes (default: 1)\n"
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 3)\n"

Usage = usage0 + usage + usage0

#default parameters:
in_file = ""
out_file = "HBSS_basic.out"
sum_file = ""
SS_detail = 1
hierarchy = 1
verbosity = 2
failmark = 0

#Secondary structure Classes based on H-bonds:
ELEMENTS = [    ["Alpha-Helix",     [], "4H", 0, 0 ], 
		["Beta-Stand-Para", [],"BSP", 0, 0 ], 
		["Beta-Strand-Anti",[],"BSA", 0, 0 ], 
		["3/10-Helix",      [], "3H", 0, 0 ], 
		["TT-Helix",        [], "5H", 0, 0 ],
		["H-bonded Turn",   [], "TU", 0, 0 ],
		["Unclassified",    [],"UNC", 0, 0 ]] 


Data = []
HB_Data = []
Input_Files = [in_file]
Output_Files = [out_file, sum_file]
Param = [SS_detail, hierarchy, HB_Data, ELEMENTS]
#Def_Args = [in_file, out_file, sum_file, SS_detail, hierarchy, failmark, verbosity, HB_Data]
Def_Args = [Input_Files, Output_Files, Param, failmark, verbosity]


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

#function to pass on process H-bond information:
def Set_HBdata(Data):
	global HB_Data
	del HB_Data[:]
	for entry in Data:
		HB_Data.append(entry)
#function to control the received H-bond information:
def Print_HBdata():
	global HB_Data
	for entry in HB_Data:
		Vprint(1,entry)

#function to remove the received H-bond information:
def Unset_HBdata():
	global HB_Data
	del HB_Data[:]

#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
	New_Args = [Input_Files, Output_Files, Param, failmark, verbosity]
	FLAGS = ["inp","write","sum","SSdet","hier","verb"]
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
				New_Args[5] = 1
		elif flag == "inp":
			New_Args[0][0] = arg
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
		elif flag == "SSdet":
			if arg in ["0","1"]:
				New_Args[2][0] = int(arg)
			else:
				Vprint(1, "@SSdet takes only '0' or '1' as arguments")
				New_Args[3] = 1		
			flag = ""	
		elif flag == "hier":
			if arg in ["0","1"]:
				New_Args[2][1] = int(arg)
			else:
				Vprint(1, "@hier takes only '0' or '1' as arguments")
				New_Args[3] = 1		
			flag = ""	
		elif flag == "verb":
			try:
				New_Args[4] = int(arg)
			except Exception:
				Vprint(1, "@verb takes on integers as arguments")
				New_Args[3] = 1
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0][0] == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1][0] in ["HBSS_basic.out",""] and acnt == 1:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args


#function to read in  H-bond trajectory:
def Read_Hbonds(in_file):
	Contacts = []
	Res_codes = []
	
	if os.path.isfile(in_file) == False:
		Vprint("Error, file not found!")
		return "None"
	
	i = open(in_file,"rb")
	flag = 0
	for line in i:
		line2 = str(line.decode("ascii")).strip("\n")
		if not line2.startswith("#") and line2 != "":
			try:
				parts = line2.split()
				frame = float(parts[0])
				resD = int(parts[1])
				resA = int(parts[2])
				data = [frame,resD,resA]
				Contacts.append(data)
				Vprint(4, data)
			except Exception:
				Vprint(2, "line could not be read:",line2)
		elif line2.startswith("#PDB  chain"):
				flag = 1
		elif line2 == "" or not line2.startswith("#"):
				flag = 0
		elif flag == 1 and line2.startswith("#"):
			try:
				parts = line2.split()
				chain_ID, pdbdata, molnum, resdata = parts[1], parts[2], parts[4],parts[-1]
				Res_codes.append([chain_ID,pdbdata,molnum,resdata])
				Vprint(4, "Residue codes for Chain %1s: %10s"%(chain_ID,resdata))
			except Exception:
				Vprint(2, "line could not be read:",line2)				
		elif line.startswith("# Number of residues"):
			parts = line.strip("\n").split()
			try:
				res_total = int(parts[-1])
			except Exception:
				Vprint(2, "residue number could  not be read:",line)


	Hbond_Data = [Contacts,Res_codes]
	return Hbond_Data

def Process_frame(Frame_data):
#	prepare container:
	Class_TMP = [   [ "4H",[] ],
			[ "BSP",[] ],
			[ "BSA",[] ],
			[ "3H",[] ],
			[ "5H",[] ],
			[ "TU",[] ],
			[ "ALL",[]] ]

	for Bond1 in Frame_data:
		res1 = Bond1[2]
		res2 = Bond1[1]
		time = Bond1[0]
		found = 0
		for Bond2 in Frame_data:
			res3 = Bond2[2]
			res4 = Bond2[1]
#Identifying SS_classes:
#			Alpha_helix:
			if (res2 == res1+4) and (res3 == res1+1) and (res4 == res3+4):
				Vprint(4, "Alpha helical segment found: %1d - %1d"%(res1,res2))
				found = 1
				cnt = res1
				while cnt <= res2:
					data = (time,cnt)	
					if not data in Class_TMP[0][1]:
						Class_TMP[0][1].append(data)
					cnt += 1

#			3/10-helix:
			elif (res2 == res1+3) and (res3 == res1+1) and (res4 == res3+3 or res4 == res3+4):
				Vprint(4, "3/10 helical segment found: %1d - %1d"%(res1,res2))
				found = 1
				cnt = res1
				while cnt <= res2:
					data = (time,cnt)	
					if not data in Class_TMP[3][1]:
						Class_TMP[3][1].append(data)
					cnt += 1

#			TT-helix:
			elif (res2 == res1+5) and (res3 == res1+1) and (res4 == res3+5 or res4 == res3+4):
				Vprint(4, "TT helical segment found: %1d - %1d"%(res1,res2))
				found = 1
				cnt = res1
				while cnt <= res2:
					data = (time,cnt)	
					if not data in Class_TMP[4][1]:
						Class_TMP[4][1].append(data)
					cnt += 1

#			Parallel strand:
			if (res2 == res3) and (res4 == res1+2) and math.fabs((res1-res2)) > 5:
				Vprint(4, "Parallel beta segment found: %3d - %3d, %3d - %3d"%(res1,res4,res2-1,res2+1))	
				found = 1
				Reslist = [res1,res1+1,res4,res2-1,res2,res2+1]
				for Res in Reslist:
					data = (time,Res)
					if not data in Class_TMP[1][1]:
						Class_TMP[1][1].append(data)
#			Anti-parallel strand:
			if (res3 == res2-2) and (res4 == res1+2) and math.fabs((res1-res2)) > 6:
				bonds = 2
				if [time,res1,res2] in Frame_data: 
					bonds += 1
				if [time,res3,res4] in Frame_data:
					bonds += 1
				if bonds >= 4:
					Vprint(4, "Anti-parallel beta segment found: %3d - %3d, %3d - %3d (%1d bonds)"%(res1,res4,res3,res2,bonds))
					found = 1
					Reslist = [res1,res1+1,res4,res3,res3+1,res2]
					for Res in Reslist:
						data = (time,Res)
						if not data in Class_TMP[2][1]:
							Class_TMP[2][1].append(data)

#		Turn structures:		
		if found != 1 and math.fabs(res1-res2) >= 2 and math.fabs(res1-res2) <= 5:
			Vprint(4, "H-bonded turn found: %3d - %3d"%(res1,res2))	
			if res1 > res2:
				cnt = res2
				stop = res1
			else:
				cnt = res1
				stop = res2
			while cnt <= stop:
				data = (time,cnt)	
				if not data in Class_TMP[5][1]:
					Class_TMP[5][1].append(data)
				cnt += 1
							
					
	
	return Class_TMP

def Get_Reslist(Rescodes):
	Res_list = []
	for Chain in Rescodes:
		try:
			parts = Chain[3].split(",")
			for Res in parts:
				Rmin,Rmax = Res.split("-")
				Rmin,Rmax = int(Rmin),int(Rmax)
				for res_num in range(Rmin,Rmax+1):
					Res_list.append(res_num)
		except Exception:
			Vprint(2, "Warning, unable to read chain coding:",Chain)
	return Res_list
	

def Process_traj(Traj_Data,Classes,param):
	Vprint(2,"\nProcessing Hbond trajectroy:")
	SS_det, hierarchy = param
	time_prev = ""
	TMP = []
	Class_new = []
	bond_num = len(Traj_Data[0])
	frame_num = 0
	res_total = 0

#	Copy class data:
	for entry in Classes:
		Class_new.append(entry)

#	determine full residue list
	All_res = []
	if Traj_Data[1] != []:
		All_res = Get_Reslist(Traj_Data[1])
		res_total = len(All_res)
		Vprint(4, All_res)
		Vprint(2, "Total residue number: %1d"%res_total)
#	turn this off (set to 0) only for debugging
#	hierarchy = 1
#	sort class series:
	def sort_byres(item):
		return item[1]

#	read in H-bonds until the end of the frame
	for i in range(bond_num):
		entry = Traj_Data[0][i]
		Vprint(4, entry)
		time,res1,res2 = entry
		if time_prev != "" and time_prev != time or i+1 == bond_num:
			if i+1 == bond_num:
				TMP.append(entry)
#		when the the time code changes, stop and identify SS elements:
			Vprint(2,"\nFrame %1d"% time_prev)
			Class_frame = Process_frame(TMP)
			Vprint(5, Class_frame)
#			Summarize and sort SS elements for the frame:
			Vprint(3, "Class Statistics:")
			ecnt = 0
			for Element in Class_frame:
				if Element != Class_frame[6]:
					resnum = len(Element[1])
#					perform hierarchical sorting:
					TMP = []
					if hierarchy == 1:
						for item in Element[1]:
							if item in Class_frame[6][1]:
								Vprint(5, "entry excluded from %1s" % Element[0])
							else:
								TMP.append(item)
								Class_frame[6][1].append(item)
						resnum_new = len(TMP)
						Vprint(3, "%10s : original residues: %1d, remaining %1d" % (Element[0],resnum,resnum_new))
#					pass results if hierarchical search is turned off
					else:
						TMP = Element[1]
						resnum_new = resnum
						Vprint(3, "%10s :  %1d residues" % (Element[0],resnum))
#					Add classification results:
					TMP.sort(key= sort_byres)
					for item in TMP:
						Class_new[ecnt][1].append(item)
						Class_new[ecnt][3] += 1

				elif Element == Class_frame[6] and hierarchy == 1:
					#List unclassified residues:
					resnum_new = 0
					for Res in All_res:
						item = (time_prev,Res)
						if not item in Element[1]:
							Class_new[ecnt][1].append(item)
							Class_new[ecnt][3] += 1
							resnum_new += 1	
					Vprint(3, "%10s : unclassifed: %1d, (from %1d)" % (Class_new[6][0],resnum_new,res_total))
				ecnt += 1
			TMP = [entry]
			frame_num += 1
#		othewise just add hbonds to the frame
		else:
			TMP.append(entry)		
		time_prev = time


#	if requested get SS element information:
	if SS_det == 1:	
		for Class in Class_new:
			seg_cnt = 0
			Res_prev = ""
			for Res in Class[1]:
				Vprint(5, Res,Res_prev)
				if  Res_prev != (Res[0],Res[1]-1):
					seg_cnt += 1
					Res_prev = Res
				else:
					Res_prev = Res
			Class[4] = seg_cnt
			Vprint(5, "%1s segments: %1d"%(Class[2],Class[4]))
		
#	determine SS statistics:		
	Class_info = []
	for Class in Class_new:
		Summary = [Class[2],0.0,0.0,0.0]
		if res_total != 0 and frame_num != 0:
			perc = float(Class[3])/(res_total*frame_num)*100
			res_avg = float(Class[3])/frame_num
			seg_avg = float(Class[4])/frame_num
			Summary[1] = perc
			Summary[2] = res_avg
			Summary[3] = seg_avg	
		Class_info.append(Summary)

#	add trajectory parameters:
	Traj_info = [frame_num, res_total, bond_num, Traj_Data[1]] 
		
	return [Class_new, Class_info, Traj_info]

def Residue_Stats(Class_Data):
	frame_num = Class_Data[2][0]
	res_total = Class_Data[2][1]
	Class_num = len(Class_Data[0])
	All_res = Get_Reslist(Class_Data[2][3])

#	prepare residue statistic block:
	Stat_data = []
	Labels = []
	cnt = 0
	for Res in All_res:
		SS_stat = [Res]
		for i in range(Class_num):
			if cnt == 0:
				Labels.append(Class_Data[0][i][2])	
			SS_class = Class_Data[0][i][1]
			scnt = 0
#			count how many times a residue appears in that class
			for entry in SS_class:
#				print entry[1],Res
				if entry[1] == Res:
					scnt += 1
			if frame_num != 0:	
				SS_perc = float(scnt)/frame_num*100
			else:
				SS_perc = 0.0
			SS_stat.append(SS_perc)			
		Stat_data.append(SS_stat)
		cnt += 1

	return [Labels,Stat_data, Class_Data[2]]
 
#function to format SS time series:
def Format_TS(Filename,Class_Data,param):
	SS_det = param[0]
	hierarchy = param[1]
	frame_num = Class_Data[2][0]
	res_total = Class_Data[2][1]
	Output = ""
#	output header:
	Header =  "#HBSS classification time series:\n"
	Header += "#H-bond file: %1s \n\n" % Filename
	Header += "#HBSS residue encoding ( %1d total residues):\n" % res_total
	Header += "#PDB  chain  residue codes ---- HBSS  mol   residue codes\n"
	chain_string =  "#      %1s     %1s        :      %6s   %1s\n"
	for Chain in Class_Data[2][3]:
		Header += chain_string % tuple(Chain)
	Output += Header+"\n"	

#	SS-class time series:
	for Class in Class_Data[0]:
		Series = "#%1s (%1s) residues\n#       time    residue\n" % (Class[0],Class[2])
		for entry in Class[1]:
			string = "   %6.1f   %6d\n" % entry
			Series += string
		Output += Series+"&\n"


#	SS-class statistics:
	Stat_block = "#Secondary structure distribution by HBSS (%):\n"
	if SS_det == 0:
		for Class in Class_Data[1]:
			string = "# %3s:   %6.1f\n" % (Class[0],Class[1])
			Stat_block += string
	else:
		Stat_block += "#name       %    residues  segments\n"
		for Class in Class_Data[1]:
			string = "# %3s:  %6.1f   %6.1f     %4.1f\n" % tuple(Class)
			Stat_block += string
	Output += Stat_block

	return Output

def Format_Sum(Filename, Sum_Data):
	Labels = Sum_Data[0]
	class_num = len(Labels)
	frame_num = Sum_Data[2][0]
	res_total = Sum_Data[2][1]
	Output = ""

#	Header:
	Header =  "#HBSS classification residue statistics:\n"
	Header += "#H-bond file: %1s \n\n" % Filename
	Header += "#HBSS residue encoding ( %1d total residues):\n" % res_total
	chain_string =  "#      %1s     %1s        :      %6s   %1s\n"
	for Chain in Sum_Data[2][3]:
		Header += chain_string % tuple(Chain)
	Output += Header+"\n"

#	Statistics Block
	Res_Block = "#residue  "
	for Label in Labels:
		Res_Block += "  %4s  " % Label
	Res_Block += "\n"
	
	Avg = []
	for j in range(class_num):
		Avg.append(0.0)

	for Res in Sum_Data[1]:
		resline = "   %4d   " % Res[0] 
		for j in range(1,class_num+1):
			resline +=  " %5.1f  " % Res[j]
			Avg[j-1] += Res[j]/res_total
		Res_Block +=  resline+"\n"
	Output += Res_Block+"\n"

#	Summary:
	Footer = "# Avg.    "
	for Aj in Avg:
		Footer += " %5.1f  " % Aj
	Output += Footer+"\n"	
	
	return Output


#Main function for script execution:
def Classify_Main(Args):
#	set run parameters:
	Input_Files, Output_Files, Param, failmark, verbosity = Args
	in_file, out_file, sum_file = Input_Files[0], Output_Files[0], Output_Files[1]
	SS_detail, hierarchy, HB_Data, Classes_0 = Param
	Set_verb(Args[4])
	Class_param = []
	Class_Data = []

#	execute main code:

#	copy SS data structure to new container:
	Class_Cont = []
	for Element in Classes_0:
		Label, code = Element[0], Element[2]
		New_Element = [Label, [], code, 0, 0]
		Class_Cont.append(New_Element)		

#	check input data
	if HB_Data != []:
		Vprint(2, "Hbond trajectory received!")
		for entry in HB_Data:
			Vprint(4, HB_Data) 
#	read in Hbond data from file
	if HB_Data == [] and in_file != "" and failmark == 0:
		Vprint(2, "\nReading Hbond trajectory:")
		HB_Data = Read_Hbonds(in_file)

#	check input data:
	if HB_Data in [ "None",[] ] and failmark == 0:
		failmark = 2

	if failmark != 0:
		Vprint(1, Usage)
		sys.exit(failmark)		

#	process trajectory frame by frame:
	Class_param = [SS_detail,hierarchy]
	Class_Data = Process_traj(HB_Data,Class_Cont,Class_param)
	
#	write output file if requested:
	if out_file != "":
		Output = Format_TS(in_file,Class_Data,Class_param)
		o = open(out_file,"wb")
		o.write(Output.encode("ascii"))
		o.close()

	if sum_file != "":
		Stat_Data = Residue_Stats(Class_Data)
		Output2 = Format_Sum(in_file,Stat_Data)
		o = open(sum_file,"wb")
		o.write(Output2.encode("ascii"))
		o.close()
#		print Stat_Data

#	reset Data and return results:
	Unset_HBdata()

	return Class_Data
	



# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Vprint(3, "\nRun parameters:\n", Custom_Args)
	Input_Files, Output_Files, Param, failmark, verbosity = Custom_Args
	in_file, out_file, sum_file = Input_Files[0], Output_Files[0], Output_Files[1]

#       executing main code
	Data_main = Classify_Main(Custom_Args)
	Vprint(2, "\nMain data:", Data_main[1])

#       print run-time message
	ftime = time.time()
	runtime = ftime-stime
	Outputs = ""
	for File in [out_file,sum_file]:
		if File != "":
			Outputs += " %1s," % File
	Outputs = Outputs[:-1]
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",Outputs)
else:
	Vprint(2,"HBSS basic classification module (HBSS_basic.py)")

