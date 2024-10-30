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
usage   = "SESCA sequence analysis module for PDB trajectories\n"
usage  += "\nusage: Pseq.py <pdb_file>  <output_file> @flag <argument>\nPossible command flags are:\n"
usage  += "   @pdb  <input_file> specify input PDB file (or trajectory)\n   @write  <outp_file> specify output file name \n"
usage  += "   @fasta <input_file> specify input sequence file (no pdb required)\n"
usage  += "   @addres  <str,str,str...> add standard residue codes for which statistics are printed\n"
usage  += "   default resiude codes: all 20 3-letter codes for natural amino acids)\n"
usage  += "   @reslist  <str,str,str,str> specify new standard resiue list (overwrites defaults) \n"
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"


Usage = usage0 + usage + usage0

#default parameters:
pdb_file = ""
fasta_file = ""
mode = "prot"
out_file = "Pseq.out"
failmark = 0
verbosity = 1

# recognized standard residue names:
RESIDUES = ['ALA','ASP','ASN','ARG','CYS','GLY','GLU','GLN','HIS','ILE',
            'LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

RES_Prot = {'A':'ALA', 'D':'ASP', 'N':'ASN', 'R':'ARG', 'C':'CYS', 'G':'GLY',
            'E':'GLU', 'Q':'GLN', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS',
            'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP',
            'Y':'TYR', 'V':'VAL', 'X':'UNK'}
RES_Nucl = {'A':'ADE', 'C':'CYT', 'G':'GUA', 'T':'THY', 'U':'URA','X':'UNK'}


Input_files = [pdb_file, fasta_file]
Output_files = [out_file]
Param = [mode, RESIDUES]
Def_Args = [Input_files, Output_files, Param, failmark, verbosity]

#function definitions:
#function to pass on defaults:
def Pass_Defaults():
        return Def_Args

#Function to pass default Residue List:
def Pass_ResList():
	return RESIDUES

#function to pass on 1-letter protein codes:
def Pass_Res_Prot():
	return RES_Prot

#function to pass on 1-letter nucleotide codes:
def Pass_Res_Nucl():
	return RES_Nucl


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
	New_Args = []
	New_Args = Pass_Defaults()
	FLAGS = ["pdb","fasta","write","mode","verb","addres","reslist"]
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
		elif flag == "fasta":
			New_Args[0][1] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "mode":
			if arg in ["prot","nucl"]:
				New_Args[2][0] = arg
			else:
				Vprint(1,"@mode takes only 'prot' and 'nucl' as arguments")				
		elif flag == "addres":
			parts = arg.split(",")
			for code in parts:
				New_Args[2][1].append(code)	
			flag = ""
		elif flag == "reslist":
			parts = arg.split(",")
			New_ResList = []
			for code in parts:
				New_ResList.append(code)	
			New_Args[2][1] = New_ResList
			flag = ""
		elif flag == "verb":
			try:
				New_Args[4] = int(arg)
			except Exception:
				Vprint(1, "Error, @verb only takes integer arguments")
				New_Args[3] = 1
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1] in ["","Pseq.out"] and acnt == 1:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args


#function to read in atoms from a PDB file:
def Read_PDB(File, trajmode):
	ATOM_LIST = []
	frame_num = 0

#	Check if file exists:
	print("File %1s"%File)
	if os.path.isfile(File) == False:
		Vprint(1, "File not found:",File)
		return "None"
	
	f = open(File,"rb")
	cnt = 0
	for line in f:
		line2 = line.decode("ascii").strip("\n")
		if line2.startswith('ATOM') == True or line2.startswith('HETATM') == True:
			try:
			
				anum = int(line2[6:11].strip(' '))
				aID = str(line2[11:16].strip(' '))
				aRI = str(line2[16:21].strip(' '))
				aCI = str(line2[21])
				aRN = int(line2[22:27].strip(' '))
				ax = float(line2[28:39].strip(' '))
				ay = float(line2[38:46].strip(' '))
				az = float(line2[46:55].strip(' '))
				data = (cnt,anum, aID, aRI, aRN, aCI, ax, ay, az)
				Vprint(4, data)
				ATOM_LIST.append(data) 
                
			except ValueError:
				Vprint(2, line,'\n could not be read in')
		elif line2.startswith("END"):
			cnt += 1
			if trajmode == 0:
				break
	f.close()
	return ATOM_LIST

#function to read fasta sequence files:
def Read_Fasta(File, Res_List):
	SEQ_LIST = []
#	Check if file exists:
	if os.path.isfile(File) == False:
		Vprint(1, "File not found:",File)
		return "None"

#	read fasta file:
	f = open(File, "rb")
	rcnt = 0
	ccnt = 0
	for line in f:
		line2 = line.decode("ascii").strip("\n")
		if line2 != "" and not line2[0] in ["#",">",";","&"]:
			line_length = len(line2)
			for i in range(line_length):
#				check each character in the line, look up 3-letter codes
				try:
					Rescode = line2[i]
					aRI = Res_List[Rescode]
					rcnt += 1
					Vprint(2,"Residue %1s read as: %3s  %1d"%(Rescode,aRI,rcnt))
						
				except Exception:
#				if no codes are available in the dictionary, dont add entry:
					Vprint(2, "Residue code: %1s not recognized"%Rescode)					
					pass
#				create a fake C-alpha atom for the new residue further processing:
				cnt, anum, aID, aRN = [0,rcnt,"CA",rcnt]
				aCI, ax, ay, az = [str(ccnt), 0.000, 0.000, 0.000]	
				data = (cnt,anum, aID, aRI, aRN, aCI, ax, ay, az)
				Vprint(4, data)
				SEQ_LIST.append(data)
#		Treat each header as a start of a new chain: 
		elif line2 != "" and line2[0] in [">"]:
			Vprint(2,"New chain: %1s"%line2)
			ccnt += 1

	return SEQ_LIST

#function to extract the sequence:
def Process_Frame(Frame,RESIDUES):
	Res_List = [[],[]]
#	set all checks and flags to 0:
	res0 = ''
	res1 = ''
	chain0 = ''
	chain1 = ''
	rcnt = 0
	hcnt = 0
	aRN = 0
	aCI = 0
#	Set Residue statistics block:
	Res_Stats = []
	for j in range(len(RESIDUES)):
		entry = [RESIDUES[j],0]
		Res_Stats.append(entry)
	Res_Stats.append(["Restotal",0])

#	cycle through the atoms:
	for Atom in Frame:
		tcnt,anum, aID, aRI, aRN, aCI, ax, ay, az = Atom
		resstop = 0
		chainstop = 0
#		check if atom belongs to a new residue, and update residue code:
		if res1 != aRN:
			res0 = res1
			res1 = aRN
			resstop = 1
#		check if the residue belongs to a new chain, and update chain code:		
		if chain1 != aCI:
			chain0 = chain1
			chain1 = aCI
			chainstop = 1
#		add a marker if we hit the end of the chain
		if chainstop == 1 and chain0 != chain1:
			mark = '\n#chain '+aCI+ ' starts:\n'
			Res_List[0].append(mark)
			rcnt = 0
#		add first residue code, when we start
		if res0 == '' and resstop == 1:
			Res_List[0].append(aRI)
			rcnt = rcnt + 1
#		add other residue codes when a new resiude number is found:
		if res0 != '' and res1 != '':
			if ((res1 == res0+1) or (res1 != (res0+1) and chainstop ==0)) and resstop == 1:
				Res_List[0].append(aRI)
				rcnt = rcnt + 1
#			Add break mark if residues are not consecutive:
			elif res1 != (res0+1) and resstop == 1 and chainstop != 1:
				mark2 = '\n# '+str(res0)+'\n#chain break: \n#'+str(res1)+'\n'
				Res_List[0].append(mark2)
				Res_List[0].append(aRI)
#			add first residue of a new chain:
			elif chainstop == 1 and resstop == 1 and res1 != res0+1:
				Res_List[0].append(aRI)
				rcnt = rcnt + 1
				

#		update residue statistics:
		if resstop == 1:
			for Res in Res_Stats:
				if aRI ==  Res[0] and Res[0] != "Restotal":
					Res[1] += 1
					Res_Stats[-1][1] += 1
#			record histidines
			if aRI == 'HIS' and resstop == 1:
				Res_List[1].append(aRN)
				hcnt = hcnt + 1
#			add a separatoe after every 20 residues (formatting):
			if rcnt == 20:
				Res_List[0].append('\n#'+str(aRN+1)+'\n')
				rcnt = 0
#	Add residue statistics:
	Res_List.append([])
	residue_num = Res_Stats[-1][1]
	Res_List[2].append(["Total residues counted",residue_num])
	for j in range(len(Res_Stats)-1):
		Res = Res_Stats[j]
		if residue_num != 0:
			Res_perc = float(Res[1])/residue_num*100
		else:
			Res_perc = 0.0
		res_data = [Res[0],Res_perc]
		Vprint(3, "%3s : %2.1f "%tuple(res_data),"%") 
		Res_List[2].append(res_data)
	
	return Res_List

#function to format output data:
def Format_Output(Filenames, Data):
	workdir, infile = Filenames
	Output = ""

#	format header
	Header = "#SESCA sequence module:\n#Workdir: %1s\n#Source file: %1s\n" % (workdir,infile)
	Output += Header

#	format Sequence data:
	Seq_data = ""
	for entry in Data[0]:
		Seq_data += entry+" "
	Output += Seq_data+"\n"

#	print Histidines:
	H_num = len(Data[1])
	Hist = "\n#List of histidines (%1d):\n" % (H_num)
	for His in Data[1]:
		Hist += "%1d\n" % His
	Output += Hist

#	format residue statistics:
	Stats = "\n#Standard Residue composition:\n"
	for entry in Data[2]:
		if entry == Data[2][0]:
			Stats += "# %1s: %1d\n" % tuple(entry)
		else:
			Stats += "# %1s: %2.1f\n" % tuple(entry)
	Output += Stats
	
	return Output	
			


#Main function for script execution:
def Seq_Main(Args):
#	set run parameters:
	Inp_files, Out_files, Param, failmark,verbosity = Args
	pdb_file, fasta_file = Inp_files
	out_file = Out_files[0]
	mode, RESIDUES = Param
	Set_verb(Args[4])

#	select residue dictionary:
	Res_Dict = {}
	if mode == "prot":
		Res_Dict = Pass_Res_Prot()
	elif mode == "nucl":
		Res_Dict = Pass_Res_Nucl()
		RESIDUES = ["ADE","CYT","GUA","THY","URA"]

#	process input data:
	Atoms = []
	in_file = ""
	if failmark == 0 and pdb_file != "":
#	parse PDB file, extract atom list from first frame:
		Vprint(1, "\nReading pdb file:",pdb_file)
		Atoms = Read_PDB(pdb_file,0)
		in_file = pdb_file
	elif failmark == 0 and fasta_file != "":
#	parse FASTA file, create a fake atom list:
		Vprint(1, "\nReading fasta file:",fasta_file)
		Atoms = Read_Fasta(fasta_file, Res_Dict)
		in_file = fasta_file	

#	check if useabledata was extracted:
	if Atoms in ["None",[]] and failmark == 0:
		Vprint(1, "Error while reading input file, script stops!")
		failmark = 2
 
	if failmark != 0:
		print(Usage)
		sys.exit(failmark)

#	determine sequence from the atom list:
	Vprint(1, "\nProcessing Residue information:")
	Res_data = Process_Frame(Atoms,RESIDUES)
	
#	write output if necessary:
	if out_file != "":
		File_path = [workdir,in_file]	
		Output_data = Format_Output(File_path, Res_data)
		Output = Output_data.encode("ascii")
		o = open(out_file,"wb")
		o.write(Output)
		o.close()

	Main_Data = Res_data

	return Main_Data
	



#executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Vprint(2, "\nRun parameters:\n", Custom_Args)
	In_files, Out_files, Param, failmark, verbosity = Custom_Args
	
#	executing main code
	Data_main = Seq_Main(Custom_Args)
	Vprint(2, '\nMain data:')
	for entry in Data_main[0]:
		Vprint(4, entry)

#	get output files:
	output = ""
	for File in Out_files:
		if File != "":
			output += " %1s,"%File
	output = output[:-1]
		

#	print run-time message
	ftime = time.time()
	runtime = ftime-stime
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",output)
else:
        Vprint(2, "SESCA sequence module (SESCA_seq.py)")

