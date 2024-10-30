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
usage   = "HBSS module for SS classification extension based on beta-sheet twist angles\n"
usage  += "usage: HBSS_extend.py <pdb_file> <hbond_file> <SS_file> <output_file> @flag <argument>\n"
usage  += "Possible command flags are:\n   @pdb specify protein trajectory in PDB format\n "
usage  += "   Pleas note that trajectory frames should be placed between lines starting with MODEL and ENDMDL respectively\n"
usage  += "   @hbond <file> specify hydrogen bond trajectory file (produced by HBSS_prep.py)\n"
usage  += "   @SS_file <file> specify time series file for basic classes (produced by HBSS_basic.py)\n"
usage  += "   @write <0, file> specify extended classification output file name ('0'- dont write, default: HBSS_ext.out)\n"
usage  += "   @twist <0, file> specify twist angle output file name ('0'- dont write, default: 0)\n"
#usage  += "   @sum   <0, file> specify summary file name ('0'- dont write,default: 0)\n"
usage  += "   @SSdet <0,1> specifies if the avg number of residues and SS segments are printed, '0' - no, SS fractions ony (default: '1')\n"
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"

usage   = "\n"
Usage = usage0 + usage + usage0

#default parameters:
pdb_file = ""
hbond_file = ""
SS_file = ""
out_file = "HBSS_ext.out"
twist_file = "HBSS_twists.out"
sum_file = ""
SS_det = 1
failmark = 0
verbosity = 1

Input_files = [pdb_file,hbond_file,SS_file]
Output_files = [out_file,twist_file, sum_file]

#default parameters that should usually not change
#backbone atom names:
Filtered_Atoms = ["N","CA","C"]
#Twist angle boundaries for a normal beta strand:
Twist_Bounds = [["para",3,23],["anti",3,23]]

#Basic and and extended class list:
Classes_Basic = [
		["Alpha-Helix",     [], "4H", 0, 0 ],
                ["Beta-Stand-Para", [],"BSP", 0, 0 ],
                ["Beta-Strand-Anti",[],"BSA", 0, 0 ],
                ["3/10-Helix",      [], "3H", 0, 0 ],
                ["TT-Helix",        [], "5H", 0, 0 ],
                ["H-bonded Turn",   [], "TU", 0, 0 ],
                ["Unclassified",    [],"UNC", 0, 0 ]]

Classes_Extended = [
		["Alpha-Helix",      [], "4H", 0, 0 ],
                ["3/10-Helix",       [], "3H", 0, 0 ],
                ["TT-Helix",         [], "5H", 0, 0 ],
                ["H-bonded Turn",    [], "TU", 0, 0 ],
                ["Unclassified",     [],"UNC", 0, 0 ],
                ["Left-Handed-Para", [],"LHP", 0, 0 ],
                ["Normal-Beta-Para", [],"NBP", 0, 0 ],
                ["Right-Handed-Para",[],"RHP", 0, 0 ],
                ["Left-Handed-Anti", [],"LHA", 0, 0 ],
                ["Normal-Beta-Anti", [],"NBA", 0, 0 ],
                ["Right-Handed-Anti",[],"RHA", 0, 0 ]]

# processed H-bond and classification data:
Data = [] 
HB_Data = []
SS_Data = []
Class_Tables = [Classes_Basic, Classes_Extended]
	
Ext_param = [Twist_Bounds, Filtered_Atoms, SS_det, HB_Data, SS_Data, Class_Tables]
Def_Args = [Input_files, Output_files, Ext_param, failmark, verbosity]


#function definitions:
#functions to pass on defaults:
def Pass_Defaults():
        return Def_Args

def Pass_Table_basic():
	return Classes_Basic

def Pass_Table_extended():
	return Classes_Extended

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

#functions to receive processed H-bond data:
def Set_HBdata(Data):
	global HB_Data
	del HB_Data[:]
	for entry in Data:
		HB_Data.append(entry)

#function to purge processed H-bond data:
def Unset_HBdata():
	global HB_Data
	del HB_Data[:]

#function to receive processed SS data:
def Set_SSdata(Data):
	global SS_Data
	del SS_Data[:]
	for entry in Data:
		SS_Data.append(entry)
	
#function to purge processed SS data:
def Unset_SSdata():
	global SS_Data
	del SS_Data[:]

#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
	New_Args = [Input_files, Output_files, Ext_param, failmark, verbosity]
	FLAGS = ["pdb","hbond","SS_file","write","twist","sum","SSdet","verb"]
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
			if arg == "0":
				New_Args[0][0] = ""
			else:
				New_Args[0][0] = arg
			flag = ""	
		elif flag == "hbond":
			if arg == "0":
				New_Args[0][1] = "" 
			else:
				New_Args[0][1] = arg
			flag = ""	
		elif flag == "SS_file":
			if arg == "0":
				New_Args[0][2] = ""
			else:
				New_Args[0][2] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "twist":
			if arg == "0":
                                New_Args[1][1] = ""
			else:
				New_Args[1][1] = arg
			flag = ""
		elif flag == "sum":
			if arg == "0":
                                New_Args[1][2] = ""
			else:
				New_Args[1][2] = arg
			flag = ""
		elif flag == "SSdet":
			if arg in ["0","1"]:
				New_Args[2][2] = int(arg)
			else:
				Vprint(1, "@det only takes '0' and '1' as arguments")
				New_Args[3] = 1
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
		elif flag == "" and New_Args[1][0] in ["", "HBSS_ext.out"] and acnt == 3:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
#		Adding global variables:
#		global HB_Data
		New_Args[2][3] = HB_Data
#		global SS_Data
		New_Args[2][4] = SS_Data

	return New_Args

#function to extract the data from a line in pdb format
def extract_pdb(frame, line2):
        try:
                anum = int(line2[6:11].strip(' '))
                aID = str(line2[11:16].strip(' '))
                aRI = str(line2[16:20].strip(' '))
                aCI = str(line2[21])
                aRN = int(line2[22:27].strip(' '))
                ax = float(line2[28:39].strip(' '))
                ay = float(line2[38:46].strip(' '))
                az = float(line2[46:55].strip(' '))
#                aen = float(line2[55:60].strip(' '))
#                abf = float(line2[61:66].strip(' '))
                aai = str(line2[76:78].strip('\n'))
                if aai.isalpha() == False:
                        if line2[12].isalpha() == True:
                                aai = str(line2[12:14])
                        else:
                                aai = str(line2[13:14])
#               acharge = float(line2[79:82].strip())
                atime = int(frame)
                Atom_Data = [anum, aID, ax, ay, az, atime, aai, aRI, aCI, aRN]
        except ValueError:
                Vprint(2,"line could not be read in:\n"+line2)
                Atom_Data = "None"
        return Atom_Data

#function to calculate the distance between points P1 and P2 
def distance(COORD):
	try:
# 3D coordinates for P1 and P2
		X1,Y1,Z1 = float(COORD[0][0]), float(COORD[0][1]), float(COORD[0][2])
		X2,Y2,Z2 = float(COORD[1][0]), float(COORD[1][1]), float(COORD[1][2])
#vector definitions for P1P2
		V1 = [(X2-X1),(Y2-Y1),(Z2-Z1)]
#length of the vector
		Distance = math.sqrt(V1[0]*V1[0] + V1[1]*V1[1] + V1[2]*V1[2])
	except ValueError:
			Vprint(2, "Input cannot be read in value set to 0")
			Vprint(2, "Coordinates  taken "+str(COORD))
			Distance = 0
	return Distance

#function to calculate the dot product of two vectors:
def dot_prod(COORD):
	try:
#getting vector and coordinates
		X1,Y1,Z1 = float(COORD[0][0]), float(COORD[0][1]), float(COORD[0][2])
		X2,Y2,Z2 = float(COORD[1][0]), float(COORD[1][1]), float(COORD[1][2])
#returning dot product:
		Dot12 =  X1*X2+Y1*Y2+ Z1*Z2
	except ValueError:	
		Vprint(2, "Input cannot be read in value set to 0")
		Vprint(2, "Coordinates  taken "+str(COORD))
		Dot12 = 0
	return Dot12

#function to calculate the cross product of two vectors:
def cross_prod(COORD):
	try:
#getting vector and coordinates
		X1,Y1,Z1 = float(COORD[0][0]), float(COORD[0][1]), float(COORD[0][2])
		X2,Y2,Z2 = float(COORD[1][0]), float(COORD[1][1]), float(COORD[1][2])
#returning dot product:
		Cross12 = [(Y1*Z2-Y2*Z1),-1*(X1*Z2-X2*Z1),(X1*Y2-X2*Y1)]
	except ValueError:	
		Vprint(2, "Input cannot be read in value set to 0")
		Vprint(2, "Coordinates  taken "+str(COORD))
		Cross12 = 0
	return Cross12

#function to initialize an SS container:
def SS_Container(Table):
	New_Container =[]
	for Element in Table:
		Label, code = Element[0], Element[2]
		New_Element = [Label, [], code, 0, 0]
		New_Container.append(New_Element)
	return New_Container

#function to read in  H-bond trajectory:
def Read_Hbonds(in_file):
	Contacts = []
	Res_codes = []
	Frames = []
	
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

	Hbond_Data = [Contacts,Res_codes,Frames]
	return Hbond_Data

#function to read HBSS_basic data:
def Read_SSdata(SS_file, Class_Table):
	Basic_Classes = Class_Table
	class_num = len(Class_Table)
	Frames = []
	Rescodes = []
	Chain_codes = []

	if os.path.isfile(SS_file) == False:
		Vprint(1, "Error, SS file nor found")
		return "None"

	i = open(SS_file,"rb")
	flag = ""
	for line in i:
		line2 = str(line.decode("ascii")).strip("\n")
#		identify which class the data belongs to:
		if line2.startswith("#") and flag == "":
			Label = line2.split()[0]
			for code in range(class_num):
				if Label == "#"+Class_Table[code][0].split()[0]:
					flag = code
#		or mark data as chain information:
			if flag == "" and line2.startswith("#PDB  chain"):
				flag = "chain"
				
		elif line2.startswith("&") or line2.strip("\n") == "":
			flag = ""
#		put chain code information into the chain list:
		elif line2.startswith("#") and flag == "chain":
			try:
				parts = line2.split()
				chain_ID, pdbdata, molnum, resdata = parts[1], parts[2], parts[4],parts[-1]
				Chain_codes.append([chain_ID,pdbdata,molnum,resdata])
				Vprint(4, "Residue codes for Chain %1s: %10s"%(chain_ID,resdata))
			except Exception:
				Vprint(2, "line could not be read:",line2)

#		put time series into the class container:
		elif not line2.startswith("#") and not flag in ["","chain"]:
			try:
				time,rescode = line2.split()
				Frames.append(float(time))
				Rescodes.append(int(rescode))
				data = (float(time),int(rescode))
				Basic_Classes[flag][1].append(data)
				Basic_Classes[flag][3] += 1
				Vprint(4, data)
			except Exception:
				Vprint(2, "line could not be read:",line2)
		Vprint(5, flag, line2)

		Frames = list(set(Frames))
		Rescodes = list(set(Rescodes))
		frame_num = len(Frames)
		res_total = len(Rescodes)
		SS_data = [frame_num, res_total, Chain_codes]

	return [Basic_Classes,SS_data]

#function to back calculated PDB information:
def Get_PDBcodes(rescode, Chain_data):
		Old_rescode = ["",0]
		for Chain in Chain_data:
#			print Chain
			chain_code = Chain[0]
			chain_min = int(Chain[1].split("-")[0])
			chain_min2 = int(Chain[3].split("-")[0])
			chain_mod = chain_min2 - chain_min
			Segments = Chain[3].split(",")
			for Seg in Segments:
				parts = Seg.split("-")
				Rmin,Rmax = int(parts[0]),int(parts[1])
				if rescode >= Rmin and rescode <= Rmax:
					res_old = rescode-chain_mod
					Old_rescode = [chain_code,res_old]
					break
			if Old_rescode[0] != "":
				break

		return Old_rescode			
				
#function to identify relevant atom codes for twist calculations:
def Find_Twist_atoms(Class_filt,HB_filt,Chain_data,mode,BB_names):
	Found_Pairs = []
	Anti_check = []

	for entry in Class_filt:
		time = entry[0]
		Pj = entry[1]
		Pi = ""
		Atom_codes = []
		Vprint(5, time, Pj)
		Bonds = [[],[]]
#		set framework for parallel strands:
		if mode == 1:
#			find Hbonds with Pi+1 as donor and Pi-1 as acceptors:
			for bond in HB_filt:
				if bond[0] == time and bond[1] == Pj:
					Bonds[0].append(bond[2])
				if bond[0] == time and bond[2] == Pj:
					Bonds[1].append(bond[1])
			Vprint(5, Bonds)
#			make sure that we found both H-bonds:
			for Acc in Bonds[0]:
				if Acc+2 in Bonds[1]:
					Pi = Acc+1
					break

#		set framework for ant-parallel strands:
		if mode == 2 and not Pj in Anti_check:
			for bond in HB_filt:
#			find Hbonds with Pj Pj+1 Pj-1:
				if bond[0] == time and bond[1] in [Pj-1,Pj,Pj+1]:
					Bonds[0].append(bond)	
				if bond[0] == time and bond[2] in [Pj-1,Pj,Pj+1]:
					Bonds[1].append(bond)
			Vprint(5, Bonds)
#			identify bond partners in bonded segment:
			for Don in Bonds[0]:
				if (Pj+1 == Don[1]) and ([time,Don[2]+2,Pj-1] in Bonds[1]) and (math.fabs(Don[1]-Don[2]) >2):
					Pi = Don[2]+1
					break 	
							
#		look up atom codes for twist calculation:
		if Pi != "":
#			print Chain_data
			ChainJ, ResJ = Get_PDBcodes(Pj, Chain_data)
			ChainI, ResI = Get_PDBcodes(Pi, Chain_data)
			Atom_codes = [
				[ChainI,ResI-1, BB_names[2]],
				[ChainI,ResI,   BB_names[0]],
				[ChainI,ResI,   BB_names[1]],
				[ChainI,ResI,   BB_names[2]],
				[ChainI,ResI+1, BB_names[0]],
				[ChainJ,ResJ-1, BB_names[2]],
				[ChainJ,ResJ,   BB_names[0]],
				[ChainJ,ResJ,   BB_names[1]],
				[ChainJ,ResJ,   BB_names[2]],
				[ChainJ,ResJ+1, BB_names[0]],	
				]
			Vprint(3, "time %1.1f res1 %1d (%1s:%1d) - res2 %1d (%1s:%1d)" % (time,Pj,ChainJ,ResJ,Pi,ChainI,ResI))	
			data = [time,Pj,Pi,Atom_codes]
			Found_Pairs.append(data)
			if mode == 2:
				Anti_check.append(Pi)

	return Found_Pairs		

def Preprocess_Frame(Class_Basic,HB_data,param):
	time, BB_names = param
	Class_update = [["Beta-Stand-Para",1],["Beta-Strand-Anti",2]]
	Chain_data = Class_Basic[1][2]
	Class_filtered = []
	Hbonds_filtered = []
	Pair_List = [[]]
	Vprint(4, "\nReprocessing frame %1d" % time)
#	Filter Hbond data by time:
	Vprint(4,"Hbonds:")
	for Hbond in HB_data[0]:
		if Hbond[0] == time:
			Hbonds_filtered.append(Hbond)
	Vprint(4, Hbonds_filtered)
#	Filter Class data by time and Label:
	for Strand in Class_update:
		Vprint(4, "Class: %1s"%Strand[0])
		mode = Strand[1]
		Class = Class_Basic[0][mode]
#		print mode,Class
		for entry in Class[1]:
			if entry[0] == time:
				Class_filtered.append(entry)
		Vprint(4, Class_filtered)

#		Find relevant residue pairs and atoms to calculate twist angles:
		Pairs_Class = Find_Twist_atoms(Class_filtered, Hbonds_filtered,Chain_data,mode,BB_names)
		for Pair in Pairs_Class:
			Pair_residues = [Pair[1]-1,Pair[1],Pair[1]+1,Pair[2]-1,Pair[2],Pair[2]+1]
			for Res in Pair_residues:
				Pair_List[0].append(Res)

		Pair_List.append(Pairs_Class)
		Pair_List[0] = list(set(Pair_List[0]))
		Pair_List[0].sort()

	return Pair_List


def Calculate_Twist(Coord,mode):
	#atom coordinates: 0-4 C(i-1),N(i),CA(i),C(i),N(i+1)
	#atom coordinates: 5-9 C(j-1),N(j),CA(j),C(j),N(j+1)
#	Calculate N-C bond centerpoints:
	B11 = [(Coord[0][0]+Coord[1][0])/2, (Coord[0][1]+Coord[1][1])/2, (Coord[0][2]+Coord[1][2])/2]	
	B12 = [(Coord[3][0]+Coord[4][0])/2, (Coord[3][1]+Coord[4][1])/2, (Coord[3][2]+Coord[4][2])/2]	
	B21 = [(Coord[5][0]+Coord[6][0])/2, (Coord[5][1]+Coord[6][1])/2, (Coord[5][2]+Coord[6][2])/2]	
	B22 = [(Coord[8][0]+Coord[9][0])/2, (Coord[8][1]+Coord[9][1])/2, (Coord[8][2]+Coord[9][2])/2]	
#	Define Beta strand vectors:
	B1  =  [(B12[0]-B11[0]),(B12[1]-B11[1]),(B12[2]-B11[2])]
	B2  =  [(B22[0]-B21[0]),(B22[1]-B21[1]),(B22[2]-B21[2])]
	D21 =  [(Coord[2][0]-Coord[7][0]),(Coord[2][1]-Coord[7][1]),(Coord[2][2]-Coord[7][2])]
#	calculating lengths, dot and cross products:
	L_B1 = distance([B11,B12])	
	L_B2 = distance([B21,B22])
	L_D21 = distance([D21,[0.0,0.0,0.0]])	
	dot_B1B2 = dot_prod([B1,B2])
	cross_B1B2 = cross_prod([B1,B2])
	Vprint(6, "Vector lengths: B1 %1.3f A, B2 %1.3f A, D21 %1.3f A" % ( L_B1, L_B2, L_D21))
	Vprint(6, "B1*B2= %1.3f A, B1xB2 = (%1.3f; %1.3f; %1.3f)" % (dot_B1B2,cross_B1B2[0],cross_B1B2[1],cross_B1B2[2]))
#	sign function and twist angle:
	try:
		angle =  math.acos(dot_B1B2 / (L_B1 * L_B2))
		sign_func = dot_prod([cross_B1B2,D21])
		if sign_func >= 0:
			sign = 1
		else:
			sign = -1
		Twist_angle = sign*math.degrees(angle)
		if mode == 2:
			Twist_angle = 180+Twist_angle
#		recalculate to  +/-180 scale
		if Twist_angle > 180:
				Twist_angle = Twist_angle-360
		Vprint(4, "Calculated Twist: %3.1f"%Twist_angle)
	except Exception:
		Vprint(2, "Could not calculate twist angle!")
		Twist_angle = 360.0

	return Twist_angle

#function to determine Twist angles in beta-strands over the trajectory:					
def Twist_Traj(pdb_file, Class_Basic, Hbonds, Filtered_Atoms):
	frame_num = Class_Basic[1][0]
	res_total = Class_Basic[1][1]
	Twist_Data = []

	def byres1(item):
        	return item[1]
	def bytype(item):
        	return item[4]
	
#	open trajectory file:
	i = open(pdb_file,"rb")
#	determine the number of lines in the file:
	line_num = 0
	for line in i:
		line_num +=1
	i.seek(0,0)
	Vprint(5, "number of lines:, %1d"%line_num)

#	read pdb trajectory:	
	Atom_TMP = []
	Twist_TMP = []
	tcnt = 1
	lcnt = 0
	Vprint(1, "\nDetermining Beta-sheet twist angles:")
	for line in i:
		line2 = str(line.decode("ascii")).strip("\n")
		lcnt += 1
		Vprint(5,line2)
#########################################################################
#		read in atoms with the right name within the frame:
		if line2.startswith("ATOM") or line2.startswith("HETATM"):
			Atom = extract_pdb(tcnt,line2)
			if Atom != "None" and Atom[1] in Filtered_Atoms:
				Atom_TMP.append(Atom)
				Vprint(5, line2, Atom)
				
#########################################################################
		elif line2.startswith("END") or (lcnt == line_num and Atom_TMP != []):
#		Pick the atoms we need for twist calculations:
			time = float(tcnt)
			mode = 0
			Frame_param = [time, Filtered_Atoms]
			Frame_data = Preprocess_Frame(Class_Basic,  Hbonds, Frame_param)
#			print "frame1:",Frame_data[1]
#			print "\nTMP:",Atom_TMP
			for Pairs in [Frame_data[1], Frame_data[2]]:
				mode += 1
#				Resolve atom coordinates 
				for Pair in Pairs:
					Pj,Pi,Atomcodes = Pair[1],Pair[2],Pair[3]
					atom_num = len(Pair[3])
					Coords = []
					Vprint(4, "\nResolving residue pair %1d -- %1d):" %(Pj,Pi))
					for k in range(atom_num):
						Coords.append([0.000,0.000,0.000])
					cnt = 0
					for Atom in Atom_TMP:
						ChainID,Resnum,AtomID = Atom[-2],Atom[-1],Atom[1]
						Label = [ChainID,Resnum,AtomID] 
						Coord = [Atom[2],Atom[3],Atom[4]]
						for k in range(atom_num):
							if Label == Atomcodes[k]:
								Coords[k] = Coord
								cnt += 1
								Vprint(5,"Atom Coordinates found!:\n",Atomcodes[k],Coord)
							if cnt >= atom_num:
								break
#				Calculate twist angles:
					if cnt >= atom_num:
						Twist_angle = Calculate_Twist(Coords,mode)
					else:
						Vprint(2,"Not all atmos were found for twist calculation!(%1d/%1d)"%(cnt,atom_num))
						Twist_angle = 0.0
					btype = "#para"
					if mode == 2:
						btype = "#anti"
					twist_data = [time,Pi,Pj,Twist_angle,btype]
					Twist_TMP.append(twist_data)
		
#			save  results for the frame:
			Twist_TMP.sort(key= byres1)
			Twist_TMP.sort(key= bytype,reverse=True)
			for entry in Twist_TMP:
				Twist_Data.append(entry)
			Twist_TMP = []
			Atom_TMP = []
			tcnt += 1	
#############################################################################################			
		
	i.close()
	return Twist_Data 

#function to twist angles for formatting:
def Format_Twists(Files, Twist_Data):
	pdb_file,hbond_file,SS_file = Files
	Output = ""

#	Header block:
	Header =  "#HBSS_extend module: Beta-sheet twist angles calculated over a PDB trajectory:\n"
	Header += "#Trajectory file: %1s\n" % Files[0]
	if Files[1] != "":
		Header += "#H-bond file: %1s\n" % Files[1]
	if Files[2] != "":
		Header += "#Basic classification file: %1s\n" % Files[2]

	Output += Header+"\n"

#	Twist angle time_series:
	Twists =  "#Twist angles are calculated between adjacent strands at res1 and res2\n"
	Twists += "#Frame    res1  -  res2  Twist (deg)   sheet-type\n"
	for entry in Twist_Data:	
		string = " %5.1f  %6d %6d     %5.1f        %5s\n" % tuple(entry)
		Twists += string
	Output += Twists

	return Output	



#function to classify beta strand residues in a segment based on their twist angle:
def Reclassify_Strands(Class_Basic, Class_Ext,Twist_data,SS_det):	
	class_num = len(Class_Basic[0])
	frame_num = Class_Basic[1][0]
	res_total = Class_Basic[1][1]
	twist_num = len(Twist_data)
	Chain_info = Class_Basic[1][-1]	
	Bounds = [["para",3,23],["anti",3,23]]

#	Updata data on each class:			
	for cnt in range(class_num):
		Class = Class_Basic[0][cnt]
#		copy data if the class is not a beta strand:
		if cnt == 0:
			for entry in Class[1]:
				Class_Ext[cnt][1].append(entry)
				Class_Ext[cnt][3] += 1
		elif not cnt in [1,2]:
			for entry in Class[1]:
				Class_Ext[cnt-2][1].append(entry)
				Class_Ext[cnt-2][3] += 1
#		reclassify residues if it is a beta strands:
		elif cnt in [1,2]:
			Stype = ["#para",5]
			if cnt == 2:
				Stype = ["#anti",8]
			for entry in Class[1]:
				time,res = entry
				Vprint(3, "Updating residue %1d classification (Frame %1d)" % (res,time))

#				determine relevant twist angles:
				Twist_filter = []
				for item in Twist_data:
					timecode, res1,res2,twist,typecode = item
					if timecode == time and (math.fabs(res-res1) <= 1 or math.fabs(res-res2) <= 1) and typecode == Stype[0]:
						Twist_filter.append(item)
					elif timecode > time:
						break	
				Vprint(5, "Relevant twist angles:",Twist_filter)
#				if twist angles with matching residue numbers were found, average over those:
				acnt = 0
				Avg_twist = 0.0
				for item in Twist_filter:
					if res == item[1] or res == item[2]:
						Avg_twist += item[3]
						acnt += 1
				if acnt == 0:
#				otherwise use neighboring twist angles:
					for item in Twist_filter:
						Avg_twist += item[3]
						acnt += 1
#				compute final twist value:
				if acnt != 0:
					Twist_final = Avg_twist/acnt
				else:	
					Vprint(2,"Residue %1d: No matching twist angles were found! Classified as non-twisted residue")
					Twist_final = 0.0

#				Add residue to its new class based on the final twist angle:
				if Twist_final < Bounds[cnt-1][1]:
					Class_code = Stype[1]+0
				if Twist_final >= Bounds[cnt-1][1] and Twist_final <= Bounds[cnt-1][2]:
					Class_code = Stype[1]+1
				if Twist_final > Bounds[cnt-1][2]:
					Class_code = Stype[1]+2
				New_abbr = Class_Ext[Class_code][2]
				Vprint (3, "%1s -> %1s (twist: %3.1f)" % (Class[2],New_abbr,Twist_final))
				Class_Ext[Class_code][1].append(entry)
				Class_Ext[Class_code][3] += 1
	
#	if requested get SS element information:
	if SS_det == 1:	
		for Class in Class_Ext:
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
	for Class in Class_Ext:
		Summary = [Class[2],0.0,0.0,0.0]
		if res_total != 0 and frame_num != 0:
			perc = float(Class[3])/(res_total*frame_num)*100
			res_avg = float(Class[3])/frame_num
			seg_avg = float(Class[4])/frame_num
			Summary[1] = perc
			Summary[2] = res_avg
			Summary[3] = seg_avg	
		Class_info.append(Summary)

#		add trajectory parameters:
		Traj_info = [frame_num, res_total, Bounds, Chain_info] 

	return [Class_Ext, Class_info, Traj_info]	

#function to format SS time series:
def Format_TS(Files,Class_Data,SS_det):
	frame_num = Class_Data[2][0]
	res_total = Class_Data[2][1]
	Output = ""
#	output header:
	Header =  "#HBSS extended classification time series:\n"
	Header += "#trajectory file: %1s \n" % Files[0]
	if Files[1] != "":
		Header += "#H-bond file: %1s \n" % Files[1]
	if Files[2] != "":
		Header += "#HBSS_basic classification: %1s \n" % Files[2]
	Header += "\n#HBSS residue encoding ( %1d total residues):\n" % res_total
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


#Main function for script execution:
def Extend_Main(Args):
#	set run parameters:
	Input_files, Output_files, Ext_param,  failmark, verbosity = Args
	pdb_file, hbond_file, SS_file = Input_files
	out_file, twist_file, sum_file = Output_files
	Twist_Bounds, Filtered_Atoms, SS_det, HB_Data, SS_Data, SS_Tables = Ext_param
	Set_verb(Args[4])

	Ext_Data = []
#	execute main code:
	
	if failmark != 0:
		Vprint(1, Usage)
		sys.exit(failmark)

#	Read hbond info:
	if hbond_file != "" and HB_Data == []:
		Vprint(2, "\nReading H-bond file:")
		HB_Data = Read_Hbonds(hbond_file)
	elif not HB_Data in [[],"None"]:
		Vprint(2, "\nH-bond data received!")
		for entry in HB_Data:
			Vprint(4, entry)
	else:
		Vprint(1, "\nNo H-bond data was obtained!")
		failmark = 2
		 

#	Read HBSS basic classification:
	if SS_file != "" and SS_Data == []:
		Vprint(2, "\nReading basic SS file:")
		Container_Basic = SS_Container(SS_Tables[0])
		SS_Data = Read_SSdata(SS_file, Container_Basic)
	elif not SS_Data in [[],"None"]:
		Vprint(2, "\nBasic SS data received!")
		for entry in SS_Data:
			Vprint(3, entry)
	else:
		Vprint(1, "\nNo basic SS data was obtained!")
		failmark = 2
		

#check if pdb file exists:
	if pdb_file == "" or os.path.isfile(pdb_file) == False:
		Vprint(1, "Error, cannot find pdb file")
		failmark = 2

	if failmark != 0:
		Vprint(1, Usage)
		sys.exit(failmark)


#	Calculate beta-strand twists for the trajectory:
	Twist_Data = Twist_Traj(pdb_file, SS_Data, HB_Data, Filtered_Atoms)
	if twist_file != "":
		Twist_output = Format_Twists(Input_files,Twist_Data)
		o = open(twist_file,"wb")
		o.write(Twist_output.encode("ascii"))
		o.close()

#	Reclassify beta strand residues based on twist angles:
	Container_Ext = SS_Container(SS_Tables[1])
	Ext_Data = Reclassify_Strands(SS_Data, Container_Ext, Twist_Data,SS_det)
	if out_file != "":
		TS_output = Format_TS(Input_files,Ext_Data,SS_det)
		o = open(out_file,"wb")
		o.write(TS_output.encode("ascii"))
		o.close()

#	purge data before returning:
	Unset_HBdata()
	Unset_SSdata()
	HB_Data = []
	SS_Data = []

	return Ext_Data
	

# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Vprint(2, "\nRun parameters:\n", Custom_Args)
	Input_file, Output_files, Ext_param, failmark,verbosity = Custom_Args
	pdb_file,hbond_file,SS_file = Input_files
	out_file, twist_file,sum_file = Output_files
	Twist_Bounds, Filtered_Atoms, SS_det, HB_Data, SS_Data = Ext_param

#       executing main code
	Data_main = Extend_Main(Custom_Args)
#	Vprint(2, "\nMain data:", Data_main)
	Vprint(2, "\nMain data:")
	for entry in Data_main[1]:
		Vprint(2,entry)

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
	Vprint(2,"HBSS extended classification module(HBSS_extend.py)")

