#!/usr/bin/env python
#

import sys
import os
import math
import time


workdir = os.getcwd()
stime = time.time()


#print workdir+"\n"

#HBSS preparation module to analyze bakcbone hydrogen bonds (and contacts) in apdb trajectory:
usage0  = "*******************************************************************\n"
usage   = "HBSS module for hydrogen bond trajectory analysis\nusage: HBSS_prep.py <pdb_file> <output_file> @flag <argument>\n"
usage  += "Possible command flags are:\n   @pdb <file> specify input structure  file (in pdb format)\n"
usage  += "   Pleas note that trajectory frames should be placed between lines starting with MODEL and ENDMDL respectively\n"
usage  += "   @write <0, file> specify output  file name (0- dont write, default: Hbond_ts.out)\n"
usage  += "   @sum   <0, file> specify summary file name (0- dont write,default: Hbond_sum.out)\n"
usage  += "   @mode <int> set work mode (defaul:t 0)\n      0 - backbone H-bond analysis based on heavy atoms only\n"
usage  += "      1 - H-bond analysis using explicit hydrogens\n      2 - contact analysis based on atom groups\n"
usage  += "   @dist <float> specify contact distance cutoff (for H-bonds or contact partners, default: d < 3.3 A )\n"
usage  += "   @angle <float> specify hydrogen bond angle cutoff (for H-bonds, default: a > 100 deg)\n"
usage  += "   @occ <float> specify occurrence cutoff (for H-bonds or contact partners, default o > 5%)\n"
usage  += "   @def <str|str|str> specify atom name filter for donor, bridge, and acceptor (default: 'CA|N|O,O1,O2,OT1,OT2')\n"
usage  += "         possible string options include 'all' and 'None' to include any or no atoms\n"
usage  += "   @groupA <str> and @groupB <str> specify atom/residue number ranges to determine H-bonds/contacts between\n"
usage  += "         strong options are 'all','None', or coma-separated range of numbers (e.g. @groupA 11-34,40-66, default: 'all' and 'all')\n"
usage  += "   @inp <0,1> specify if groups are defined by atom (0) or residue (1) numbers (default: 1 for residue numbers)\n"
usage  += "   @det <0,1> specify if atom details are printed in the summary file, '0' means only residue information is printed (default: 1)\n"
usage  += "   @verb <int> set verbosity level from 0 to 5 (default: 1)\n"
Usage = usage0 + usage + usage0

#default parameters:
pdb_file = ""
out_file = "Hbond_ts.out"
sum_file = "Hbond_sum.out"
input_type = 1
details = 1
GROUPS = ["",""]
IO_param = [pdb_file, out_file, sum_file, details, input_type, GROUPS] 

mode = 0
occ_cutoff = 5.0
dist_cutoff = 3.3
angle_cutoff = 100.0
Hbond_atoms = [["CA"],["N"],["O","O1","O2","OT1","OT2"]]
MODIFIED = []
EXCLUDED = ["H2O","SOL","HOH","SOLV","WAT"]
HB_param = [mode, dist_cutoff, angle_cutoff, occ_cutoff, Hbond_atoms, MODIFIED]

failmark = 0
verbosity = 1


Def_Args = [IO_param, HB_param, failmark, verbosity]


#function definitions:
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

#function to read in arguments:
def Read_Args(Args):
	argnum = len(Args)
	Vprint(4, "Reading in %1d arguments:\n" % argnum, Args)
	IO_param, HB_param, failmark, verbosity = Pass_Defaults()

	New_Args = [IO_param, HB_param, failmark, verbosity]
	FLAGS = ["pdb","write","sum","det","def","inp","groupA","groupB","mode","dist","angle","occ","verb"]
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
				New_Args[2] = 1
#		basic I/O flags:
		elif flag == "pdb":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
				New_Args[0][1] = ""
			else:
				New_Args[0][1] = arg
			flag = ""
		elif flag == "sum":
			if arg == "0":
				New_Args[0][2] = ""
			else:
				New_Args[0][2] = arg
			flag = ""
		elif flag == "det":
			if arg in ["0","1"]:
				New_Args[0][3] = int(arg)
			else:
				Vprint(1, "@det takes only '0' and '1' as arguments!")
				New_Args[2] = 1
			flag = ""
		elif flag == "inp":
			if arg in ["0","1"]:
				New_Args[0][4] = int(arg)
			else:
				Vprint(1, "@inp takes only '0' and '1' as arguments!")
				New_Args[2] = 1
			flag = ""
		elif flag == "groupA":
				New_Args[0][5][0] = arg
				New_Args[1][5].append("A")
				flag = ""
		elif flag == "groupB":
				New_Args[0][5][1] = arg
				New_Args[1][5].append("B")

#		flags for hydrogen-bond parameters:
		elif flag == "mode":
			if arg in ["0","1","2"]:
				New_Args[1][0] = int(arg)
			else:
				Vprint(1, "@det takes only <0,1,2> as arguments!")
				New_Args[2] = 1
			flag = ""
		elif flag == "dist":
			try:
				New_Args[1][1] = float(arg)
				New_Args[1][5].append("dist")
			except Exception:
				Vprint(1, "@dist takes only floats as arguments!")
				New_Args[2] = 1
			flag = ""
		elif flag == "angle":
			try:
				New_Args[1][2] = float(arg)
				New_Args[1][5].append("ang")
			except Exception:
				Vprint(1, "@angle takes only floats as arguments!")
				New_Args[2] = 1
			flag = ""
		elif flag == "occ":
			try:
				New_Args[1][3] = float(arg)
			except Exception:
				Vprint(1, "@occ takes only floats as argums!")
				New_Args[2] = 1
			flag = ""
		elif flag == "def":
#			try:
				parts = arg.split("|")
				New_Atomdef = []
				for part in parts:
					if part == "None":
						New_Atomdef.append([])
					else:
						New_Atomdef.append(part.split(","))
				Vprint(4, New_Atomdef)
				if len(New_Atomdef) != 3:
					Vprint(1, "@def takes three '|'-separated strings")
					New_Args[2] = 1
				else:
					New_Args[1][4] = New_Atomdef
					New_Args[1][5].append("def")
#			except Exception:
#				Vprint(1, "@def takes three '|'-separated strings")
#				New_Args[2] = 1
		elif flag == "verb":
			New_Args[3] = int(arg)
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0] == "" and acnt == 0:
			New_Args[0] = arg
			flag = ""
		elif flag == "" and New_Args[1] == "" and acnt == 1:
			New_Args[1] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args


#function to set Hbond and atom filtering parameters based on the operation mode:
def Set_Filters(HB_param, input_type, GROUPS,failmark):
	mode, dist_cutoff, angle_cutoff, occ_cutoff, Hbond_atoms, MODIFIED = HB_param
	FILTER = [4,[],[],[]]
	Mode_info = ""

#	default hydrogen bond parameters:
#	heavy atom only H-bonds:
	if mode == 0:
		Mode_info += "\nHbond analysis based on backbone heavy atoms\n"
#		use N---O  distance < 3.3 A and CA-N---O angle > 100 degrees 
		if not "def" in MODIFIED:
			Hbond_atoms = [["CA"],["N"],["O","O1","O2","OT1","OT2"]]
		if not "dist" in MODIFIED:
			dist_cutoff = 3.3
		if not "ang" in MODIFIED:
			angle_cutoff = 100.0
		Mode_info += "Hbond definition:  N---O  distance < %1.1f A and CA-N---O angle > %2.1f degrees" % (dist_cutoff,angle_cutoff)
#	full atom H-bond parameters:
	if mode == 1:
		Mode_info += "\nHbond analysis based on backbone hydrogen-bonds\n"
		if not "def" in MODIFIED:
			Hbond_atoms = [["N"],["H","HN","H1","H2","H3"],["O","O1","O2","OT1","OT2"]]
		if not "dist" in MODIFIED:
			dist_cutoff = 2.3
		if not "ang" in MODIFIED:
			angle_cutoff = 135.0
		Mode_info += "Hbond definition:  H---O  distance < %1.1f A and N-H---O angle > %2.1f degrees" % (dist_cutoff,angle_cutoff)
#	contact analysis mode:
	if mode == 2:
		Mode_info += "\nCustom contact analysis based on atom distance between contact groups\n"
		if not "def" in MODIFIED:
			Hbond_atoms = [[],["all"],["all"]]
		if not "dist" in MODIFIED:
			dist_cutoff = 3.5
		if not "ang" in MODIFIED:
			angle_cutoff = 0
		Mode_info += "contact definition:  A1---A2  distance < %1.1f A" % (dist_cutoff)

	Vprint(2, Mode_info)
#	add residue codes to the filter:
	for Codes in Hbond_atoms:
		for Code in Codes:
			FILTER[3].append(Code)

# 	set up default group filters:
	if mode == 2:
		if GROUPS[0] == "" or GROUPS[1] == "":
			Vprint(1, "Incomplete atom/residue list: Please define two groups (A and B)!")
			failmark = 1
			return [HB_param,FILTER,failmark]
	else:
		if not "A" in MODIFIED and not "B" in MODIFIED:
			GROUPS = ["all","all"] 
		elif "A" in MODIFIED and not "B" in MODIFIED:
			GROUPS[1] = GROUPS[0]	
		elif "B" in MODIFIED and not "A" in MODIFIED:
			GROUPS[0] = GROUPS[1]

#	add atom/residue numbers to the filter:
	gcnt = 1
	for Group in GROUPS:
		if Group == "all":
			FILTER[gcnt].append("all")
		else:
			parts = Group.split(",")
			for entry in parts:
				if entry.find("-") != -1:
					smaller_parts = entry.split("-")
					pair = (int(smaller_parts[0]),int(smaller_parts[1]))
				else:
					number = int(entry)
					pair = (number,number)
				FILTER[gcnt].append(pair)			
		gcnt += 1
		
# 	determine if the numbers stand for atoms or residues:
	if input_type == 0:
		Vprint(2, "Atom numbers will be used for filter groups")
		FILTER[0] = 0
	else:
		Vprint(2, "Residue numbers will be used for filter groups")
		FILTER[0] = -1
	
	Vprint(3, "Filter parameters:",FILTER)
#	return final Hbond definitions and Filters:
	New_HB = [mode, dist_cutoff, angle_cutoff, occ_cutoff, Hbond_atoms, MODIFIED]
	New_Filters = [New_HB,FILTER,failmark]	
	
	return New_Filters	
	
#function to filter desired atoms:
def Filter_Atoms(Atom_data,Filters):
	AtomID = Atom_data[1]
	Fnum = Filters[0]
	Atom_num = Atom_data[Fnum] 
	
	Passed_Filters = 0
#	Check if the Atom/Residue number fits in group1:
	found = 0
	if "all" in Filters[1]:
		found = 1
	else:
		for group in Filters[1]:
			lower,upper = group
			if Atom_num >= lower and Atom_num <= upper:
#				Vprint(5,"G1",lower,upper,Atom_num)
				found = 1
				break
#	Check if Atom/Residue number fits in group2:	
	if "all" in Filters[2]:
		Passed_Filters = 1+found
		found += 1
	else:
		for group in Filters[2]:
			lower,upper = group
			if Atom_num >= lower and Atom_num <= upper:
				Passed_Filters = 1+found
#				Vprint(5,"G2",lower,upper,Atom_num)
				found += 1
				break
#	Set filter to -1 if no match was found:
	if found == 0:
		Passed_Filters = -1
			
#	Check if the Atom ID-s match:
	if not AtomID in Filters[3] and not "all" in Filters[3]:  
		Passed_Filters = -1
	if Passed_Filters != -1:
		Vprint(4, AtomID,Atom_num,Passed_Filters)
	return Passed_Filters
		
#class definition for atoms in pdb format
class ATOMS:
	def __init__(self,data):
		number, ID, X, Y, Z, time, PID, RID, CID, rnum   = data
		self.num = int(number)
		self.ID = str(ID)
		self.x = float(X)
		self.y = float(Y)
		self.z = float(Z)
		self.t = float(time)
		self.pi = str(PID)
		self.ri = str(RID)
		self.ci = str(CID)
		self.rn = int(rnum)
#return all raw data on the atom:
	def raw(self):
		raw_data =  (self.t, self.num, self.ID, self.x, self.y, self.z, self.rn, self.ri, self.ci, self.pi)
		return raw_data


#function to extract the data from a line in pdb format
def extract_pdb(frame, line2):
        try:
                anum = int(line2[6:11].strip(' '))
                aID = str(line2[11:16].strip(' '))
                aRI = str(line2[16:20].strip(' '))
                aCI = str(line2[21])
                aRN = int(line2[22:27].strip(' '))
                ax = float(line2[30:38].strip(' '))
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

#retrun data in pdb formatted string:
def print_pdb(data):
	time, num, aid, ax, ay, az, arn, rid, cid, el = data
	alt = " "
#	el = " "
	ch = " "
	temp = 0.0
	occ = 1.0
	raw_data = ("ATOM",num,aid,alt,rid,cid,arn," ",ax,ay,az,occ,temp,el,ch)
	fstring = '%-6s%5d  %-3s%1s%-3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%1s\n'
#	pdb format:'ATOM', atomcounter, atomtype, altLoc, resName, chainID, resNum, iCode, x, y, z, occupancy, 
#		    tempFactor, element, charge
	pdb_string = fstring % raw_data
	return pdb_string

#distance function:
def dist_square(COORD):
	try:
#	 3D coordinates for P1 and P2
		X1,Y1,Z1 = float(COORD[0][0]),float(COORD[0][1]),float(COORD[0][2])
		X2,Y2,Z2 = float(COORD[1][0]),float(COORD[1][1]),float(COORD[1][2])
        #length of the vector, squared
		Distance_square = ((X2-X1)**2 + (Y2-Y1)**2 + (Z2-Z1)**2)
	except ValueError:
		Vprint(1, "Input cannot be read in value set to 0")
		Vprint(1, "Coordinates  taken "+str(COORD))
		Distance_square = 0
	return Distance_square

#function to calculate the angle between vectors defined by P1-P2 and P2-P3
def angle(COORD):
	try:
		X1,Y1,Z1 = float(COORD[0][0]),float(COORD[0][1]),float(COORD[0][2])
		X2,Y2,Z2 = float(COORD[1][0]),float(COORD[1][1]),float(COORD[1][2])
		X3,Y3,Z3 = float(COORD[2][0]),float(COORD[2][1]),float(COORD[2][2])
#		vector definitions for bonds
		V1 = [(X2-X1),(Y2-Y1),(Z2-Z1)]
		V2 = [(X3-X2),(Y3-Y2),(Z3-Z2)]
#		length of the vectors:         
		L1 = math.sqrt(V1[0]*V1[0] + V1[1]*V1[1] + V1[2]*V1[2])
		L2 = math.sqrt(V2[0]*V2[0] + V2[1]*V2[1] + V2[2]*V2[2])
#		getting dot product:
		Dot12 =  V1[0]*V2[0]+ V1[1]*V2[1]+ V1[2]*V2[2]

		if L1 == 0.0 or L2 == 0.0:
			Vprint(1, "Zero denominator error, value set to 0")
			Angle = 0
		else:
#			print Dot12,L1,L2
			Angle1 = math.acos(Dot12/(L1*L2))
			Angle = 180.0 - math.degrees(Angle1)
#		print "angle in radians: "+ str(Angle1)
	except ValueError:
		Vprint(1, "Input cannot be read in value set to 0")
		Vprint(1, "Coordinates  taken "+str(COORD))
		Angle = 0
	return Angle

#function to find Hydrogen bonds:
def Find_Contacts(GroupD,GroupA,HB_param):
	mode, dist_cutoff, angle_cutoff, occ_cutoff, Hbond_atoms, MODIFIED = HB_param
	dist_cut2 = dist_cutoff**2
	CONTACTS = []
	Lstring = "%1s : %4d %4s : %4s"
	for Atom1 in GroupD:
#		find the bridge atom first:
		if ("all" in Hbond_atoms[1] or Atom1.ID in Hbond_atoms[1]):
			Coords1 = (Atom1.x, Atom1.y, Atom1.z)
			Label1 = (Atom1.ci, Atom1.rn, Atom1.ri, Atom1.ID)
#			find Acceptor atom if we really look for H-bonds:
			if mode < 2:
				for AtomD in GroupD:
					if ("all" in Hbond_atoms[0] or AtomD.ID in Hbond_atoms[0]) and (AtomD.ci,AtomD.rn,AtomD.ri) == (Atom1.ci,Atom1.rn,Atom1.ri):
						Atom0 = AtomD
						Coords0 = (Atom0.x, Atom0.y, Atom0.z)
						Label0 = (Atom0.ci, Atom0.rn, Atom0.ri, Atom0.ID)
						break
	
#			find suitable acceptors/contacts:
			for Atom2 in GroupA:
				Coords2 = (Atom2.x, Atom2.y, Atom2.z)
				Label2  = (Atom2.ci, Atom2.rn, Atom2.ri, Atom2.ID)

#				first calculate if acceptor is within distance:
				squared_distance = dist_square([Coords1,Coords2]) 	
				if ("all" in Hbond_atoms[2] or Atom2.ID in Hbond_atoms[2]) and squared_distance <= dist_cut2 and squared_distance > 0.0:
					if mode < 2:
						hbond_angle = angle([Coords0, Coords1, Coords2])
						if hbond_angle >= angle_cutoff:
							hbond_distance = math.sqrt(squared_distance)
							Vprint(3, "H-bond found (with dist %2.1f A and angle %2.1f deg):" % (hbond_distance,hbond_angle))
							Vprint(3, "",Lstring%Label0,"\n", Lstring%Label1,"\n", Lstring%Label2,"\n")
							HB_props = [Atom1.t,Atom1.num, Atom2.num, hbond_distance, hbond_angle]
							Labels = [Label0,Label1,Label2]
							Atoms = [Atom0, Atom1, Atom2]
							Contact = (HB_props,Labels,Atoms)
							CONTACTS.append(Contact)	
					else: 
						atom_distance = math.sqrt(squared_distance)
						Vprint(3, "Contact found (with dist %2.1f A):" % (atom_distance))
						Vprint(3, Lstring%Label1,"\n",Lstring%Label2,"\n")
						HB_props = [Atom1.t, Atom1.num, Atom2.num, atom_distance]
						Labels = [Label1,Label2]
						Atoms = [Atom1, Atom2]
						Contact = (HB_props,Labels,Atoms)
						CONTACTS.append(Contact)	
								
	return CONTACTS 


#function to read a trajectory:
def Do_traj(pdb_file, GROUPS,Filter,HB_param):
#	run parameters:
	mode, dist_cutoff, angle_cutoff, occ_cutoff, Hbond_atoms, MODIFIED = HB_param
	Vprint(1, "Processing PDB file:",pdb_file)
	f = open(pdb_file, "rb")
	tcnt = 0
	TMP_Atoms = [[],[]]
	TMP_Chains = []
	RES_LIST = []
	HBONDS = []
	res_cnt = 0
	res_prev = 0
	chan_ID = ""
	chain_prev = ""
	Chain = ["",[[0,0]],0]
########################################################################
	for line in f:
		line2 = str(line.decode("ascii")).strip("\n")
		Vprint(5, line2)
#		read in atoms within the frame:
		if line2.startswith("ATOM"):
#		if line2.startswith("ATOM") or line2.startswith("HETATM"):
			Atom = extract_pdb(tcnt,line2)
			belongs = Filter_Atoms(Atom,Filter)
			if not Atom == "None" and belongs != -1:
#		check if atom passes through the filters:
				Chosen_Atom = ATOMS(Atom)
				if belongs == 2:
					TMP_Atoms[0].append(Chosen_Atom)
					TMP_Atoms[1].append(Chosen_Atom)
				else:
					TMP_Atoms[belongs].append(Chosen_Atom)
#		Follow residue and chain codes:
			if Atom != "None":
				res_num = Atom[-1]
				chain_ID = Atom[-2]
#		track chain codes for residue renumbering:
				if chain_ID != chain_prev:
#					apend the previous chain entry
					if Chain[0] != "":
						Chain[1][-1][1] = res_prev
						TMP_Chains.append(Chain)
						Vprint(5,"Adding Chain:",Chain)
#					determine the largerst residue number so far
					res_last = res_prev + Chain[2]
					res_shift = max(res_last,res_cnt)
					if res_num < res_shift:
						mod = res_shift
					else:
						mod = 0
#					start new chain entry
					Chain = [chain_ID,[[res_num,0]],mod]
#		track residue numbers for trajectory check:
				if res_num != res_prev:
#					add a new entry to chain data if there is a sequence gap
					if res_num != res_prev+1 and res_prev != 0 and chain_ID == chain_prev:
						Chain[1][-1][1] = res_prev
						Chain[1].append([res_num,0])
						Vprint(5, "Chain_gap",Chain)
	
					res_cnt += 1
				res_prev = res_num
				chain_prev = chain_ID
									
##############################################################									
#		at the end of the frame process filtered atoms:
		elif line2.startswith("END"):
			for Atom in TMP_Atoms[0]:
				Vprint(5, Atom.raw())

#			find hydrogen bonds by checking bridge atoms:
			tcnt += 1
			hbond_cnt = 0
			New_Contacts = Find_Contacts(TMP_Atoms[0],TMP_Atoms[1],HB_param)
			for Contact in New_Contacts:
				HBONDS.append(Contact)
				hbond_cnt += 1

#			if the grouping is asymmetric, consider the other Hbonds from the other side:
			if TMP_Atoms[0] != TMP_Atoms[1] and mode != 2:
				New_Contacts2 = Find_Contacts(TMP_Atoms[1],TMP_Atoms[0],HB_param)
				for Contact in New_Contacts2:
					HBONDS.append(Contact)
					hbond_cnt += 1
			Vprint(2, "Frame %1d: %1d contacts found" % (tcnt,hbond_cnt))
			 							
#			add residue and chain codes:
			Chain[1][-1][1] = res_prev
			TMP_Chains.append(Chain)
			TMP_Chains.append(res_cnt)
			RES_LIST.append(TMP_Chains)
			Vprint(3,TMP_Chains)
			if res_cnt != RES_LIST[0][-1]:
				Vprint(2, "Warining: frame %1d contained %3d residues,(%3d expected)"%(tcnt,res_cnt,RES_LIST[0][-1]))
#			reset frame variables:
			Vprint(2, "\n\n")
			res_cnt = 0
			res_prev = 0
			TMP_Atoms = [[],[]] 
			TMP_Chains = []
			chain_ID = ""
			chain_prev = ""
			Chain = ["",[[0,0]],0]
#			break

#	Collect data if no END was provided in the file:
	if TMP_Atoms != [[],[]]:
#		find hydrogen bonds by checking bridge atoms:
		tcnt += 1
		hbond_cnt = 0
		New_Contacts = Find_Contacts(TMP_Atoms[0],TMP_Atoms[1],HB_param)
		for Contact in New_Contacts:
			HBONDS.append(Contact)
			hbond_cnt += 1
#		if the grouping is asymmetric, consider the other Hbonds from the other side:
		if TMP_Atoms[0] != TMP_Atoms[1] and mode != 2:
			New_Contacts2 = Find_Contacts(TMP_Atoms[1],TMP_Atoms[0],HB_param)
			for Contact in New_Contacts2:
				HBONDS.append(Contact)
				hbond_cnt += 1
#		add residue and chain info:
		Chain[1][-1][1] = res_prev
		TMP_Chains.append(Chain)
		TMP_Chains.append(res_cnt)
		RES_LIST.append(TMP_Chains)

		Vprint(2, "Frame %1d: %1d contacts found\n\n" % (tcnt,hbond_cnt))
		TMP_Atoms = [[],[]] 

	f.close()
#	print(HBONDS[0])
	Stats = RES_LIST[0]
	Stats.append(tcnt)
		
	return [Stats,HBONDS]

def Renum_Contact(Contact, Resinfo):
	New_Contact = [0,0,0,0]
	res_total, frame_num  = Resinfo[-2],Resinfo[-1]
	chain_num = len(Resinfo)-2
#	get the old chain and residue codes:
	Chain1, ResN1 = Contact[1][-2][0],Contact[1][-2][1]
	Chain2, ResN2 = Contact[1][-1][0],Contact[1][-1][1]
#	modify code based on Resinfo:
	found = 0
	for i in range(chain_num):
		if Chain1 ==  Resinfo[i][0]:
			New_Contact[0] = i
			New_Contact[1] = ResN1 + Resinfo[i][2]
			found += 1
		if Chain2 ==  Resinfo[i][0]:
			New_Contact[2] = i
			New_Contact[3] = ResN2 + Resinfo[i][2]
			found += 1
	if found != 2:
		New_Contact = []
		Vprint(1, "Could not relabel contact:",Contact) 

	return New_Contact

#function to format chain encoding:
def Format_Chain(Chain_data):
	entry_num = len(Chain_data)-2
	res_total = Chain_data[-2]
	Chain_info  =   "#HBSS residue encoding (%4d total residues):\n" % res_total
	Chain_info +=   "#PDB  chain  residue codes ---- HBSS  mol   residue codes\n"
	chain_string =  "#      %1s     %1s        :      %6d   %1s\n"
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
		Chain_info += chain_string % (chain_code,old_codes,i,new_codes[:-1])
	Chain_info = Chain_info + "\n\n"	
	return Chain_info

#function to summarize Hbond statistics:
def Hbond_Summary(Data,Stats,mode):
	contact_total = 0
	res_total, frame_num  = Stats[-2],Stats[-1]
	SUM_RESULTS =[]
	Vprint(1, "\nCalculating Contact Summary:")
	for Contact in Data:
#		loop through the identified Hbonds/contacts
		Label1 = Contact[1][-2]
		Label2 = Contact[1][-1]
		L1_main	= (Contact[2][-2].num,Contact[2][-2].ID)
		L2_main = (Contact[2][-1].num,Contact[2][-1].ID)
		New_Labels = Renum_Contact(Contact,Stats)

		if mode < 2:
#			add donor atom if Hbonds are searched
			Label0 = Contact[1][0]
			L0_main = (Contact[2][0].num,Contact[2][0].ID)
			found = 0
#			check if H-bond was seen before, increase occurrence if yes
			for Result in SUM_RESULTS:
				if Result[0] == L1_main and Result[1] == L2_main:
					found = 1
					Result[2] += 1 
					contact_total += 1
					break
#			if H-bond was not seen yet, add it to the list
			if found == 0:
				new_entry = [L1_main, L2_main, 1, [Label0, Label1, Label2], New_Labels]
				SUM_RESULTS.append(new_entry)
				contact_total += 1
		if mode == 2:
			# if we are looking for contacts organize data according to atom number:
			if L1_main[0] <= L2_main[0]:
				Pos1_main = L1_main
				Pos2_main = L2_main
			else:
				Pos1_main = L2_main
				Pos2_main = L1_main
#			check if Contact was seen before, increase occurrence if yes
			found = 0
			for Result in SUM_RESULTS:
				if Result[0] == Pos1_main and Result[1] == Pos2_main:
					found = 1
					Result[2] += 1 
					contact_total += 1
					break
#			if Contact was not seen yet, add it to the list
			if found == 0:
				new_entry = [Pos1_main, Pos2_main, 1,[Label1,Label2], New_Labels]
				SUM_RESULTS.append(new_entry)
				contact_total += 1

#	Second sorting based on residue identity:
	SUM_RESULTS2 = []
	for Contact in SUM_RESULTS:
		if Contact[4] != []:
#			define residue and atom pair labels:
			pair_occ = 0.0
			if frame_num != 0:
				pair_occ = float(Contact[2])/frame_num*100
			Atom_Pair =[Contact[0],Contact[1],pair_occ]
			Res_Label = (Contact[4][0],Contact[4][1],Contact[3][-2][2], Contact[4][2],Contact[4][3],Contact[3][-1][2])

			found  = 0
			for Result in SUM_RESULTS2:
#			add atom pair if the residue label was already seen:
				if Res_Label == Result[0]:
					found = 1
					Result[1] += pair_occ
					Result[2].append(Atom_Pair)
			if found == 0:
				New_entry = [Res_Label,pair_occ,[Atom_Pair]]
				SUM_RESULTS2.append(New_entry)

				
	return SUM_RESULTS2

#function to format H-bond summary:
def Format_Sum(Run_param,Sum_data, Stats):
#	get raw data:
	workdir, pdb_file, dist_cut, angle_cut, occ_cut, details, mode = Run_param
	res_total, frame_num  = Stats[-2],Stats[-1]
	chain_num = len(Stats)-2
	Output = ""

#	header: input parameter summary
	Header =  "#HBSS Hydrogen bond analysis module\n#Work directory: %1s\n#Source file: %1s\n" % (workdir,pdb_file)
	if occ_cut > 0.0:
		Header += "#Occurrence cutoff at: %1.1f" % occ_cut + "%\n" 
	Header += "#Distance cutoff: %2.1f A\n" % dist_cut
	Output += Header
	
#	footer: frame and residue numbers
	Footer = "# Number of processed frames: %6d\n# Number of residues (frame1): %6d\n" %(frame_num, res_total)
	if frame_num == 0 or res_total == 0:
		Output += "\n\n"+Footer
		return Output

#	chain encoding:
	Chain_codes = Format_Chain(Stats)

#	main data block:
	Cont_stat = ""
	Cont_count = [0,0]
	if mode < 2:
		Cont_stat += "#Angle cutoff: %3.1f deg\n\n" %(angle_cut)
		Cont_stat += Chain_codes
		Cont_stat += "#Searching for backbone H-bonds:\n#        Donor------------Acceptor\n"
	else:
		Cont_stat += "\n#Searching for contact atoms:\n#       ID1  atom1    ID2  atom2   %\n"

	Cont_stat += "# mol1  res1  resID   mol2 res2 resID      %\n"
	for Contact in Sum_data:
		Cont_count[1] += Contact[1]/100
		if Contact[1] >= occ_cut:
			Residue_string = (2*"   %1d  %6d %4s   ") % Contact[0] + "%2.1f\n" % Contact[1]
			Cont_stat += Residue_string
			Cont_count[0] += Contact[1]/100
			if details == 1:
				for Pair in Contact[2]:
					Pair_string = "#  "+"%6d    %4s" % Pair[0] + 4*" " + "%6d    %4s" % Pair[1]+ 3*" " + "%2.1f\n" % Pair[2]
					Cont_stat += Pair_string
	Output += Cont_stat
	
	Footer += "# Average number of Contacts above cutoff: %2.1f / %2.1f " % (tuple(Cont_count)) 	
	Output += "\n"+Footer

	return Output

#function generate H-bond time series:
def Hbond_TS(Data,Stats):
	Hbond_Ts = []
	for Contact in Data:
#		extract H-bond data:
#		print(Contact)
		time,dist = Contact[0][0]+1.0, Contact[0][3]
		New_ResN = Renum_Contact(Contact,Stats)
		if New_ResN != []:
			Res1,Res2 = New_ResN[1],New_ResN[3]
		else:
			Res1,Res2 = Contact[1][-2][1], Contact[1][-1][1]
		
		Data_point = [time,Res1,Res2,dist]
		if len(Contact[0]) >= 5:
			angle = Contact[0][4]
			Data_point.append(angle)
		
		Hbond_Ts.append(Data_point)
	return Hbond_Ts

#function to format H-bond time series:
def Format_TS(TS_param, TS_data, Stats):
#	get raw data:
	workdir, pdb_file, dist_cut, angle_cut, mode = TS_param
	res_total, frame_num  = Stats[-2],Stats[-1]
	Output = ""

#	header: input parameter summary
	Header =  "#HBSS Hydrogen bond analysis module\n#Work directory: %1s\n#Source file: %1s\n" % (workdir,pdb_file)
	Header += "#Distance cutoff: %2.1f A\n" % dist_cut
	Output += Header

#	footer: frame and residue numbers
	Footer = "# Number of processed frames: %6d\n# Number of residues (frame1): %6d\n" %(frame_num, res_total)
	if frame_num == 0 or res_total == 0:
		Output += "\n\n"+Footer
		return Output

#	chain encoding:
	Chain_codes = Format_Chain(Stats)

#	main data block:
	TS_stat = ""
	if mode < 2:
		TS_stat +="#Angle cutoff: %3.1f deg\n\n" %(angle_cut)
		TS_stat += Chain_codes
		TS_stat +="#Searching for backbone H-bonds:\n# frame   res(D)---res(A)  dist(nm) angle(deg)\n"
		TS_string = "%6.1f     %4d    %4d      %1.2f     %2.1f \n"
	else:
		TS_stat += "\n#Searching for contact atoms:\n# frame  atom1 --- atom2   dist (nm)\n"
		TS_string = "%6.1f     %4d    %4d      %1.2f\n"	
	for entry in TS_data:
		TS_stat += TS_string % tuple(entry)	
	Output += TS_stat
			
	Output += "\n"+Footer
	return Output

#Main function for script execution:
def Prep_Main(Args):
#	set run parameters:
	IO_param, HB_param, failmark,verbosity = Args
	pdb_file, out_file, sum_file, details, input_type, GROUPS = IO_param	
	Set_verb(Args[3])

	Main_Data = []
#	execute main code:

#	first set atom filters for collecting atom data:
	FILTER = []
	Get_Filters = Set_Filters(HB_param, input_type, GROUPS, failmark)
	FILTER = Get_Filters[1]
	HB_param = Get_Filters[0]
	mode, dist_cutoff, angle_cutoff, occ_cutoff, Hbond_atoms, MODIFIED = HB_param
	failmark = Get_Filters[2]

#	check if no errors occured so far:
	if pdb_file == "" or os.path.isfile(pdb_file) == False:
		failmark = 2
		Vprint(1, "Error: input file not found:",pdb_file)
	if failmark != 0:
		print(Usage)
		sys.exit(failmark)

#	Find Hbonds/Contacts in the trajectory:
	Main_Data = Do_traj(pdb_file,GROUPS,FILTER,HB_param)

#	Summarize data:
	Sum_Data = Hbond_Summary(Main_Data[1],Main_Data[0],mode)
	for entry in Sum_Data:
		Vprint(3,entry)
	if sum_file != "":
#		write summary
		Sum_param = [workdir, pdb_file, dist_cutoff, angle_cutoff, occ_cutoff, details, mode]
		Sum_output = Format_Sum(Sum_param, Sum_Data, Main_Data[0])
		s = open(sum_file,"wb")
		s.write(Sum_output.encode("ascii"))
		s.close()
				
#	Obtain renumbered time series:
	TS_Data = Hbond_TS(Main_Data[1],Main_Data[0])
	for entry in TS_Data:
		Vprint(4,entry)
	if out_file != "":
#		write time series:		
		TS_param = [workdir, pdb_file, dist_cutoff, angle_cutoff, mode]
		TS_output = Format_TS(TS_param,TS_Data,Main_Data[0])
		t = open(out_file,"wb")
		t.write(TS_output.encode("ascii"))
		t.close()	

	return [TS_Data,Main_Data[0]]
	



# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	IO_param, HB_param, failmark, verbosity = Custom_Args
	Set_verb(Custom_Args[3])
	Vprint(2, "\nRun parameters:\n", Custom_Args)

#       executing main code
	Data_main = Prep_Main(Custom_Args)
	Vprint(3, "\nMain data:")
	for entry in Data_main[1]:
		Vprint(3, entry)

#       print run-time message
	ftime = time.time()
	runtime = ftime-stime
	Out_files = ""
	for i in range(1,3):
		if IO_param[i] != "":
			Out_files += " %1s,"%IO_param[i]
	Out_files = Out_files[:-1]
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",Out_files)
else:
	Vprint(2,"HBSS preparation module (HbSS_prep.py)")

