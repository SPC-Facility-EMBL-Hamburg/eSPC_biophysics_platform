#!/usr/bin/env python
#

#License information:
#    DISICL dihedral angle calculation program for structure analysis.
#    Copyright (C) 2014 University of Natural Resources and Life Sciences, Vienna (Universitaet fuer Bodenkultur)
# 
#    Contact Information:
#    written by: Gabor Nagy
#    supervisor: Chris Oostenbrink
#    University of Natural Resources and Life Sciences (Universitaet fuer Bodenkultur)
#    Institute for Molecular Modeling and Simulation
#    Muthgasse 18,A-1190 Vienna, Austria
#    e-mails:
#    disicl@boku.ac.at
#    gabor.nagy@boku.ac.at
#    chris.oostenbrink@boku.ac.at
#    	
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License version 3
#    along with this program (see license.txt in the program directory).
#    If not, see <http://www.gnu.org/licenses/>.
#
#    If you found DISICL useful in your research, please site the following publications:
#
#    Nagy, G. & Oostenbrink, C. Dihedral-based segment identification and classification
#    of biopolymers I: Proteins. J. Chem. Inf. Model (2014). doi:10.1021/ci400541d
#
#    Nagy, G. & Oostenbrink, C. Dihedral-based segment identification and classification
#    of biopolymers II: Polynucleotides. J. Chem. Inf. Model (2014). doi:10.1021/ci400542n
#

import os
import sys
import math
workdir = os.getcwd()
stime = os.times()

print("DISICL_dihed program starts")
print(workdir+"\n")

#usage:
usage0 = '********************************************\n'
usage1 = 'DISICL dihedral analysis program\nUsage to prepare dihedral-angle time series:\n\n'
usage2 = 'DISICL_dihed.py <input file> <options>\n\n'
usage3 = 'Input file options:\n   <file-name>\n\t   Single trajectory file in pdb 1.0 compatible format (default)\n'
usage4 = '   @file <file-name>\n\t   Give argument file (no other arguments read from command line)\n\n'
usage5 = 'Additional options:\n   @res <res,res0-resmax>\n\t   Define residue range for analyisis (comma separated)\n\t'
usage6 = '   arguments are integers, or ranges (- separated)\n\t   if no range is given all valid residues are analysed\n\n'
usage7 = '   @fi <str,str,str,str>\n   @psi <str,str,str,str>\n   @gamma <str,str,str,str>\n   @xi <str,str,str,str>\n'
usage8 = '\t   Modify hard coded dehidral names,\n\t   strings are the atom identifiers (or a single 0 for none)\n\t   Alternative names can be given (/ separated)\n\n'
usage9 = '\t  Default dihedrals for proteins:\n\t     fi C(-1),N,CA,C\n\t     psi N,C,CA,N(+1)\n\t     gamma (none)\n\t     xi (none)\n'
usage10 = "\t  Default dihedrals for nucleotides:\n\t     fi (none)\n\t     psi C4',C3',O3',P(+1)\n\t     gamma C3',O3',P(+1),O5'(+1)\n\t     xi O4',C1',N1/N9,C2/C4\n\n"
usage11 = '   @type <t>\n\t   Set standard dihedral types:\n\t   1 - proteins (default)\n\t   2 - DNA (pdb naming)\n\t   3 - DNA (GROMOS naming)\n\t   0 - custom (no hardcoded vaules)\n\n'
usage12 = '   @time <time,timestep>\n\t   Define timestamps for trajectory (comma separated)\n\t'
usage13 = '   use floats for starting time and timestep (default 0.0,1.0)\n\n   @stride <n>\n\t   Calculate only every nth frame (default 1)\n\n'
usage14 = '   @write <w>\n\t   Write out results in different formats, w can take values:\n\t   1 - Ramachandran plot format\n\t'
usage15 = '   2 - DISICL time series (default)\n\t   3 - prepare density map (basic format)\n\t   4 - prepare density map (xmgrace heat maps)\n\n'
usage16 = '   @ext <.str>\n\t   Give input file format (if different):\n\t   g96,cnf,trc,trc.gz [GROMOS format]\n\t   pdb [ protein databank format] (default)\n\n'
usage17 = '   @path <path>\n\t   Set different source directory for structure files (default ./)\n\n   @debug\n\t   Use debug mode (write out calculation specs at every frame)\n\n'
usage18 = '   @lib <library file>\n\t   Read in angle definitions from custom library\n\t   @lib help for format information\n\n'
usage19 = '   @traj\n\t   Accelerate trajectory analysis with\n\t   (safety checks tunred off, list based atom parsing from frame 2) \n\n'
usage20 = '   @command\n\t   Flag marks the argument file identical with the input file\n\t   arguments are read until the first line starting with an "&"\n\t'
usage21 = '   For argument files ony, do not use @command in the command line!\n\n   @help\n\t   Print this message\n'


usage = usage0+usage1+usage2+usage3+usage4+usage5+usage6+usage7+usage8+usage9+usage10+usage11
usage = usage+usage12+usage13+usage14+usage15+usage16+usage17+usage18+usage19+usage20+usage21+usage0

# input definitions

path = ""
filename = ""
namebase = ""
ext = "pdb"

type = 1
fi = ["C","N","CA","C"]
#fi_dna = []
fi2 = [-1,0,0,0]

xi = []
#xi_dna2 = ["O4*","C1*","N1/N9","C2/C4"]

xi2 = [0,0,0,0]

psi = ["N","CA","C","N"]
#psi_dna2 = ["C4*","C3*","O3*","P"]
psi2 = [0,0,0,1]

gamma = []
#gamma_dna2  = ["C3*","O3*","P","O5*"]
gamma2 = [0,0,1,1]

ARGUMENTS = []
resdata = ''
RESDATA = []

TIME = [0,1]

library = ''
ANGLES = []

write = 2
debug = 0
stride = 1
traj = 0

#plotting parameters:
min1 = -180
max1 = 180
step1 = 36
min2 = -180
max2 = 180
step2 = 36
hstep = 10
emflag = 1

#if @file flag is provided read in arguments from there:
infofile = 0
argfile = ''
for arg in sys.argv:
    if arg.find('@file') != -1:
        infofile = 1
    elif infofile == 1:
        argfile = arg
        infofile = 0


if argfile == '':
    for arg in sys.argv:
        ARGUMENTS.append(arg)
elif os.path.isfile(argfile) == True:
	ARGUMENTS.append(argfile)
	i = open(argfile,"rb")
	for line0 in i:
		line = str(line0.decode("ascii"))
		if line.startswith("&"):
			break
		if not line.startswith('#'):
			parts = line.strip("\n").split()
			for part in parts:
				if not part == '':
					ARGUMENTS.append(part)
	i.close()
else:
	print("\nError: Argument file not found")
	failmark = 1
#command line interface:
flag = ''
flaglist = ['traj','debug','write','path','file','ext','stride','time','res','ang','type','lib','command','help']
failmark = 0
emark = 0
argnum = len(ARGUMENTS)-1
print('number of arguments: '+str(argnum))
print('flags found:')
cnt = 1
while cnt <= argnum:
	arg = ARGUMENTS[cnt]
	if arg == '@traj':
		traj = 1
		print("traj")
		flag = ''
	elif arg == '@debug':
		debug = 1
		print("debug")
		flag = ''
	elif arg == '@help':
		print("help")
		failmark = 1
		flag = ''
	elif arg == '@command':
		filename = ARGUMENTS[0]
	elif arg.startswith('@'):
		flag = arg.strip('@')
		print(flag)
		if not flag in flaglist:
			print('\nERROR: Unkown flag, program stops!')
			failmark = 1
			flag = ''
	elif flag == '' and cnt == 1:
		filename = arg
	elif flag == 'ext':
		ext = arg
		emark = 1
		flag = ''
	elif flag == 'res':
		resdata = arg
#		print resdata
		flag = ''
	elif flag == 'path':
		path = arg
		flag = ''
	elif flag == 'write':
		if arg == '1' or arg == '2' or arg == '3' or arg == '4':
			write = int(arg)
		else:
			print("\nError: wrong option for write: "+arg)
			failmark == 1
		flag = ''
	elif flag == 'fi':
		fi = arg.split(',')
		if fi == ['0']:
			fi = []
		elif not len(fi) == 4:
			print("\nError: need 4 atoms for dihedral calculation (fi)")
			failmark = 1
		flag = ''
	elif flag == 'psi':
		psi = arg.split(',')
		if psi == ['0']:
			psi = []
		elif not len(psi) == 4:
			print("nError: need 4 atoms for dihedral calculation (psi)")
			failmark = 1
		flag = ''
	elif flag == 'gamma':
		gamma = arg.split(',')
		if gamma == ['0']:
			gamma = []
		if not len(gamma) == 4:
			print("Error: need 4 atoms for dihedral calculation (gamma)")
		flag = ''
	elif flag == 'xi':
		xi = arg.split(',')
		if xi == ['0']:
			xi = []
		if not len(xi) == 4:
			print("\nError: need 4 atoms for dihedral calculation (xi)")
			failmark = 1
		flag = ''
	elif flag == 'time':
		try:
			t1, t2 = arg.split(',')
			t1 = float(t1)
			t2 = float(t2)
			TIME = [t1,t2]
		except ValueError:
			print("\nError: give start time and time step for @time flag (floats)")
			failmark = 1
		flag = ''
#	elif flag == 'write':
#		print 'not yet implemented'
	elif flag == 'type':
		try:
			type = int(arg)
			if not (type >= 0 and type <= 3):
				failmark = 1
				print("\nError: unknown type option: "+str(type))
		except ValueError:
			print("\nError: @type flag takes only integers")
			failmark = 1
		flag = ''
	elif flag == 'stride':
		try:
			stride = int(arg)
		except ValueError:
			print("\nError: @stide flag takes only integers")
			failmark = 1
		flag = ''
	elif flag == 'ang':
		ANGLES.append(arg)
		flag = ''
	elif flag == 'lib':
		library = arg
		flag = ''
	cnt = cnt + 1

#checking work directory and file path validity

zipmark = 0

if path == "":
	path = workdir
else:
	path1 = os.path.join(workdir,path)
	if os.path.isdir(path1) == True:
		path = path1
		print("\npath to work directory:\n"+path)
		os.chdir(path)
	elif os.path.isdir(path) == True:
		print("\npath to work directory:\n"+path)
		os.chdir(path)
	else:
		print("\nError: Could not locate path to work directory!")
		print(path)
		failmark = 1
		
#reading in data files
infile = filename
if filename == '':
	print('\nError: No input file was given program stops')
	failmark = 1
#	sys.exit(0)
elif os.path.isfile(infile) == False:
	parts = filename.split('.')
	namebase = parts[0]
	newinfile = os.path.join(workdir,infile)
	newinfile2 = os.path.join(path,infile)
	if os.path.isfile(newinfile) == True:
		print("\nInput file detected:",newinfile)
		fullpath = newinfile
	elif os.path.isfile(newinfile2) == True:
		print("Input file detected:",newinfile2)
		fullpath = newinfile2
	else:
		print("\nError: failed to locate input file!")
		print(infile)
		failmark = 1
else:
	print("\nInput file detected:",infile)
	parts = filename.split('.')
	namebase = parts[0]
	fullpath = infile
	
if failmark != 1:
	if parts[1] == 'cnf' or parts[1] == 'g96' or parts[1] == 'trc':
		ext == 'g96'
	elif parts[1] == 'trc.gz':
		zipmark = 1


if ext == 'pdb':
        print('using deafult pdb format')
elif ext == 'g96':
	print('using gromos format based on ext flag ')
	print('\nError: not yet implemented. Sorry for the trouble!')
	sys.exit(0)

#sorting residue information
data = ''
if not resdata == '':
	print("\n@res flag detected, residues found are:")
	data = resdata.split(',')
elif failmark != 1:
	print("\nNo @res flag was given! Attempting to read DISICL remarks from file:")
	i = open(fullpath,'rb')
	remarks = ''
	for line0 in i:
		line = str(line0.decode("ascii"))
		if line.startswith("DISICL"):
			remarks = line.strip('\n')
	if remarks != '':
		protlist = ''
		dnalist = ''
		try:
			print(remarks)
			parts = remarks.split()
			files = int(parts[1])
			frames = int(parts[2])
			resnumber = int(parts[3])
			protlist = parts[4].strip('[').strip(']')
			dnalist = parts[5].strip('[').strip(']')
			otherlist = parts[6].strip('[').strip(']')
		except ValueError:
			print("Error: Failed to read remark information! Program stops")
			failmark = 1
		if protlist != '' and type == 1:
			data = protlist.split(',')
			print("protein remark information found:",protlist)
		elif (type == 2 or type == 3) and dnalist != '':
			data = dnalist.split(',')
			print("dna remark information found:",dnalist)
		elif type == 0 and otherlist != '':
			data = otherlist.split(',')
			print("custom remark information found:",otherlist)
		else:
			print("Error: no applicable resides found in remarks, check @type!")
			failmark == 1
	else:
		print("    No remarks detected, program stops!")
		failmark = 1
		
reslist = ''
if not data == '':	
	for d in data:
		try:
			if reslist == '':
				reslist = d
			else:
				reslist = reslist+','+d
			if not d.find('-') == -1:
				dmin, dmax = d.split('-')
				dmin = int(dmin)
				dmax = int(dmax)
				cnt = dmin
				while cnt <= dmax:
					RESDATA.append(cnt)
					cnt = cnt + 1
			else:
				dnum = int(d)
				RESDATA.append(dnum)
		except ValueError:
			print("\nError: could not read in residue information ("+d+")")
	
#	print RESDATA
	print("\nResidue list updated:\n",reslist) 	
if not len(RESDATA) == 0:
	res0 = RESDATA[0]
	resmax = len(RESDATA)-1
	resn = len(RESDATA)
	resf = RESDATA[resmax]
	reslist = "["+reslist+"]"
	print("\nres0, resmax, resnum\n",res0,resf,(resmax+1))
else:
	print("\nError: no residues are given, program stops")
	failmark = 1

if TIME != [0,1]:
	print("\nTimestamps udated to:\ntime0= "+str(TIME[0])+" timestep= "+str(TIME[1]))

if type == 2:
		print('\nSetting dihedrals for default DNA analysis')
		fi = []
		psi = ["C4'","C3'","O3'","P"]
		gamma = ["C3'","O3'","P","O5'"]
		xi = ["O4'","C1'","N1/N9","C2/C4"]
if type == 3:
		print('\nSetting dihedrals for default GROMOS DNA analysis')
		fi = []
		psi = ["C4*","C3*","O3*","P"]
		gamma = ["C3*","O3*","P","O5*"]
		xi = ["O4*","C1*","N1/N9","C2/C4"]
custom_lib = 0
libpath = ''
if library != '':	
	libpath1 = os.path.join(workdir,library)
	libpath2 = os.path.join(path,library)
	if os.path.isfile(library) == True:
		print("\nLibrary file found:",library)
		libpath = library
	elif os.path.isfile(libpath1) == True:
		print("\nLibrary file found:",libpath1)
		libpath = libpath1
	elif os.path.isfile(libpath2) == True:
		print("\nLibrary file found:",libpath2)
		libpath = libpath2
	elif library == "help":
		failmark = 2
	else:
		print("\nError: Unable to locate library file!\n",libpath)
		failmark = 1

if libpath != '':
	l = open(libpath,'rb')
	for line0 in l:
		line = str(line0.decode("ascii"))
		if line.startswith("@ang"):
			parts =line.strip('\n').split()
			Angle = parts[1]
			ANGLES.append(Angle)
			custom_lib = 1


lib0 = "**********************************************************************\n"
lib1 = "Angle definition format for DISICL_dihed program:\n"
lib2 = "Give custom defined dihedral angles with:\n   @ang <name,atom1,atom2,atom3,atom4,mod1,mod2,mod3,mod4'>\n"
lib3 = "\t   name is a string used for the dihedral angle\n\n"
lib4 = "\t   atom1-atom4 are atom name strings\n\n\t   mod1-mod4 are residue number modifiers (integers)\n"

libinfo = lib0+lib1+lib2+lib3+lib4+lib0

CUSTOM_DIHEDRALS = []
for Angle in ANGLES:
	try:
		parts = Angle.split(',')
		name = parts[0]
		atom1 = parts[1]
		atom2 = parts[2]
		atom3 = parts[3]
		atom4 = parts[4]
		mod1 = int(parts[5])
		mod2 = int(parts[6])
		mod3 = int(parts[7])
		mod4 = int(parts[8])
		atoms = [atom1,atom2,atom3,atom4]
		mods = [mod1,mod2,mod3,mod4]
		dihedral = [name,atom1,atom2,atom3,atom4,mod1,mod2,mod3,mod4]
		CUSTOM_DIHEDRALS.append(dihedral)
		print('Angle '+name+' read in by atom atoms:',atoms,mods)
	except Exception:
		print('\nError: Angle definition could not be read in!\n',Angle)
		failmark = 2
if failmark == 1:
	print("\n"+usage)
	sys.exit(0)	
elif failmark == 2:
	print("\n"+usage)
	print("\n"+libinfo)
	sys.exit(1)
#mathematical function to calculate dihedral angle

#phi = arccos(|Ua*Ub|/(|Ua|*|Ub|)) = arcsin(|Ua X Ub|/(|Ua|*|Ub|))

#ohter words: AXB= a*b*sin(phi) A*B = a*b*cos(phi)
# sin(phi)= AXB/a*b cos(phi) = A*B/a*b
# A*B = Xa*Xb+Ya*Yb+Za*Zb
# AXB = [(Ya*Zb-Yb*Za),-(Xa*Zb-Xb*Za),(Xa*Yb-Xb*Ya)]
# |A| = a = (Xa^2+Ya^2+Za^2)^0.5

def dihedral(COORD):
#coordinate definitions for atoms
	try:
		X1 = float(COORD[0])
		Y1 = float(COORD[1])
		Z1 = float(COORD[2])
		X2 = float(COORD[3])
		Y2 = float(COORD[4])
		Z2 = float(COORD[5])
		X3 = float(COORD[6])
		Y3 = float(COORD[7])
		Z3 = float(COORD[8])
		X4 = float(COORD[9])
		Y4 = float(COORD[10])
		Z4 = float(COORD[11])
#		print ("atoms detected:")	
#		print ("atom1("+str(X1)+","+str(Y1)+","+str(Z1)+")")
#		print ("atom2("+str(X2)+","+str(Y2)+","+str(Z2)+")")
#		print ("atom3("+str(X3)+","+str(Y3)+","+str(Z3)+")")
#		print ("atom4("+str(X4)+","+str(Y4)+","+str(Z4)+")\n")
#vector definitions for bonds
		V1 = [(X2-X1),(Y2-Y1),(Z2-Z1)]
		V2 = [(X3-X2),(Y3-Y2),(Z3-Z2)]
		V3 = [(X4-X3),(Y4-Y3),(Z4-Z3)]
#		print "vectors defined as :\nV1("+str(V1)+")\nV2("+str(V2)+")\nV3("+str(V3)+")\n"
#normal vector definition for 2 planes
#UA = V1 X V2
#UB = V2 X V3
		UA = [(V2[1]*V1[2]-V2[2]*V1[1]),-(V2[0]*V1[2]-V2[2]*V1[0]),(V2[0]*V1[1]-V2[1]*V1[0])]
		UB = [(V3[1]*V2[2]-V3[2]*V2[1]),-(V3[0]*V2[2]-V3[2]*V2[0]),(V3[0]*V2[1]-V3[1]*V2[0])]
		
#		print "plane norms defined as :\nUA("+str(UA)+")\nUB("+str(UB)+")\n"
#needed units for angle calculation
#	       	AXB = [(UA[1]*UB[2]-UA[2]*UB[1]),-(UA[0]*UB[2]-UA[2]*UB[0]),(UA[0]*UB[1]-UA[1]*UB[0])]
		AB = (UA[0]*UB[0]+UA[1]*UB[1]+ UA[2]*UB[2])

#sign = used for the convention equals the sign of V1*UB 
		signer = (V1[0]*UB[0]+V1[1]*UB[1]+ V1[2]*UB[2])
		if signer >= 0:
			sign = -1
		else:
			sign = 1
		if AB >= 0:
			ABS = AB
		else:
			ABS = -1*AB		 		
#		LX = math.sqrt(AXB[0]*AXB[0]+AXB[1]*AXB[1]+AXB[2]*AXB[2])
		LA = math.sqrt(UA[0]*UA[0]+UA[1]*UA[1]+UA[2]*UA[2])
		LB = math.sqrt(UB[0]*UB[0]+UB[1]*UB[1]+UB[2]*UB[2])
#		print "Calculating dihedral angle:"
#		print	"AXB vectors defined as :\nAxB("+str(AXB)+")"

#		print "A*B skalar = "+str(AB)+", sign skalar V1*B = "+str(signer)
#		print "vectors vector lenghts :\nUA("+str(LA)+")\nUB("+str(LB)
#		print "AXB length:" +str(LX)+"\n"
		if LA == 0 or LB == 0:
			print("Zero denominator error, value set to 0")
			Angle = 0
		else:
			Dihed1 = math.acos(AB/(LA*LB))
#			Dihed2 = math.asin(LX/(LA*LB))
			Angle = math.degrees(Dihed1*sign)
#			print "dihedal in radians: "+ str(Dihed1*sign)
#			print str(Dihed2)1f = ""

	except ValueError:
			print("Input cannot be read in value set to 0")
			print("Coordinates  taken "+str(COORD))
			Angle = 0
	
	return Angle

#reading in pdb informations for dihedral calculation:

#class dfined for pdb code
class Inputs:
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

# reading in file contents
	

outfile = namebase+".out"

i = open(fullpath,'rb')
o = open(outfile,'wb')

if zipmark == 1:
	print("decompressing input file")
	os.system('gunzip '+fullpath)
	infile = infile.strip('.gz')


ATOM = []
FULL = []
INPS = []
MOD = []
LIST = []
D = []
#adding hardcoded dihedrals if no library is used
if custom_lib == 0 and type >= 1 and type <= 3:
	print("\nNo library definition is specified using hardcoded libararies")
	FI = []
	GAMMA = []
	PSI = []
	XI= []
	INPS = [fi,psi,gamma,xi]
	MOD = [fi2,psi2,gamma2,xi2]
	LIST = [FI,PSI,GAMMA,XI]
	if type == 1:
		D = ["PHI ","PSI ","GAMMA","XI  "]
	if type >= 2:
		D = ["PHI","EPS","ZET","CHI"]
#adding custom dihedrals:
if CUSTOM_DIHEDRALS != '':
	for Dih in CUSTOM_DIHEDRALS:
		length = len(Dih[0])
		if length >= 5:
			name = Dih[0]
		else:
			name = Dih[0]+(5-length)*' '
		D.append(name)
		Atoms = [Dih[1],Dih[2],Dih[3],Dih[4]]
		INPS.append(Atoms)
		Mods = [Dih[5],Dih[6],Dih[7],Dih[8]]
		MOD.append(Mods)
		LIST.append([])
	
inum = len(LIST)
if inum == 0:
	print("\nError: no available dihedrals are detected! Program stops")
	print(usage)
	sys.exit()
elif inum > 2 and write >= 3:
	print("\nWarning: more than two dihedrals are defined, but output format supports only two!")
	print("Program continues, first two dihedrals will be used to create deinsity maps!")
	
	
tcnt = 0
step = (stride-1)

vec = 0
torsnum = 0
while vec < inum:
	if not INPS[vec] == []:
		print(D[vec]+" dihedral defined by"+str(INPS[vec])+","+str(MOD[vec]))
		torsnum = torsnum + 1
	vec = vec + 1

if debug == 1:
	print("\nDebug mode is on, detailed messages are printed")

#reading in atomic information from input
print("\nStarting Dihedral angle calculations:")
for line0 in i: 
	line = str(line0.decode("ascii"))
	if line.startswith('ATOM') or line.startswith('HETATM'):
		line2 = line.strip().strip('\n')
		try:
			anum = int(line2[6:11].strip(' '))
			aID = str(line2[11:16].strip(' '))
			aRI = str(line2[16:20].strip(' '))
			aCI = str(line2[21])
			aRN = int(line2[22:27].strip(' '))
			ax = float(line2[30:38].strip(' '))
			ay = float(line2[38:46].strip(' '))
			az = float(line2[46:54].strip(' '))
#			aen = float(line2[55:60].strip(' '))
#			abf = float(line2[61:66].strip(' '))
			aai = str(line2[76:78].strip('\n'))
			if aai.isalpha() == False:
				if line2[12].isalpha() == True:
					aai = str(line2[12:14])
				else:
					aai = str(line2[13:14])
#			acharge = float(line2[79:82].strip())
			atime = TIME[0] + tcnt*TIME[1]
#			print( str(anum)+"\t"+aID+"\t"+aRI+"\t"+aCI+"\t"+str(aRN)+"\t"+str(ax)\
#			+"\t"+str(ay)+"\t"+str(az)+"\t"+str(aen)+"\t"+str(abf)+"\t"+aai)

			data = [anum, aID, ax, ay, az, atime, aai, aRI, aCI, aRN]

			Atom = Inputs(data)
			ATOM.append(Atom)
		except ValueError:
			print("line could not be read in:\n"+line2)
			print(anum, aID, ax, ay, az, atime, aai, aRI, aCI, aRN)
#colecting atomic infromation at the end of each frame:
	if line.startswith('END'):
		tcnt = tcnt + 1
		step = step + 1
		print("frame "+str(tcnt)+" found")
		cnt = 0
		Atoms = []
		TMP = []
		vec = 0
		inpvec = []
		dcnt = []
		if step == stride:
			while vec < inum:
#			if not INPS[vec] == []:
				inpvec.append(0)
				avec = ["","","",""]
				Atoms.append(avec)
				TMP.append([])
				vec = vec + 1
			inpvec0 = inpvec
			secnum = len(inpvec)
#atom selection for dihedrals
			icnt = 0
			while icnt < inum:
				inp = INPS[icnt]
				if not inp == []:
#defining atom identifiers for parsing:					
					if (traj == 1 and tcnt == 1) or (traj == 0):
						xcnt = 0
						labels = []
						while xcnt < 4:
							if inp[xcnt].find("/") == -1:
								labels.append([inp[xcnt]])
							else:
								parts = inp[xcnt].split("/")
								labels.append(parts)
							xcnt = xcnt +1
					cnt = 0
#				print str(labels)
					while cnt <= resmax:
						A11 = ""
						A12 = ""
						A13 = ""
						A14 = ""
						num11 = ""
						num12 = ""
						num13 = ""
						num14 = ""
						vec = 0
						dcnt = []
						check = 0
						dcnt = [0,0,0,0]
						while vec < inum:
							vec =vec + 1
#parsing pdb and seariching appropriate atoms:
						if (traj == 1 and tcnt == 1) or (traj == 0):	
							acnt = 0	
							for A in ATOM:
								seqnum = (RESDATA[cnt]+MOD[icnt][0])
								for label in labels[0]:
									if (A.ID == label and A.rn == seqnum and dcnt[0] == 0):
										A11 = A
										num11 = acnt
										dcnt[0] = dcnt[0] + 1
	#								print str(icnt+1),"1",str(A.num),str(A.ri),str(A.rn),str(A.ID)
	#								print label,seqnum
										
								seqnum = (RESDATA[cnt]+MOD[icnt][1])
								for label in labels[1]:
									if (A.ID == label and A.rn == seqnum and dcnt[1] == 0):
										A12 = A
										num12 = acnt
										dcnt[1] = dcnt[1] + 1
	#								print str(icnt+1),"2",str(A.num),str(A.ri),str(A.rn),str(A.ID)
	#								print label,seqnum
										
								seqnum = (RESDATA[cnt]+MOD[icnt][2])
								for label in labels[2]:
									if (A.ID == label and A.rn == seqnum and dcnt[2] == 0):
										A13 = A
										num13 = acnt
										dcnt[2] = dcnt[2] + 1
	#								print str(icnt+1),"3",str(A.num),str(A.ri),str(A.rn),str(A.ID)
	#								print label,seqnum
										
								seqnum = (RESDATA[cnt]+MOD[icnt][3])
								for label in labels[3]:
									if (A.ID == label and A.rn == seqnum and dcnt[3] == 0):
										A14 = A
										num14 = acnt
										dcnt[3] = dcnt[3] + 1
	#								print str(icnt+1),"4",str(A.num),str(A.ri),str(A.rn),str(A.ID)
	#								print label,seqnum
								acnt = acnt + 1
#using obtained atomlists (if available):
						else:
							num1 = LIST[icnt][cnt][0]
							num2 = LIST[icnt][cnt][1]
							num3 = LIST[icnt][cnt][2]
							num4 = LIST[icnt][cnt][3]
							if not num1 == "X":
								A11 = ATOM[num1]
								dcnt[0] = 1
								A12 = ATOM[num2]
								dcnt[1] = 1
								A13 = ATOM[num3]
								dcnt[2] = 1
								A14 = ATOM[num4]
								dcnt[3] = 1
#calculate dihedral					
						check = dcnt[0]*dcnt[1]*dcnt[2]*dcnt[3]
						if check == 1:
							COORD = [A11.x,A11.y,A11.z,A12.x,A12.y,A12.z,A13.x,A13.y,A13.z,A14.x,A14.y,A14.z]
							Angle = dihedral(COORD)
							inpvec[icnt] = inpvec[icnt] + 1
#check residue definitions on first frame or in debug mode:
							if (debug == 0 and tcnt == 1) or (debug == 1):
								infin = (D[icnt],RESDATA[cnt],Angle)
								print("%5s for residue %-3d is %3.1f deg, defined by atoms:" % infin)
								print(A11.num,A12.num,A13.num,A14.num)
								print(A11.num,A11.ri,A11.rn,A11.ID)
								print(A12.num,A12.ri,A12.rn,A12.ID)
								print(A13.num,A13.ri,A13.rn,A13.ID)
								print(A14.num,A14.ri,A14.rn,A14.ID)
							COORD = []							
							if (traj == 1 and tcnt == 1):
								dihnums = [num11,num12,num13,num14]
								
#handle missing atoms or exceptions
						elif check == 0:
							Angle = "X"
							dihnums = ["X","X","X","X"]
							infin = (D[icnt],RESDATA[cnt])
							print("Warning! %5s for residue %-3d is missing at least one atom, dihedral left out!" % infin)
						else:
							Angle = "X"
							dihnums = ["X","X","X","X"]
							infin = (D[icnt],RESDATA[cnt])
							print("Warning! %5s for residue %-3d has than 4 atom meeting the search criterion, dihedral left out!" % infin)
						TMP[icnt].append(Angle)
						if (tcnt == 1 and traj == 1):
							LIST[icnt].append(dihnums)
						dihnums = []
						cnt = cnt + 1
				icnt = icnt + 1
# collecting final data
			cnt = 0		
			while cnt <= resmax:
				time = TIME[0]+(tcnt-1)*TIME[1]
				resnum = RESDATA[cnt]
				diheds = [resnum,time]
				numbers = [tcnt]
				icnt = 0
				while icnt < inum:
					if INPS[icnt] != []:
						diheds.append(TMP[icnt][cnt])
						numbers.append(inpvec[icnt])
						numbers.append(D[icnt])
					icnt = icnt + 1
				if not "X" in diheds:
					FULL.append(diheds)
				cnt = cnt + 1
			numbers.append(resn)
			infofr = tuple(numbers)
#			print infofr
			stringfr = "frame %-3d processed with "+torsnum*"%3d %3s,"+" for requested %-3d residues"
			print(stringfr % infofr)
			step = 0
		ATOM = []
		
		
Fnumber = tcnt
print("file read in, "+str(Fnumber)+" frame(s) found\n")

	
#writing out dihedrals:
o.write(('#Ramachandran plot based on '+infile+'\n').encode("ascii"))
icnt = 0
List = []
while icnt < inum:
	if not INPS[icnt] == []:
		List.append(icnt)
	icnt = icnt + 1
torsnum = len(List)

if write == 1:
	rcnt = 0
	while rcnt <= resmax:
		header = '#Ramachandran plot for residue '+str(RESDATA[rcnt])+'\n#'
		for tors in List:
			header = header +'   '+D[tors]+'(deg)'
		header = header + '\n'
		o.write(header.encode("ascii"))
		for Res in FULL:
			if Res[0] == RESDATA[rcnt]:
				In = []
				cnt = 0
				while cnt < torsnum:
					In.append(Res[(cnt+2)])
					cnt = cnt + 1
				In = tuple(In)
				form = torsnum*'%10.2f   '+'\n'
				o.write((form%In).encode("ascii"))
		o.write("&\n".encode("ascii"))
		rcnt = rcnt + 1
elif write == 2:
	tcnt = 0
	Dinfo1 = 'DISICL'+3*' '+'time0'+3*' '+'timestep'+3*' '+'frames'+3*' '+'residues'+3*' '+'residue list\n'
	Dinfo2 = 'DISICL'+4*' '+str(TIME[0])+6*' '+str(TIME[1]*stride)+7*' '+str(int(math.floor(float(Fnumber)/stride)))+8*' '+str(resn)+7*' '+str(reslist)+'\n'
	Dinfo = '#'+Dinfo1+'#'+Dinfo2
	print(Dinfo1,Dinfo2)
	header = '#DISICL time series for dihedral angles\n'+Dinfo+'#time  Residue'
	for tors in List:
		header = header+('     '+D[tors]+'(deg)')
	header = header+"\n"
	o.write(header.encode("ascii"))
	while tcnt <= Fnumber:
		time = TIME[0] + tcnt*TIME[1]
		for Res in FULL:
			if Res[1] == time: 
				cnt = 0
				In1 = [Res[1],Res[0]]
				In2 = []
				while cnt < torsnum:
					In1.append(Res[(cnt+2)])
					cnt = cnt + 1
				In1 = tuple(In1)
				fstring1 = '%6.6f %4d   '+torsnum*'%10.2f   '+'\n'				
				o.write((fstring1%In1).encode("ascii"))
				cnt = cnt + 1
#		o.write('&\n')
		tcnt = tcnt + 1
if write == 3 or write == 4:
	rcnt = 0
	while rcnt <= resmax:
#making the grid
		o.write(('#Ramachandran density map for residue '+str(RESDATA[rcnt])+'\n').encode("ascii"))
		o.write(('#grid number	density (%)\n').encode("ascii"))
		GRID = []
		size1 = float((max1-min1)/step1)
		size2 = float((max2-min2)/step2)
		print('grid parameters:\nmin,max,n step, stepsize')
		print(str(min1),str(max1),str(step1),str(size1))
		print(str(min2),str(max2),str(step2),str(size2))
		print('number of gridpoints:'+str(step1*step2))
		n = 0
		print('making density map for '+str(RESDATA[rcnt]))
		while n < step1:
			m = 0
			while m < step2:
				binum = n*step2+m
				GRID.append(0)
				m = m + 1
			n = n + 1
#filling the grid
		pcnt = 0
		for Res in FULL:
			if Res[0] == RESDATA[rcnt]:
				fi = Res[2]
				psi = Res[3]
				print('('+str(fi)+','+str(psi)+')')
				a = int(math.floor((fi-min1)/size1))
				b = int(math.floor((psi-min2)/size2))
				print('truncated ('+str(a)+','+str(b)+')')
				gridnum = int(a*step2+b)
				print('bin:'+str(gridnum))
				GRID[gridnum] = (GRID[gridnum] + 1)
				pcnt = pcnt + 1
		pointnum = pcnt
		o.write(("#number of points: "+str(pointnum)+"\n").encode("ascii"))
#plotting the grid
		if write == 3:
			print("using output method 3")
			cnt = 0
			for bin in GRID:
				posx = math.floor(cnt/step2)
				posy = cnt-posx*step2
				X = min1+(posx*size1)
				Y = min2+(posy*size2)
				dat = (float(bin)/pointnum*100)
#				print str(dat)
				In = (X, Y, dat)
				form = '%10.1f   %10.1f   %10.2f\n'
				o.write((form%In).encode("ascii"))
				cnt = cnt + 1
			o.write("&\n".encode("ascii"))
			o.close()
			os.system('mv '+outfile+' dih_'+str(RESDATA[rcnt])+'.map\n' )
			o = open(outfile,'wb')
			rcnt = rcnt + 1
		if write == 4:
			print("using output method 4")
			binmin = pointnum
			binmax = 0
			for bin in GRID:
#				print bin
				bin = int(bin)
				if bin > binmax:
					binmax = bin
				if bin < binmin and bin != 0:
					binmin = bin
#			print(binmin,binmax)
			binlength = float(binmax-binmin)
			size3 = float(binlength/hstep)
#			print str(size3)
			bcnt = 0
			while bcnt <= hstep:
				cnt = 0
				scnt = 0
				down = float(binmin + (bcnt)*size3)
				up = float(binmin + (bcnt+1)*size3)
#				print str(hstep),str(up),str(down)
#				o.write("#bins with "+str(down/pointnum)+"-"+str(up/pointnum)+" abundance\n")
				o.write(("#bins with "+str(down)+"-"+str(up)+" abundance\n").encode("ascii"))
				for bin in GRID:
					bin = int(bin)
					if bin >= down and bin < up:
						posx = math.floor(cnt/step2)
						posy = cnt-posx*step2
						X = min1+(posx*size1)
						Y = min2+(posy*size2)
						In4 = (X, Y)	
						form4 = '%10.1f   %10.1f\n'
						o.write((form4%In4).encode("ascii"))
						scnt = scnt + 1
					cnt = cnt + 1
				if scnt == 0:
					X = 180
					Y = -180
					In4 = (X, Y)	
					form4 = '%10.1f   %10.1f\n'
					o.write((form4 % In4).encode("ascii"))	
				o.write("&\n".encode("ascii"))
				bcnt = bcnt + 1
				cnt = 0
			if emflag == 1:	
				o.write("#empty bins\n".encode("ascii"))	
				for bin in GRID:
					if bin < binmin:
						posx = math.floor(cnt/step2)
						posy = cnt-posx*step2
						X = min1+(posx*size1)
						Y = min2+(posy*size2)
						In4 = (X, Y)	
						form4 = '%10.1f   %10.1f\n'
						o.write(form4 % In4)
					cnt = cnt + 1
				o.write("&\n".encode("ascii"))
			o.close()
			os.system('mv '+outfile+' dih_'+str(RESDATA[rcnt])+'.map\n' )
			o = open(outfile,'wb')
			
# output formatting for density map
i.close()
o.close()

typemark = ''	
if type == 0:
	typemark = "custom"
if type == 1:
	typemark = "protein"
if type == 2 or type == 3:
	typemark = "nucleotide"
rtime = os.times()
#print rtime, sttime
runtime = rtime[4]-stime[4]
if runtime == 0.0:
	runtime = rtime[0]+rtime[1]-(stime[0]+stime[1])
print("program runtime: %1.2f  seconds (%1.2f h)" % (runtime,(float(runtime)/3600)))
if write == 1 or write == 2:
	print("Results are written in: "+outfile)
print("Dihedral program ("+typemark+") finished successfully! ("+infile+")")
