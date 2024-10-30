#!/usr/bin/env python
#

#License information:
#    DISICL visualization script to show structure classifications on the 3D models they are based on.
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


import sys
import os
import math
import __main__
import threading
import time

#usage
usage0 = '**********************************************************************\n'
usage1 = 'DISICL visualization module\n'
usage2 = 'Usage for visualization script:\n\nDISICL_vis.py [input file] [options]\n\n'
usage3 = 'Input file options:\n   <name>\n\t   Filename of the DISICL input file  (without extension)\n\n'
usage4 = '   @file <file-name>\n\t   Give argument file (no other arguments read from command line)\n\n'
usage5 = '   @pdb <file1> @class <file2> [<file3> ...]\n\t   Specify pdb trajectory file (file1) and \n\t'
usage6 = '   DISICL classification file(s) (file2, file3, ...)\n\t   (extensions are required)\n'
usage7 = '\nAdditional options:\n'
usage8 = '   @type <t>\n\t   Set visualized molecule types, t takes the value of:\n\t   1 - proteins (default)\n\t'
usage9 = '   2 - nucleotides (DNA/RNA)\n\t   3 - proteins and nucleotides\n\t   0 - custom (no hardcoded parameters)\n\n'
usage10 = '   @dest <path>\n\t   Set different goal directory for results (relative, default ./)\n\n'
usage11 = '   @time <time,timestep>\n\t   Define timestamps for classification file (comma separated)\n\n'
usage12 = '   @stride <n>\n\t   Use only every nth frame from the .pdb file (default 1)\n\n'
usage13 = '   @write <w>\n\t   Write out results in different formats, w can take values:\n\t   1 - PyMol session (default)\n\t'
usage14 = '   2 - PyMol session + .png images (1/frame)\n\t   3 - Animated PyMol session\n\n\t'
usage15 = '   Option 3 requires high amounts of memory and time,\n\t   recommended for small molecules and trajectories only!\n\n'
usage16 = '   @col <name,color,type> [<name,color,type> ...]\n\t   Provide custom color definition(s) (see @lib help)\n\n'
usage17 = '   @lib <library file1> [<library file2> ...]\n\t   Read in coloring schemes from library file(s)\n\t'
usage18 = '   Use @lib help for further info on library definitions\n\n   @help\n\t   Print this message\n\n'
usage19 = '   @debug\n\t   Use debug mode (write out colored residues)\n\n'
usage20 = 'WARNING: DISICL_vis.py requires PyMol, and the env. variable PYMOL_PATH !!!\n\n'
usage21 = 'NOTE: DNA/RNA coloring is shifted +1 residue\n      to better represent the nucleotide segments \n\n'
usage22 = 'NOTE: Unclassfied central residues that belong to a classified segment\n      are added automatically for proteins and nucleotides\n'

usage = usage0+usage1+usage2+usage3+usage4+usage5+usage6+usage7+usage8+usage9+usage10+usage11+usage12
usage = usage+usage13+usage14+usage15+usage16+usage17+usage18+usage19+usage20+usage21+usage22+usage0

lib0 = "**********************************************************************\n"
lib1 = "Format specifications for the DISICL_vis color schemes:\n\n"
lib2 = "Define a color definition using:\n   @col <name,color,type>\n\t"
lib3 = "   name is the collector code used for strucutre element\n\t   (put beween brackets in the class header)\n\n\t"
lib4 = "   color is the PyMol reference name for the color associated\n\t   with the structure (see PyMol for details)\n\n\t"
lib5 = "   type is the molecule type code\n\t   1 - protein\n\t   2 - nucleotide\n\t   0 - custom\n"

libinfo = lib0+lib1+lib2+lib3+lib4+lib5+lib0


#setting basic variables
ARGUMENTS = []
namebase = ''
workdir = os.getcwd()
Pymol_Path = ''
stime = os.times()
write = 1
CUSTOM_INFO = []
mode = 1
path = ''
infile = ''
namebase = ''
LIBS = []
timecodes = ''
framenumber = 1
TIME = [0.0,1.0]
COLDEF = []
Dest = ''
stride = 1
debug = 0

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
else:
	ARGUMENTS.append(argfile)
	if os.path.isfile(argfile) == True:
		i = open(argfile,"rb")
		for line in i:
			if line.startswith("&"):
				break
			if not line.startswith('#'):
				parts = line.strip("\n").split()
				for part in parts:
					if not part == '':
						ARGUMENTS.append(part)
		i.close()
	else:
		print "\nERROR: Argument file not found!"
		failmark = 1
#command line interface:

failmark = 0
argnum = len(ARGUMENTS)-1
print('number of arguments: '+str(argnum))
print 'flags found:'
flaglist = ['lib','dest','file','time','type','command','help','name','pdb','class','write','stride','col','debug']
flag = ''
cnt = 1
if argnum > 0:
	while cnt <= argnum:
		arg = ARGUMENTS[cnt]
		if arg == "@command":
			infile = ARGUMENTS[0]
		elif arg == "@help":
			failmark = 1
			print "help"
			flag = ''
		elif arg == "@debug":
			debug = 1
			print "debug"
			flag = ''
		elif arg.startswith("@"):
			flag = arg.strip('@')
			print flag
			if not flag in flaglist:
				print '\nERROR: Unkown flag, program stops!'
				failmark = 1
				flag = ''
		elif (flag == '' and cnt == 1) or flag == 'name':
			namebase = arg
			flag = ''
		elif flag == 'pdb':
			infile = arg
			flag = ''
		elif flag == 'class':
			CUSTOM_INFO.append(arg)
		elif flag == 'lib':
		    LIBS.append(arg)
		elif flag == 'col':
		    COLDEF.append(arg)    
		elif flag == 'time':
			timecodes = arg
			flag = ''
		elif flag == 'dest':
			Dest = arg
			flag = ''
		elif flag == 'type':
		    try:
			mode = int(arg)
			if not (mode >= 0 and mode <= 3):
				print "\nError: type not recognized!"
				failmark = 1
		    except Exception:
			print "\nError: type not recognized!"
			failmark = 1
		    flag = ''
		elif flag == 'stride':
		    try:
			stride = int(arg)
		    except Exception:
			print "\nError: The @stride flag takes integers only!"
			failmark = 1
		    flag = ''
		elif flag == 'write':
		    try:
			write = int(arg)
		    except Exception:
			print "\nError: format type not recognized!"
			failmark = 1
		    flag = ''
		cnt = cnt + 1
else:
	print "Error: no arguments are given, program stops!"
	failmark = 1

Pymol_Path = os.getenv("PYMOL_PATH")

if (Pymol_Path == '' or Pymol_Path == None):
    print "\nERROR:PyMol environmental variable PYMOL_PATH not found, script stops!"
    failmark = 1
else:
    print "\nPyMol is found in:",Pymol_Path
    sys.path.append(Pymol_Path)
    import pymol
    from pymol import cmd

#reading in color definitions form library files:
if LIBS != [] and LIBS != ['help']:
    for libfile in LIBS:
	lib = ''
	fullpath = os.path.join(workdir,libfile)
	if os.path.isfile(libfile) == True:
	    lib = libfile
	elif os.path.isfile(fullpath) == True:
	    lib = fullpath
	if lib != '':
	    l = open(lib,"rb")
	    print "\nReading in coloring information from file:"
	    print lib
	    for line in l:
		if not line.startswith("#") or not line.strip() == "":
			    line2 = line.strip().strip("\n")
			    if line2.startswith("@col"):
				    line3 = line2.strip("@col").strip()
				    COLDEF.append(line3)
elif LIBS == ["help"]:
    print usage
    print libinfo
    sys.exit()


#reading in coloring information
COLOR_DNA = []
COLOR_PROT = []
COLOR_CUST = []

if COLDEF != []:
    print "\nDefining Color schemes:"
    for entry in COLDEF:
	parts = entry.split(",")
	partnum = len(parts)
	if partnum == 3:
	    try:
		clname = parts[0]
		color = parts[1]
		cltype = int(parts[2])
		data = clname+","+color
		if cltype == 0:
		    COLOR_CUST.append(data)
		    print "custom class:",clname," is colored as:",color
    		elif cltype == 1:
		    COLOR_PROT.append(data)
		    print "protein class:",clname," is colored as:",color
		elif cltype == 2:
		    COLOR_DNA.append(data)
		    print "nucleotide class:",clname," is colored as:",color
	    except Exception:
		print "Color information entry:",entry
		print "Could not be read in"
elif failmark == 1:
    pass
else:
    print "No coloring scheme was defined, using hard-coded schemes!"
    COLORS1 = ["BI,green","BII,yellow","BIII,limon","AH,red","ZH,skyblue","AB,deepolive","AB2,sand","AZ,magenta","QL,gray40"]
    COLORS2 = ["ZB,purpleblue","AD,hotpink","BD,slate","ZD,lightblue","AL,orange","BL,teal","TL,yelloworange","ST,cyan"]
    COLORS3 = ["BH,green","IB,limon","IA,yelloworange","TR,sand","TA,sand"]
    COLORS4 = []
    COLORS5 = []
    C1 = [COLORS1,COLORS2,COLORS3,COLORS4,COLORS5]
    print "Hard-coded nucleotide colors:"
    for cols in C1:
	    for color in cols:
		    COLOR_DNA.append(color)
		    print color
    
    
    CPROT1 = ["ALH,forest","3H,limon","PIH,cyan","EBS,red","NBS,orange","PP,sand","TI,lime","TII,magenta","TVIII,hotpink","GXT,slate"]
    CPROT2 = ["SCH,marine","HP,yelloworange","HC,gray40","BC,deepolive","TC,teal","BU,wheat","LHII,purpleblue","LHH,density"]
    CPROT3 = ["HEL,green","3HT,cyan","BS,red","IRB,sand","BT,magenta","OTT,teal","LHT,purpleblue"]
    C2 = [CPROT1,CPROT2,CPROT3]
    print "Hard-coded protein colors:"
    for cols in C2:
	    for color in cols:
		    COLOR_PROT.append(color)
		    print color

	
#looking up classification files:

protinfo  = ["DISICL_p",".out"]
dnainfo = ["DISICL_n",".out"]
custinfo = ["DISICL_c",".out"]
INFOS = [protinfo,dnainfo,custinfo]
INFOFILES = []
if CUSTOM_INFO == []:
	FILES = os.listdir(workdir)
	for File in FILES:
		for info in INFOS:
			if File.find(info[0]) != -1 and File.find(info[1]) != -1 and File.find(namebase) != -1:
				INFOFILES.append(File)
else:
	for File in CUSTOM_INFO:
		INFOFILES.append(File)
	
if timecodes != '':
    parts = timecodes.split(',')
    t0 = float(parts[0])
    ts = float(parts[1])
    TIME = [t0,ts]
    print "Timestamps reset to : ",TIME
	
# adding residues to classes and getting the color for them
print "\nReading in info files:\n"
INFLIST = []
if INFOFILES != []:
	pinfo = 0
	dinfo = 0
	cinfo = 0
	for i in  range(len(INFOFILES)):
	    File = INFOFILES[i]
	    print File,i
	    if CUSTOM_INFO == []:
		if pinfo == 0 and File.find(protinfo[0]) != -1 and (mode == 1 or mode == 3):
		    pinfo = 1
		    INFLIST.append(i)
		    print File,"added as prot." 
		if dinfo == 0 and File.find(dnainfo[0]) != -1 and (mode == 2 or mode == 3):
		    dinfo = 1
		    INFLIST.append(i)
		    print File,"added as nucl."
		if cinfo == 0 and File.find(custinfo[0]) != -1 and mode == 0:
		    cinfo = 1
		    INFLIST.append(i)
		    print File,"added as cust."
	    else:
		if os.path.isfile(File) == True:
		    INFLIST.append(i)
		    print "User defined classification added:",File
		else:
		    print "Warning: User defined classification not found:",File
if INFLIST == []:
    print "\nERROR: No valid classification file was detected, script stops!"
    failmark = 1
    
if infile == '' and namebase != '':
	infile = namebase+"_DISICL.pdb"
if os.path.isfile(infile) == True:
    print "\nCoordniate trajectory file found:",infile
else:
    print "\nERROR: No valid coordinate file was detected, script stops!"
    failmark = 1

if failmark == 1:
	print usage
	sys.exit()
	
DATA = []
DATA2 = []
DATA3 = []
remark = ''

for i in INFLIST:    
    infofile = INFOFILES[i]
    info = open(infofile,'rb')
    TMP = []
    print "classification file:",infofile
    for line in info:
	if line.startswith("#DISICL") and line.find("time") == -1:
	    remark = line
	    
    if remark != '':
	try:
	    parts = remark.strip("\n").split()
	    t0 = float(parts[1])
	    ts = float(parts[2])
	    fn = int(parts[3])
	    seg = int(parts[4])
	    Rlist = parts[5]
	    print "\nDISICL remarks read in:\nnumber of expected frames:",fn,"\nnumber of expected segments:",seg
	    framenumber = fn
	    if timecodes == '' and i == INFLIST[0]:
		TIME = [t0,ts]
		print "Timestamps read from remarks:\ntime0 = ",t0,"timestep = ",ts
	except Exception:
	    print "\nDISICL Remark could not be read in:",remark
    else:
	print "\nWarning, No DISICL remarks were found in data series!"
    info.close()
    info = open(infofile,'rb')
    for line in info:
	    if line.startswith('#'):
		    if line.find('time') == -1:
			parts = line.strip('#').split()
			header = line   
	    elif line.startswith('&'):
		    tone = "gray90"
		    name = 'UC'
		    dna = 0
		    prot = 0
		    for code in COLOR_DNA:
			abbr, color = code.split(",")
			if header.find("("+abbr+")") != -1:
				name = abbr
				tone = color
				dna = 1
		    if dna == 0:
			    for code in COLOR_PROT:
				abbr, color = code.split(",")
				if header.find("("+abbr+")") != -1:
					name = abbr
					tone = color
					prot = 1
		    if prot == 0 and dna == 0:
			for code in COLOR_CUST:
			    abbr, color = code.split(",")
			    if header.find("("+abbr+")") != -1:
				    name = abbr
				    tone = color

		    listdata = (name,tone,TMP)
		    if (dna == 1) or (prot == 0 and dna == 0 and mode == 2):
			DATA.append(listdata)
		    elif (prot == 1) or (prot == 0 and dna == 0 and (mode == 1 or mode == 3)):
			DATA2.append(listdata)
		    elif prot == 0 and dna == 0 and mode == 0:
			DATA3.append(listdata)
		    TMP = []
	    else:
		    try:	
			    parts = line.split()
			    time = float(parts[0])
			    res = int(parts[1])
			    frame = int((time-TIME[0])/TIME[1])+1
			    resdata = (frame,res)
			    TMP.append(resdata)
		    except Exception:
			    print line,"could not be read"
    info.close()


#adding unmarked residues that belong to a class:
print "\nAdding class infromation:"
enum = len(DATA)
TMP = []
cnt = 0
while cnt < enum:
	cont = []
	TMP.append(cont)
	cnt = cnt + 1

ecnt = 0
noclass = "UC"
print "Adding nucleic acid codes:" 
for e in DATA:
	if e[0].find(noclass) == -1 and e[2] != []:
		rcnt = 0
		for res in e[2]:
			rcnt = rcnt + 1
			timep = res[0]
			resp = res[1]-1
			res1 = (timep,resp)
			if res1 not in e[2]:
				cnt = 0
				found = 0
				while cnt < enum:
					if (res1 in DATA[cnt][2]) and (DATA[cnt][0].find(noclass) == -1):
						found = 1
						if debug == 1:
						    print("residue "+str(res)+" found in "+e[0])
					cnt = cnt + 1
				if found == 0:
					TMP[ecnt].append(res1)
					rcnt = rcnt + 1
					if debug == 1:
					    print("residue "+str(res)+" added to "+e[0])
		print e[0]+" done with "+str(rcnt)+" residues"
	ecnt = ecnt + 1
ecnt = 0
while ecnt < enum:
	for add in TMP[ecnt]:
		DATA[ecnt][2].append(add)	
	ecnt = ecnt + 1

 
enum2 = len(DATA2)
TMP = []
cnt = 0
while cnt < enum2:
	cont = []
	TMP.append(cont)
	cnt = cnt + 1
ecnt = 0
print "\nAdding protein codes:"
if enum2 != 0:
	for e in DATA2:
		if e[0].find(noclass) == -1 and e[2] != []:
			rcnt = 0
			for res in e[2]:
				if debug == 1:
				    print("residue "+str(res)+" found in "+e[0])
				resn = res[1]+1
				timen = res[0]
				res1 = (timen,resn)
				rcnt += 1
				if res1 not in e[2]:
					if res1 in DATA2[enum2-1][2] and DATA2[enum2-1][0].find(noclass) != -1:
						TMP[ecnt].append(res1)
						rcnt += 1
						if debug == 1:
						    print("residue "+str(res1)+" added to "+e[0])
			print e[0]+" done with "+str(rcnt)+" residues"
		ecnt = ecnt + 1					
ecnt = 0
while ecnt < enum2:
	for add in TMP[ecnt]:
		DATA2[ecnt][2].append(add)	
	ecnt = ecnt + 1

enum3 = len(DATA3)
ecnt = 0
print "\nAdding custom codes:" 
for e in DATA3:
	if e[0].find(noclass) == -1 and e[2] != []:
		rcnt = 0
		for res in e[2]:
			rcnt = rcnt + 1
			if debug == 1:
			    print("residue "+str(res)+" found in "+e[0])
		print e[0]+" done with "+str(rcnt)+" residues"


print "\nReading in class information finished.\n"

#start pymol in a new thread
exitFlag = 0
import time
class Mythread (threading.Thread):
	def __init__(self,name,ID,command):
		threading.Thread.__init__(self)
		self.ID = int(ID)
		self.name = str(name)
		self.command = command
	def run(self):
		print "Starting  "+self.name
		pymol.pymolargv = ["pymol","-pc",infile]
#		pymol.start_pymol()
		pymol.finish_launching()		
		print "Exiting "+self.name
	
if Dest != '':
    fulldest = os.path.join(workdir,Dest)
    if os.path.isdir(fulldest) == True:
	print "Destination directory found!, changing path"
    else:
	print "Creating destination directory:",fulldest
	os.mkdir(fulldest)
    os.chdir(fulldest)
    
pytime = os.times()	
pymolthread = Mythread("pymol-1",1,"pymol.start_pymol()")
pymolthread.start()
print "\nStarinting PyMol thread, formatting results!\n"

if namebase != '':
	name_format= namebase+'_DISICL'
else:
	parts = infile.split(".")
	ext_length = len(parts[-1])+1
	name_format =infile[:-ext_length]
print "loading session: "+name_format
from pymol import stored
time.sleep(10)
if Dest == '':
    cmd.load(infile)
else:
    infile2 = os.path.join(workdir,infile)
    cmd.load(infile2)
cmd.color("gray90","all")
cmd.hide('lines','all')
#cmd.set('depth_cue',0)
#cmd.set('fog',0)
#myspace = {"cnt" : 0}

obtained_framenum = cmd.count_states(selection='(all)')
if framenumber != 1:
	if (1+(framenumber-1)*stride) <= obtained_framenum and (1+framenumber*stride) > obtained_framenum:
		print "\nExpected number of frames : "+str(framenumber)+" Obtained: "+str(obtained_framenum)
		if stride != 1:
		    print "Striding was set to: "+str(stride)+", Culling snapshots"
	else:
	    print "\nError: Obtained ("+str(obtained_framenum)+") and expected ("+str(framenumber*stride)+") frame numbers do not match!"
	    suggest = obtained_framenum/framenumber
	    if suggest > 1:
		print "Use the @stride "+str(suggest)+" command if you want cull pdb coordinates!"
	    else:
		print "Check @stride parameters (current:",stride,")"
	    cmd.quit()
	    sys.exit()
else:
	framenumber = obtained_framenum
	
fcnt = 1
while fcnt <= framenumber:
	cmd.create("frame_"+str(fcnt), name_format, (1+(fcnt-1)*stride), 1)
	if fcnt != 1:
	    cmd.align("frame_"+str(fcnt), "frame_1")
	elif fcnt  == 1:
	    cmd.orient("frame_"+str(fcnt))
	fcnt = fcnt + 1

cmd.delete(name_format)

cmd.disable('all')

acc = 30

DATAS = [DATA,DATA2,DATA3]
for Data in DATAS:
	print len(Data)
	if Data == DATA:
		modnum = 1
	else:
		modnum = 0
	for e in Data:
		if  e[2] != [] and e[0].find(noclass) == -1:
			name = e[0]
			tone = e[1]
			mark = 0
			current_frame = 0
			added_res = 0
#			print e[2][-1], e[2]
			print "Loading "+name+" class elments"
			for res in e[2]:
				rescode = res[1]
#				framenum = int((res[0]-TIME[0])/TIME[1])+1
				framenum = res[0]
				added = 0
				if mark == 0:
					cmdstring = ' /frame_'+str(framenum)+'///'+str(rescode+modnum)
					current_frame = framenum
					cmd.select(name,cmdstring)
					added_res = added_res + 1
					added = 1
					mark = 1
				elif mark == 1:
					cmdstring = name+' + '+'/frame_'+str(framenum)+'///'+str(rescode+modnum)
					current_frame = framenum
					added_res = added_res+1
					added = 1
					mark = mark + 1
				elif (1 < mark < acc) and (framenum == current_frame):
						cmdstring = cmdstring+"+"+str(rescode+modnum)
						mark = mark + 1
						added_res = added_res + 1
						added = 1
				if (mark == acc) or (framenum != current_frame) or (res == e[2][-1]):
					cmd.select(name,cmdstring)
					mark = 1
#					if debug == 1:
#					    print cmdstring
#					    print str(added_res)+" residues added"
					if added == 0:
					    cmdstring = name+' + '+'/frame_'+str(framenum)+'///'+str(rescode+modnum)
					    added_res = added_res+1
					    current_frame = framenum
					    mark = mark + 1
			print "   "+str(added_res)+" residues added"
			cmd.color(tone,name)
	
cmd.disable('all')
cmd.enable('frame_1')

if write == 3:
    fcnt = 1
    print "\nPreparing scenes....."
    while fcnt <= framenumber: 
	    if fcnt != 1:
		    cmd.disable('frame_'+str(fcnt-1))
	    cmd.enable('frame_'+str(fcnt))
	    cmd.show('cartoon','frame_'+str(fcnt))
	    cmd.scene('F'+str(fcnt), action='store', active=1, color=1, rep=1,view=0)
	    fcnt=fcnt+1
    print "done"
    
    print "Setting up movie parameters for session..."
    cmd.mset('1x1')
    cmd.set('scene_loop')
    cmd.set('movie_fps', 1.0)
    cmd.mdo(1, 'scene auto, next')
    cmd.scene('F'+str(1), action='recall')
    cmd.mplay()
    
elif write == 2:
    print "\nRendering images....."
    cmd.show('cartoon','all')
    fcnt = 1
    while fcnt <= framenumber:
	print "makig image",fcnt
	if fcnt != 1:
		cmd.disable('frame_'+str(fcnt-1))
	cmd.enable('frame_'+str(fcnt))
	cmd.center('frame_'+str(fcnt),'0','1')
#	cmd.ray("1600", "1200")
#	time.sleep(2)
	cmd.save('frame_'+str(fcnt)+'.png')
	print "image ",fcnt,"saved"
	time.sleep(2)
	fcnt=fcnt+1
    print "done"
else:
    cmd.enable('frame_1')
    fcnt = 2
    while fcnt <= framenumber: 
	cmd.disable('frame_'+str(fcnt))
	fcnt = fcnt + 1
    cmd.show('cartoon','all')
    print "done"
    	
print "\nSaving Session"	
cmd.save(name_format+".pse")
cmd.quit()
if debug == 1:
    print("Pymol session written! Script finished successfully!\n")
os.chdir(workdir)

rtime = os.times()
runtime = rtime[4]-stime[4]
if runtime == 0.0:
    runtime = rtime[0]+rtime[1]-(stime[0]+stime[1])
print "Script fininshed in %3.2f seconds (%1.2f h)" % (runtime,runtime/3600)
print "Visualization script ("+name_format+") finished successfully!" 
