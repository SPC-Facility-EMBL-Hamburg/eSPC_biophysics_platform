#!/usr/bin/env python
#

#License information:
#    DISICL main module to scynchronize DISICL program workflow.
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
print(workdir)

#usage
usage0 = '************************************\n'
usage1 = 'DISICL: Dihedral based Segment Identification and CLassification\nUsage for DISICL_main.py:\n\n'
usage2 = 'DISICL_main.py [input file] [library specification] [options]....\n\nInput file options:\n'
usage3 = '   <file-name>\n\t   Single trajectory file in pdb 1.0 compatible format (default)\n'
usage4 = '   @file <file-name>\n\t   Give argument file (no other arguments read from command line)\n'
usage5 = '   @traj <file1> <file2> ...\n\t   Give multiple files as a single trajectory\n'
usage6 = '   @single <file1> <file2> ...\n\t   Give multiple files as a database of independent structures\n'
usage7 = '   <file-name>\n\t   Single dihedral time series file (.dih)\n\t   (if DISICL_prep and DISICL_dihed are turned off)\n\n'
usage8 = 'Library specification options:\n'
usage9 = '   @protlib [library options] (default)\n\t   Define protein library for classification\n'
usage10 = '   @dnalib [library options]\n\t   Define nucleotide library for classifcaion\n'
usage11 = '   @custlib [library options]\n\t   Define custom biopolymer library for classifcation\n\n'
usage12 = 'Library options:\n   1 - Use DISICL detailed library (default)\n   2 - Use DISICL simple library\n'
usage13 = '   3 - Use both DISICL detailed and simple libraries\n   <file-name> - Use custom library file\n'
usage14 = '   \tThere are no DISICL libraries for the @custlib type\n\n'
usage15 = 'Additional options:\n   @prep  <0/1>\n   @dihed <0/1>\n   @dbssp <0/1>\n'
usage16 = '   \tTurn on(1) or off(0) DISICL preparation, dihedral or dbssp modules\n   \tTurning off modules is advised only for advanced users!\n\n'
usage17 = '   @vis <v>\n\t   Visualize results with DISICL_vis module, v can take values:\n\t   0 - no visualization (default)\n\t'
usage18 = '   1 - PyMol session (.pse)\n\t   2 - PyMol session (.pse) + images (1/frame, .png) \n\t   3 - Animated PyMol session (.pse) \n\n\t'
usage19 = '   Option 3 requires high amounts of memory and time,\n\t   recommended for small molecules and trajectories only!\n\n'
usage20 = '   @res <res,res0-resmax>\n\t   Define residue range for analyisis (comma separated)\n\t   arguments are integers, or ranges (- separated)\n\t'
usage21 = '   if no range is given all valid residues are analysed\n\n   @time <time,timestep>\n\t   Define timestamps for trajectory (comma separated)\n\t'
usage22 = '   use floats for starting time and timestep (default 0.0,1.0)\n\n   @stride <n>\n\t   Calculate only every nth frame (default 1)\n\n'
usage23 = '   @path <path>\n\t   Set different source directory for structure files (default ./)\n\n'
usage24 = '   @dest <path>\n\t   Set different goal directory for results (relative, default ./)\n\n'
usage25 = '   @protres <str1,str2,...>\n   @nuclres <str1,str2,...>\n   @custres <str1,str2,...>\n'
usage26 = '\t   Add standard protein, nucleotide, or custom residue codes\n\t   (comma separated)\n\n'
usage27 = '   @rename <0/1>\n\t   Control DNA/RNA renaming to standard format\n\t   1 - on (default)\n\t   0 - off\n\n'
usage28 = '   @exclude <str1,str2,...>\n\t   Exclude residue codes from analysis (comma separated)\n\n'
usage29 = '   @ext <.str>\n\t   Give input file extension (. followed by a string)\n\t   default extensions:\n'
usage30 = '\t   .pdb (protein databank format)\n\t   .dih (DISICL dihedral time series)\n\n'
usage31 = '    @filter <f>\n\t   Set filter level to restrict residues, f can take values:\n\t   0 - no filter\n\t   1 - only excluded\n\t'
usage32 = '   2 - all nonstandard residues\n\t   3 - all nonstandard residues and alternative conformations (default)\n\n'
usage33 = '   @libhelp\n\t   Library format guide\n\n   @debug\n\t   Print additonal messages in the logfile\n\n'
usage34 = '   @command\n\t   Flag marks the argument file identical with the input file\n\t   arguments are read until the first line starting with an "&"\n\t'
usage35 = '   For argument files ony, do not use @command in the command line!\n\n   @help\n\t   Print this message\n\n'
usage36 = 'NOTE: Visualization requires PyMol, see DISICL_vis.py usage for more details\n\n'

#usage14 = '@out to define database output filename (def: analysed.list)\n'
#usage15 = '@grom to use GROMOS residue ID-s for dssp_dna:\n 1-yes\n 2-no (def)\n'
#usage27 = '   @refine [int]\n\t   Write refined pdb files:\n\t   1-yes (def)\n\t   0-no\n'

usage = usage0+usage1+usage2+usage3+usage4+usage5+usage6+usage7+usage8+usage9+usage10+usage11+usage12
usage = usage+usage13+usage14+usage15+usage16+usage17+usage18+usage19+usage20+usage21+usage22+usage23+usage24
usage = usage+usage25+usage26+usage27+usage28+usage29+usage30+usage31+usage32+usage33+usage34+usage35+usage36+usage0

lib0 = "**********************************************************************\n"
lib1 = "Format specifications for the DISICL library files:\n\n"
lib2 = "Give custom defined dihedral angles with:\n   @ang <name,atom1,atom2,atom3,atom4,mod1,mod2,mod3,mod4'>\n"
lib3 = "\t   name is a string used for the dihedral angle\n\n"
lib4 = "\t   atom1-atom4 are atom name strings\n\n\t   mod1-mod4 are residue number modifiers (integers)\n\n"
lib5 = "Regions are defined as a rectangular/cubic parts of the dihedral angle space\n"
lib6 = "Define Ramachandran regions using:\n   @reg <name,xmin,xmax,ymin,ymax,zmin,zmax,...>\n"
lib7 = "\t   name is a string used for strucutre definitions\n\n\t   xmin,xmax,... are pairs of floats, one for each dihedral angle\n\n"
lib8 = "Define segment for a secondary structure element using:\n   @sec <name,reg1,reg2,....,bin>\n"
lib9 = "\t   name is the full name of secondary structure element\n\n\t   reg1, reg2, ... are the names of sebsequent regions\n"
lib10 = "\t   the residues of the segment should fall in\n\n\t   bin is the collector code for the secondary structure element\n\t   (if you have multiple areas).\n\n"
lib11 = "\t   It is recommended that name and bin is always changed together!\n"
lib12 = "Define a color definition using:\n   @col <name,color,type>\n\t"
lib13 = "   name is the collector code used for strucutre element\n\t   (put beween brackets in the class header)\n\n\t"
lib14 = "   color is the PyMol reference name for the color associated\n\t   with the structure (see PyMol for details)\n\n\t"
lib15 = "   type is the molecule type code\n\t   1 - protein\n\t   2 - nucleotide\n\t   0 - custom\n"

libhelp = lib0+lib1+lib2+lib3+lib4+lib5+lib6+lib7+lib8+lib9+lib10+lib11+lib2+lib3+lib4+lib5+lib0



######System variables, should be edited upon installtion#######
#local definitions for Windows / Linux commands:
win = 0
# 0 = use linux format, 1 = assume windows format
#Set main DISICL path and library path for DISICL
DISICL_dir =  "/home/os/spc_shiny_servers/circularDichroismApp/appFiles/ChiraKit/SESCA_v097/DISICL"
libs = os.path.join(DISICL_dir,"libs")

######Nothing should be modified from this point, apply changes at your own risk####### 

#program definitions and paramteres:
if win == 1:
    rm = "del"
    cp = "copy"
    mv = "move"
else:
    rm = "rm"
    cp = "cp"
    mv = "mv"

protein = 1
dna = 0
traj = 1
single = 0
debug = ''
prep = 1
dihed = 1
dbssp = 1
vis = 0
prot_det = 1
prot_sim = 0
prot_cust = 0
dna_det = 0
dna_sim = 0
dna_cust = 0
other = 0

prep_prog = os.path.join(DISICL_dir,"DISICL_prep.py")
dihed_prog = os.path.join(DISICL_dir,"DISICL_dihed.py")
dbssp_prog = os.path.join(DISICL_dir,"DISICL_dbssp.py")
vis_prog = os.path.join(DISICL_dir,"DISICL_vis.py")
dna_det_lib = os.path.join(libs,"DISICL_nucl_det.lib")
dna_sim_lib = os.path.join(libs,"DISICL_nucl_sim.lib")
dna_cust_lib = ''
prot_det_lib = os.path.join(libs,"DISICL_prot_det.lib")
prot_sim_lib = os.path.join(libs,"DISICL_prot_sim.lib")
prot_cust_lib = ''
custom_lib = ''
listfile = 'analysed.list'
LOGFORMAT = ['DISICL_','.log']
ext = '.pdb'

ARGUMENTS = []
TRAJFILES = []
JOBLIST = []
MODLIST = []
TIME = [0,1]
RES = []
resdata = ''
timecodes = ''
stride = ''
source = ''
dest = ''
rename = ''
ext = '.pdb'
refine = ''
resfilt = ''
customprot = ''
customnucl = ''
customres = ''
customex = ''

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
            print("\nERROR: Argument file not found!")
            failmark = 1

#command line interface:
flaglist = []
flaglist1 = ['traj','single','filter','ext','path','dest','file','help','command','debug','prep','dihed','dbssp','vis']
flaglist2 = ['protres','nuclres','custres','exclude','rename','stride','time','res','protlib','dnalib','custlib','libhelp']
for flag in flaglist1:
    flaglist.append(flag)
for flag in flaglist2:
    flaglist.append(flag)

failmark = 0
flag =''    
argnum = len(ARGUMENTS)-1
#print ARGUMENTS
print('number of arguments: '+str(argnum))
if argnum >= 1:
    print('flags found:')
    cnt = 1
    while cnt <= argnum:
        arg = ARGUMENTS[cnt]
#        print(arg,flag)
        if arg == '@help':
            print(usage)
            sys.exit()
        elif arg == '@libhelp':
            print(libhelp)
            sys.exit()
        elif arg == '@debug':
            debug = 1
            print("debug")
            flag = ''
        elif arg == '@command':
            TRAJFILES.append(ARGUMENTS[0])
            flag = ''
        elif arg == '@traj':
            traj = 1
            flag = "traj"
            single = 0
            print("traj")
        elif arg == '@single':
            flag = 'single'
            traj = 0
            single = 1
            print("single")    
        elif arg.startswith('@'):
            flag = arg.strip('@')
            print(flag)
            if not flag in flaglist:
                  print('\nERROR: Unkown flag, program stops!')
                  failmark = 1
                  flag = ''
        elif (flag == '' and cnt == 1) or flag == 'traj' or flag == 'single':
            TRAJFILES.append(arg)
        elif flag == 'path':
            source = arg
            flag = ''
        elif flag == 'dest':
            dest = arg
            flag = ''
        elif flag == 'refine':
            try:
                refine = int(arg)    
                if refflag < 0 or refflag > 1:
                    print('\nERROR: refine code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: refine code not recognized!')
            flag = ''
        elif flag == 'prep':
            try:
                prep = int(arg)
                if prep < 0 or prep > 1:
                    print('\nERROR: prep code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: prep code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'prep':
            try:
                prep = int(arg)
                if prep < 0 or prep > 1:
                    print('\nERROR: prep code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: prep code not recognized!')
                falimark = 1
            flag == ''
        elif flag == 'dihed':
            try:
                dihed = int(arg)
                if dihed < 0 or dihed > 1:
                    print('\nERROR: dihed code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: dihed code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'dbssp':
            try:
                dbssp = int(arg)
                if dbssp < 0 or dbssp > 1:
                    print('\nERROR: dbssp code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: dbssp code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'vis':
            try:
                vis = int(arg)
                if vis < 0 or vis > 3:
                    print('\nERROR: vis code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: vis code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'filter':
            try:
                resfilt = int(arg)
                if resfilt < 0 or resfilt > 3:
                    print('\nERROR: filter code not recognized!')
                    falimark = 1
            except Exception:
                print('\nERROR: filter code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'ext':
                ext = arg
                flag = ''
        elif flag == "res":
            resdata = arg
            flag = ''
        elif flag == "time":
            timecodes = arg
            try:
                    parts = timecodes.split(',')
                    timezero = float(parts[0])
                    timestep = float(parts[1])
                    TIME = [timezero,timestep]
            except Exception:
                    print("\nError: could not read in time information!")
                    failmark = 1
            flag = ''
        elif flag == 'protres':
            customprot = arg
            flag = ''
        elif flag == 'nuclres':
            customnucl = arg
            flag = ''
        elif flag == 'custres':
            customres = arg
            flag = ''
        elif flag == 'exclude':
            customex = arg
            flag = ''
        elif flag == 'rename':
            try:
                rename = int(arg)
                if not(rename == 0 or rename == 1):
                    print("\nError: @rename code not recognized")
                    failmark = 1
            except ValueError:
                print("\nError: @rename code not recognized")
                failmark = 1
            flag = ''
        elif flag == 'stride':
            try:
                    stride = int(arg)
            except ValueError:
                    print("\nError: @stide flag takes only integers")
                    failmark = 1
            flag = ''
        elif flag == 'protlib':
            protein = 1
            if arg == "0":
                protein = 0
            elif arg == "1":
                prot_det = 1
                prot_sim = 0
            elif arg == "2":
                prot_det = 0
                prot_sim = 1
            elif arg == "3":
                prot_det = 1
                prot_sim = 1
            else:
                prot_det = 0
                prot_sim = 0
                prot_cust = 1
                prot_cust_lib = arg
        elif flag == 'dnalib':
            dna = 1
            if arg == 0:
                dna = 0
            elif arg == "1":
                dna_det = 1
                dna_sim = 0
            elif arg == "2":
                dna_det = 0
                dna_sim = 1
            elif arg == "3":
                dna_det = 1
                dna_sim = 1
            else:
                dna_det = 0
                dna_sim = 0
                dna_cust = 1
                dna_cust_lib = arg
        elif flag == 'custlib':
            other = 1
            if arg == "0":
                other = 0
            else:
                protein = 0
                dna = 0
                custom_lib = arg
        cnt =cnt + 1

else:
    print('no arguments provided, program stops!')
    failmark = 1

#checking module information:
if dihed == 0 and prep == 1 and dbssp == 1:
    print("\nWarning: dihed module is turnod off, while prep and dbssp module is on!\n")
    print("DISICL result will not be consistent with input files!!!!")
elif prep == 0 and dihed == 0 and dbssp == 1:
    print('\nWarning: Modules prep and dihed are turned off, setting file extension to ".dih"')
    ext = '.dih'
elif prep == 0 and dihed == 0 and dbssp == 0:
    print("\nERROR: All modules are turned off, program stops!")
    failmark = 1
if dihed == 0 and protein+dna+other > 1:
    print("Warning: ")
#checking path and file validity
if single == 1:
    print("@single flag detected, Database mode is used\n")
    for arg2 in ARGUMENTS:
                if arg2 == "@traj":    
                    print("\nERROR: @single and @traj flags cannot be provided at the same time!")
                    failmark = 1
if other == 1:
    print("\nCustom library definition detected,\nall hardcoded libraries and values are disabled!")
    if customres == '':
        print("Warning: custom library provided, but no custom definiton for residues")
                    
if TRAJFILES != []:    
    print("\nFiles to process:")
    for filename in TRAJFILES:
        print(filename)
else:
    print("\nERROR: no input files detected!")
    print(usage)
    sys.exit()

Datadir1 = os.path.join(workdir,source)
print("Source is: %1s"%source)
if source == '':
    DataDir = workdir
elif os.path.isdir(Datadir1) == True: 
    DataDir = Datadir1
elif os.path.isdir(source) == True:
    Datadir = source
else:
    print("Error: Could not find source directory!")
    print(usage)
    sys.exit()

print("Source Directory found:")
print(DataDir)
os.chdir(DataDir)
Filelist = os.listdir(DataDir)

print('\nValid Files found:')
cnt = 0
FirstFile = ''
for File in Filelist:
    if (File in TRAJFILES) and (os.path.isfile(File) == True):
        ext_length = len(ext)
        namebase  = File[0:(-1*ext_length)]
        file_extension = File[-1*ext_length:]
        if cnt == 0:
            FirstFile = namebase 
        if file_extension != ext:
            print('\nERROR:',File,"\nexpected extension: "+ext)
            print("found: "+file_extension)
            print("use @ext if extension is correct")
            print(usage)
            sys.exit()
        print(namebase)
        line = [namebase]
        JOBLIST.append(line)
        cnt = cnt + 1
if cnt == 0:
    print('\nERROR: no valid files are detcted in source directory !!!')
    failmark = 1

Libs = []
cnt = 0
if prot_cust_lib != '':
    Libs.append(prot_cust_lib)
if dna_cust_lib != '':
    Libs.append(dna_cust_lib)
if custom_lib != '':
    Libs.append(custom_lib)
for lib in Libs:
    newlib = ''
    libpath1 = os.path.join(workdir,lib)
    libpath2 = os.path.join(DataDir,lib)
#    libpath3 = os.path.join(GoalDir,lib)
    if os.path.isfile(libpath1) == True:
        newlib = libpath1
    elif os.path.isfile(libpath2) == True:
        newlib = libpath2
    elif os.path.isfile(lib) == True:
        newlib = lib
    else:
        print("defined library not found!")
        print(newlib)
        failmark = 1
    if failmark == 0:
        print("custom library found:\n"+lib)
        if lib == prot_cust_lib:
            prot_cust_lib = newlib
        if lib == dna_cust_lib:
            dna_cust_lib = newlib
        if lib == custom_lib:
            custom_lib = newlib
            

os.chdir(workdir)
    
if dest == '':
    GoalDir = workdir
else:
    GoalDir = os.path.join(workdir,dest)

if refine == 1 or refine == '':
    if os.path.isdir(GoalDir) == True:
        print('\ndestination directory found:')
        print(GoalDir)
    else:
        print('\ndestination directroy not found, attemp to create it')
        print(GoalDir)
        os.system('mkdir '+dest)
        if os.path.isdir(GoalDir) == True:
            print('\ndestination directory created:')
            print(GoalDir)
        else:
            print('Error, could not create destination directory, program stops!')
            print(usage)
            sys.exit()

if resdata != '':
    print("Custom residue range defined:",resdata,"\nDISICL analysis is carried out on these residues\n")

if timecodes != '':
    print("Custom timestamp labels defined:","\ntime0= ",str(TIME[0]),"timestep=",str(TIME[1]),"\n")

    
if failmark == 1:
    print(usage)
    sys.exit()
    
#preparing information for classification protocols:

separator = '\n************************************************\n'
if single == 1:
    logfile = listfile
if traj == 1:
    logfile = LOGFORMAT[0]+FirstFile+LOGFORMAT[1]
    info = ''
module_times = 0.0
log = os.path.join(GoalDir,logfile)
logheader = "# DISICL_main.py logfile"
if win == 0: 
    os.system('echo "'+logheader+'\n" > '+log)
else:
    os.system('echo '+logheader+' > '+log)
    os.system("echo ''  >> "+log)
#gathering custom command flags for DISICL_prep

if prep == 1:
    preptype = 0
    if protein == 1 and dna == 1:
        preptype = 3
    elif dna == 1:
        preptype = 2
    elif protein == 1:
        preptype = 1
    elif other == 1:
        preptype = 0
    print(separator+"\nDISICL_prep module started, using command flags:\n"+prep_prog)
    prepheader = 10*'*'+' DISICL_prep  '+10*'*'
    if win == 0:
        os.system("echo '"+prepheader+"\n\n' >> "+log)
    else:
        os.system("echo "+prepheader+" >> "+log)
        os.system("echo ''  >> "+ log)
    Prep_flags = [refine,resfilt,ext,source,dest,customprot,customnucl,customres,customex,rename,traj,single,]
    prepadd = ' '
    add = "@type "+str(preptype)
    prepadd = prepadd+add+' '
    print(add)
    if refine != '':
        add = "@refine "+str(refine)
        prepadd = prepadd+add+' '
        print(add)
    if resfilt != '':
        add = "@filter "+str(resfilt)
        prepadd = prepadd+add+' '
        print(add)
    if ext != '.pdb':
        add = "@ext "+str(ext)
        prepadd = prepadd+add+' '
        print(add)
    if source != '':
        add = "@path "+str(source)
        prepadd = prepadd+add+' '
        print(add)
    if dest != '':
        add = "@dest "+str(dest)
        prepadd = prepadd+add+' '
        print(add)
    if customprot != '':
        add = "@protres "+str(customprot)
        prepadd = prepadd+add+' '
        print(add)
    if customnucl != '':
        add = "@nuclres "+str(customnucl)
        prepadd = prepadd+add+' '
        print(add)
    if customres != '':
        add = "@custres "+str(customres)
        prepadd = prepadd+add+' '
        print(add)
    if customex != '':
        add = "@exclude "+str(customex)
        prepadd = prepadd+add+' '
        print(add)
    if rename != '':
        add = "@rename "+str(rename)
        prepadd = prepadd+add+' '
        print(add)
    if traj == 1:
        add = "@traj "
        prepadd = prepadd+add
        print(add)
    elif single == 1:
        add = "@single "
        prepadd = prepadd+add
        print(add)    
    prep_command =  prep_prog+prepadd   
    for File in JOBLIST:
        print(File[0]+ext)
        prep_command = prep_command+File[0]+ext+' '
#    print "\nDISICL_prep command:\n",prep_command
    print(" >> "+log+"\n")
    if win == 0:
        os.system("echo '"+prep_command+"\n' >> "+log)
    else:
        os.system("echo "+prep_command+" >> "+log)
        os.system("echo '' >> "+log)
    os.system(prep_command+' >> '+log)
    l = open(log,'rb')
#checking module success and results
    timeline = ''
    modtime = 0
    for line0 in l:
        line = str(line0.decode("ascii"))
        message = line
        if traj == 1:
          if line.startswith("DISICL"):
            info = line
        if line.startswith("program runtime"):
            timeline = line
    if timeline != '':
        try:
            parts = timeline.split()
            modtime = float(parts[2]) 
            module_times = module_times+modtime
        except Exception:
            print("\nWarning: runtime on module was not found")
    if message.find("Preparation script finished successfully!") != -1:
        print("\nDISICL_prep finished!")
    else:
        print("\nError during data preparation, see logfile for details!")
        failmark = 3    

os.chdir(GoalDir)
    
#getting proper joblist
JOBLIST = []
if traj == 1:
    if prep == 0:
        f = open(FirstFile+ext,"rb")
        info = ''
        for line0 in f:
            line = str(line0.decode("ascii"))
            if line.startswith("DISICL") and dihed == 1:
                info = line
            if line.startswith("#DISICL") and dihed == 0:
                info = line
        f.close()
        if info == '':
            print(FirstFile)
            print("\nERROR: DISICL_prep turned off and no remark in input file! Program stops!")
            print(usage)
            sys.exit()
    if prep == 1 or (prep == 0 and dihed == 1):
        try:
            parts = info.strip('\n').split()
            filename = FirstFile
            modnum = int(parts[2])
            resnum = int(parts[3])
            resstring = parts[4]+parts[5]+parts[6]
        except Exception:
            print(FirstFile)
            print("ERROR: Unable to read in remark information in file, program stops!")
            print(usage)
            sys.exit()
    if prep == 0 and dihed == 0:
        try:
            parts = info.strip('\n').split()
            filename = FirstFile
            modnum = int(parts[3])
            resnum = int(parts[4])
            resstring = parts[5]
        except Exception:
            print(FirstFile)
            print("ERROR: Unable to read in remark information in file, program stops!")
            print(usage)
            sys.exit()
    data = (filename,modnum,resnum,resstring)
    JOBLIST.append(data)
elif single == 1:
    infile = 'database.list'
    i = open(infile,'rb')
    for line0 in i:
        line = str(line0.decode("ascii"))
        if not line.startswith('#'):
            mark = 0
            try:
                print('\n'+line.strip('\n'))
#                part = line.split('[')
                parts = line.split()
                filename = parts[0]
                modnum = int(parts[1])
                resnum = int(parts[2])
                exclude = int(parts[3])
                resstring = parts[4]+parts[5]+parts[6]
                if os.path.isfile(filename+'_1'+ext) == False:
                    mark = 1
                    print('file not found for entry: '+filename+'\n entry excluded!')
                else:
                    print('name: '+filename+'\nmodels: '+str(modnum)+'\nresnum:'+str(resnum)+'\nreslist:\n'+resstring)
                data = [filename,modnum,resnum,resstring]
                if resnum != 0 and resstring != '' and mark == 0:
                    JOBLIST.append(data)
            except IndexError:
                print('line ignored: '+line)
    i.close()

entrynum = len(JOBLIST)
#print entrynun
#for job in JOBLIST:
 #   print job
 
cnt = 0
fcnt = 0
rcnt = 0
scnt = 0
RMAX = ['',0]
RMIN = ['',0]
MMAX = ['',0]
if single == 1:
    o = open(listfile,'wb')  
    o.write('#list of database entries analysed\n#Name '+10*' '+'   resnum    success    reason to fail\n')
for File in JOBLIST:
#getting data statistics
    print(' ')
#    print File
    name = File[0]
    mnum = File[1]
    rnum = File[2]
    reslist = File[3]
    cnt = cnt + 1
    fcnt = fcnt + mnum
    rcnt = rcnt + mnum*rnum
    if single == 1:
        logfile = LOGFORMAT[0]+name+LOGFORMAT[1]
        print('*****************************')
        print('Database entry ('+str(cnt)+'/'+str(entrynum)+')')
        os.system('echo "Logfile for database entry: '+name+'.pdb" > '+logfile)
    ecnt = 1
    while ecnt <= mnum:
        failmark = 0
        bad = ''
        protok = 1
        dnaok = 1
        custok = 1
        if single == 1:
            name2 = name+'_'+str(ecnt)
            print('model '+name2+' ( from '+str(mnum)+')')
        else:
            name2 = name+'_DISICL'
#gathering information for DISICL_dihed
        if dihed == 1:
            if prep == 0:
                name2 = name
            dihedfile = name2+ext
            print(separator)
            print("DISICL_dihed module started, using command flags:")
            print(dihed_prog+"\n"+dihedfile)
            dihed_header = 10*'*'+' Dihed program '+10*'*'
            if win == 0:
                os.system("echo '"+dihed_header+"\n' >> "+logfile)
            else:
                os.system("echo "+dihed_header+" >> "+logfile)
                os.system("echo ''  >> "+logfile)
            DIHED_MODS= []
            if protein == 1:
                DIHED_MODS.append("prot")
            if dna == 1:
                DIHED_MODS.append("nucl")
            if other == 1:
                DIHED_MODS.append("cust")
                
#relevant command flags
            dihadd = ' '
            dihed_flags = [resdata,timecodes,stride,traj,debug,prot_cust_lib,dna_cust_lib,custom_lib]
            if resdata != '':
                add = '@res '+resdata+' '
                print(add)
                dihadd = dihadd + add
            if timecodes != '':
                add = '@time '+str(TIME[0])+','+str(TIME[1])+' '
                print(add)
                dihadd = dihadd + add
            if stride != '':
                add = '@stride '+str(stride)+' '
                print(add)
                dihadd = dihadd + add
            if traj == 1:
                add = '@traj '
                print(add)
                dihadd = dihadd + add
            if debug == 1:
                add = '@debug '
                print(add)
                dihadd = dihadd + add
            if protein == 1 and prot_cust_lib != '':
                add = '@lib '+prot_cust_lib+' '
                print(add+" (protein)")
            if dna == 1 and dna_cust_lib != '':
                add = '@lib '+dna_cust_lib+' '
                print(add+" (nucleotide)")
            if other == 1 and custom_lib != '':
                add = '@lib '+prot_cust_lib+' '
                print(add+" (custom)")
            print(" >> "+logfile+"\n")
            
#running dihedral calculations:
            for Type in DIHED_MODS:
                lib = ''
                typevalue = ''
                if Type == 'prot':
                    typename = "proteins"
                    lib = prot_cust_lib
                    typevalue = 1
                elif Type == 'nucl':
                    typename = "nucleotides"
                    lib = dna_cust_lib
                    typevalue = 2
                elif Type == 'cust':
                    typename = "custom residues"
                    lib = custom_lib
                    typevalue = 0
                dihedtype = ' @type '+str(typevalue)
                if lib != '':
                    dihedtype = dihedtype+' @lib '+lib
                    
                if traj == 1:
                    dihed_out = name+'_'+Type+'.dih'
                elif single == 1:
                    dihed_out = name2+'_'+Type+'.dih'
                    if resdata == '':
                        parts = reslist.split("[")
                        if Type == 'prot':
                            add = parts[1].strip("]")
                        if Type == 'nucl':
                            add = parts[2].strip("]")
                        if Type == 'cust':
                            add = parts[3].strip("]")
                        dihadd = dihadd+"@res "+add+" "
                print('Running DISICL_dihed program for '+typename)
                dihed_command = dihed_prog+' '+dihedfile+dihedtype+dihadd+' >> '+logfile
                if win == 0:
                    os.system('echo "DISICL_dihed for '+typename+'" >> '+logfile)
                    os.system('echo "'+dihed_command+'\n" >>' +logfile)
                else:
                    os.system("echo DISICL_dihed for "+typename+" >> "+logfile)
                    os.system('echo '+dihed_command+' >> '+logfile)
                os.system(dihed_command+' >> '+logfile)
#checking module success
                l = open(logfile,'rb')
                modtime = 0
                timeline = ''
                for line0 in l:
                    line = str(line0.decode("ascii"))
                    last = line
                    if line.startswith("program runtime"):
                        timeline = line 
                l.close()
                if not (last.find('Dihedral') != -1 and last.find(name2) != -1 and last.find(Type) != -1):
                    failmark = 4
                    if Type == 'prot':
                        protok = 0
                    elif Type == 'nucl':
                        dnaok = 0
                    elif Type == 'cust':
                        custok = 0
                    print(30*' '+'failed')   
                    bad = bad+' dihed_'+typename
                else:
                    os.system(mv+' '+name2+'.out '+dihed_out)
                    rename = 'renaming '+name2+'.out to '+dihed_out
                    if win == 0:
                        os.system('echo "'+rename+'" >> '+logfile)
                    else:    
                        os.system('echo '+rename+' >> '+logfile)
                    if timeline != '':    
                        try:
                                parts = timeline.split()
                                modtime = float(parts[2]) 
                                module_times = module_times+modtime
                        except Exception:
                            print("\nWarning: runtime on module was not found")
                    print(30*' '+'done')
                print('')
#gathering information for classification:
        if dbssp == 1:
            print(separator)
            print("DISICL_dbssp module started, using command flags:")
            print(dbssp_prog)
            dbssp_header = 7*'*'+' DISICL classification program '+7*'*'
            if win == 0:
                os.system("echo '"+dbssp_header+"\n' >> "+logfile)
            else:
                os.system("echo "+dbssp_header+" >> "+logfile)
                os.system("echo ''  >> "+logfile)
            if traj == 1:
                base = name
            else:
                base = name2
            LIBS= []
            SEGMENTS = [0,0,0]
            prot_dihed_file = ''
            dna_dihed_file = ''
            custom_dihed_file = ''
            if protein == 1 and protok == 1:
                if dihed == 1:
                    prot_dihed_file = base+'_prot.dih'
                else:
                    prot_dihed_file = base+ext
                print(prot_dihed_file+" (for proteins)")
                if prot_det == 1:
                    LIBS.append("pdet")
                if prot_sim == 1:
                    LIBS.append("psim")
                if prot_cust == 1:
                    LIBS.append("pcus")
            if dna == 1 and dnaok == 1:
                if dihed == 1:
                    dna_dihed_file = base+'_nucl.dih'
                else:
                    dna_dihed_file = base+ext
                print(dna_dihed_file+" (for nucleotides)")
                if dna_det == 1:
                    LIBS.append("ndet")
                if dna_sim == 1:
                    LIBS.append("nsim")
                if dna_cust == 1:
                    LIBS.append("ncus")
            if other == 1 and custok == 1:
                if dihed == 1:
                    custom_dihed_file = base+'_cust.dih'
                else:
                    custom_dihed_file = base+ext
                LIBS.append("cust")
#relevant command flags:
            dbsspadd = ' '
            dbssp_flags = [resdata,timecodes,debug,prot_cust_lib,dna_cust_lib,custom_lib]
            if resdata != '':
                add = '@res '+resdata+' '
                print(add)
                dbsspadd = dbsspadd + add
            if timecodes != '':
                add = '@time '+str(TIME[0])+','+str(TIME[1])+' '
                print(add)
                dbsspdd = dbsspadd + add
            if debug == 1:
                add = '@debug '
                print(add)
                dbsspadd = dbsspadd + add
            if protein == 1 and prot_cust_lib != '':
                add = '@lib '+prot_cust_lib+' '
                print(add+" (protein)")
            if dna == 1 and dna_cust_lib != '':
                add = '@lib '+dna_cust_lib+' '
                print(add+" (nucleotide)")
            if other == 1 and custom_lib != '':
                add = '@lib '+prot_cust_lib+' '
                print(add+" (custom)")
            print(" >> "+logfile+"\n")
#running classifications claculations:
#            print LIBS
            for job in LIBS:
                if job == 'pdet':
                    dfile = prot_dihed_file
                    lib = prot_det_lib
                    jobname = "protein - detailed"
                elif job == 'psim':
                    dfile = prot_dihed_file
                    lib = prot_sim_lib
                    jobname = "protein - simple"
                elif job == 'pcus':
                    dfile = prot_dihed_file
                    lib = prot_cust_lib
                    jobname = "protein - custom"
                elif job == 'ndet':
                    dfile = dna_dihed_file
                    lib = dna_det_lib
                    jobname = "nucleotide - detailed"
                elif job == 'nsim':
                    dfile = dna_dihed_file
                    lib = dna_sim_lib
                    jobname = "nucleotide - simple"
                elif job == 'ncus':
                    dfile = dna_dihed_file
                    lib = dna_cust_lib
                    jobname = "nucleotide - custom"    
                elif job == 'cust':
                    dfile = custom_dihed_file
                    lib = custom_lib
                    jobname = "custom library"
                if traj == 1:    
                    outformat = "DISICL_"+job+'_'+name
                elif single == 1:
                    outformat = "DISICL_"+job+'_'+name2
                    
                print('Running DISICL_dbssp program ('+jobname+')... ')
                job_header = 'DISICL_dbssp for '+jobname
                if win == 0:
                    os.system('echo "\n'+job_header+'\n" >> '+logfile)
                else:
                    os.system('echo '+job_header+'  >> '+logfile)
                    os.system("echo '' >> "+logfile)
                dbssp_command = dbssp_prog+' '+dfile+' @lib '+lib+dbsspadd+' >> '+logfile
                if win == 0:
                    os.system('echo "'+dbssp_command+'" >>' +logfile)
                else:
                    os.system('echo '+dbssp_command+' >>' +logfile)
                os.system(dbssp_command)
# checking module success and results
                l = open(logfile,'rb')
                last = ''
                remark = ''
                frames = 0
                segments = 0
                segnum = 0
                dbout = ''
                dbstat = ''
                timeline = ''
                modtime = 0
                for line0 in l:
                    line = str(line0.decode("ascii"))
                    last = line
                    if line.startswith("DISICL"):
                        remark = line
                    if line.startswith("writing classification time series"):
                        parts = line.strip('\n').split()
                        dbout = parts[-1]
                    if line.startswith("writing statistics to"):
                        parts = line.strip('\n').split()
                        dbstat = parts[-1]
                    if line.startswith("program runtime"):
                        timeline = line
                l.close()
                if not (last.find('Classification program') != -1 and last.find(dfile) != -1) and dbout != '' and dbstat != '':
                    failmark = 5
                if failmark != 5:                        
                    os.system(mv+' '+dbout+' '+outformat+'.out')
                    os.system(mv+' '+dbstat+' '+outformat+'.stat')
                    rename1 = 'renaming '+dbout+' to '+outformat+'.out'
                    rename2 = 'renaming '+dbstat+' to '+outformat+'.stat'
                    if win == 0:
                        os.system('echo "'+rename1+'" >> '+logfile)
                        os.system('echo "'+rename2+'\n" >> '+logfile)
                    else:
                        os.system('echo '+rename1+' >> '+logfile)
                        os.system('echo '+rename2+' >> '+logfile)
                        os.system("echo '' >> "+logfile)
                    if timeline != '':
                        try:
                                parts = timeline.split()
                                modtime = float(parts[2])
                                module_times = module_times+modtime    
                        except Exception:
                            print("\nWarning: runtime on module was not found")
                    print(30*' '+'done')
                else:
                    print(30*' '+'failed')
                    bad = bad+' dbssp'
                if remark != '':
                    try:
                        parts = remark.strip('\n').split()
                        frames = int(parts[3])
                        segments = int(parts[4])
                        segnum = segments*frames
                    except Exception:
                        print(remark)
                        print("\nWarning: no statistics found for "+outformat)
                    if job == "pdet" or job == "psim" or job == "pcus":
                        if SEGMENTS[0] == 0 and segnum != 0:
                            SEGMENTS[0] = segnum
                    if job == "ndet" or job == "nsim" or job == "ncus":
                        if SEGMENTS[1] == 0 and segnum != 0:
                            SEGMENTS[1] = segnum
                    if job == "cust":
                        if SEGMENTS[2] == 0 and segnum != 0:
                            SEGMENTS[2] = segnum
            segtotal = SEGMENTS[0]+SEGMENTS[1]+SEGMENTS[2]            
            scnt = scnt + segtotal
            print('')
 
 #gathering information for visualization:       
        if vis != 0:
            print(separator)
            print("DISICL_vis module started, using command flags:")
            print(vis_prog)
            vis_header = 6*'*'+' DISICL visualization script '+6*'*'
            if win == 0:
                os.system("echo '"+vis_header+"\n' >> "+logfile)
            else:
                os.system("echo "+vis_header+" >> "+logfile)
                os.system("echo ''  >> "+logfile)
            if traj == 1:
                base = name
            else:
                base = "@pdb "+name2+".pdb "
                
            if protein == 1 and dna == 1:
                vistype = 3
            elif dna == 1:
                vistype = 2
            elif protein == 1:
                vistype = 1
            else:
                vistype = 0
#relevant command flags:
            visadd = ' '
            vis_flags = [vis,timecodes,stride,prot_cust_lib,dna_cust_lib,custom_lib,debug]

            LIBS = []
            CLS = []
            if prot_cust_lib != '':
                LIBS.append(prot_cust_lib)
                CLS.append('pcus')
            elif prot_det == 1:
                LIBS.append(prot_det_lib)
                CLS.append('pdet')
            elif prot_sim == 1 and prot_det == 0:
                LIBS.append(prot_sim_lib)
                CLS.append('psim')
                
            if dna_cust_lib != '':
                LIBS.append(dna_cust_lib)
                CLS.append('ncus')
            elif dna_det == 1:
                LIBS.append(dna_det_lib)
                CLS.append('ndet')
            elif dna_sim == 1 and dna_det == 0:
                LIBS.append(dna_sim_lib)
                CLS.append('nsim')
                
            if custom_lib != '':
                LIBS.append(custom_lib)
                CLS.append('cust')
                
            if single == 1:
                base = base+' @class '
                for clfile in CLS: 
                    add = "DISICL_"+clfile+"_"+name2+".out "
                    base = base + add
            print(base)
                    
            if LIBS != []:
                add = '@lib'
                for libs in LIBS:
                    add = add+' '+libs
                add = add+' '
                print(add)
                visadd = visadd + add
                
            if vis != 1:
                add = '@write '+str(vis)+' '
                print(add)
                visadd = visadd + add
            if vistype != 1:
                add = '@type '+str(vistype)+' '
                print(add)
                visadd = visadd + add
            if stride != '':
                add = '@stride '+str(stride)+' '
                print(add)
                visadd = visadd + add
            if debug == 1:
                add = '@debug '
                print(add)
                visadd = visadd + add
#            if timecodes != '':
#                add = '@time '+str(TIME[0])+','+str(TIME[1])+' '
#                print add
#                visadd = visadd + add
            print("Starting visualization in PyMol")
            vis_command = vis_prog+' '+base+visadd+' >> '+logfile
            if win == 0:
                os.system('echo "'+vis_command+'" >>' +logfile)
            else:
                os.system('echo '+vis_command+' >>' +logfile)
            os.system(vis_command)
#checking module success
            l = open(logfile,'rb')
            modtime = 0
            timeline = ''
            if single == 1:
                base = name2
            else:
                base = name
            for line0 in l:
                line = str(line0.decode("ascii"))
                last = line
                if line.startswith("Script finished"):
                    timeline = line 
            l.close()
            if not (last.find('Visual') != -1 and last.find(base) != -1):
                failmark = 6
                print(30*' '+'failed')   
                bad = bad+' vis'
            else:
                if single == 1:
                    os.system(mv+' '+'frame_1.png '+base+'.png')
                    rename = 'renaming frame_1.png to '+base+'.png'
                    if win == 0:
                        os.system('echo "'+rename+'" >> '+logfile)
                    else:    
                        os.system('echo '+rename+' >> '+logfile)
                if timeline != '':    
                    try:
                            parts = timeline.split()
                            modtime = float(parts[2]) 
                            module_times = module_times+modtime
                    except Exception:
                        print("\nWarning: runtime on module was not found")
                print(30*' '+'done')
            print('')
#closing database entry, and moving on                
        if single == 1:
            if failmark == 0:
                data = (name2,rnum,'yes','  ')
            else:
                data = (name2,rnum,'no ',bad)
            dstring = '%-20s'+4*' '+'%5d'+4*' '+'%3s'+4*' '+'%6s\n'
            o.write(dstring % data)
            ecnt = ecnt + 1
        if traj == 1:
            ecnt = ecnt + mnum
        

print(separator)
print('number of files: '+str(cnt))
print('number of models: '+str(fcnt))
print('number of valid residues: '+str(rcnt))
print('number of classified segments:'+str(scnt))
ftime = os.times()
runtime = ftime[4]-stime[4]
if runtime == 0.0:
    runtime = (ftime[0]+ftime[1])-(stime[0]+stime[1])+module_times
print("program runtime: %1.2f  seconds (%1.2f h)" % (runtime,(float(runtime)/3600)))

#i.close()
if single == 1:
    print('DISICL program finished sucessfully! Results listed in analysed.list ')
    o.close()
else:
    print('DISICL program finished successfully! Results in '+GoalDir)
