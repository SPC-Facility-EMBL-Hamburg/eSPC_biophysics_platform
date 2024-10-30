#!/usr/bin/env python
#

#License information:
#    DISICL structure preparation script for pdb-like molecular trajectories .
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

print("DISICL_prep script starts:")
print(workdir+"\n")

#usage
usage0 = '************************************\n'
usage1 = 'DISICL structure preparation script\nUsage for sorting pdb-like trajectories:\n\n'
usage2 = 'DISICL_prep.py <trajectory file> <options>\n\n'
usage3 = 'Trajectory file options:\n   <file-name>\n\t   Single trajectory file in pdb 1.0 compatible format (default)\n'
usage4 = '   @file <file-name>\n\t   Give argument file (no other arguments read from command line)\n'
usage5 = '   @traj <file1> <file2> ...\n\t   Give multiple files as a single trajectory\n'
usage6 = '   @single <file1> <file2> ...\n\t   Give multiple files as a database of independent structures\n\n'
usage7 = 'Additional options:\n   @refine <0/1>\n\t   Write refined .pdb file:\n\t   1 - yes (default)\n\t   2 - no\n\n'
usage8 = '   @type <t>\n\t   switch between molecule types, t can take values:\n\t   0 - custom\n\t   1 - proteins (default)\n'
usage9 = '\t   2 - nucleotides (DNA-RNA)\n\t   3 - both proteins and nucleotides\n\n   @ext <.str>\n\t   Give file extension (default .pdb)\n\n'
usage10 = '   @path <path>\n\t   Set different source directory for structure files (default ./)\n\n'
usage11 = '   @dest <path>\n\t   Set different goal directory for results (relative, default ./)\n\n'
usage12 = '   @protres <str1,str2,...>\n   @nuclres <str1,str2,...>\n   @custres <str1,str2,...>\n'
usage13 = '\t   Add standard protein, nucleotide, or custom residue codes\n\t   (comma separated)\n\n'
usage14 = '   @rename <0/1>\n\t   Control DNA/RNA renaming to standard format\n\t   1 - on (default)\n\t   0 - off\n\n'
usage15 = '   @exclude <str1,str2,...>\n\t   Exclude residue codes from analysis (comma separated)\n\n'
usage16 = '   @filter <f>\n\t   Set filter level to restrict residues, f can take values:\n\t   0 - no filter\n\t   1 - only excluded\n\t'
usage17 = '   2 - all nonstandard residues\n\t   3 - all nonstandard residues and alternative conformations (default)\n\t'
usage18 = '   See script exclusion section for filter lists\n\n   @list <str>\n\t   Set library list name\n\n   @help\n\t   Print this message\n'

usage = usage0+usage1+usage2+usage3+usage4+usage5+usage6+usage7+usage8+usage9+usage10
usage = usage+usage11+usage12+usage13+usage14+usage15+usage16+usage17+usage18+usage0

#######System variables, should be edited upon installation#####
#local definitions for Windows / Linux commands:
win = 0
# 0 = use linux format, 1 = assume windows format

######Nothing should be modified from this point, apply changes at your own risk#######

if win == 1:
    rm = "del"
    cp = "copy"
    mv = "move"
else:
    rm = "rm"
    cp = "cp"
    mv = "mv"
    
#hardcoded residue names:
#for proteins
PROTRES1 = ['ALA','ASP','ASN','ARG','CYS','GLY','GLU','GLN','HIS','ILE']
PROTRES2 = ['LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
CUSTOM_PROT = []
#for DNA
DNARES1 = ['DA','DC','DG','DT','DU']
DNARES2 = ['DADE','DCYT','DGUA','DTHY','DURA']
CUSTOM_NUCL = []
#for RNA
RNARES1 = ['A','C','G','T','U']
RNARES2 = ['ADE','CYT','GUA','THY','URA']
#other residues
CUSTOM_RES = []

#hardcoded exclusions:
EXCLUDE1 = ['SOLV','H2O','WAT','SPC','TIP3','HEM','FAD','NA','CL','K']
EXCLUDE2 = ['Ca','MG','ZN','NA+','CL-','K+','CA2+','MG2+','ZN+','ZN2+']
CUSTOM_EXCLUDE = []
   
    
#command flags for modules
refflag = 1
resfilt = 3
ext = '.pdb'
protein = 1
dna = 0
Dest = ''
Path = ''
listname = 'database.list'
customprot = ''
customnuc = ''
customres = ''
customex = ''
TRAJFILES = []
ARGUMENTS = []
single = 0
rename = 1

#command line interface:
failmark = 0
flaglist = ['type','traj','single','refine','filter','ext','path','dest','list','file','help','protres','nuclres','custres','exclude','rename']
flag = ''

#if @file flag is provided read in arguments from there:
infofile = 0
infile = ''
for arg in sys.argv:
    if arg.find('@file') != -1:
        infofile = 1
    elif infofile == 1:
        infile = arg
        infofile = 0


if infile == '':
    for arg in sys.argv:
        ARGUMENTS.append(arg)
else:
    ARGUMENTS.append(infile)
    if os.path.isfile(infile) == True:
        i = open(infile,"rb")
        for line0 in i:
            line = str(line0.decode("ascii"))
            if not line.startswith('#'):
                parts = line.strip("\n").split()
                for part in parts:
                    if not part == '':
                        ARGUMENTS.append(part)
        i.close()
    else:
            print("\nERROR: Argument file not found!")
            failmark = 1
    
        
#read in and interpret arguments:
argnum = len(ARGUMENTS)-1
print('number of arguments: '+str(argnum))
if not argnum == 0:
    print('flags found:')
    cnt = 1
    while cnt <= argnum:
        arg = ARGUMENTS[cnt]
        if arg == '@help':
            failmark = 1
        if arg.startswith('@'):
            flag = arg.strip('@')
            print(flag)
            if not flag in flaglist:
                print('\nERROR: Unkown flag, program stops!')
                failmark = 1
                flag = ''
        elif (flag == '' and cnt == 1) or flag == 'traj':
            TRAJFILES.append(arg)
        elif flag == 'single':
            single = 1
            TRAJFILES.append(arg)
        elif flag == 'type':
            arg = int(arg)
            flag = ''
            if arg  == 1:
                protein = 1
                dna = 0
            elif arg == 2:
                dna = 1
                protein = 0
            elif arg == 3:
                protein = 1
                dna = 1
            elif arg == 0:
                protein = 0
                dna =   0
            else:
                print('\nERROR: type code not recognized!')
                falimark = 1
        elif flag == 'refine':
            refflag = int(arg)
            if refflag < 0 or refflag > 1:
                print('\nERROR: refine code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'filter':
            resfilt = int(arg)
            if resfilt < 0 or resfilt > 3:
                print('\nERROR: filter code not recognized!')
                falimark = 1
            flag = ''
        elif flag == 'ext':
            ext = arg
            flag = ''
        elif flag == 'path':
            Path = arg
            flag = ''
        elif flag == 'dest':
            Dest = arg
            flag = ''
        elif flag == 'list':
            listname = arg
            flag = ''
        elif flag == 'protres':
            customprot = arg
            flag = ''
        elif flag == 'nuclres':
            customnuc = arg
            flag = ''
        elif flag == 'custres':
            customres = arg
            flag = ''
        elif flag == 'exclude':
            customex = arg
            flag = ''
        elif flag == 'rename':
            rename = int(arg)
            flag = ''
        cnt = cnt + 1
        
if single == 1:
    for arg2 in ARGUMENTS:
                if arg2 == "@traj":    
                    print("\nERROR: @single and @traj flags cannot be provided at the same time!")
                    failmark = 1

if customprot != '':
    print("\ncustom protein residue codes added:")
    parts = customprot.split(',')
    for code in parts:
        print(code)
        CUSTOM_PROT.append(code)
        
if customnuc != '':
    print("\ncustom nucleotide residue codes added:")
    parts = customnuc.split(',')
    for code in parts:
        print(code)
        CUSTOM_NUCL.append(code)
        
if customres != '':
    print("\ncustom residue codes added:")
    parts = customres.split(',')
    for code in parts:
        print(code)
        CUSTOM_RES.append(code)
        
if customex != '':
    parts = customex.split(',')
    print("\ncustom exclusions added:")
    for code in parts:
        print(code)
        CUSTOM_EXCLUDE.append(code)
 
print(CUSTOM_RES,protein,dna)  
if CUSTOM_RES == [] and protein == 0 and dna == 0:     
    print("\nError: custom classification with @type 0 invoked, but no custom residues defined!")
    print("Use the @custres flag to give custom residue names!")
    failmark = 1
        
if failmark == 1:
    print(usage)
    sys.exit()

if not TRAJFILES == []:    
    print("\nFiles to process:")
    for filename in TRAJFILES:
        print(filename)
else:
    print("\nERROR: no input files detected!")
    print(usage)
    sys.exit()

Datadir1 = os.path.join(workdir,Path)
if Path == '':
    DataDir = workdir
elif os.path.isdir(Datadir1) == True: 
    DataDir = Datadir1
elif os.path.isdir(Path) == True:
    Datadir = Path
else:
    print("Error: Could not find source directory!")
    print(usage)
    sys.exit()

os.chdir(DataDir)
Filelist = os.listdir(DataDir)
    
if Dest == '':
    GoalDir = workdir
else:
    GoalDir = os.path.join(workdir,Dest)


# finalizing exclusion and selected residue lists:
RES = []
PROTRES = [PROTRES1,PROTRES2,CUSTOM_PROT]
NUCLRES = [DNARES2,RNARES2,CUSTOM_NUCL]
if protein  == 1:
    for Reslist in PROTRES:
        for res in Reslist:
            RES.append(res)
if dna == 1:
    for Reslist in NUCLRES:
        for res in Reslist:
            RES.append(res)
if CUSTOM_RES != []:
    for res in CUSTOM_RES:
        RES.append(res)
    

EXCLUDE = []
EXCLUSIONS = [EXCLUDE1,EXCLUDE2,CUSTOM_EXCLUDE]
for Exlist in EXCLUSIONS:
    for ex in Exlist:
        EXCLUDE.append(ex)
    

#checking directory and file list
LIST = []
EXLIST = []

if refflag == 1:
    if os.path.isdir(GoalDir) == True:
        print('\ndestination directory found:')
        print(GoalDir)
    else:
        print('\ndestination directroy not found, attemp to create it')
        print(GoalDir)
        os.system('mkdir '+Dest)
        if os.path.isdir(GoalDir) == True:
            print('\ndestination directory created:')
            print(GoalDir)
        else:
            print('Error, could not create destination directory, program stops!')
            print(usage)
            sys.exit()
        
print('\nValid Files found:')
cnt = 0
totframe = 0
for File in Filelist:
    if (File in TRAJFILES) and (os.path.isfile(File) == True):
        ext_length = len(ext)
        namebase  = File[0:(-1*ext_length)]
        file_extension = File[-1*ext_length:]
        if file_extension != ext:
            print('\nERROR:',File,"\nexpected extension: "+ext)
            print("found: "+file_extension)
            print("use @ext if extension is correct")
            print(usage)
            sys.exit()
        print(namebase)
        line = [namebase]
        LIST.append(line)
        cnt = cnt + 1
if cnt == 0:
    print('\nERROR: no valid files are detcted in source directory !!!')
    print(usage)
    sys.exit()

filenum = cnt
print('\nTotal number of files: '+str(filenum))

#reading in file information for residues:
filenum = len(LIST)
if refflag == 1:
    outfile = os.path.join(GoalDir,LIST[0][0]+'_'+"DISICL"+ext)
    o = open(outfile,'wb')

rnum = 0
enum = 0
Rlist = ''
cnt = 0
for File in LIST:
    infile = File[0]+ext
    i = open(infile,'rb')
    fcnt = 1
    if refflag == 1:
        if single == 1:
            outfile = os.path.join(GoalDir,File[0]+'_'+str(fcnt)+ext)
            o = open(outfile,'wb')
    R0 = -100
    R1 = -99
    R2 = 0
    L2 = 0
    Reslist = []
    TMPLIST = []
    First_data = [0,0,"",[]]
    rmin = 100000
    rmax = 0
    rcnt = 0
    ecnt = 0
    endflag = 1
    oflag = 0
    oadd = 0
    chainflag = 0
    writeflag = 0
    newflag = 0
    for line0 in i:
        line = str(line0.decode("ascii"))
#handling model separations
        if line.startswith('END'):
            if refflag == 1:
                enum = ecnt
                rnum = rcnt
                fcnt = fcnt + 1
                if rnum != 0:
                    o.write('ENDMDL\n'.encode("ascii"))
                if single == 1:
                    o.close()
                    outfile = os.path.join(GoalDir,File[0]+'_'+str(fcnt)+ext)
                    o = open(outfile,'wb')
#                o.write('#filtered model for '+File[0]+ext+'\nMODEL'+8*' '+str(fcnt)+'\n')
                R0 = -100
                R1 = -99
                R2 = 0
                L2 = 0
                Reslist = []
                rmin = 100000
                rmax = 0
                rcnt = 0
                ecnt = 0
                endflag = 1
                oflag = 0
                chainflag = 0
                oadd = 0
                newflag = 0
                Reslist  = []
#writing out new residue lists
                chaininfo = ('c',0)
                TMPLIST.append(chaininfo)
                cdmin = ''
                cdmax = ''
                cpmin = ''
                cpmax = ''
                comin = ''
                comax = ''
                protlist = ''
                dnalist = ''
                otherlist = ''
                for tmp in TMPLIST:
#                    print tmp
#if residue marked as protein
                    if tmp[0] == 'p':
                        if cpmin == '':
                            cpmin = tmp[1]
                        elif cpmin > tmp[1]:
                            cpmin = tmp[1]
                        if cpmax == '':
                            cpmax = tmp[1]
                        elif cpmax < tmp[1]:
                            cpmax = tmp[1]
#residue marked as dna
                    elif tmp[0] == 'd':
                        if cdmin == '':
                            cdmin = tmp[1]
                        elif cdmin > tmp[1]:
                            cdmin = tmp[1]
                        if cdmax == '':
                            cdmax = tmp[1]
                        elif cdmax < tmp[1]:
                            cdmax = tmp[1]
#residue marked as other
                    elif tmp[0] == 'o':
                        if comin == '':
                            comin = tmp[1]
                        elif comin > tmp[1]:
                            comin = tmp[1]
                        if comax == '':
                            comax = tmp[1]
                        elif comax < tmp[1]:
                            comax = tmp[1]
#when chain ends write out minimal and maximal values
                    elif tmp[0] == 'c':
                        if cpmax != '' and cpmin != '' and cpmax > cpmin+1:
                            if protlist == '':
                                pdata = str(cpmin+1)+'-'+str(cpmax-1)
                            else:
                                pdata = ','+str(cpmin+1)+'-'+str(cpmax-1)
                            protlist = protlist + pdata
                        cpmin = ''
                        cpmax = ''
                        if cdmax != '' and cdmin != '' and cdmax > cdmin:
                            if dnalist == '':
                                ddata = str(cdmin)+'-'+str(cdmax-1)
                            else:
                                ddata = ','+str(cdmin)+'-'+str(cdmax-1)
                            dnalist = dnalist + ddata
                        cdmin = ''
                        cdmax = ''
                        if comax != '' and comin != '' and comax >= comin:
                            if otherlist == '':
                                odata = str(comin)+'-'+str(comax)
                            else:
                                odata = ','+str(comin)+'-'+str(comax)
                            otherlist = otherlist + odata
                        comin = ''
                        comax = ''
#               print Newreslist
#		print TMPLIST
                Newreslist = [protlist,dnalist,otherlist]
                Rlist = '['+str(Newreslist[0])+'] ['+str(Newreslist[1])+'] ['+str(Newreslist[2])+']'
                print(Rlist)
#		save the first residue list in the file 
                if fcnt == 2:
                        TMPLIST0 = []
                        for rcode in TMPLIST:
                              TMPLIST0.append(rcode)
                        First_data = [rnum, enum, Rlist, TMPLIST0]
#			print "data on first frame",First_data 
                TMPLIST = []
#handling chain ends
        elif line.startswith('TER'):
            if refflag == 1 and chainflag != 1:
                o.write(line.encode("ascii"))
            chainflag = 1
            chaininfo = ('c',0)
            TMPLIST.append(chaininfo)
#            print R2,aRN,rmax,rmin
            
#handling chain ends
        elif line.startswith('TER'):
            if refflag == 1 and chainflag != 1:
                o.write(line.encode("ascii"))
            chainflag = 1
            chaininfo = ('c',0)
            TMPLIST.append(chaininfo)
#            print R2,aRN,rmax,rmin
            
            if rmax != 0 and (aRN+oadd) <= R2 :
                writeflag = 2
#                print "chain ends:",R2,rmin,rmax
#reading in line information
        elif (line.startswith('ATOM') or line.startswith('HETATM')):
            banflag = 0
            try:
                anum = int(line[6:11].strip(' '))
                aID = str(line[11:16].strip(' '))
                aRI = str(line[16:21].strip(' '))
                aCI = str(line[21])
                try:
                    aRN = int(line[22:27].strip(' '))
                    lflag = 0
#                    print aRN
                except ValueError:
                    RN1 = line[22:27].strip(' ')
                    last = RN1[-1]
                    RN2 = RN1.strip(last)
#handling RNA loops                    
                    if last.isalpha() == True and RN2.isdigit() == True:
                        RN2 = int(RN2)
                        Alpha = ['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
                        abc = len(Alpha)
                        cnt = 1
                        while cnt < abc:
                            code = Alpha[cnt-1]    
                            if last == code:
                                nadd = cnt
                            cnt = cnt + 1
                        aRN = RN2+nadd
                        lflag = 1
                        oflag = 1
#                        print RN1, last, aRN
                    else:
                        lflag = 1
                        aRN = int(line[22:27].strip(' '))
                ax = float(line[28:39].strip(' '))
                ay = float(line[38:46].strip(' '))
                az = float(line[46:55].strip(' '))
#filtering nonstandard residues                
                if ( (aRI in EXCLUDE or line.startswith('HETATM')) and resfilt > 0 ):
                    banflag = 1
                if resfilt > 1:
                    banflag = 1
                    for res in RES:
                        if aRI.find(res) != -1:
                            banflag = 0
                    if dna == 1 and banflag == 1:
                        for res in DNARES1:
                            if res == aRI:
                                banflag = 0
                        for res in RNARES1:
                            if res == aRI:
                                banflag = 0
                        
                    if resfilt > 2:
                        if banflag == 0:
                            altpos1 = line[16]
                            if not(altpos1 == ' ' or altpos1 == 'A' or (altpos1 == 'D' and dna == 1)):
#                                print 'line excluded:'
#                                print line.strip('\n')
#                                print anum,',',aID,',',aRI,',',aCI,',',aRN,',',altpos1,'\n'
                                banflag = 1
#handling DNA uniform codes (A,C,T,G,U)
                if banflag == 0 and dna == 1 and rename == 1:
                    bcnt = 0
                    bnum = len(RNARES1)
                    while bcnt < bnum:
                        if (aRI == RNARES2[bcnt]) or (aRI == DNARES1[bcnt]) or (aRI == DNARES2[bcnt]):
                            first = line[0:16]
                            last = line[21:]
                            ncode = 2*' '+RNARES1[bcnt]+2*' '
                            newline = first+ncode+last
                            line = newline
                        bcnt = bcnt + 1
                    if resfilt == 3 and aID.find('*') != -1:
                        first = line[0:11]
                        last = line[16:]
                        nID = aID.strip('*')+"'"
                        if len(nID) < 5:
                            ins = '%4s' % nID
                            newline = first+' '+ins+last
                        else:
                            ins = '%5s' % nID
                            newline = first+ins+last
                        line = newline
#handling residue code changes
                if banflag == 0:
                    if (R2 != aRN and oflag == 0) or (R2 != (aRN+oadd) and oflag == 1) or (R2 == (aRN+oadd) and L2 != lflag):
                        if chainflag == 1:
#                            print "chain change",aRN,rmax
                            if aRN < rmax:
                                oadd = (rmax+1)-aRN
                                oflag = 1
                                if refflag == 0 and oflag == 1:
                                    break
                                    o.close()
                            chaininfo = ("c",0)
                            TMPLIST.append(chaininfo)
                            chainflag = 0
                        if aRN <= 0 and oflag == 0:
                            oadd = (1-aRN)
                            oflag = 1
                        if aRN < R2 or (aRN == R2 and L2 > lflag):
                            oadd = (R2+1-aRN)
                            oflag = 1
#updating neighbor list                            
                        if newflag == 2:
                            R0 = -100
                            R1 = -99
                            newflag = 0
                        elif newflag == 1:
                            R0 = -100
                            R1 = R2
                            newflag = 0
                        else:
                            R0 = R1
                            R1 = R2
                        if oflag == 0:
                            R2 = aRN
                            L2 = lflag
                        elif oflag == 1:
                            R2 = (aRN+oadd)
                            L2 = lflag
#writing residues to TEMPLIST:                        
                        NewresID = str(line[16:21].strip(' '))
                        Newresnum = R2
                        OldresID = aRI
                        Oldresnum = aRN
                        Oldchain = aCI
                        newcode = (5-len(NewresID))*' '+NewresID+2*' '+str(Newresnum)+(6-len(str(Newresnum)))*' '+Oldchain
                        oldcode = (5-len(OldresID))*' '+OldresID+2*' '+str(Oldresnum)+(6-len(str(Oldresnum)))*' '+Oldchain
                        if newcode != oldcode:
                            print('residue:'+oldcode+' renamed to:'+newcode)
                        if R2 != (R1+1):
                            chaininfo = ("c",0)
                            TMPLIST.append(chaininfo)
                        restype = ''
                        if protein == 1 and restype == '':
                            for Codes in PROTRES:
                                for amino in Codes:
                                    if NewresID.find(amino) != -1:
                                        restype = 'p'
                        if dna == 1 and restype == '':
                            DNACODES = [DNARES1,RNARES1]
                            for Codes in DNACODES:
                                for base in Codes:
                                    if NewresID == base:
                                        restype = 'd'
                            for Codes in NUCLRES:
                                for base in Codes:
                                    if NewresID.find(base) != -1:
                                        restype = 'd'
                        if CUSTOM_RES != [] and restype == '':
                                for other in CUSTOM_RES:
                                    if NewresID == other:
                                        restype = 'o'
                        if restype == '':
                            restype = 'c'
                        resinfo = (restype,R2)
                        TMPLIST.append(resinfo)
#                        print R0,R1,R2
#getting chain ends
                        if R1 <= rmin and (R0 > 0 and (R0+1) == R1):
                            rmin = R1
                        if R1 >= rmax and (R2 == R1+1):
                            rmax = R1
                        if R2 != (R1+1):
#                            print R1,R2
                            if rmax != 0 and rmax != R2 and (R1 > 0 and R0 > 0) and (rmax >= rmin):
                                writeflag =1
#                                print 'Residue break:',R2,R1,R0,rmin,rmax
#                                print aRN,oadd
#                                print line
                            else:
                                rmin = 10000
                                rmax = R2
                        rcnt = rcnt + 1
#writing out line if needed                    
                    if refflag == 1:
                        if endflag == 1:
                            if single == 1:
                                modnum = fcnt
                            else:
                                modnum = fcnt + totframe
                            new_line = '#filtered model for '+File[0]+ext+'\nMODEL'+8*' '+str(modnum)+'\n'
                            o.write(new_line.encode("ascii"))
                            print('\nModel '+str(fcnt)+' from '+File[0])
                            endflag = 0
                        if oflag == 0:
                            o.write(line.encode("ascii"))
                            
                        elif oflag == 1:
                            rescode = line[22:27]
                            parts = line.split(rescode)
                            if len(parts) != 2:
                                print('Warning ! in file '+File[0]+' line:\n',line)
                                print('trying alternative method')
                                first = line[0:22]
                                last = line[27:]
                                digit = len(str(R2))
                                if len <= 3:
                                    text = (4-digit)*' '+str(R2)
                                    midle = ('%4s ' % text)
                                else:
                                    text = str(R2)
                                    midle = ('%4s ' % text)
                                newline = first+midle+last
                                print(newline)
                                print(line)
                                if len(newline) == len(line):
                                    o.write(newline.encode("ascii"))
                                else:
                                    
                                    print('Error!! could not modify residue ID!!')  
                                    o.write(('#'+line).encode("ascii"))
                            else:
                                first = parts[0]
                                last = parts[1]
                                midle = ('%5d' % R2)
                                newline = first+midle+last
                                o.write(newline.encode("ascii"))
#listing exclusions
                else:   
#                    print line
#                    print ' line was excluded'
                    ecnt = ecnt + 1
                    if not aRI in EXLIST:
                        EXLIST.append(aRI)
            except ValueError:
                print(line,'\n could not be read in')
#                print anum,',',aID,',',aRI,',',aCI,',',aRN,',',ax,',',ay,',',az,'\n'
#writing residue list data
        if writeflag != 0:
            if ( rmin < rmax) and rmin != 100000 and rmax != 0:
                data = str(rmin)+'-'+str(rmax)
            elif rmax != 0:
                data = str(rmax)
            if rmax != 0 and rmax != R2 and (R1 > 0 and R0 > 0) and (rmax >= rmin):
                newflag = writeflag
                rmin = 100000
                rmax = R2
            writeflag = 0
         
#finish up things if no ENDS were added to the file!   
    if endflag == 0:
        enum = ecnt
        rnum = rcnt
        chaininfo = ('c',0)
        TMPLIST.append(chaininfo)
        cdmin = ''
        cdmax = ''
        cpmin = ''
        cpmax = ''
        comin = ''
        comax = ''
        protlist = ''
        dnalist = ''
        otherlist = ''
        for tmp in TMPLIST:
#                    print tmp
#if residue marked as protein
            if tmp[0] == 'p':
                if cpmin == '':
                    cpmin = tmp[1]
                elif cpmin > tmp[1]:
                    cpmin = tmp[1]
                if cpmax == '':
                    cpmax = tmp[1]
                elif cpmax < tmp[1]:
                    cpmax = tmp[1]
#residue marked as dna
            elif tmp[0] == 'd':
                if cdmin == '':
                    cdmin = tmp[1]
                elif cdmin > tmp[1]:
                    cdmin = tmp[1]
                if cdmax == '':
                    cdmax = tmp[1]
                elif cdmax < tmp[1]:
                    cdmax = tmp[1]
#residue marked as other
            elif tmp[0] == 'o':
                if comin == '':
                    comin = tmp[1]
                elif comin > tmp[1]:
                    comin = tmp[1]
                if comax == '':
                    comax = tmp[1]
                elif comax < tmp[1]:
                    comax = tmp[1]
#when chain ends write out minimal and maximal values
            elif tmp[0] == 'c':
                if cpmax != '' and cpmin != '' and cpmax > cpmin+1:
                    if protlist == '':
                        pdata = str(cpmin+1)+'-'+str(cpmax-1)
                    else:
                        pdata = ','+str(cpmin+1)+'-'+str(cpmax-1)
                    protlist = protlist + pdata
                cpmin = ''
                cpmax = ''
                if cdmax != '' and cdmin != '' and cdmax > cdmin:
                    if dnalist == '':
                        ddata = str(cdmin)+'-'+str(cdmax-1)
                    else:
                        ddata = ','+str(cdmin)+'-'+str(cdmax-1)
                    dnalist = dnalist + ddata
                cdmin = ''
                cdmax = ''
                if comax != '' and comin != '' and comax >= comin:
                    if otherlist == '':
                        odata = str(comin)+'-'+str(comax)
                    else:
                        odata = ','+str(comin)+'-'+str(comax)
                    otherlist = otherlist + odata
                comin = ''
                comax = ''
#                print Newreslist
        Newreslist = [protlist,dnalist,otherlist]
        Rlist = '['+str(Newreslist[0])+'] ['+str(Newreslist[1])+'] ['+str(Newreslist[2])+']'
        TMPLIST = []
#    Rlist = Reslist
#    Rlist = str(Newreslist[0])+' '+str(Newreslist[1])

#Collect final file statistics:
#check number of frames: 
    if fcnt > 1:
        framenum = fcnt -1
    else:
        framenum = 1
# check if first/last frame included residues:	
    if rnum  != 0:
        note = "(last frame)"
    elif rnum == 0 and First_data[0] != 0:
        rnum, enum, Rlist, note = First_data[0],First_data[1],First_data[2], "(first frame)"
        framenum -= 1
    else:
        print("Warning, no valid residues were detected for the first and last frames!")
        rnum, enum, Rlist, note = First_data[0],First_data[1],First_data[2], "first frame"
        framenum -= 1
        print(First_data)

    print('\n'+File[0]+ext+' ('+str(cnt+1)+'/'+str(filenum)+'), number of counted residues: %1d '%rnum + note) 
    print(Rlist,framenum,rnum,enum,'\n')
    File.append(Rlist)
    File.append(framenum)
    File.append(rnum)
    File.append(enum)
    totframe = totframe + framenum
    if refflag == 1:
        if endflag == 0:
            o.write('ENDMDL\n')
        if single == 1:
            o.close()
    i.close()
#removing last file if empty
    if single == 1:
        lastfile = os.path.join(GoalDir,File[0]+'_'+str(fcnt)+ext)
        u = open(lastfile,'rb')
        uflag = 0
        for line0 in u:
            line = str(line0.decode("ascii"))
            if line.startswith('ATOM') or line.startswith('HETATM'):
                uflag = 1
                break
        u.close()
        if uflag != 1:
            os.system(rm+' '+lastfile)
    cnt = cnt + 1
    
# summarize trajectory information:
if single == 0:
    totframe = 0
    prevlist = ''
    cnt = 0
#    print "\nList:\n",LIST
    while cnt < filenum:
        File = LIST[cnt]
        name = File[0]
        reslist = str(File[1])
        length = len(reslist)
        framenum = File[2]
        rnum =  File[3]
        enum = File[4]
        if cnt == 0:
            reslist0 = reslist
            resnum0 = rnum
            print('\n',resnum0,"residues found in frame 1:",reslist0)
        if (prevlist != '') and reslist != '':
            if reslist != prevlist:
                print("WARNING: residue list for frame "+str(cnt)+" is different from previous frame!")
                print(prevlist," not equals to", reslist)
        totframe = totframe + framenum
        cnt = cnt + 1
    data = (cnt,totframe,reslist0)
    summary1 = 'DISICL'+5*' '+'files'+5*' '+'frames'+5*' '+'residues'+5*' '+'residue list\n'
    summary2 = 'DISICL'+7*' '+str(cnt)+10*' '+str(totframe)+10*' '+str(resnum0)+5*' '+reslist0
    summary = summary1+summary2
    print('\n'+summary)
    o.write((summary+'\n').encode("ascii"))
    o.close()
#writing out file list if individual snapshots are required
if single == 1:
    os.chdir(GoalDir)
    listfile = listname    
    l = open(listfile,'wb')
    l.write('#Database list file for '+str(filenum)+' analysed pdb files:\n')
    l.write('# pdb name'+5*' '+'models'+5*' '+'residues'+4*' '+'exclded atoms'+15*' '+'residue list'+'\n')
    for File in LIST:
        name = File[0]
        reslist = str(File[1])
        length = len(reslist)
        framenum = File[2]
        rnum =  File[3]
        enum = File[4]
        data = (name,framenum,rnum,enum,reslist)
        lstring = '%-15s'+1*' '+'%3d'+8*' '+'%5d'+6*' '+'%7d'+3*' '+'%30s\n'
        l.write(lstring % data)
    
    l.write('\n#list of excluded residues:\n')
    for ex in EXLIST:
        l.write ('#'+ex+'\n')
    l.close()
rtime = os.times()
runtime = rtime[4]-stime[4]
if runtime == 0.0:
    runtime = rtime[0]+rtime[1]-(stime[0]+stime[1])
print("\nprogram runtime: %1.2f  seconds (%1.2f h)" % (runtime,(float(runtime)/3600)))
print("File list and/or refined files written in "+GoalDir)
print("Preparation script finished successfully! ("+DataDir+")")
