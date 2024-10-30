#!/usr/bin/env python
#

#License information:
#    DISICL secondary structure classification program for biopolymers.
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

workdir = os.getcwd()
stime = os.times()

print("DISICL_dbssp program starts")
print(workdir+"\n")

#usage:
usage0 = '********************************************\n'
usage1 = 'DISICL Dihedral Based Secondary Structure Partioning\n'
usage2 = 'Usage for classification program is:\n\nDISICL_dbssp.py [input file] [library specification] [options]\n\n'
usage3 = 'Input file options:\n   <file-name>\n\t   Single dihedral time series file (.dih)\n'
usage4 = '   @file <file-name>\n\t   Give argument file (no other arguments read from command line)\n\nLibrary specification:\n'
usage5 = '   @lib <library file>\n\t   Read in region and segment definitions from library file\n\t    Use @lib help for further info on library definitions\n\n'
usage6 = 'Additional options:\n   @res <res,res0-resmax>\n\t   Define residue range for analyisis (comma separated)\n\t'
usage7 = '   arguments are integers, or ranges (- separated)\n\t   if no range is given all valid residues are analysed\n\n'
usage8 = '   @time <time,timestep>\n\t   Define timestamps for trajectory (comma separated)\n\t   use floats for starting time and timestep (default 0.0,1.0)\n\n'
usage9 = '   @path <path>\n\t   Set different source directory for structure files (default ./)\n\n'
usage10 = '   @debug\n\t   Use debug mode (write out calculation specs at every frame)\n\n'
usage11 = '   @command\n\t   Flag marks the argument file identical with the input file\n\t   arguments are read until the first line starting with an "&"\n\t'
usage12 = '   For argument files ony, do not use @command in the command line!\n\n   @help\n\t   Print this message\n'

usage = usage0+usage1+usage2+usage3+usage4+usage5+usage6+usage7+usage8+usage9+usage10+usage11+usage12+usage0

# default variables:
filebase = ["Dih_",".out"]
resdata = ''
ARGUMENTS = []
outbase =  ["DISICL",".out"]
outbase2 = ["DISICL",".stat"]
defflag = "prot"
path = ''
infile = ''
TIME = [1.0,1.0]
timecodes = ''
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

failmark = 0
argnum = len(ARGUMENTS)-1
print('number of arguments: '+str(argnum))
print('flags found:')
flaglist = ['lib','path','res','file','time','debug','command','help']
flag = ''
cnt = 1
if argnum > 0:
	while cnt <= argnum:
		arg = ARGUMENTS[cnt]
		if arg == "@command":
			infile = ARGUMENTS[0]
		elif arg == "help":
			failmark = 1
			print("help")
			flag = ''
		elif arg == "@debug":
			debug = 1
			print("debug")
			flag = ''
		elif arg.startswith("@"):
			flag = arg.strip('@')
			print(flag)
			if not flag in flaglist:
				print('\nERROR: Unkown flag, program stops!')
				failmark = 1
				flag = ''
		elif flag == '':
			infile = arg
		elif flag == "res":
			resdata = arg
			flag = ''
		elif flag == "path":
			path = arg
#			os.chdir(path)
			flag = ''
		elif flag == "lib":
			defflag = arg
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
		cnt = cnt + 1	
elif argnum == 0 or infile == '':		
	print("no input file was given")
	print(usage)
	sys.exit(0)

print("\nChecking input file and path:")
if path != '':
#checking file and path validity
	newpath = os.path.join(workdir,path)
	if os.path.isdir(newpath) == True:
		path = newpath
		os.chdir(path)
		print("\nChanging directroy to path:\n"+str(path))
	elif os.path.isdir(path) == True:
		os.chdir(path)
		print("\nChanging directroy to path:\n"+str(path))
	else:
		print("\nError: could not find given path!")
		print(path)
		failmark = 1
if os.path.isfile(infile) == False:
	newinfile = os.path.join(workdir,infile)
	newinfile2 = os.path.join(path,infile)
	if os.path.isfile(newinfile) == True:
		print("Input file detected:")
		print(newinfile)
		infile = newinfile
	else:
		print("\nError: failed to locate input file!")
		print(infile)
		failmark = 1
else:
	print("\nInput file detected:",infile)
	
if not resdata == '':
	print("\n@res flag detected, residues found are:")
	data = resdata.split(',')
elif resdata == '' and infile != '':
	print("\nNo @res flag was given! Attempting to read DISICL remarks from file:")
	i = open(infile,'rb')
	remarks = ''
	for line0 in i:
		line = str(line0.decode("ascii"))
		if line.startswith("#DISICL"):
			remarks = line.strip('\n')
	if remarks != '':
		try:
			print(remarks)
			parts = remarks.split()
			tzero = float(parts[1])
			tstep = float(parts[2])
			framenum = int(parts[3])
			resnum = int(parts[4]) 
			reslist = parts[5].strip('[').strip(']')
			data = reslist.split(',')
			if timecodes == '':
				TIME[0] = tzero
				TIME[1] = tstep
				
			if data == '':
				print("Error: residue list in remarks is empty! Program stops")
				failmark = 1
		except ValueError:
			print("Error: Failed to read remark information! Program stops")
			failmark = 1
	else:
		"    No remarks detected, program stops!"
		failmark = 1
else:
	data = ''
	failmark = 1
RES = []
Newreslist = ''
for Res in data:
	try:
		if Newreslist == '':
			Newreslist = str(Res)
		else:
			Newreslist = Newreslist+','+Res
		if Res.find("-") == -1:
			RES.append(int(Res))
		else:
			resrange = Res.split("-")
			low = int(resrange[0])
			high = int(resrange[1])
			rcnt = low
			while rcnt <= high:
				RES.append(rcnt)
				rcnt = rcnt + 1
	
	except Exception:
		print("\nError: could not read in residue information ("+Res+")")
if not RES == []:
	print("\nResidue list updated to: ",Newreslist)
else:
	print("\nError: no valid residues were found, Program stops!")
	failmrak = 1
if TIME != [1.0,1.0] or timecodes != '':
	print("Timestamp information updated to: time0= "+str(TIME[0])+" timestep= "+str(TIME[1]))
else:
	print("Standard timecodes are used: time0= "+str(TIME[0])+" timestep= "+str(TIME[1]))
		
if (defflag != "prot") and (defflag != "help"):
	if os.path.isfile(defflag) == False:
		libpath1 = os.path.join(workdir,defflag)
		libpath2 = os.path.join(path,defflag)
		if os.path.isfile(libpath1) == True:
			print("\nLibrary file found:\n"+str(libpath1))
			defflag = libpath1
		elif os.path.isfile(libpath2) == True:
			print("\nLibrary file found:\n"+str(libpath2))
			defflag = libpath2
		else:
			print("\nError: definition library not found!")
			failmark = 1
	else:
		print("\nLibrary file found:\n"+str(defflag))
#boundariy definitions:
#Class difinitions:
class Regions:
	def __init__(self, data):
		dim = len(data)-1
		Boundaries = []
		cnt = 1
		while cnt <= dim:
			Boundaries.append(float(data[cnt]))
			cnt = cnt + 1
		self.name = str(data[0])	
		self.dim = int(dim)
		self.bound = Boundaries
class Segments:
	def __init__(self, data):
		length = len(data)-2
		Residues = []
		cnt = 1
		while cnt <=  length:
			Residues.append(str(data[cnt]))
			cnt = cnt + 1
		self.name = str(data[0])
		self.bin = str(data[-1])	
		self.length = int(length)
		self.res = Residues
	
def addreg (rawdata, Class):
	try:
		In = Class(rawdata.split(","))
		print(str(In.name)+" added to class: "+str(Class))
	except ValueError:
		print("inputs could not be read in")
	return In

REGIONS = []
SECS = []
regdim = ''
segdim = ''

lib0 = "**********************************************************************\n"
lib1 = "Format specifications for the DISICL_dbssp calssification libraries:\n\n"
lib2 = "Regions are defined as a rectangular/cubic parts of the dihedral angle space\n"
lib3 = "Define Ramachandran regions using:\n   @reg <name,xmin,xmax,ymin,ymax,zmin,zmax,...>\n"
lib4 = "\t   name is a string used for strucutre definitions\n\n\t   xmin,xmax,... are pairs of floats, one for each dihedral angle\n\n"
lib5 = "Define segment for a secondary structure element using:\n   @sec <name,reg1,reg2,....,bin>\n"
lib6 = "\t   name is the full name of secondary structure element\n\n\t   reg1, reg2, ... are the names of sebsequent regions\n"
lib7 = "\t   the residues of the segment should fall in\n\n\t   bin is the collector code for the secondary structure element\n\t   (if you have multiple areas).\n\n"
lib8 = "\t   It is recommended that name and bin is always changed together!\n"

libinfo = lib0+lib1+lib2+lib3+lib4+lib5+lib6+lib7+lib8+lib0
#default protein regions:
if defflag == "prot":
	print("Error: No definition library provided!")
	if debug == 1:
		print("Debug mode is called, program runs with hardcoded examples")
	else:
		print("Reading in hardcoded library entries, then program stops!")
		print(usage)
		failmark = 2
	
#Early library example! for test purposes only!!!!! 
#region definitions
#xh2 = [ 3-10 helix -113,-48,-22,22]
	regdim = 2
	xh1 = "delta,-113,-48,-22,22"
	REGIONS.append(addreg(xh1,Regions))
#xh2 = [ alfa helix -58,-43]
	xh2 = "alfa,-83,-43,-62,-22"
	REGIONS.append(addreg(xh2,Regions))
#ebs = [ extended beta -147,158]
	ebs = "beta2,-167,-127,136,179"
	REGIONS.append(addreg(ebs,Regions))
	
#secondary structure dfeinitions 
#H1 = ["#3-10 helix residues\n#       time\tresidue\n"]
	segdim = 2
	hel31 = "3-10 helix,alfa,delta,H1"
	SECS.append(addreg(hel31,Segments))
	hel32 = "3-10 helix,delta,delta,H1"
	SECS.append(addreg(hel32,Segments))
	hel33 = "3-10 helix,delta,alfa,H1"
	SECS.append(addreg(hel33,Segments))
	ahelix = "alfa helix,alfa,alfa,H2"
	SECS.append(addreg(ahelix,Segments))
	ebs = "extended beta sheet,beta2,beta2,EBS"
	SECS.append(addreg(ebs,Segments))


	
elif defflag == "help":
	failmark = 2
	
else:
#custom definition library
	l = open(defflag,"rb")
	print("\nReading in library definitions:")
	for line0 in l:
		line = str(line0.decode("ascii"))
		if not line.startswith("#") or not line.strip() == "":
			line2 = line.strip().strip("\n")
#adding region definitons
			if line2.startswith("@reg"):
				line3 = line2.strip("@reg").strip()
				parts = line3.split(",")
				partnum = len(parts)
				print("\nregion definition with "+str(partnum)+" arguments:")
				try:
					name = str(parts[0])
					cnt = 1
					while cnt < partnum:
						boundary = float(parts[cnt])
						cnt = cnt + 1
					if partnum % 2 == 1:
						dimension = (partnum-1)/2
						print(name,"region with",str(dimension),"dimensions detected")
						if regdim == '':
							regdim = dimension
						elif regdim != dimension:
							print("\nError: based on previous definitions "+str(regdim)+" dimensions are expected, found:",str(dimension))
							print(line)
							failmark = 2
						cnt = 0
						while cnt < dimension:
							boundmin = float(parts[2*cnt+1])
							boundmax = float(parts[2*cnt+2])
							if boundmax <= boundmin:
								print("\nError: maximal boundary for dimension "+str(cnt+1)+" is smaller than minimal value")
								print(line)
								failmark = 2
							cnt = cnt + 1
					else:
						print("\nERROR: wrong number of arguments for region object")
						print("name + min/max boundaries per dimension is required divided by coma")
						print(line)
						failmark = 2
				except ValueError:
					print("library error at line:")
					print(line)
					failmark = 2
				REGIONS.append(addreg(line3,Regions))
#adding segment definitions:
			if line2.startswith("@sec"):
				line3 = line2.strip("@sec").strip()
				parts = line3.split(",")
				partnum = len(parts)
				try:
					name = str(parts[0])
					code =	str(parts[-1])
					length = partnum-2
					if length >= 1:
						if segdim == '':
							segdim = length
						elif segdim != length:
							print("\nError: based on previous definitions "+str(segdim)+" length is expected, found:",str(length))
							print(line)
							failmark = 2
						string = "("
						cnt = 1
						while cnt <= length:
							if not cnt == 1:
								string = string+","
							string = string+str(parts[cnt])
							cnt = cnt + 1
						string = string+")"
						print("\nsegment definition with",str(length)," residue length detected")
						print(name,"(",code,") class received segment:",string)
					else:
						print(line)
						print("\nError: segment definition does not have the minimal amount of arguments")
						print("Segment definiton requires a class name, at least one region code, and a class code")
						failmark = 2
				except Exception:
					print("\nError: segment definition could not be read in")
					print(line)
					failmark = 2
				SECS.append(addreg(line3,Segments))
			
if failmark == 0:		
	empty= 'unclassified,"X","X",UC'
	SECS.append(addreg(empty,Segments))
if failmark == 1:
	print(usage)
	sys.exit(0)
if failmark == 2:
	print(libinfo)
	sys.exit(0)
#stating operations
if debug == 1:
	print("\nprinting regions:")
	for Region in REGIONS:
		print(Region.name,Region.dim,Region.bound,str(math.fabs(Region.bound[0]-Region.bound[1])))
	print("\nprinting segments:")
	for Segment in SECS:
		print(Segment.name,Segment.bin,Segment.length,Segment.res)
	
#print "lets stop here!"
#sys.exit()
#getting residue list:


#creating secondary structure lists:
fin = []
FIN =[]
ssvector = []
for struct in SECS:
	if struct.bin not in fin:
		fin.append(struct.bin)
		header = "#"+struct.name+" ("+struct.bin+") residues\n#       time\tresidue\n"
		list = [header]
		FIN.append(list)
		ssvector.append(0)
ssnum = len(fin)
print("\nSecondary structure elements:")
print(str(fin),str(ssnum))
#print str(ssvector)
#print str(FIN)

#preparing classification lists and time series:
outfile = outbase[0]+"_"+str(RES[0])+"-"+str(RES[-1])+outbase[1]
o = open(outfile,"wb")
o.write("# DISICL calssification of the biopolymer conformational states\n".encode("ascii"))
if argnum == 0:
	infile = filebase[0]+str(RES[0])+"-"+str(RES[-1])+filebase[1]
i = open(infile,"rb")
cnt0 = 0
Segment = []
cnt = 0
while cnt < segdim:
	Segment.append([0.0,0,""])
	cnt = cnt + 1
n = segdim -1

if debug == 1:
	print("\nDebug mode is on, detailed messages are printed")

print("\nStarting segment classification based on dihedrals:")
fcnt = 0
FRAME = []
time1 = ''
time2 = ''
for line0 in i:
#getting residue information from input file
	line = str(line0.decode("ascii"))
	if not (line.startswith("#") or line.startswith("&") or line.startswith("@")):
		line2 = line.strip().strip("\n")
		ok = 0
		wmark = 0
		time0 = Segment[0][0]
		resold = Segment[0][1]
		try:
			
			parts = line2.split()
			partnum = len(parts)
			time = float(parts[0])
			res = int(parts[1])
			Angles = []
			if partnum > 2:
				dcnt = 2
				while dcnt < partnum:
					angle = float(parts[dcnt])
					Angles.append(angle)
					dcnt = dcnt + 1
			
			datain = (time0,resold)
			dstring = "   %10.1f\t%4d\n"
			data = dstring % datain
			ok = 1
		except Exception:
			print(line+"\nline could not be read")
			wmark = 2
#assign regions to residues in the Segment
		if ok == 1:
			regnew = ""
			for reg in REGIONS:
				reglength = reg.dim
				dimension = reg.dim/2
				if (partnum-2)  == dimension:
					dcnt = 0
					fit = 0
					while dcnt < dimension:
						dmin = reg.bound[2*dcnt+0]
						dmax = reg.bound[2*dcnt+1]
						angle = Angles[dcnt]
						if (dmin < angle <= dmax):
							fit = fit + 1
						dcnt = dcnt + 1
					if fit == dimension:
						regnew = reg.name

				elif (partnum-2) > dimension:
					if wmark == 0:
						print("Warning, library for "+str(dimension)+" variables, and "+str(partnum-2)+" are present, first ones are taken!")
					wmark = 1
#					if ( reg.xl < fi <= reg.xu ) and ( reg.yl < psi <= reg.yu):
#						res1 = reg.name
				else:
					if wmark == 0:
						print("Warning, library requires "+str(dimension)+" varibles, "+str(partnum-2)+" found,skipping line")
					wmark = 2
			if time0 != time and dimension > 1:
				Segment[0][2] = "X"
			if regnew == "":
				regnew = "X"
			Segment[n][0] = time
			Segment[n][1] = res
			Segment[n][2] = regnew
			if (resold in RES or res in RES) and (time == TIME[0] or debug == 1):
				headline = ''
				showline = '  '
				cnt = 0
				while cnt < segdim:
					headline = headline+"time"+str(cnt)+","+1*" "+"res"+str(cnt)+","+1*" "+"region"+str(cnt)+" | "
					showline = showline+str(Segment[cnt][0])+3*" "+str(Segment[cnt][1])+3*" "+str(Segment[cnt][2])+1*"  |  "
					cnt = cnt + 1
				headline  = headline+6*" "+"[ANGLES for res"+str(n)+"]"
				showline = showline+3*" "+str(Angles)
				print(headline)
				print(showline)
				
#			if ( reg.xl < fi <= reg.xu ) and ( reg.yl < psi <= reg.yu) and ( reg.zl < xi <= reg.zu):
#						res1 = reg.name
#		except Exception:
#			print line
#			print "could not be read in"
# calassify segment based on residue assignment
		ID = ""
		cont = 1
		cnt = 0
		while cnt < n:
			Res0 = Segment[cnt][1]
			Res1 = Segment[cnt+1][1]
			if Res1 != (Res0+1):
				cont = 0
			cnt = cnt + 1
		
		if cont == 1  and wmark != 2 and resold in RES:
			for struct in SECS:
				fit = 0
				dcnt = 0
				if struct.length == segdim:
					while dcnt < segdim:
						regcode = Segment[dcnt][2]
#						print dcnt,struct.res
						refcode = struct.res[dcnt]
						if regcode == refcode:
							fit = fit + 1
						dcnt = dcnt + 1
					if fit == struct.length:
						ID = struct.bin
			if ID == "":
				ID = "UC"
			if time == TIME[0] or debug == 1:
				print("Segment ",str(resold)," assigned to class:",ID,"\n")
			cnt = 0
			while cnt < ssnum:
				if ID == fin[cnt]:
					ssvector[cnt] = ssvector[cnt] + 1
					FIN[cnt].append(data)
					if time1 == '':
						time1 = datain[0]
					if time2 == '' and time1 != '' and time1 != datain[0]:
						time2 = datain[0]
					frame = int((datain[0]-TIME[0])/TIME[1]+1)
					resid = datain[1]
					info = (frame,resid)
					FRAME.append(info)
				cnt = cnt + 1
			cnt0 = cnt0 + 1	
#storing residue info for next segment
		if segdim != 1 and wmark != 2:
			cnt = n
			while cnt > 0:
				Segment[cnt-1][0] = Segment[cnt][0]
				Segment[cnt-1][1] = Segment[cnt][1]
				Segment[cnt-1][2] = Segment[cnt][2]
				cnt = cnt -1
#		res0 = res1 
#		resold = res
#		time0 = time
		
i.close()

#checking basic results: 
print("\n")
fmin = ''
fmax = ''
wmark = 0
tst = 0
if time2 != '' and time1 != '': 
	tst = time2-time1
if time1 != TIME[0]:
	print("Warning! first timecode in input file ("+str(time1)+") does not match predefined starting time ("+str(TIME[0])+")\n")
	wmark = 1
if tst != 0 and tst != TIME[1]:
	print("Warning! timestep in input file ("+str(tst)+") does not match predefined starting time ("+str(TIME[1])+")\n")
	wmark = 1 
 
for entry in FRAME:
	if fmin == '':
		fmin = entry[0]
	elif fmin > entry[0]:
		fmin = entry[0]
	if fmax == '':
		fmax = entry[0]
	elif fmax < entry[0]:
		fmax = entry[0]

if fmin == '':
	fmin = 0
if fmax == '':
	fmax = 0
fcnt = fmin
framenum = 0
resnum0 = 0

while fcnt <= fmax:
	resnum = 0
	for entry in FRAME:
		if entry[0] == fcnt:
			resnum = resnum + 1
	print("Frame "+str(fcnt)+" was found with "+str(resnum)+" classified residues")
	if resnum != 0:
		framenum = framenum + 1
	if fcnt == fmin:
		resnum0 = resnum
	fcnt = fcnt + 1
print("Total number of frames found: "+str(framenum)+"\n"+str(resnum0)+" applicable segments in the first frame")

if wmark == 1:
	print("\nDue to Warnings about the timestamps, summary should be treated with caution!!!\n")	
summary1 = "DISICL"+2*" "+"time0"+2*" "+"timestep"+2*" "+"frames"+2*" "+"segment0"+3*" "+"reslist"
summary2 = "DISICL"+3*" "+str(TIME[0])+5*" "+str(TIME[1])+7*" "+str(framenum)+7*" "+str(resnum0)+6*" "+"["+Newreslist+"]"
print("\n"+summary1+"\n"+summary2+"\n")

#writing out dssp output information
o.write(("# based on dihedral angle file: "+infile+"\n# and library file: "+defflag+"\n").encode("ascii"))
if wmark != 1:
	o.write(("#\n#"+summary1+"\n#"+summary2+"\n#\n").encode("ascii"))
for sec in FIN:
	print(" writing DISICL input for:\n"+str(sec[0]))
	for data in sec:
		o.write(data.encode("ascii"))
	o.write("&\n".encode("ascii"))


o.write("#Secondary structure distribution of requested segments (%):\n".encode("ascii"))
cnt = 0
while cnt < ssnum:
	if cnt0 != 0:
		value = float(ssvector[cnt]*100)/float(cnt0)
	else:
		value = 0
	o.write(("# "+str(fin[cnt])+":  %2.1f \n" % (value)).encode("ascii"))
	cnt = cnt + 1
o.close()
print("Segment analysis complete")
print("writing classification time series to "+outfile)

#making statistics

print("calculating residue statistics")
i2 = open(outfile,"rb")
outfile2 = outbase2[0]+"_"+str(RES[0])+"-"+str(RES[-1])+outbase2[1]
o2 = open(outfile2,"wb")

MATRIX = []
for sec in FIN:
	line = []
	for res in RES:
		line.append(0)
	MATRIX.append(line)
cntl = 0 
for line0 in i2:
	line = str(line0.decode("ascii"))
	if not line.startswith("#"):
		if line.startswith("&"):
			cntl = cntl + 1
		if not line.startswith("&"):
			line2 = line.strip().strip("\n")
			parts = line2.split()
			res = int(parts[1])
#			print line2,res
			rnum = len(RES)
			cntr = 0
			while cntr < rnum:
				if RES[cntr] == res:
#					print (cntl,cntr)
					MATRIX[cntl][cntr] = (MATRIX[cntl][cntr]) + 1
				cntr = cntr + 1
i2.close()

print("writing statistics to "+outfile2)
linenum = cntl-1
o2.write("#Distribution of DISICL secondary structure elements\n".encode("ascii"))
o2.write(("#based on dihedral file: "+infile+"\n#and library file: "+defflag+"\n").encode("ascii"))
head = "#residue"+7*" "
average = []
for sec in fin:
	sp = len(str(sec))
	head = head+(str(sec)+(16-sp)*" ")

head = head + "\n"
o2.write(head.encode("ascii"))	
rcnt = 0
lcnt = 0

rnum = len(RES)

while rcnt < rnum:
	info = [RES[rcnt]]
	tot = 0
	for mat in MATRIX:
		tot = tot + mat[rcnt]
	while lcnt <= linenum:
		number = MATRIX[lcnt][rcnt]
		if not tot == 0:
			perc = float(number)/float(tot)*100
		else:
			perc = 0.0
		info.append(number)
		info.append(perc)
		lcnt = lcnt + 1
#	ostring = " %4d"+(linenum+1)*("\t %8d  %-2.1f ")+"\n"
	ostring = " %4d"+(linenum+1)*(" %8d  %4.1f ")+"\n"
	

#	print info
#	print ostring
	o2.write((ostring%tuple(info)).encode("ascii"))
	rcnt = rcnt + 1
	lcnt = 0

LAST=[] 
for l in ssvector:
	if cnt0 != 0:
		perc = float(l)/float(cnt0)*100
	else:
		perc = 0
	LAST.append(l)
	LAST.append(perc)
info = tuple(LAST)
#astring = "\n#av:"+(linenum+1)*("\t %8d  %-2.1f ")+"\n"
astring = "\n#av: "+(linenum+1)*(" %8d  %4.1f ")+"\n"

o2.write((astring%info).encode("ascii"))
i2.close()
o2.close()
ftime = os.times()
runtime = ftime[4]-stime[4]
if runtime == 0.0:
	runtime = ftime[0]+ftime[1]-(stime[0]+stime[1])
print("program runtime: %1.2f  seconds (%1.2f h)" % (runtime,(float(runtime)/3600)))
print("Classification program finished successfully! ("+infile+")")
