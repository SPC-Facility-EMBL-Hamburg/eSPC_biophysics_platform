#!/usr/bin/env python
#

import sys
import os
import math
import time

#script to quickly install SESCA on a given machine:
workdir = os.getcwd()
s_time = time.time()

usage0 = "**********************************************\n"
usage  = "SESCA setup script. Usage:\n"
usage += "\nsetup.py @flag <argument>\n\nPossible command flags are:\n"
usage += "   @prefix <path>  install SESCA to specified directory (default '.')\n"
usage += "   @source <path>  use SESCA source files from specified directory (default '.')\n"
usage += "   @win <0,1> install for windows, 1 - True , 0 - False (default 0)\n"
usage += "   @help (or -h or --help) print this message\n"
Usage = usage0 + usage + usage0


#############################################################
#default variables:
prefix = ""
win = False
source = ""
DIRS = [["DSSP","libs", "examples"], ["DISICL","HBSS","scripts"]]
DSSP_ver = ["dssp-2.0.4-linux-amd64","dssp-2.0.4-win32.exe"]
#############################################################


#############################################################
#handle supplied arguments:
argnum = len(sys.argv)-1
flag = ""
for r in range(argnum):
	arg = sys.argv[r+1]
	if arg.startswith("@"):
		flag = arg.strip("@")
	elif flag == "prefix":
		prefix = arg
		flag = ""
	elif flag == "win":
		if arg in ["1","True"]:
			win = True
		elif arg in ["0","False"]:
			win = False
		else:
			print("invalid argument")

		flag = ""
	elif flag == "source":
		source = arg
		flag = ""
	if arg in ["-h", "--help", "@help"]:
		print(Usage)
		sys.exit()


#set up basic command line commands:
if win == False:
	copy = "cp"
	copy_dir = "cp -r"
	delete = "rm"
	move = "mv"
	show = "ls"
else:
	copy = "copy"
	copy_dir = "Xcopy /E /I"
	delete = "del"
	move = "move"
	show = "dir"

########################################################	
#defining a helper function here: 

#function to change lines according to keywords:
def Modify_File(File, Mods):
	Output = ""
	mod_num = len(Mods)

#	make sure we have the file:
	if os.path.isfile(File) == False:
		print("Unable to modifiy file: no file found!")
		return Output

	f = open(File,"rb")
	cnt = 0
	for line0 in f:
		line = str(line0.decode("ascii"))
		cnt += 1
		modified = 0
		for Mod in Mods:
			keyword = Mod[0]
			if line.startswith(keyword) and Mod[2] == 0:
				Output +=  Mod[1]
				modified += 1
				Mod[2] += 1
				print("line %3d modified to:\n %1s"%(cnt,Mod[1]))
		if modified == 0:
			Output += line
			
	f.close()
				
	return Output

#function to control the necessary modifications:
def Set_Mods(Prefix, Subdir, File, win, Dssp_files):
	Modify = []
	if Subdir == "scripts":
#		modifications on  core SESCA modules:								
		dssp_version = ""
		Modify.append(["SESCA_Dir",'SESCA_Dir =     "%1s"\n'%Prefix, 0])
		if win == False:
			Modify.append(["win = ", "win = 0\n", 0])
			dssp_version = Dssp_files[0]
		if win == True:
			Modify.append(["win = ", "win = 1\n", 0])
			dssp_version = Dssp_files[1]
		if File == "SESCA_dssp.py":
			DSSP_dir = os.path.join(Prefix,"DSSP")
			dssp_main = os.path.join(DSSP_dir,dssp_version)					
			Modify.append(["dssp_main", 'dssp_main = "%1s"\n'%dssp_main, 0])

	if Subdir == "HBSS":
		HBSS_dir = os.path.join(Prefix, "HBSS")
		Modify.append(["HBSS_dir", 'HBSS_dir =  "%1s"\n'%HBSS_dir, 0])

	if Subdir == "DISICL":
		DISICL_dir = os.path.join(Prefix, "DISICL")
		Modify.append(["DISICL_dir", 'DISICL_dir =  "%1s"\n'%DISICL_dir, 0])
		if win == True:
			Modify.append(["win = ", "win = 1\n", 0])
		else:
			Modify.append(["win = ", "win = 0\n", 0])				
#	print(Modify)
	return Modify

#helper function to hande backslash error:
def Fix_prefix(prefix):
	prefix_fixed = prefix
#	Specials = ["\\n","\\b","\\t","\\r","\\f"]
	Specials = ["\\"]

#	make sure the prefix is absolute:
	if prefix_fixed.startswith("."):
		curr_dir = os.getcwd()
		prefix_joined = os.path.join(curr_dir, prefix_fixed)
		prefix_fixed = os.path.normpath(prefix_joined)

#	cycle through all special characters:
	for Special in Specials:
#		split the text if they are found:
		parts = prefix_fixed.split(Special)
		p_num = len(parts)
		print(Special,p_num)
		prefix_fixed = ""
#		insert double slash:
		for p in range(p_num-1):
			prefix_fixed += parts[p]+1*"\\"+Special
			
		prefix_fixed += parts[(p_num-1)]


	print("\nFinal prefix: %1s"%prefix_fixed)
	
	return prefix_fixed

############# Run setup script from here:

#determine source directory if was not determined:
if source == "":
	source = os.getcwd()
source = os.path.normpath(source)
print("\nUsing source directory:\n%1s"%source)

#determine target directory:
if prefix == "" and source != "":
	prefix = source
	
#fix backslash errors in the prefix:
prefix_path = Fix_prefix(prefix)

#create target directory it does not exist yet
if prefix != "" and os.path.isdir(prefix) == False:
	if prefix.startswith("."):
		prefix = os.path.join(workdir, prefix)
	prefix = os.path.normpath(prefix)
	dir_command = "mkdir %1s"% prefix
	print("\nCreating target directory:\n%1s"%prefix) 
	os.system(dir_command)
elif prefix != "" and prefix != source:
	if prefix.startswith("."):
		prefix = os.path.join(workdir, prefix)
	prefix = os.path.normpath(prefix)
	print("\nTarget directory found:\n%1s"%prefix) 
#go to source directory:
os.chdir(source)

#modify local files if no prefix was given:
if prefix == source:
	Modify = []
#	cycle through the marked sub-directories
	for Subdir in DIRS[1]:
		os.chdir(Subdir)
#		find python modules:
		Subdir_Files = os.listdir(".")
		print("\nUpdating variables in: %1s"%Subdir)
		for File in Subdir_Files:
			extension = File.split(".")[-1]
#			if the file is a python file, check for modifications:
			if os.path.isfile(File)== True and extension in ["py"]:
				print("modifying file: %1s"%File)
				Modify = Set_Mods(prefix_path, Subdir, File, win, DSSP_ver)
				Content = Modify_File(File, Modify)
#				overwrite old file:
				f = open(File,"wb")
				f.write(Content.encode("ascii"))
				f.close()
#				give run permissions on linux:
				if win == False:
					os.system("chmod %1d %1s "%(744, File))	

		os.chdir("..")		
	

#copy files if a prefix was given:
if prefix != source:
	Files = os.listdir(".")
	for File in Files:
#		copy basic files:
		if os.path.isfile(File) == True:
			extension = File.split(".")[-1]
			if extension in ["py","txt"]:
				copy_command = "%1s %1s %1s" %(copy,File,prefix)
				os.system(copy_command)
				print(copy_command)
#		copy basic directories:
		elif os.path.isdir(File) == True and not File.startswith("__"):
			if File in DIRS[0]:
#				copy nonproblematic sub-directories:
				new_dir = os.path.join(prefix,File)
				copy_command = "%1s %1s %1s" % (copy_dir, File, new_dir)
				os.system(copy_command)
				print(copy_command)
			elif File in DIRS[1]:
#				make sub-directories where module files are stored:
				new_dir = os.path.join(prefix,File)
				dir_command = "mkdir %1s"%new_dir
				print(dir_command)
				os.system(dir_command)
				Modify =[]

#				enter subdir and copy files:				
				os.chdir(File)
				Subdir_Files = os.listdir(".")
				Subdir_path = os.path.join(prefix, File)
				for Sub_file in Subdir_Files:
#					copy sub-directories:
					if os.path.isdir(Sub_file) == True and not Sub_file.startswith("__"):
						new_dir2 = os.path.join(Subdir_path, Sub_file)
						copy_command = "%1s %1s %1s" % (copy_dir, Sub_file, new_dir2)
						os.system(copy_command)
						print(copy_command)
#					copy non-module files:
					elif os.path.isfile(Sub_file) == True:
						extension = Sub_file.split(".")[-1]
						if not extension in ["py","pyc"]:
							copy_command = "%1s %1s %1s" %(copy,Sub_file,Subdir_path)
							os.system(copy_command)
							print(copy_command)
#						modify python modules:
						elif extension in ["py"]:
							print("modifying file: %1s"%Sub_file)
							Modify = Set_Mods(prefix_path, File, Sub_file, win, DSSP_ver)
							Content = Modify_File(Sub_file, Modify)
							new_File = os.path.join(Subdir_path,Sub_file)
							f = open(new_File,"wb")
							f.write(Content.encode("ascii"))
							f.close()
#							give run permissions on linux:
							if win == False:
								os.system("chmod %1d %1s "%(744, new_File))	
				os.chdir("..")
f_time = time.time()
runtime = f_time - s_time
print("\nSetup script finshed!")
print("Script runtime was %2.2f seconds" % runtime)
