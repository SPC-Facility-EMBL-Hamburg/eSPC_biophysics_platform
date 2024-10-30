#!/usr/bin/env python
#

#Step 1: import basic modules
import sys
import os

#check if the os module works:
Workdir = os.getcwd()
print("Current work directory: %10s"%Workdir)



#Step 2: import SESCA modules
#define SESCA_dir based on current directory:
SESCA_dir = os.path.join(Workdir,"../..")
SESCA_dir = os.path.normpath(SESCA_dir)
Script_dir = os.path.join(SESCA_dir,"scripts")
#make sure that the SESCA scripts directory is in the system path:
print("\nImporting SESCA from directory: %10s"%Script_dir)
sys.path.append(Script_dir)

#import SESCA_main.py and SESCA_pred.py from the scripts directory:
import SESCA_main as Main
import SESCA_pred as Pred
import SESCA_bayes as Bayes

#get all imports for SESCA main as well:
Imports = Main.Import_Custom()

SESCA_args     = " @spect ../5-Bayesian_SS/CD_igg2.dat @lib DS-dTSC3 @corr ../5-Bayesian_SS/SC_calc.out @write Bayes_est_1.out @proj Bayes_map_1.out    @iter 250  @verb 3"
Processed_Args = Bayes.Read_Args(SESCA_args.split())
print(Processed_Args)
Data           = Bayes.SSbayes_Main(Processed_Args)

#print(len(Data))
sys.exit()

#Step 3: predict  CD spectra from three structures:
Files = ["1scd.pdb","9pap.pdb", "1hnn.pdb"]
cnt = 1
CD_Spectra = []
SS_Comp = []
Data = []
for File in Files:
	print("\nNow processing: %1s"%File)
#	define SESCA arguments:
	SESCA_args = " @pdb %1s @write Spectrum_%1d.out @lib HBSS-3 @verb 1" % (File,cnt)
	print("Arguments: %1s"%SESCA_args)
	Processed_Args = Main.Read_Args(SESCA_args.split())
#	call SESCA_Main function to execute the program:
	Data = Main.SESCA_Main(Processed_Args)
#	Data array 0: SS composition, CD spectra
#	Data array 1: comparison results (now empty)
#	save the SS composotion:
	SS_Comp.append(Data[0][0])
#	save the calculated CD spectrum:
	CD_Spectra.append(Data[0][1])

#	print secondary strucutre summary:
	print("\nSecondary structure:")
	for Coeff in Data[0][0]:
		print("%10s: %1.3f"%tuple(Coeff))

	print("\nPredicted CD spectrum:")
	for Wave in Data[0][1]:
		print("%2.1f   %1.3f"%tuple(Wave))

	cnt += 1



#Step 4: Compare how similar the CD spectra are to one another:
RMSD_Mat = []
f_num = len(CD_Spectra)
for i in range(f_num):
#	get first spectrum i:
	Spectrum_i = CD_Spectra[i]
	print(Spectrum_i)
	Row = [(i+1)]
	for k in range(f_num):
#		get second spectrum k:
		Spectrum_k = CD_Spectra[k]
#		use the Spectrum comparison function from SESCA_pred
		Comparison_ik = Pred.Compare_Spectra(Spectrum_i, Spectrum_k)
		print(Comparison_ik[1])
#		Add the RMSD to the comparison:
		match,RMSD,MUSE,NRMSD,Miss = Comparison_ik[1] 	
		Row.append(RMSD)
	RMSD_Mat.append(Row)



#Step 5: Display spectrum devations:
#generate output header:
Output = "#Spectrum Deviations (kMRE):\n"+5*" "
for i in range(f_num):
	Output += "   %1d  " % (i+1)
Output += "\n"
#add RMSD matrix to output:
for Row in RMSD_Mat:
	Output += "   %1d   %2.2f  %2.2f  %2.2f\n"%tuple(Row)  
#print and write output: 
print("\n"+Output)
o = open("Spectrum_comparison.out","wb")
o.write(Output.encode("ascii"))
o.close()
