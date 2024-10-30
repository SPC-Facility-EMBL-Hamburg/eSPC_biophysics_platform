#Init module for the SESCA package
#written by:     Gabor Nagy in 2019
#affiliation:    Max Planck Institute for Biophysical Chesmistry
#		 Department of Theoretical and Computational Biophysics
#		 Am Fassberg 11, 37077 - Gottingen, Germany
#email:		 gabor.nagy@mpibpc.mpg.de
#
#	This is a python package for validating 3D structural models of proteins (in PDB format)
#       based on measured Circular Dichroism spectra


__all__ = [
		'SESCA_bayes.py'
		'SESCA_deconv.py',
		'SESCA_dssp.py',
		'SESCA_main.py',
		'SESCA_min.py',
		'SESCA_pred.py',
		'SESCA_process',
		'SESCA_scale.py',
		'SESCA_seq.py',
		'SESCA_solver.py',
		
	  ]    

# a test function to set print Verbose messages:
verb_init = 2
def Vprint(level,*Messages):
	global verb_init
	if level <= verb_init:
		string = ""
		for message in Messages:
			string += str(message)+" "
			print(string)
	else:
		pass
 
def Set_verb_init(int):
	global verb_init
	verb_init = int

def Print_verb_init():
	print(verb_init)

#set global verbosity:
Set_verb_init(verb_init)
Vprint(2, "Loading SESCA_modules:")
import SESCA_main
import SESCA_pred
import SESCA_scale
import SESCA_bayes
import SESCA_process
