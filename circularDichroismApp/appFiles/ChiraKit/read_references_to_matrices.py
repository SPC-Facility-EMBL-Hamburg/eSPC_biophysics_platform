import numpy  as np 
import pandas as pd

'''
Helper functions to read input files containing the F and C matrices 

C (also called 'A' in the literature) is a m × n matrix containing the CD spectra in 
delta epsilon units (mean residue molar extinction).
m and n are respectively the number of wavelengths and proteins.
The lower wavelength must be lower or equal than 190.

F is a l × n matrix containing the known secondary 
structure elements. l and n are respectively the number of 
secondary structure elements and proteins. 
To prevent misuse, if the C matrix wavelength lower limit is between 
180 and 190, only upto four secondary structure elements
are accepted. Similarly, if the C matrix wavelength lower limit
is between 170 and 180, only upto six secondary structure elements are accepted.

'''

def check_line_is_valid_for_matrix(splitted_line):

	try:

		abs(float(splitted_line[0])) < 1e2 
		abs(float(splitted_line[1])) < 1e2 
		abs(float(splitted_line[2])) < 1e2 

		return len(splitted_line) >  5

	except:

		return False

# heuristic algorithm to guess the delimiter and decimal characters of the F and C data files
def find_delimiter_and_decimalChar(file):

	with open(file,'r') as f:

		firstLine = f.read().splitlines()[0]

		delimiter_list = [',','\t',' ',';']

		sel_delimiter  = None

		for delimiter in delimiter_list:

			fl_split = firstLine.split(delimiter)

			if check_line_is_valid_for_matrix(fl_split):

				sel_delimiter = delimiter
				decimal_char  = '.'
				break

		# Find if we have ',' instead of '.' as decimal
		if sel_delimiter is None:

			for delimiter in delimiter_list[1:]:

				firstLine = firstLine.replace(',', '.')
				fl_split  = firstLine.split(delimiter)

				if check_line_is_valid_for_matrix(fl_split):

					sel_delimiter = delimiter
					decimal_char  = ','
					break

	return sel_delimiter, decimal_char

def read_matrix_file(file):
	
	sel_delimiter, decimal_char = find_delimiter_and_decimalChar(file)

	df = pd.read_csv(file,sep = sel_delimiter, decimal = decimal_char,header=None)

	return df.to_numpy()
