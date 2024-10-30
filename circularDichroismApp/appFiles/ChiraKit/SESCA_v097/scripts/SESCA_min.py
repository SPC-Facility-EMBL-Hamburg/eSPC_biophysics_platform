#!/usr/bin/env python
#

import sys
import os
import math
import time

workdir = os.getcwd()
stime = time.time()


#print workdir+"\n"

usage0  = "*******************************************************************\n"
usage   = "SESCA minimization module, based on the Nelder-Mead Simplex method \n"
Usage = usage0 + usage + usage0

#default parameters:
in_file = ""
out_file = ""
test = 1
#nelder-mead paramteters for reflection, expansion, contraction, and shrink:
adaptive = 0
NM_param = [1.0, 2.0, 0.5, 0.5]
#minimization parameters:
tolerance = 0.00000001
Max_iter = 100
display = 1
Min_param = [adaptive,tolerance, Max_iter,display]
failmark = 0
verbosity = 2

Input_files = [in_file]
Output_files = [out_file]
Param = [NM_param, Min_param, test]

Def_Args = [Input_files, Output_files, Param, failmark, verbosity]


#function definitions:
#function to pass on defaults:
def Pass_Defaults():
        return Def_Args

def Pass_NMparam():
	return NM_param

def Pass_ANMparam(dim):
	ANM_param = [1.0, (1.0+2.0/dim), (0.75-(1.0/(2*dim))), (1.0-(1.0/dim))]
	return ANM_param

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
	New_Args = [Input_files, Output_files, Param, failmark, verbosity]
	FLAGS = ["inp","write","verb"]
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
				New_Args[3] = 1
		elif flag == "inp":
			New_Args[0][0] = arg
			flag = ""	
		elif flag == "write":
			if arg == "0":
                                New_Args[1][0] = ""
			else:
				New_Args[1][0] = arg
			flag = ""
		elif flag == "verb":
			New_Args[4] = int(arg)
			flag = ""
#		setting default files if no flags are provided:
		elif flag == "" and New_Args[0] == "" and acnt == 0:
			New_Args[0][0] = arg
			flag = ""
		elif flag == "" and New_Args[1] == "" and acnt == 1:
			New_Args[1][0] = arg
			flag = ""
		else:
			Vprint(1,"unknown argument:",arg)
			New_Args[3] = 1

		acnt += 1
	return New_Args

def value_key(Point):
	return Point[0]	

#####################################################################
#Class definition for a Simplex:
class Simplex:
#	Define a simplex for an N dimensional function, it should contain N+1 points
#   All point should have a function value and and array of N coordinates
	def __init__(self,Points):
		dimension = len(Points)-1
		self.dim = int(dimension)
		self.points = []
#		check that the Simplex has the right amount of points
		for i in range(dimension+1):
			Point =  Points[i]
			P_dim = len(Point[1])
			if P_dim != dimension:
				Vprint(2, "Error: Point has unexpected dimensions (%1d, expected %1d)"%(Pdim,dimension))
				raise TypeError()
			else:
				Vprint(4, "New Simplex point:",Point)
				self.points.append(Point)
				
#	set points in the simplex:
	def Set_Pi(self,num, new_Pi):
		if num <= self.dim:
			if self.dim == len(new_Pi[1]):
				self.points[num] = new_Pi
			else:
				Vprint(1, "Simplex points need %1d coordinates" % self.dim)
		else:
			Vprint(1, "Simplex has %1d points!" % self.dim)

#	get a given point in the simplex:
	def Get_Pi(self,num):
		if num <= self.dim:
			return self.points[num]
		else:
			Vprint(1, "Simplex has %1d points!" % self.dim)
			
#	list the value of all points in the simplex:			
	def values(self):
		Simp_values = []
		for Point in self.points:
			Simp_values.append(Point[0])
		return Simp_values

#	a function to sort points by their value:
	def sort_points(self):
		self.points.sort(key= value_key, reverse=True)
		Vprint(5, "Sorted simplex:\n", self.points) 
			
#	function to print simplex
	def print_pi(self):
		Info = "Simplex dimension: %1d\n" % self.dim
		Info += "Current Points:\n"
		for Pi in self.points:
			Info += " with F(xi) = %1.3f, Xi={" % Pi[0]
			for coord in Pi[1]:
				Info += "%1.3f,"%coord
			Info = Info[:-1]+"}\n"
		return Info
			
#	function to calculate the center for the best n simplex vertices:
	def Xncent(self):
		Center_coords = []
		for j in range(self.dim):
			Sumj = 	0.0
			for i in range(1,self.dim+1):
					Pi = self.points[i][1]
					Sumj += Pi[j]/self.dim
			Center_coords.append(Sumj)
		return Center_coords
#####################################################################

#function reflect the worst point of simplex:
def Reflect_worst(Simplex, NM_param, Func, Fparam, test):
	Xref = [0.0, []]
	dim = Simplex.dim
	alpha = NM_param[0]
	Xworst = Simplex.Get_Pi(0)
	Xcent = Simplex.Xncent()
	shrink_flag = 0
	Vprint(3, "Reflecting worst point:")
#	reflect worst vertex through the center of the rest:
	for j in range(dim):
		Xj = Xcent[j] + alpha*(Xcent[j]-Xworst[1][j])
		Xref[1].append(Xj)
#	calculate function value at reflection point:
	try:
		New_Val = Func(Xref[1],Fparam)
		if not New_Val in ["",None]:
			Xref[0] = New_Val
			Vprint(4, "Reflection point:",Xref)
		else:
			raise ValueError
	except Exception:
		Vprint(2, "Could not calculate function value at reflection point:",Xref)
		Xref[0] = 10000.0
		return Xref
	
#  try expanding the simplex in the direction of Xref:
	Xbest = Simplex.Get_Pi(dim)
	beta = NM_param[1]
	Xexp = ["", []]
	Xcont =["",[]]
	if (Xbest[0] > Xref[0]) or test == 1:
#		if Xref looks like the best point, try expanding by beta
		Vprint(4, "expanding simplex towards reflection")
		for j in range(dim):
			Xj = Xcent[j] + beta*(Xref[1][j]-Xcent[j])
			Xexp[1].append(Xj)
#		calculate function value at expansion point:
		try:
			New_Val2 = Func(Xexp[1],Fparam)
			if not New_Val2 in ["",None]:
				Xexp[0] = New_Val2
				Vprint(4, "Expansion point:",Xexp)
			else:
				raise ValueError
		except Exception:
			Vprint(2, "Could not calculate function value at expansion point:",Xexp)
		
#	if Xref is bad, try value at contraction point:
	if dim >= 2:
		Xsecond = Simplex.Get_Pi(1)
		gamma = NM_param[2]
		if Xref[0] >= Xworst[0]:
			# if the reflection point worse than Xworst, contract even more
			gamma = -1*gamma
		if (Xsecond[0] <= Xref[0] ) or test == 1:
#			if Xref is better than the second worst, try contracting back towards Xcent
			Vprint(4, "contracting reflection towards centroid")
			for j in range(dim):
				Xj = Xcent[j] + gamma*(Xref[1][j]-Xcent[j])
				Xcont[1].append(Xj)
		#	calculate function value at contraction point:
			try:
				New_Val3 = Func(Xcont[1],Fparam)
				if not New_Val3 in ["",None]:
					Xcont[0] = New_Val3
					Vprint(4, "Contraction point:",Xcont)
#					determine if we need a shrinking step instead a contraction:
					if gamma >= 0 and Xcont[0] > Xref[0]:
						shrink_flag = 1
					if gamma < 0 and Xcont[0] > Xworst[0]:
						shrink_flag = 1
				else:
					raise ValueError
			except Exception:
				Vprint(2, "Could not calculate function value at contraction point:",Xcont)
		
# determine which point is the best:	
	if Xexp[0] != "" and Xexp[0] < Xref[0]:
		Xnew = Xexp
	elif Xcont[0] != "" and Xcont[0] <= Xref[0] and shrink_flag == 0: 
		Xnew = Xcont
	elif Xcont[0] != "" and shrink_flag == 1:
		Xnew = Xworst
	else:
		Xnew = Xref
		
	return [Xnew, shrink_flag]

#function to shrink the whole the simplex:
def Shrink_Simplex(Simplex,NM_param, Func, Fparam):
	dim = Simplex.dim
	delta = NM_param[3]
	Xbest = Simplex.Get_Pi(dim)
	# all vertices except the best, shrink towards the center:
	for i in range(dim):
#		determine new coordinates for the shrunken point
		Xold = Simplex.Get_Pi(i)
		Xnew = ["",[]]
		Vprint(4, "Shrinking vertex:", Xold)
		for j in range(dim):
			Xj = Xbest[1][j] + delta*(Xold[1][j]-Xbest[1][j])
			Xnew[1].append(Xj)
#		Calculate vertex value at new point:
		try:
			New_Val = Func(Xnew[1],Fparam)
			if not New_Val in ["",None]:
				Xnew[0] = New_Val
				Vprint(4, "shrunken vertex:",Xnew)
			else:
				raise ValueError
		except Exception:
			Xnew[0] = 10000.0+Xold[0]
			Vprint(2, "Could not calculate function value at shrunken point:",Xnew)
#replace vertex and sort at the end:
		Simplex.Set_Pi(i,Xnew)
	Simplex.sort_points()
	
	return Simplex

#function to check if the simplex converged:
def Check_conv(Simplex,tolerance):
	Vprint(4, "Checking convergence")
	Diff = [tolerance*10,tolerance*10]

#	get the best vetrex, its value, and simplex dimensions:
	Svalues = Simplex.values()
	Best_val = Svalues[-1]
	dim = len(Svalues)-1
	Pb = Simplex.Get_Pi(dim)

#	calculate the difference in vertex values:
	Func_diff = []
	for i in range(dim):
		Fconv = math.fabs(Svalues[i] - Best_val)
		Func_diff.append(Fconv)
	Vprint(4,"Func_diff:", Func_diff)
	Diff[0] = max(Func_diff)

#	calcualte the largest distance between the best point and other points:
	Point_diff = []
	for i in range(dim):
		Pi = Simplex.Get_Pi(i)
		dist = 0.0
		for j in range(dim):
			dist += (Pb[1][j] - Pi[1][j])**2
		dist = math.sqrt(dist)
		Point_diff.append(dist)
	Vprint(4,"Poin_diff:", Point_diff)
	Diff[1] = max(Point_diff)
				
	return Diff
	

def Iterate_Min(Func , X0=[1.0], Fparam= [-1.0,[0.0]], Min_par= Min_param, NM_par= NM_param, *args):
	adaptive,tolerance, Max_iter,display = Min_par
	dim = len(X0)
	Final = []
	Log = []
	
#	stop if the point has 0 dimensions, or the function at X0 is not calculable:
	if dim == 0:
		Vprint("Error: X0 has zero dimensions!")
		return None
	F_x0 = ""
	try:
		F_x0 = Func(X0,Fparam)
	except Exception:
		Vprint(1, "Error: Function value cannot be determined at X0'")
		return None
# modify NM parameters if adaptive NM is requested:
	if adaptive == 1:
		NM_par = Pass_ANMparam(dim)
		Vprint(2, "Using adaptive parameters:",NM_par)
	else:
		Vprint(2, "Using simplex move parameters:",NM_par)	
	
#	iteratively move the simplex until convergence is reached:
	Vprint(2, "\nStarting Simplex minimization:")
	step = 0
	Current = []
	Prev_step = []
	step_info = ""
	while step <= Max_iter:
		if step == 0:
			Xi = Generate_Simp(Func, X0, Fparam, *args)
			Simp_iter = Simplex(Xi)
			Simp_iter.sort_points()
			Vprint(3, "\nGenerated initial simplex for %1d variables:" % dim)
			Vprint(3, Simp_iter.print_pi())
			
		else:
			Worst_vertex = Simp_iter.values()
			New_vertex,shrink_flag = Reflect_worst(Simp_iter, NM_par, Func, Fparam,0)
			if New_vertex[0] < Worst_vertex[0]:
				Simp_iter.Set_Pi(0, New_vertex)
				Simp_iter.sort_points()
			elif New_vertex[0] > Worst_vertex[0] or shrink_flag == 1:
				Simp_iter = Shrink_Simplex(Simp_iter, NM_par, Func, Fparam)
				
#		check results for this step:
		Best_point = Simp_iter.Get_Pi(dim)
		Worst_point = Simp_iter.Get_Pi(0)
		Current = [step, Best_point[0],Best_point[1], Worst_point[0],""]
		step_info = "Step: %1d, Best:  F(Xi)= %1.3f @ Xi ={" % (Current[0],Current[1])
		for coord in Current[2]:
			step_info += " %1.3f," % coord
		step_info = step_info[:-1]+"}"
		if step != 0:
			Diffs = Check_conv(Simp_iter,tolerance)	
#			diff = math.fabs(Current[3] - Prev_step[3])
			Current[4] = Diffs[0]
			step_info += " Worst: F(Xi)= %1.3f conv: %1.3e" % (Current[3],Current[4])
			
			if min(Diffs) < tolerance:
#			stop iterating if the change is smaller than the tolerance value
				Vprint(3,step_info)
				Vprint(2, "\nSimplex convergence reached: %1.4E (tol: %1.4E)" % (Current[4],tolerance))
				if Diffs[0] <= Diffs[1]:
					Vprint(2, "Reason: Function value tolerace (%1.1e)" % Diffs[0])
				else: 	
					Vprint(2, "Reason: Vertex distance tolerace (%1.1e)" % Diffs[1]) 	
				break
		Vprint(3,step_info)
#		prepare for the next step
		Log.append(Current)
		Prev_step = Current
		step += 1

#finish minization:
	Final_step = Current
	if step >= Max_iter:
		Vprint(2, "Maximum allowed iteration steps reached!")
	if display == 1:
		Display_string = "\nTotal number of iterations: %1d\nFinal simplex:\n" % step
		Display_string += Simp_iter.print_pi()
		Display_string += step_info
		Vprint(1, Display_string)
	
	
	return Final_step

	
	
#define harmonic test function:
def Test_func(Point,param):
	force_k = param[0]
	P0 = param[1]
	Diffs = []
	if len(Point) != len(P0):
		Vprint(2, "Mismatching, coordinate dimensions!")
		return None

	for j in range(len(Point)):
		Diffs.append((Point[j]-P0[j]))
		
	Value = 0.0
	for diff in Diffs:
		Value += diff**2
	Value = force_k * Value
	
	return Value

#function to generate an inital simplex around a point:
def Generate_Simp(Func, P0, Fparam, *args):
	dim_P = len(P0)
	argnum = len(args)
	random = 0
	step = 0.1
	step0 = 0.005
	if argnum >= 1:
		step = param[0]
	if argnum >= 2:
		step0 = param[1]
	if argnum >= 3:
		random = param[2]

#determine the value P0
	P0_value = Func(P0, Fparam)
	New_P0 = [P0_value,P0]
	New_Points = [New_P0]
# generate N points near the point:
	for i in range(dim_P):
		Pi = [0.0,[]]
		for j in range(dim_P):
			if j == i and P0[j] != 0.0:
				Pi[1].append((1.0+step)*P0[j])
			elif j == i and P0[j] == 0.0:
				Pi[1].append(step0)
			else:
				Pi[1].append(P0[j])
#		determine function value at Pi:
		Pi_value = Func(Pi[1],Fparam)
		if Pi_value != None:
			Pi[0] = Pi_value
#		append point:
		New_Points.append(Pi)
	return New_Points

	

#define function for unit tests:
def Run_Tests(NM_param, Min_param):
	
		Failed_Checks = 0
#		run basic unit tests:
		Vprint(2, "Starting unit test:")
		P0 = [0.0,0.0,0.0]
		P1 = [2.0,10.0,1.0]
		test_k = 1.0
		U_exp = 105.0
		test_Fparam = [test_k,P0]
		test_iter_param = [0, 0.0000001, 200, 1]
#		check test function:
		U_check = Test_func(P1, test_Fparam)
		if U_check == U_exp:
			Vprint(2, "Test on basic function passed: (Ucheck= %2.3f)" % U_check)
		else:
			Vprint(1, "Error, test on basic function failed (Ucheck= %2.3f, expected %2.3f)" % (U_check, U_exp))
			Failed_Checks += 1

#		generate initial points for a simplex:
		Test_points = Generate_Simp(Test_func,P1,test_Fparam)
		Vprint(2, "Generated simplex points:")
		for Pi in Test_points:
			Vprint(3, Pi)
		point_num = len(Test_points)
		point_exp = 4
		if len(Test_points) != point_exp:
				Failed_Checks += 1
				Vprint(1, "Initial Simplex generation failed, generated points: %1d (expected: %1d)" % (point_num,point_exp))
				
#		check basic simplex functions:				
		Simp0 = Simplex(Test_points)
		Simp0.sort_points()
		Vprint(2, Simp0.print_pi())
		
#		replace 1 point
		Lowpoint = Simp0.Get_Pi(0)
		Newpoint = [0.0, P0]
		Vprint(2, "Replacing Simplex vertex",Lowpoint," with ",Newpoint)
		Simp0.Set_Pi(0,Newpoint)
		Simp0.sort_points()
		Vprint(2, Simp0.print_pi())
		
#		check if everything is OK:
		Min_val = Simp0.points[-1][0]
		Vals = Simp0.values()
		Vprint(2, "Simplex values:",Vals)
		if Vals[-1] == Min_val:
			Vprint(2, "Simplex class test passed (Min value= %1.3f)" % Min_val)
		else:
			Vprint(1, "Simplex calss test failed: Min value= %1.3f, (expected %1.3f" % (Min_val, Vals[-1]))
			Failed_Checks += 1

#		check simplex operations:
		Vprint(2, "\nTesting Simplex moves:")
		Simp0.Set_Pi(3, Lowpoint)
		Simp0.sort_points()
		Vprint(2, Simp0.print_pi())
#		get centroid
		Xcent = Simp0.Xncent()
		Vprint(2, "Centroid:",Xcent)
#		reflect worst point:
		X_ref = Reflect_worst(Simp0, NM_param, Test_func, test_Fparam,1)
		Vprint(2, "\nRefelction point:",X_ref)
#		shrink simplex:
		Vprint(2, "\nShrinking Simplex:")
		Simp1 = Shrink_Simplex(Simp0,NM_param,Test_func,test_Fparam)
		Vprint(2, Simp1.print_pi())
#		iterate simplex:
		Test_iter = Iterate_Min(Test_func,P1,test_Fparam,test_iter_param, NM_param)
		
		return Failed_Checks



# Main function for testing module functions:
def Mintest_Main(Args):
#	set run parameters:
	Input_files, Out_files, Param, failmark,verbosity = Args
	Set_verb(Args[4])
	NM_param, Min_param, test = Param
	
	Main_Data = []
#	execute main code:

	if test == 1:
		Test_results = Run_Tests(NM_param, Min_param)
			
			
		
				
	return Main_Data
	



# executing standalone script:
if __name__ == '__main__':
#       handling command line arguments:        
	Arguments = []
	for arg in sys.argv:
		if not arg == sys.argv[0]:
			Arguments.append(arg)
	Custom_Args = Read_Args(Arguments)
	Vprint(2, "\nRun parameters:\n", Custom_Args)
	Input_files, Output_files, Param, failmark, verbosity = Custom_Args
	out_file = Output_files[0]

#       executing main code
	Data_main = Mintest_Main(Custom_Args)
	Vprint(2, "\nMain data:", Data_main)

#       print run-time messages
	ftime = time.time()
	runtime = ftime-stime
	outfiles = ""
	for File in Output_files:
		if File != "":
			outfiles += " %1s," % File
	outfile = outfiles[:-1]	
	Vprint(1, "\nScript runtime was %2.2f seconds" % runtime)
	Vprint(1, "Script finished sucessfully! Output written to:",outfiles)
else:
	Vprint(1,"Simplex minimzation module (SESCA_min)")

