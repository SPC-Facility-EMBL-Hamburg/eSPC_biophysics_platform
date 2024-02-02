import numpy as np

def normalised_euclidean_distance(v1,v2):
	
	'''
	Takes two np arrays and computes the normalised euclidean distance

	'''

	v1_diff = v1 - np.mean(v1)
	v2_diff = v2 - np.mean(v2)

	num = np.linalg.norm(  v1_diff - v2_diff )**2
	den = np.linalg.norm(  v1_diff )**2 + np.linalg.norm(  v2_diff )**2

	ned = num / den * 0.5

	return ned