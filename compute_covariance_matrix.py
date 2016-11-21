'''
Program Purpose:  compute covariance matrices and inverses with respect to vQE and HKLE bases

Status: Working
'''

import numpy as np 
from numpy import linalg as LA
import compute_vQ_to_HKL_basis_change_matrix as HKL_basis


def get_sigma_vQE(J, M):
	J_T = np.transpose(J)

	A = np.dot(J,M)
	B = np.dot(A,J_T)
	Sigma = B

	return Sigma


def get_sigma_inv_vQE(J, M):
	Sigma = get_sigma_vQE(J,M)
	SigmaInv = LA.inv(Sigma)
	return SigmaInv


'''
Description:
	Arguments:

		Sigma_vQE (numpy 2D array): covariance matrix in vQ,E basis

		lattice_param_vectors (list):  [a1, a2, a3], where ai is a lattice vector (np array), in the crystal's own local Cartesian coordinate system, which depends on the lattice paramters a, b, c

		uList, vList (3-entry lists):  vectors in 3-dimensional HKL basis.  u is a vector in HKL which is parallel to the instrument z-axis (beam axis), and v is any vector in the xz plane which is linearly independent of u

		angle (float):  sample rotation angle
'''
def get_sigma_HKLE(Sigma_vQE, lattice_param_vectors, uList, vList, angle, debugMode=0):

	# get change of basis matrix
	P = HKL_basis.HKL_basis_change_matrix(lattice_param_vectors, uList, vList, angle,debugMode)

	# compute inverse matrix
	Pinv = LA.inv(P)

	# compute new sigma
	A = np.dot(P,Sigma_vQE)
	SigmaHKLE = np.dot(A, Pinv)

	# return covariance matrix in HKLE basis
	return SigmaHKLE




def get_sigma_inv_HKLE(SigmaHKLE):
	SigmaInv = LA.inv(SigmaHKLE)
	return SigmaInv



'''
Description:
	Summary:  This function takes two given independent vectors, x and y, in HKLE or vQE space, and then computes a 2-D covariance matrix of x and y
'''
def get_2D_sliced_covariance(Sigma, x, y=np.array([0.0, 0.0, 0.0, 1.0])):

	x = x / np.linalg.norm(x)
	y = y / np.linalg.norm(y)

	# make 2x4 matrix to map 4-D HKLE or vQE space into a 2D subspace
	transform = np.stack((x,y))

	Reduced_Sigma = np.dot(transform, Sigma)
	Reduced_Sigma = np.dot(Reduced_Sigma, np.transpose(transform))

	return Reduced_Sigma


def get_2D_sigma_inv(Sigma_Reduced):
	SigmaInv = LA.inv(Sigma_Reduced)
	return SigmaInv

##############################

if __name__ == "__main__":

	# test case:
	a = 3.0
	b = 3.0
	c = 3.0
	a1 = np.array([a, 0, 0])
	a2 = np.array([0, b, 0])
	a3 = np.array([0, 0, c])
	lattice = [a1, a2, a3]

	u = [1, 0, 2]
	v = [1, 0, 0]
	angle = 0.0

	basis_change = HKL_basis.HKL_basis_change_matrix(lattice, u, v, angle,1)

	print ("change of basis matrix  = \n")
	print (basis_change)

	import ARCS_error_analysis as arcs

	# test parameters
	t = 3900.9				# microseconds
	tof = t*10**-6 			# seconds
	theta = 2.24506			# radians (polar)
	phi = 0.0 - 0.56543		# radians (azimuthal)
	
	Ei_meV = 150.0

	JM = arcs.get_jacobian_and_params_matrices_from_event_data(tof, theta, phi, Ei_meV)

	J = JM[0]
	M = JM[1]

	print ("J = ")
	print (J)
	print ("\n")

	print ("M = ")
	print (M)
	print ("\n")

	Sigma_vQE = get_sigma_vQE(J, M)
	print ("Sigma_vQE = ")
	print (Sigma_vQE)
	print ("\n")

	Sigma_HKLE = get_sigma_HKLE(Sigma_vQE, lattice, u, v, angle)
	print ("Sigma_HKLE = ")
	print (Sigma_HKLE)
	print ("\n")

	SigmaInv_HKLE = get_sigma_inv_HKLE(Sigma_HKLE)
	print ("SigmaInv_HKLE = ")
	print (SigmaInv_HKLE)
	print ("\n")


	# test 4D -> 2D covariance reduction
	Reduced_Sigma = get_2D_sliced_covariance(Sigma_HKLE, np.array([2.0, 2.0, 0.0, 0.0]))
	print("Reduced_Sigma = ")
	print(Reduced_Sigma)


	# test display
	import plot_covariance_ellipse as plot_cov
	x1 = 2
	x2 = 3
	A = plot_cov.get_2D_subcovariance(Sigma_vQE, x1, x2)
	print ("A = ")
	print (A)
	print ("\n")
	k = 2
	alpha = 0.5
	chi2 = plot_cov.get_critical_chi_squared(k, alpha)

	plot_cov.plot_quadform(A, x1, x2, chi2, "Qz deviation (A^-1)", "E deviation (meV)")
