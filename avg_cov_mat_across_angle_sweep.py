'''
Program Purpose:  Compute average covariance matrix for a particular (H,K,L,E) point, as the sample is rotated through a range of angles

Status:  INCOMPLETE / NON-FUNCTIONAL
'''


import numpy as np
from numpy import linalg as LA 
import compute_vQ_to_HKL_basis_change_matrix as HKL_basis
import ARCS_error_analysis as arcs


'''
Description:
	Returns:
		instrument_coords:  [theta, phi, t12, tof], where theta is the spherical polar angle of the detector pixel, phi is the azimuthal angle, t12 is the time between detector pixel 1 and 2, and tof is the time-of-flight (seconds)
'''
def HKLE_to_instrumentCoords(lattice, u, v, sample_angle, HKLE, Ei_meV, L12, Lms, Lsp, debugMode=1):

	# extract H,K,L,E values (left in inverse Angstrom units)
	# H = HKLE[0]
	# K = HKLE[1]
	# L = HKLE[2]
	# E = HKLE[3]

	# compute vQE to HKLE change of basis matrix
	vQE_to_HKLE_mat = HKL_basis.HKL_basis_change_matrix(lattice, u, v, sample_angle, debugMode)

	# compute inverse basis change matrix (HKLE to vQE)
	HKLE_to_vQE_mat = LA.inv(vQE_to_HKLE_mat)

	# compute vQE at specified HKLE point
	vQE = np.dot(HKLE_to_vQE_mat, HKLE)

	# check if E is correct -- if not, change it:

	# first, compute scalar Q (in inverse Angstroms)
	Q_magnitude = np.sqrt(vQE[0]**2 + vQE[1]**2 + vQE[2]**2)
	# change units to inverse meters
	#Q_inverse_meters = Q_magnitude * 10**-10
	Q_inverse_meters = Q_magnitude * 10**10

	# compute final energy
	E = HKLE[3]  # E = Ei - Ef
	Ef_meV = Ei_meV - E

	# compute initial and final energies in joules
	Ei = arcs.meV_to_joules(Ei_meV)
	Ef = arcs.meV_to_joules(Ef_meV)

	# compute initial and final speeds (m/s)
	vi = arcs.get_v_from_E(Ei)
	vf = arcs.get_v_from_E(Ef)

	# compute instrument coordinates corresponding to the specified vQE point, for the given sample rotation angle and Ei

	Qx = vQE[0]
	Qy = vQE[1]
	Qz = vQE[2]
	Q = np.sqrt(Qx**2 + Qy**2 + Qz**2)

	# polar angle of detector pixel
	theta = np.arccos(Qz / Q)

	# Experimental test:  numpy only ouputs in [0,pi], but I believe I want the range [-pi/2, pi/2]:
	#if theta > np.pi / 2.0:
	#	theta = theta - np.pi

	Qx_tmp1 = Qx / Q
	Qx_tmp2 = Qx_tmp1 / np.sin(theta)

	Qy_tmp1 = Qy / Q
	Qy_tmp2 = Qy_tmp1 / np.sin(theta)

	quadrants = set([1,2,3,4])
	final_quadrant = 0

	if np.sign(Qx) > 0:
		quadrants.discard(2)
		quadrants.discard(3)
	else:
		quadrants.discard(1)
		quadrants.discard(4)
	if np.sign(Qy) > 0:
		quadrants.discard(3)
		quadrants.discard(4)
	else:
		quadrants.discard(1)
		quadrants.discard(2)
	if 1 in quadrants:
		final_quadrant = 1
		phi = np.arccos(Qx_tmp2)
	elif 2 in quadrants:
		final_quadrant = 2
		phi = np.arccos(Qx_tmp2)
	elif 3 in quadrants:
		final_quadrant = 3
		phi_part1 = np.arccos(Qx_tmp2)
		angle_to_x_axis = np.pi - phi_part1
		phi = phi_part1 + 2.0 * angle_to_x_axis
	elif 4 in quadrants:
		final_quadrant = 4
		phi = np.arcsin(Qy_tmp2)

	# azimuthal angle of detector pixel
	# phi = np.arccos(Qx_tmp2)

	# phi2 = np.arcsin(Qy_tmp2) #- np.pi / 2.0
	# if phi2 < 0.0:
	# 	phi2 = phi2 + np.pi
	# if abs(phi - phi2) > 0.001:
	# 	return "Error:  unable to compute azimuthal angle phi:  phi = " + str(phi) + ", phi2 = " + str(phi2)

	# compute t12 (seconds)
	t12 = L12 / vi

	# now, compute tof (seconds)
	tof = (Lms / vi) + (Lsp / vf)

	instrument_coords = [theta, phi, t12, tof]

	return instrument_coords






##########################

#if __name__ == "__main__":

	# do something to test






