#  Test file to check scan lines

import numpy as np
from matplotlib import pyplot as plt
import avg_cov_mat_across_angle_sweep as avg_cov






if __name__ == '__main__':

	# Jiao's test case for silicon:  

	scale = -5. -1./3
	H = 1.0*scale
	K = .5*scale
	L = -.5*scale
	E = 35.0

	# lattice parameters
	a = 5.4907
	b = 5.4907
	c = 5.4907

	# lattice vectors
	a1 = np.array([a, 0, 0])
	a2 = np.array([0, b, 0])
	a3 = np.array([0, 0, c])

	# create a list of lattice vectors
	lattice = [a1, a2, a3]


	# sample rotation angle
	#angle = 44.3427
	start_angle = -5.0
	stop_angle = 90.0
	num_angles = 100
	psi_range = np.linspace(start_angle, stop_angle, num=num_angles)
	theta_range = np.zeros(num_angles)
	phi_range = np.zeros(num_angles)

	# a bit of hacking:
	# define distances (eventually, this will be done automatically)
	L12 = 18.50 - 11.83  # distance from beam monitor 1 to 2
	Lms = 13.60  # distance from moderator to sample
	Lsp = 3.45  # distance from sample to detector pixel -- this will also need to be changed eventually to take a particular values


	# incident beam energy
	Ei_meV = 100.0  # meV

	u = [-1, 1, -1]
	v = [2, 1, -1]

	# now, convert to real units (inverse angstroms)
	h = H * 2*np.pi / a
	k = K * 2*np.pi / b
	l = L * 2*np.pi / c
	x = np.array([-1., 1., -1., 0.0])


	# make array of hklE values
	HKLE = np.array([h, k, l, E])


	for i in range(num_angles):
		psi = psi_range[i]

		# compute instrument coordinates
		instr_coords = avg_cov.HKLE_to_instrumentCoords(lattice, u, v, psi, HKLE, Ei_meV, L12, Lms, Lsp, 1)

		theta_range[i] = instr_coords[0]
		phi_range[i] = instr_coords[1]

	plt.ylabel("theta")
	plt.xlabel("phi (azimuthal)")
	plt.plot(phi_range, theta_range)
	plt.show()
