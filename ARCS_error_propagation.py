'''
Program Purpose:  perform error propagation for ARCS instrument; compute Jacobian with respect to instrument/experiment parameters

Status: Working
'''

import numpy as np 

###################################
# define instrumental parameters:

# Lms = length along beam from moderator to sample (meters)
Lms = 13.60

# Lsp = length along beam from sample to pixel
# TODO: get this from pixel data
Lsp = 3.0

# distance between beam monitors
L12 = 18.50 - 11.83


####################################
# define physical constants:

# m = mass of neutron (kg)
m = 1.674929*(10**-27)

# hbar = reduced Planck's constant (J s)
hbar = 1.0545718*(10**-34)


#####################################
# auxiliary functions to compute quantities from given info:

def meV_to_joules(Ei_meV):
	Ei = Ei_meV * 1.602177 * 10**(-22)
	return Ei

def joules_to_meV(E):
	E_meV = E * 6.242 * 10**21
	return E_meV

def vi_from_Ei(Ei):
	vi = np.sqrt(2.0*Ei/m)
	return vi

def t12_from_vi(vi):
	t12 = L12 / vi
	return t12

def vf_from_tof_and_t12(tof, t12, Lsp, Lms, L12):
	D = tof - t12 * (Lms / L12)
	N = Lsp
	vf = N / D
	return vf

def magnitude_Q(vf, vi):
	Q = (m/hbar) * np.sqrt((vf - vi)**2)
	return Q

def get_E(Ei, Ef):
	return Ei - Ef   # yes, it is supposed to be backwards

def magnitude_E(E):
	magE = np.sqrt(E**2)
	return magE

#####################################

# functions to compute partial derivatives

def vi_partial_t12(L12, t12):
	partial = 0.0 - L12 / (t12**2)
	return partial

def vf_partial_t12(Lsp, Lms, L12, tof, t12):
	A = Lsp * Lms
	B = L12 * (tof - t12*(Lms / L12))**2
	partial = A / B
	return partial

# Qz partials
def Qz_partial_t12(theta, L12, t12, Lsp, Lms, tof):
	vf_t12 = vf_partial_t12(Lsp, Lms, L12, tof, t12)
	A = np.cos(theta) * vf_t12 + L12 / (t12**2)
	partial = (m / hbar) * A
	return partial

def Qz_partial_theta(vf, theta):
	partial = 0.0 - (m/hbar) * vf * np.sin(theta)
	return partial

def Qz_partial_phi():
	return 0.0

# Qx partials
def Qx_partial_t12(Lsp, Lms, L12, tof, t12, theta, phi):
	vf_t12 = vf_partial_t12(Lsp, Lms, L12, tof, t12)
	partial = (m/hbar)*np.sin(theta)*np.cos(phi) * vf_t12
	return partial

def Qx_partial_theta(vf, theta, phi):
	partial = (m/hbar) * vf * np.cos(theta) * np.cos(phi)
	return partial

def Qx_partial_phi(vf, theta, phi):
	partial = 0.0 - (m/hbar) * vf * np.sin(theta) * np.sin(phi)
	return partial

# Qy partials
def Qy_partial_t12(Lsp, Lms, L12, tof, t12, theta, phi):
	vf_t12 = vf_partial_t12(Lsp, Lms, L12, tof, t12)
	partial = (m/hbar)*np.sin(theta)*np.sin(phi) * vf_t12
	return partial

def Qy_partial_theta(vf, theta, phi):
	partial = (m/hbar) * vf * np.cos(theta) * np.sin(phi)
	return partial

def Qy_partial_phi(vf, theta, phi):
	partial = (m/hbar) * vf * np.sin(theta) * np.cos(phi)
	return partial

# E partials
def E_partial_t12(vi, vf, vi_t12, vf_t12):
	partial = m*(vi*vi_t12 - vf*vf_t12)
	return partial


##################################

# Set up matrix (Jacobian)
def setup_jacobian(vi, vf, Ei, Ef, theta, phi, L12, Lms, Lsp, tof, t12, debugMode=1):
	
	magQ = magnitude_Q(vf, vi)
	E = get_E(Ei, Ef)
	magE = magnitude_E(E)

	vi_t12 = vi_partial_t12(L12, t12)
	vf_t12 = vf_partial_t12(Lsp, Lms, L12, tof, t12)

	Qz_t12 = Qz_partial_t12(theta, L12, t12, Lsp, Lms, tof)
	Qz_theta = Qz_partial_theta(vf, theta)
	Qz_phi = Qz_partial_phi()

	Qx_t12 = Qx_partial_t12(Lsp, Lms, L12, tof, t12, theta, phi)
	Qx_theta = Qx_partial_theta(vf, theta, phi)
	Qx_phi = Qx_partial_phi(vf, theta, phi)

	Qy_t12 = Qy_partial_t12(Lsp, Lms, L12, tof, t12, theta, phi)
	Qy_theta = Qy_partial_theta(vf, theta, phi)
	Qy_phi = Qy_partial_phi(vf, theta, phi)

	E_t12 = E_partial_t12(vi, vf, vi_t12, vf_t12)
	if debugMode == 1:
		print "E_partial_t12 = " + str(E_t12) + "\n"

	E_theta = 0.0
	E_phi = 0.0

	# define matrix entries (converted to inverse angstroms)
	J_11 = 10**-10 * Qx_t12 #/ magQ 
	J_12 = 10**-10 * Qx_theta #/ magQ 
	J_13 = 10**-10 * Qx_phi #/ magQ 

	J_21 = 10**-10 * Qy_t12 #/ magQ 
	J_22 = 10**-10 * Qy_theta #/ magQ 
	J_23 = 10**-10 * Qy_phi #/ magQ 

	J_31 = 10**-10 * Qz_t12 #/ magQ 
	J_32 = 10**-10 * Qz_theta #/ magQ 
	J_33 = 10**-10 * Qz_phi #/ magQ 

	# converted to meV
	J_41 = joules_to_meV(E_t12) #/ magE
	J_42 = joules_to_meV(E_theta) #/ magE 
	J_43 = joules_to_meV(E_phi) #/ magE 

	J = np.array([[J_11, J_12, J_13], [J_21, J_22, J_23], [J_31, J_32, J_33], [J_41, J_42, J_43]])

	return J



#####################################3

if __name__ == "__main__":

	# test case:
	tof = 3900.9			# microseconds
	t = tof*10**-6 			# seconds
	Lsp = 3.44735 			# meters
	theta = 2.24506			# radians (polar)
	phi = 0.0 - 0.56543		# radians (azimuthal)
	Ei_meV_expected = 100.0
	Ei = meV_to_joules(Ei_meV_expected)
	vi = vi_from_Ei(Ei)
	t12 = t12_from_vi(vi)
	vf = vf_from_tof_and_t12(t, t12, Lsp, Lms, L12)
	Ef = 0.5*m*vf**2

	# (guessed) uncertainties:
	sigma_t12 = 10.0*10**-6
	sigma_theta = 1.5 / 360.0 * 2.0*np.pi
	sigma_phi = 1.5 / 360.0 * 2.0*np.pi

	var_t12 = sigma_t12**2
	var_theta = sigma_theta**2
	var_phi = sigma_phi**2

	J = setup_jacobian(vi, vf, Ei, Ef, theta, phi, L12, Lms, Lsp, tof, t12)

	print "J = "
	print J
	print "\n"