'''
Program Purposes: 
	- define instrument-dependent constants + parameters with uncertainties
	- create parameter covariance matrix
	- perform error propagation for ARCS instrument; 
	- compute Jacobian with respect to instrument/experiment parameters

Notes:
	Currently, this program considers the quantity "Lsp" (distance, or length of beam travel, from the sample to the pixel) to be fixed by default, with an uncertainty which is also set to some default value.
		Eventually, I would like to automatically look up pixel location (radial spherical coordinate = Lsp) to more accurately estimate this value

Status: Tentatively complete
'''

# import required modules
import numpy as np 

####################################
# define physical constants:

# m = mass of neutron (kg)
m = 1.674929*(10**-27)

# hbar = reduced Planck's constant (J s)
hbar = 1.0545718*(10**-34)



#####################################
# define instrumental parameters:

# Lms = length along beam from moderator to sample
Lms = 13.60				# meters

# Lsp = length along beam from sample to pixel
# TODO: get this from pixel data
default_Lsp = 3.0		# meters

# distance between beam monitors
L12 = 18.50 - 11.83		# meters


######################################
# define instrumental parameter uncertainties
sigma_t12 = 2.5*10**-6 				# seconds
sigma_theta = 0.25 / 360.0 * 2.0*np.pi 	# radians
sigma_phi = 0.25 / 360.0 * 2.0*np.pi 	# radians
sigma_Lsp = 0.015						# meters

var_t12 = sigma_t12**2
var_theta = sigma_theta**2
var_phi = sigma_phi**2
var_Lsp = sigma_Lsp**2

######################################
# unit conversion assistant functions

def meV_to_joules(E_meV):
	E = E_meV * 1.602177 * 10**(-22)
	return E

def joules_to_meV(E):
	E_meV = E * 6.242 * 10**21
	return E_meV

######################################
# intermediate variable functions

# time from moderator to sample from initial speed
def get_tms_from_vi(vi):
	tms = Lms / vi
	return tms

def get_tms_from_t12(t12):
	tms = t12 * Lms / L12
	return tms

# sample-to-pixel time from time-of-flight and initial speed
def get_tsp_from_tof_vi(tof, vi):
	t_ms = tms_from_vi(vi)
	t_sp = tof - t_ms
	return t_sp

# convert energy to speed (*note: this E is NOT the energy transfer, which is also confusingly called "E" by neutron scientists)
def get_v_from_E(E):
	v = np.sqrt(2.0*E/m)
	return v

# time to travel from beam monitor 1 to 2 from initial speed
def get_t12_from_vi(vi):
	t12 = L12 / vi
	return t12

# final speed from tof and t12
def get_vf_from_tof_t12_Lsp(tof, t12, Lsp):
	D = tof - t12 * (Lms / L12)
	N = Lsp
	vf = N / D
	return vf

# alternatively, final speed from tms
def get_vf_from_tms_Lsp(tms, Lsp):
	vf = Lsp / tms
	return vf

# TODO: handle the case where this becomes infinite (near vertical scattering)
def get_Lsp_from_pixel_angles(theta, phi):
	s = np.sin(theta) * np.cos(phi - np.pi/2.0)
	Lsp = default_Lsp / np.sqrt(1.0 - s**2)
	return Lsp

######################################
# compute primary vQE variables

def get_Qz_inv_meters(vf, vi, theta):
	Qz = (m / hbar) * (vf*np.cos(theta) - vi)
	return Qz

def get_Qx_inv_meters(vf, theta, phi):
	Qx = (m / hbar) * vf * np.sin(theta) * np.cos(phi)
	return Qx

def get_Qy_inv_meters(vf, theta, phi):
	Qy = (m / hbar) * vf * np.sin(theta) * np.sin(phi)
	return Qy

# energy transfer
def get_E_J(vf, vi):
	E = (m / 2.0)*(vi**2 - vf**2) 
	return E

######################################
# partial derivative functions (in seconds, meters, joules, radians)

# speed partials
def vi_partial_t12(t12):
	partial = 0.0 - L12 / (t12**2)
	return partial


def vf_partial_t12(tof, t12, Lsp):
	A = Lsp * Lms
	B = L12 * (tof - t12*(Lms / L12))**2
	partial = A / B
	return partial

def vf_partial_Lsp(tof, tms):
	partial = 1.0 / (tof - tms)
	return partial


# Qz partials
def Qz_partial_t12(theta, t12, tof, Lsp):
	vf_t12 = vf_partial_t12(Lsp, tof, t12)
	A = np.cos(theta) * vf_t12 + L12 / (t12**2)
	partial = (m / hbar) * A
	return partial

def Qz_partial_Lsp(vf_Lsp, theta):
	partial = (m / hbar) * vf_Lsp * np.cos(theta)
	return partial 

def Qz_partial_theta(vf, theta):
	partial = 0.0 - (m/hbar) * vf * np.sin(theta)
	return partial

def Qz_partial_phi():
	return 0.0


# Qx partials
def Qx_partial_t12(tof, t12, theta, phi, Lsp):
	vf_t12 = vf_partial_t12(Lsp, tof, t12)
	partial = (m/hbar)*np.sin(theta)*np.cos(phi) * vf_t12
	return partial

def Qx_partial_Lsp(vf_Lsp, theta, phi):
	partial = (m / hbar) * vf_Lsp * np.sin(theta) * np.cos(phi)
	return partial

def Qx_partial_theta(vf, theta, phi):
	partial = (m/hbar) * vf * np.cos(theta) * np.cos(phi)
	return partial

def Qx_partial_phi(vf, theta, phi):
	partial = 0.0 - (m/hbar) * vf * np.sin(theta) * np.sin(phi)
	return partial

# Qy partials
def Qy_partial_t12(tof, t12, theta, phi, Lsp):
	vf_t12 = vf_partial_t12(Lsp, tof, t12)
	partial = (m/hbar)*np.sin(theta)*np.sin(phi) * vf_t12
	return partial

def Qy_partial_Lsp(vf_Lsp, theta, phi):
	partial = (m / hbar) * vf_Lsp * np.sin(theta) * np.sin(phi)
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

def E_partial_Lsp(vf_Lsp, vf):
	partial = 0.0 - m * vf * vf_Lsp
	return partial




#######################################

# Set up matrix (Jacobian)
def setup_jacobian(vi, vf, theta, phi, tof, t12, Lsp=default_Lsp, debugMode=1):

	Lsp = get_Lsp_from_pixel_angles(theta, phi)
	
	E = get_E_J(vf, vi)
	tms = get_tms_from_vi(vi)

	vi_t12 = vi_partial_t12(t12)
	vf_t12 = vf_partial_t12(tof, t12, Lsp)
	vf_Lsp = vf_partial_Lsp(tof, tms)

	Qz_t12 = Qz_partial_t12(theta, t12, Lsp, tof)
	Qz_Lsp = Qz_partial_Lsp(vf_Lsp, theta)
	Qz_theta = Qz_partial_theta(vf, theta)
	Qz_phi = Qz_partial_phi()

	Qx_t12 = Qx_partial_t12(tof, t12, theta, phi, Lsp)
	Qx_Lsp = Qx_partial_Lsp(vf_Lsp, theta, phi)
	Qx_theta = Qx_partial_theta(vf, theta, phi)
	Qx_phi = Qx_partial_phi(vf, theta, phi)

	Qy_t12 = Qy_partial_t12(tof, t12, theta, phi, Lsp)
	Qy_Lsp = Qy_partial_Lsp(vf_Lsp, theta, phi)
	Qy_theta = Qy_partial_theta(vf, theta, phi)
	Qy_phi = Qy_partial_phi(vf, theta, phi)

	E_t12 = E_partial_t12(vi, vf, vi_t12, vf_t12)
	E_Lsp = E_partial_Lsp(vf_Lsp, vf)

	if debugMode == 1:
		print ("E_partial_t12 = " + str(E_t12) + "\n")

	E_theta = 0.0
	E_phi = 0.0

	# define matrix entries (converted to inverse angstroms)
	J_11 = 10**-10 * Qx_t12
	J_12 = 10**-10 * Qx_theta 
	J_13 = 10**-10 * Qx_phi 
	J_14 = 10**-10 * Qx_Lsp 

	J_21 = 10**-10 * Qy_t12 
	J_22 = 10**-10 * Qy_theta
	J_23 = 10**-10 * Qy_phi 
	J_24 = 10**-10 * Qy_Lsp

	J_31 = 10**-10 * Qz_t12 
	J_32 = 10**-10 * Qz_theta 
	J_33 = 10**-10 * Qz_phi 
	J_34 = 10**-10 * Qz_Lsp

	# converted to meV
	J_41 = joules_to_meV(E_t12) #/ magE
	J_42 = joules_to_meV(E_theta) #/ magE 
	J_43 = joules_to_meV(E_phi) #/ magE 
	J_44 = joules_to_meV(E_Lsp)

	J = np.array([[J_11, J_12, J_13, J_14], [J_21, J_22, J_23, J_24], [J_31, J_32, J_33, J_34], [J_41, J_42, J_43, J_44]])

	return J


def setup_params_matrix():
	M = np.array([[var_t12, 0.0, 0.0, 0.0], [0.0, var_theta, 0.0, 0.0], [0.0, 0.0, var_phi, 0.0], [0.0, 0.0, 0.0, var_Lsp]])
	return M



def get_jacobian_and_params_matrices_from_event_data(tof, theta, phi, Ei_meV, Lsp=default_Lsp, debugMode=1):

	Lsp = get_Lsp_from_pixel_angles(theta, phi)

	Ei = meV_to_joules(Ei_meV)
	vi = get_v_from_E(Ei)
	t12 = get_t12_from_vi(vi)
	vf = get_vf_from_tof_t12_Lsp(tof, t12, Lsp)

	if debugMode == 1:
		print("t12 = " + str(t12))

	J = setup_jacobian(vi, vf, theta, phi, tof, t12)

	M = setup_params_matrix()

	return [J, M]



#######################################

if __name__ == "__main__":

	# test case:
	t = 3900.9				# microseconds
	tof = t*10**-6 			# seconds
	theta = 2.24506			# radians (polar)
	phi = 0.0 - 0.56543		# radians (azimuthal)
	
	Ei_meV = 150.0

	JM = get_jacobian_and_params_matrices_from_event_data(tof, theta, phi, Ei_meV)

	J = JM[0]
	M = JM[1]

	#J = setup_jacobian(vi, vf, Ei, Ef, theta, phi, L12, Lms, Lsp, tof, t12)

	print ("J = ")
	print (J)
	print ("\n")

	Sigma = np.transpose(J)
	Sigma = np.dot(M, Sigma)
	Sigma = np.dot(J, Sigma)

	print ("Sigma = ")
	print (Sigma)