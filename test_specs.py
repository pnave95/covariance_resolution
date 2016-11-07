# test case specification file:

import numpy as np
import compute_covariance_matrix as cov
import ARCS_error_propagation as arcs
import plot_covariance_ellipse as plot_cov
import compute_vQ_to_HKL_basis_change_matrix as HKL_basis 

# test case:
# a = 3.0
# b = 3.0
# c = 3.0
a = 5.4907
b = 5.4907
c = 5.4907
# a1 = np.array([a, 0, 0])
# a2 = np.array([0, b, 0])
# a3 = np.array([0, 0, c])
a1 = np.array([a, 0, 0])
a2 = np.array([0, b, 0])
a3 = np.array([0, 0, c])
# a1 = np.array([0, 2.715, 2.715])
# a2 = np.array([2.715, 0, 2.715])
# a3 = np.array([2.715, 2.715, 0])
lattice = [a1, a2, a3]

#u = [1, 0, 2]
#v = [1, 0, 0]
u = np.array([-1, 1, -1])
v = np.array([2, 1, -1])
#angle = 0.0
angle = 0.773926

basis_change = HKL_basis.HKL_basis_change_matrix(lattice, u, v, angle,1)

print "change of basis matrix  = \n"
print basis_change

Lms = 13.60
L12 = 18.50 - 11.83

# m = mass of neutron (kg)
m = 1.674929*(10**-27)

# hbar = reduced Planck's constant (J s)
hbar = 1.0545718*(10**-34)

# test parameters
# tof = 3900.9			# microseconds
# t = tof*10**-6
# Lsp = 3.44735 			# meters
# theta = 2.24506			# radians (polar)
# phi = 0.0 - 0.56543		# radians (azimuthal)
Ei_meV_expected = 100.0
Ei = arcs.meV_to_joules(Ei_meV_expected)
vi = arcs.vi_from_Ei(Ei)
t12 = arcs.t12_from_vi(vi)
# vf = arcs.vf_from_tof_and_t12(t, t12, Lsp, Lms, L12)
# Ef = 0.5*m*vf**2

# # (guessed) uncertainties:
# sigma_t12 = 10.0*10**-6
sigma_theta = 1.5 / 360.0 * 2.0*np.pi
sigma_phi = 1.5 / 360.0 * 2.0*np.pi

# test silicon parameters
tof = 3968.2
t = tof*10**-6
Lsp = np.sqrt(1.3915**2 + 2.65775**2)
r_x = 1.39152986 # meters
r_y = 0.0
r_z = 2.6577518
r = Lsp
theta = np.arctan(r_x / r_z )
phi = np.pi / 2.0
vf = arcs.vf_from_tof_and_t12(t, t12, Lsp, Lms, L12)
Ef = 0.5*m*vf**2

sigma_t12 = 8.25 * 10**-6 




var_t12 = sigma_t12**2
var_theta = sigma_theta**2
var_phi = sigma_phi**2

J = arcs.setup_jacobian(vi, vf, Ei, Ef, theta, phi, L12, Lms, Lsp, tof, t12)

print "J = "
print J
print "\n"

M = cov.setup_params_matrix(var_t12, var_theta, var_phi)
print "M = "
print M
print "\n"

Sigma_vQE = cov.get_sigma_vQE(J, M)
print "Sigma_vQE = "
print Sigma_vQE
print "\n"

Sigma_HKLE = cov.get_sigma_HKLE(Sigma_vQE, lattice, u, v, angle)
print "Sigma_HKLE = "
print Sigma_HKLE
print "\n"

SigmaInv_HKLE = cov.get_sigma_inv_HKLE(Sigma_HKLE)
print "SigmaInv_HKLE = "
print SigmaInv_HKLE
print "\n"


# test display
x1 = 0
x2 = 2
A = plot_cov.get_2D_subcovariance(Sigma_vQE, x1, x2)
print "A = "
print A
print "\n"
plot_cov.plot_quadform(A, x1, x2, "Qx deviation (A^-1)", "E deviation (meV)")
