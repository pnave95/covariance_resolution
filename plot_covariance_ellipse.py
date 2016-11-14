'''
Program Purpose:  plotting and visualization of resolution ellipsoids

Status: INCOMPLETE
'''

import numpy as np 
from numpy import linalg as LA 
from matplotlib import pyplot as plt



# Compute critical chi-squared value for given confidence (alpha) level and number of degrees of freedom
def get_critical_chi_squared(k, alpha):
	import scipy.special as sp
	chi2 = sp.chdtri(k,alpha)
	return chi2

# Extract 2x2 submatrix from covariance matrix for plotting
def get_2D_subcovariance(C, x1, x2):

	# extract relevant entries
	A_11 = C[x1,x1]
	A_12 = C[x1,x2]
	A_21 = C[x2,x1]
	A_22 = C[x2,x2]
	A = np.array([[A_11, A_12], [A_21, A_22]])

	return A


# plot 2x2 quadratic form matrix A
def plot_quadform(A, x1, x2, chi2, x1_title, x2_title, debugMode=1):

	# compute eigenvalues and eigenvectors of the covariance submatrix A
	L, v = LA.eigh(A)

	# compute angle between first eigenvector and "x1" axis
	theta = np.arctan(v[0,1] / v[0,0])

	# compute coefficients for ellipse equation
	a = np.sqrt(L[0]) * np.sqrt(chi2)
	b = np.sqrt(L[1]) * np.sqrt(chi2)

	# debugging stuff
	if debugMode==1:
		print ("\n")
		print ("eigenvalues: ")
		print (L)
		print ("\n")
		print ("eigenvectors: \n")
		print (v )
		print ("\n")
		print ("angle = " + str(theta) + "\n")
		print ("(a,b) = " + str(a) + ", " + str(b) + " \n")

	# create a vector of angles from 0 to 2pi for graphing
	phi = np.linspace(0.0, 2.0*np.pi, num=200)

	# compute x-prime and y-prime values (ellipse in eigencoordinates)
	xp = a*np.cos(phi)
	yp = b*np.sin(phi)

	# stack x-prime, y-prime values into a matrix
	Xp = np.stack((xp, yp))
	
	if debugMode==1:
		print ("np.shape(Xp) = " + str(np.shape(Xp)))

	# create rotation matrix to "move" the ellipse back into x,y (actually: x1,x2) coordinates instead of eigencoordinates
	R = np.array([[np.cos(theta), - np.sin(theta)], [np.sin(theta), np.cos(theta)]])

	# compute the coordinates for each ellipse point in the original x1,x2 coordinates
	X = np.dot(R,Xp)
	x = X[0,:]
	y = X[1,:]

	# plotting
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel(x1_title)
	ax.set_ylabel(x2_title)
	plt.plot(x,y)
	plt.savefig("testing")
	plt.show()


def plot_quadform_method2(A, x1, x2, chi2, x1_title, x2_title, debugMode=1):

	#alternative version (this produces ellipses which appear more segmented for some reason)
	Ainv = LA.inv(A)
	print ("A^-1 = \n")
	print (Ainv)

	# create vector of angles from 0 to 2pi for graphing
	phi = np.linspace(0.0, 2.0*np.pi, num=200)

	xx = np.cos(phi)
	yy = np.sin(phi)
	xy = np.stack((xx,yy))
	xyT = np.transpose(xy)
	print ("np.shape(xyT) = " + str(np.shape(xyT)))

	for i in range(len(xx)):
		XY = xy[:,i]
		XYT = np.transpose(XY)
		r1 = np.dot(XYT, Ainv)
		r2 = np.dot(r1,XY)
		r = np.sqrt(chi2 / r2)
		xy[:,i] = r*xy[:,i]

	xx = xy[0,:]
	yy = xy[1,:]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel(x1_title)
	ax.set_ylabel(x2_title)
	plt.plot(xx,yy)
	plt.savefig("testing2")
	plt.show()


'''
INCOMPLETE
'''
def plot_3D_confidence_ellipsoid(A, x1, x2, x3, chi2, x1_title, x2_title, x3_title, debugMode=1):

	# compute inverse of 3x3 matrix A
	Ainv = LA.inv(A)

	# create vector of polar and azimuthal angles
	phi = np.linspace(0.0, 2.0*np.pi, num=200)
	theta = np.linspace(0.0, np.pi, num=100)
	#spherical_grid = np.meshgrid(phi, theta)

	x = np.zeros((100,200))
	y = np.zeros((100,200))
	z = np.zeros((100,200))

	for i in range(100):
		for j in range(200):
			x[i][j] = np.sin(theta[i]) * np.cos(phi[i])
			y[i][j] = np.sin(theta[i]) * np.sin(phi[i])
			z[i][j] = np.cos(theta[i])

	xx = np.flat(x)
	yy = np.flat(y)
	zz = np.flat(z)

	# xyz = np.stack((xx, yy, zz))
	# xyT = np.transpose(xy)
	# print ("np.shape(xyT) = " + str(np.shape(xyT)))

	# x = np.zeros(spherical_grid)
	# y = np.zeros(spherical_grid)
	# z = np.zeros(spherical_grid)


###############################

if __name__ == "__main__":

	# test 
	A = np.array([[1,2,3], [2,4.5, 1.4], [3, 1.4, 2.2]])

	x1 = 0
	x1_title = "x"
	x2 = 1
	x2_title = "y"
	x3 = 2
	x3_title = "z"
	k = 3
	alpha = 0.5
	chi2 = get_critical_chi_squared(k, alpha)
	plot_3D_confidence_ellipsoid(A, x1, x2, x3, chi2, x1_title, x2_title, x3_title)