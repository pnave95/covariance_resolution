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
def plot_quadform(A, x1, x2, x1_title, x2_title, debugMode=1):
	from matplotlib import pyplot as plt
	from numpy import linalg as LA


	L, v = LA.eigh(A)
	theta = np.arctan(v[0,1] / v[0,0])
	a = np.sqrt(L[0]) * np.sqrt(1.386)
	b = np.sqrt(L[1]) * np.sqrt(1.386)

	if debugMode==1:
		print "\n"
		print "eigenvalues: "
		print L
		print "\n"
		print "eigenvectors: \n"
		print v 
		print "\n"
		print "angle = " + str(theta) + "\n"
		print "(a,b) = " + str(a) + ", " + str(b) + " \n"

	phi = np.linspace(0.0, 2.0*np.pi, num=200)
	xp = a*np.cos(phi)
	yp = b*np.sin(phi)
	#X = np.array([[a*np.cos(phi)], [b*np.sin(phi)]])
	Xp = np.stack((xp, yp)) #.reshape(2,50)
	
	print "np.shape(Xp) = " + str(np.shape(Xp))
	R = np.array([[np.cos(theta), - np.sin(theta)], [np.sin(theta), np.cos(theta)]])

	X = np.dot(R,Xp)
	x = X[0,:]
	y = X[1,:]
	fig = plt.figure()
	ax = fig.add_subplot(111)
	ax.set_xlabel(x1_title)
	ax.set_ylabel(x2_title)
	plt.plot(x,y)
	plt.savefig("testing")
	plt.show()

	#alternative version (this produces ellipses which appear more segmented for some reason)
	Ainv = LA.inv(A)
	print "A^-1 = \n"
	print Ainv

	xx = np.cos(phi)
	yy = np.sin(phi)
	xy = np.stack((xx,yy))
	xyT = np.transpose(xy)
	print "np.shape(xyT) = " + str(np.shape(xyT))

	for i in range(len(xx)):
		XY = xy[:,i]
		XYT = np.transpose(XY)
		r1 = np.dot(XYT, Ainv)
		r2 = np.dot(r1,XY)
		#r = np.sqrt(1.386 / abs(r2))
		r = np.sqrt(1.386 / r2)
		#rp = 
		#r = np.sqrt(abs(r2))
		#print np.shape(r)
		#r = r[0,0]
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



###############################

if __name__ == "__main__":

	# test 
	A = np.array([[1,2], [2,1.0]])