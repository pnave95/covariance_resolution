# test for python angle quadrant detection (assumes x,y != 0)
import numpy as np

x = -0.2
y = 1.2

quadrants = set([1,2,3,4])
final_quadrant = 0

if np.sign(x) > 0:
	quadrants.discard(2)
	quadrants.discard(3)
else:
	quadrants.discard(1)
	quadrants.discard(4)
if np.sign(y) > 0:
	quadrants.discard(3)
	quadrants.discard(4)
else:
	quadrants.discard(1)
	quadrants.discard(2)
if 1 in quadrants:
	final_quadrant = 1
elif 2 in quadrants:
	final_quadrant = 2
elif 3 in quadrants:
	final_quadrant = 3
elif 4 in quadrants:
	final_quadrant = 4


if final_quadrant == 0:
	print "Error"
else:
	print final_quadrant
