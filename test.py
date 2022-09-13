'''
Created on 2020.09.15

@author: dd
'''
import numpy as np
fd = np.array([[0,2,0],[1,0,0],[0,0,30]])

p1 = fd[0]
p2 = fd[1]
p3 = fd[2]
# These two vectors are in the plane
v1 = p3 - p1
v2 = p2 - p1
# the cross product is a vector normal to the plane
cp = np.cross(v1, v2)
a, b, c = cp
d = np.dot(cp, p3)
print('The equation is {0}x + {1}y + {2}z = {3}'.format(a, b, c, d))