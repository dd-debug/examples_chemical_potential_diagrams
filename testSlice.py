'''
Created on 2020.06.11

@author: dd
'''
import numpy as np
from scipy.spatial import ConvexHull, convex_hull_plot_2d
points = np.array([[0.2, 0.2],
                        [0.4, 0.4],
                        [0.2, 0.4],
                        [0.4, 0.2],
                        [0.3, 0.6]])  # 30 random points in 2-D
hull = ConvexHull(points)
print(points)


import matplotlib.pyplot as plt
plt.plot(points[:,0], points[:,1], 'o')
for simplex in hull.simplices:
    print(simplex)
    print(points[simplex])
    print(points[simplex, 0])
    print(points[simplex][0])
    plt.plot(points[simplex, 0], points[simplex, 1], 'k-')
plt.show()

# We could also have directly used the vertices of the hull,
#  which for 2-D are guaranteed to be in counterclockwise order:

print(hull.vertices)
print(points[hull.vertices])
print(points[hull.vertices][0])
print(points[hull.vertices,0])
plt.plot(points[hull.vertices,0], points[hull.vertices,1], 'r--', lw=2)
plt.plot(points[hull.vertices[0],0], points[hull.vertices[0],1], 'ro')


print(hull.equations)
plt.show()