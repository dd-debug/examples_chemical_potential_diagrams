'''
Created on 2021.09.03

@author: dd
'''
import numpy as np

import matplotlib.pyplot as plt
import matplotlib
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

fig, ax = plt.subplots()
patches = []
num_polygons = 1
num_sides = 5

for i in range(num_polygons):
    aa = np.random.rand(num_sides ,2)
    print(aa)
    polygon = Polygon(aa, True)
    patches.append(polygon)

p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)

colors = 100*np.random.rand(len(patches))
p.set_array(np.array(colors))

ax.add_collection(p)
plt.show()