import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pymatgen.core.composition import Composition
from pymatgen.ext.matproj import MPRester
MPR = MPRester("2d5wyVmhDCpPMAkq")
entry = MPR.get_entry_by_material_id("mp-1029775")
formE = 0.8564851286899993
# sum(ux) - G = 0
Ba = entry.composition.get_atomic_fraction("Ba")
Mn = entry.composition.get_atomic_fraction("Mn")
N = entry.composition.get_atomic_fraction("N")
a,b,c,d = Ba,Mn,N,formE
#BaMnN2
x = np.linspace(-10,0,10)
y = np.linspace(-10,0,10)

X,Y = np.meshgrid(x,y)
Z = (-d - a*X - b*Y) / c

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlim([-10,0])
ax.set_ylim([-10,0])
ax.set_zlim([-10,0])
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

surf = ax.plot_surface(X, Y, Z)
plt.show()