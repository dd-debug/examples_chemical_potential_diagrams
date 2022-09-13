# -*- coding: UTF-8 -*-
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram,PDPlotter
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
MPR = MPRester("2d5wyVmhDCpPMAkq")
import numpy as np
import mpl_toolkits.mplot3d as a3
# elsList = ['Mn','O'] 
# PDentries = getOrigStableEntriesList(elsList)
# # PDentries = MPR.get_entries_in_chemsys(elsList)
# CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
# # limits = []
# # for i in range(len(elsList)):
# #     limits += [[-10,0]]
# cp = ChemPotDiagram(CPentries,elsList)
# ChemPotPlotter(cp).show()
# # ChemPotPlotter(cp).show(limits = limits)
 
elsList = ['Mn','O'] 
PDentries = getOrigStableEntriesList(elsList)
pd = PhaseDiagram(PDentries)
D = {}
arraylist = []
xx =[]
yy = []
zz= []
PDentries = sorted(PDentries, key = lambda e: e.composition.get_atomic_fraction('Mn'))
for entry in PDentries:
    alist = []
    a = pd.get_form_energy_per_atom(entry)
    zz.append(a)
    for el in elsList:
        frac = entry.composition.get_atomic_fraction(el)
        alist.append(frac)
        if el == elsList[0]:
            xx.append(frac)
        if el == elsList[1]:
            yy.append(frac)
    alist.append(a)
    D[entry.name] = alist
    arraylist.append(alist)
print(D)
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
for e in D:
    ax.plot(xx,yy,zz)
    
    print(xx)
    ax.scatter(D[e][0],D[e][1],D[e][2], marker = "o", color = "r")
    ax.text(D[e][0],D[e][1],D[e][2]-0.1, e,ha='center',
            va='center')
ax.set_xlabel('Mn Label')
ax.set_ylabel('O Label')
ax.set_zlabel('G Label')

plt.show()
# jet= plt.get_cmap('gist_rainbow')
# n = len(PDentries)
# color=iter(jet(np.linspace(0,1,n)))
# fig = plt.figure(figsize=(9.2, 7))
# ax = a3.Axes3D(fig) 
# 
# for e in D:
#     
#     hull = ConvexHull(np.array(arraylist))
#     print(hull)
#     simplices = hull.simplices  
#     org_triangles = [CP_domain_vertices[e][s] for s in simplices]
#     c=next(color)
#     pc = a3.art3d.Poly3DCollection(org_triangles, alpha = alpha, facecolor=c,edgecolor = edc)
#     ax.add_collection3d(pc)
#     center = np.average(CP_domain_vertices[e], axis=0)
# #                 print(org_triangles)
# #                 print(e.name,"vertices",CP_domain_vertices[e])
#     print("center",center)
#     text_font = {'fontname':'Consolas', 'size':'15', 'color':c, 'weight':'normal'}

# ax.dist=10
# ax.azim=30
# ax.elev=10
# ax.set_xlim(limits[0])
# ax.set_ylim(limits[1])
# ax.set_zlim(limits[2])
# ax.set_xlabel('Chem Pot '+self.elementList[0],fontname='Consolas',fontsize = 12)
# ax.set_ylabel('Chem Pot '+self.elementList[1],fontname='Consolas',fontsize = 12)
# ax.set_zlabel('Chem Pot '+self.el