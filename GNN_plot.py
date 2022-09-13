'''
Created on 2021.10.07

@author: jiadongc
'''

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.composition import Composition
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as a3
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
# import mpl_toolkits.mplot3d as a3

def plot_G_N1_N2(elsList, colordict=None):
    PDentries = getOrigStableEntriesList(elsList)
    pd = PhaseDiagram(PDentries)
    D = {}
    PDentries = sorted(PDentries, key = lambda e: e.composition.get_atomic_fraction(elsList[0]))
    for entry in PDentries:
        alist = []
        a = pd.get_form_energy_per_atom(entry)
        for el in elsList:
            print(entry.composition)
            frac = entry.composition.get_atomic_fraction(el)
            alist.append(frac)
        alist.append(a)
        D[entry.name] = alist
    
    print(D)
    fig = plt.figure(figsize=(9.2, 7))
    ax = fig.add_subplot(111, projection='3d')
    
    '''plot convex hull, it on surface N1+N2=1'''
    xx = [D[e][0] for e in D]
    yy = [D[e][1] for e in D]
    zz = [D[e][2] for e in D]
    colorlist = [colordict[e] for e in D]
    
    for i in range(len(xx)-1):
        xxx = np.linspace(xx[i],xx[i+1],9)
        yyy = np.linspace(yy[i],yy[i+1],9)
        zzz = np.linspace(zz[i],zz[i+1],9)
        scolor = [colorlist[i],colorlist[i+1]]
        for j in range(len(xxx)-1):
#             ax.plot(xxx[j:j+2],yyy[j:j+2],zzz[j:j+2],color=scolor[j%2],linewidth = 3)
            ax.plot(xxx[j:j+2],yyy[j:j+2],zzz[j:j+2],color="black",linewidth = 3)

    '''plot surface N1+N2=1'''
    x = np.array([[1, 0], [1, 0]])
    y = np.array([[0, 1], [0, 1]])
    z = np.array([[1, 1], [min(zz)*2.5, min(zz)*2.5]])
    ax.plot_surface(x, y, z, color = "#94b8b8",alpha=0.2)
    

    '''scatter and text each stable material'''
    for e in D:
#         if colordict != None:
#             c = colordict[e]
#         else:
#             c = "r"
        c="green"
        ax.scatter3D(D[e][0],D[e][1],D[e][2], marker = "o",s = 80, color = c)
#         ax.text(D[e][0],D[e][1],D[e][2]+min(zz)*0.15, e,fontsize=10)
    ax.grid(b=None)
#     ax.set_xlabel('x_'+elsList[0],fontname='Arial',fontsize = 12)
#     ax.set_ylabel('x_'+elsList[1],fontname='Arial',fontsize = 12)
#     ax.set_zlabel('G',fontname='Arial',fontsize = 12)
#     ax.set_xticklabels([])
#     ax.set_yticklabels([])
#     ax.set_zticklabels([])
    ax.set_xlim([0,1.2])
    ax.set_ylim([0,1.2])
    ax.set_zlim([min(zz)*2.5,2])
    fig.tight_layout()
    return D,ax
colordict = {'O2': [1.0, 0.0, 0.16, 1.0], 'MnO2': [1.0, 0.918918918918919, 0.0, 1.0], 'Mn2O3': [0.0, 1.0, 0.0, 1.0], 'Mn3O4': [0.0, 0.9239130434782604, 1.0, 1.0], 'MnO': [0.16304347826086973, 0.0, 1.0, 1.0], 'Mn': [1.0, 0.0, 0.75, 1.0]}

elsList = ['Mn','O']
D,ax = plot_G_N1_N2(elsList,colordict)

values = [-4.05734142,-3.94265858,1.92447197]
mu1=values[0]
mu2=values[1]
phi=values[-1]
G=-3
formula = "Mn3O4"
coords = [[0,0,phi],[(phi-G)/(-mu1),0,G],[0,(phi-G)/(-mu2),G]]
pc = a3.art3d.Poly3DCollection(coords, alpha = 0.2, 
                               facecolor=colordict[formula],edgecolor = None)
ax.add_collection3d(pc)
ax.plot([1,0],[0,1],[phi+mu1,phi+mu2],color=colordict[formula],linewidth = 2)
ax.scatter3D(1,0,phi+mu1,marker = "o",s = 60, color = colordict[formula])
ax.scatter3D(0,1,phi+mu2,marker = "o",s = 60, color = colordict[formula])
# 
# ax.plot([1,0],[0,1],[mu1,mu2],color="r",linewidth = 2)
# ax.scatter3D(1,0,mu1,marker = "o",s = 60, color = "r")
# ax.scatter3D(0,1,mu2,marker = "o",s = 60, color = "r")
# 
# ax.plot([Composition(formula).get_atomic_fraction("Mn"),Composition(formula).get_atomic_fraction("Mn")],
#         [Composition(formula).get_atomic_fraction("O"),Composition(formula).get_atomic_fraction("O")],
#         [D[formula][2],D[formula][2]-abs(phi)],
#         color="r",linewidth = 2,linestyle='dashed')
# ax.scatter3D(0,0,phi,marker = "o",s = 60, color = "k")
plt.show()

























