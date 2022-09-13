'''
Created on 2020.09.11

@author: dd
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
import numpy as np
import matplotlib.pyplot as plt
# from adjustText import adjust_text
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection

label_font = {'fontname':'Calibri', 'size':'20', 'color':'black', 'weight':'normal'}

text_font = {'fontname':'Calibri', 'size':'15', 'weight':'normal'}

fig = plt.figure(figsize=(9.2, 7))
ax = plt.gca()
ax2 = ax.twiny()
ax.tick_params(direction='out',labelsize= 15, length=2, width=2, colors='k',
       grid_color='k', grid_alpha=0.5)
ax2.tick_params(direction='out',labelsize= 15, length=2, width=2, colors='k',
       grid_color='k', grid_alpha=0.5)
ax.set_yticks(np.arange(-2.0, 0.05, 0.5))
ax.set_xticklabels([])
ax.set_yticklabels([])
# plt.xlabel("Composition "+els[1] + "(%)",**label_font)
# plt.ylabel('Formation energy' +"(eV/atom)",**label_font)
'''get from get_interfacial function by jiadong'''
molfrac=[1.0, 0.8333333333333335, 0.714285714285714, 0.6363636363636364, 0.5789473684210525, 0.4999999999999997, 0.33333333333333326, 0.0]
molfrac=[1.0, 0.5000000000000001, 0.333333333333333, 0.09523809523809497, 0.0]
# molfrac = [1.0, 0.8823529411764706, 0.7142857142857145, 0.6923076923076924, 0.68, 0.6666666666666669, 0.49999999999999983, 0.0]

molfrac = [round(i,2) for i in molfrac]
molfrac2 = ["" for i in molfrac]
x=[1,0.6716666666666665, 0.5049999999999976, 0.4067647058823529, 0.34983870967741904, 0.2807142857142855, 0.16166666666666663,0]
y=[0,-0.49595025225332884, -0.6321929041899971, -0.6243806900976465, -0.5661817263283752, -0.49442007353713147, -0.32262756900777767,0]
x=[1,0.5050000000000001, 0.3383333333333299, 0.09023809523809512,0]
y=[0,-0.6010791749700026, -0.7033746549600007, -0.2825588115528507,0]
# x=[1,0.7550000000000001, 0.5050000000000002, 0.4786842105263159, 0.4544594594594595, 0.4394444444444445, 0.2807142857142857,0]
# y=[0,-0.5042596677849996, -0.7850543871899972, -0.790769769811575, -0.7790130073097467, -0.7663201143911063, -0.5495044040133353,0]

plt.setp(ax.spines.values(), color='#cc33ff',linewidth=2)
ax.set_ylim([-2.0,0.1])
ax.set_xlim([-0.03,1.02])
ax.plot(x,y,linewidth = 3,c= 'k')


ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(x)
ax2.set_xticklabels(molfrac2)
# ax2.set_xlabel("mol fraction")
# plt.show()
# fig.savefig("-".join(els)+"_convex_hull",dpi = 200)
colors = [[0.8852941176470588, 0.5953026209214163, 0.34162975414777047, 1.0],
[0.49999999999999994, 0.9110226492460883, 0.6911987140551711, 1.0],
[0.36470588235294116, 0.4555113246230441, 0.9201720358189464, 1.0],
[0.6147058823529412, 0.7239446234451342, 0.5599633321173183, 1.0],
[0.6147058823529412, 0.45551132462304417, 0.4201720358189464, 1.0],
[0.7950980392156862, 0.49499010664035353, 0.3276419250067269, 1.0]]
cbamnn2 = [0,1,0.78009699,1]
cbamnn2 = "#84F8DF"
for i in range(len(x)-2):
    cc=colors[i]
    if i == 1:
        cc = cbamnn2
        x1 = x[i+2]
        x2 = x[i+1]
        y1 = y[i+2]
        y2 = y[i+1]
        x0 = x[i]
        y0 = y[i]
        x3 = x[i+3]
        y3 = y[i+3]
        print(y0,y2,y1)
        triangle = [[1,(y2-y1)/(x2-x1) * (1-x2) + y2],[x[i+1],y[i+1]],[1,(y0-y2)/(x0-x2) * (1-x2) + y2]]
        polygon = Polygon(triangle, True)
        print(triangle)
        pc=PatchCollection([polygon],alpha = 0.5, facecolor = cc,edgecolor = cc,linewidths=1)
        ax.add_collection(pc)
        ax.plot([1,1],[(y2-y1)/(x2-x1) * (1-x2) + y2,(y0-y2)/(x0-x2) * (1-x2) + y2],c = cc,linewidth = 3)
        
#         triangle = [[0,(y2-y1)/(x2-x1) * (1-x2) + y2],[x[i+1],y[i+1]],[0,(y0-y2)/(x0-x2) * (1-x2) + y2]]
        triangle = [[0,-(y2-y1)/(x2-x1) * (x2) + y2],[x[i+1],y[i+1]],[0,-(y0-y2)/(x0-x2) * (x2) + y2]]
        polygon = Polygon(triangle, True)
        print(triangle)
        pc=PatchCollection([polygon],alpha = 0.5, facecolor = cc,edgecolor = cc,linewidths=1)
        ax.add_collection(pc)
        ax.plot([0,0],[-(y2-y1)/(x2-x1) * (x2) + y2,-(y0-y2)/(x0-x2) * (x2) + y2],c = cc,linewidth = 3)

for comp, form_E in zip(x,y):
    ax.plot(comp,form_E,c = '#bf00ff',marker = 'o',markersize = 12)

plt.show()




































