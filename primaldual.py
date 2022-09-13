'''
Created on Nov 15, 2021

@author: jiadongc
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from pymatgen.ext.matproj import MPRester
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
import numpy as np
import matplotlib.pyplot as plt
# from adjustText import adjust_text
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries,\
    ChemPotEntry
from pymatgen.entries.computed_entries import ComputedEntry
from copy import deepcopy


jet= plt.get_cmap('rainbow')
n = 10
colors=['r','#ff9900',"#00ccff","#ff66ff","#ffcc66","#8c8c8c","#0040ff","#00cc44","#bf8040","#ac00e6"]
# colors=['r','orange',"lightblue","pink","yellow","gray","blue","green","brown","purple"]



els = ["Ba",'N']
elsE = [Element(i) for i in els]
entries = getOrigStableEntriesList(els)
entries = sorted(entries,key = lambda e: e.composition.get_atomic_fraction(els[1]))
complist = [e.composition.get_atomic_fraction(els[1]) for e in entries]
formElist = [PhaseDiagram(entries).get_form_energy_per_atom(e)-0.01 for e in entries]

# fig = plt.figure(figsize=(9.2, 7))
# ax = plt.gca()
# fig1 = plt.figure(figsize=(9.2, 7))
# ax1 = plt.gca()
fig, ax = plt.subplots()
fig1,ax1 = plt.subplots()
ax.plot(complist,formElist,linewidth = 1.5,c= 'gray',zorder=1) 
ax.plot([1,0],[-0.75,-0.8],c="#1a1a1a",linewidth = 1.5,zorder=2) #0040ff
for dashx,dashy in zip(complist,formElist):
    ax.plot([dashx,dashx], [dashy,0.05*dashx-0.8], color="lightgray", linestyle='dashed',zorder=1)


with MPRester("SZXJWLvi8njBGvA4sT") as MPR:
    entries = MPR.get_entries_in_chemsys(els)
# pd = PhaseDiagram(entries)
# PDPlotter(pd).show()
print(len(entries))
for i in entries:
    print(i.name)
old_CPentries = trans_PD_to_ChemPot_entries(entries,els)
CPentries =[]
for e in old_CPentries:
    if e.form_E>0:
        if len(e.composition.elements) == 1:
            continue
    if e.name == "BaN2":
        ban2 = deepcopy(e)
    CPentries.append(e)
ban2.form_E=-0.1
CPentries.append(ban2)
label_font = {'fontname':'Calibri', 'size':'20', 'color':'black', 'weight':'normal'}

text_font = {'fontname':'Calibri', 'size':'15', 'weight':'normal'}
# need to sort based on composition
CPentries = sorted(CPentries,key = lambda e: e.entry.composition.get_atomic_fraction(els[1]))
complist = []
formElist = []
namelist = []


for e in CPentries:
    comp = e.entry.composition.get_atomic_fraction(els[1])
    complist.append(comp)
    formElist.append(e.form_E)
    print(e.name,e.form_E)
    namelist.append(e.name)
    print(comp,e.form_E)

colorDict = {}
for e,c in zip(CPentries,colors):
    colorDict[e.name]=c
    mu1s = np.linspace(-10, 1, 1100)
    zz=2
    if e.name == "Ba2N":
        mu1s=np.linspace(-0.5, 0.5, 100)
    if e.name == "Ba3N2":
        mu1s=np.linspace(-1,0.25,200)
    if e.name == "Ba3N":
        mu1s=np.linspace(-0.2,0.2,100)
        zz=10
    if e.name == "N2":
        mu1s=np.linspace(-4.5,1,500)
    if e.name == "BaN6":
        mu1s=np.linspace(-4.2,0.5,500)
    if e.name == "BaN2":
        mu1s=np.linspace(-3,0.75,500)
    if e.name=="BaN":
        mu1s=np.linspace(-1,0.25,300)
    mu2s = []
    print(e.composition.get_atomic_fraction(els[1]))
    fe = e.form_E
    mu2s=(fe-mu1s*e.composition.get_atomic_fraction(els[0]))/e.composition.get_atomic_fraction(els[1])
    print(len(mu1s),len(mu2s))
    ax1.plot(mu1s,mu2s,linewidth = 1.5,c=c,zorder=zz)
ax1.plot([0,0],[-5,1],c="r",linewidth = 1.5,zorder=1)
colorDict["BaN2"]="#0040ff"

namelist = ['N2', 'Ba2N', 'BaN2', 'BaN6', 'Ba3N', 'Ba']
vertices = [[[[-4.5], [-3.569602414279995]], [[0.0], [0.0]]], [[[-0.28486814500000257], [-0.010507995000000214]], [[-1.0047829123799974], [-1.5535032123800003]]], [[[-0.28486814500000257], [-1.6568497474999968]], [[-1.0047829123799974], [-0.3187921111299996]]], [[[-3.569602414279995], [-1.6568497474999968]], [[0.0], [-0.3187921111299996]]], [[[0.0], [-0.010507995000000214]], [[-1.585027197380001], [-1.5535032123800003]]], [[[0.0], [0.0]], [[-10.0], [-1.585027197380001]]]]
stable_entries = getOrigStableEntriesList(els)
for e in stable_entries:
    ind = namelist.index(e.name)
    ax1.plot(vertices[ind][0],vertices[ind][1],c=colorDict[e.name],linewidth=4,zorder=20,linestyle="dashed")

ax1.scatter(-0.86,-0.71,s=150,c="#333333",zorder=22)
ax1.plot([-0.86,-0.86],[-0.71,1.8],c="gray",linestyle="dashed",linewidth=1.5)

for comp, form_E, color in zip(complist,formElist,colors):
    ax.scatter(comp,form_E,marker = 'o',c=color,s = 250,zorder=10)

# plt.xlabel("Composition "+els[1] + "(%)",**label_font)
# plt.ylabel('Formation energy' +"(eV/atom)",**label_font)
# fig.savefig("-".join(els)+"_convex_hull",dpi = 200)

ax.set_xlim([-0.1,1.1])
ax.set_ylim([-0.9,0.2])
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_xlabel(None)
ax.set_ylabel(None)
ax.tick_params(direction='out',labelsize= 15, length=2, width=2, colors='k',
       grid_color='k', grid_alpha=0.5)
ax.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)

ax1.set_xlim([-5,1.5])
ax1.set_ylim([-3.5,2])
ax1.set_xticklabels([])
ax1.set_yticklabels([])
# ax1.set_xlabel("muBa")
# ax1.set_ylabel("muN")
ax1.tick_params(direction='out',labelsize= 15, length=2, width=2, colors='k',
       grid_color='k', grid_alpha=0.5)
ax1.tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False)
plt.show()




























































