'''
Created on 2020.06.03

@author: dd
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
import numpy as np
import matplotlib.pyplot as plt
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
els = ["O","Mn"]
elsE = [Element(i) for i in els]
entries = getOrigStableEntriesList(els)
CPentries = trans_PD_to_ChemPot_entries(entries,els)
for a in CPentries:
    if a.name == "Mn3O4":
        print("Mn3O4",a.form_E,a.entry.composition.get_atomic_fraction(els[1]))
# need to sort based on composition
CPentries = sorted(CPentries,key = lambda e: e.entry.composition.get_atomic_fraction(els[1]))
mu0 = []
mu1 = [-7]
complist = []
complist2 = []
namelist = []
for e in CPentries:
    formE = e.form_E
    ind = CPentries.index(e)
    if ind != len(CPentries)-1:
        formEnext = CPentries[ind+1].form_E
        compnext = CPentries[ind+1].entry.composition.get_atomic_fraction(els[1])
        comp = e.entry.composition.get_atomic_fraction(els[1])
        slope = (formEnext-formE)/(compnext-comp)
        intercept = formE-slope*e.entry.composition.get_atomic_fraction(els[1])
        mu0.append(intercept)
        mu0.append(intercept)
        mu1.append(intercept+slope)
        mu1.append(intercept+slope)
        complist.append(comp)
        complist.append(compnext)
        complist2.append(comp)
        complist2.append(comp)
    namelist.append(e.name)
mu0.append(-7)
complist.append(1)
complist2.append(1)
print(np.array(mu0))
print(mu1)
# gaga
print(np.array(complist))
print(namelist)
# mu0comp = np.column_stack((complist,mu0))
label_font = {'fontname':'Calibri', 'size':'20', 'color':'black', 'weight':'normal'}

text_font = {'fontname':'Calibri', 'size':'15', 'weight':'normal'}
fig = plt.figure(figsize=(9.2, 7))
ax = plt.gca()
# index = 0
# for e in mu0comp[]
ax.tick_params(direction='out',labelsize= 15, length=2, width=2, colors='k',
       grid_color='k', grid_alpha=0.5)
ax.plot(complist,mu0,label = "mu"+els[0],linewidth = 3)
print(len(complist2),len(mu1))
ax.plot(complist2,mu1,label = "mu"+els[1],linewidth = 3)
for e in complist:
    inx = complist.index(e)
    if inx != len(complist)-1:
        if e == complist[inx+1]:
            print()
#             ax.text(e,(mu0[inx]+mu0[inx+1])/2,namelist[int((inx+1)/2)],text_font)
# for e in complist2:
#     inx = complist2.index(e)
#     if inx != len(complist2)-1:
#         if e == complist2[inx+1]:
#             ax.text(e,(mu1[inx]+mu1[inx+1])/2,namelist[int((inx+1)/2)])
# ax.text(0,0,namelist[0],text_font)
# ax.text(1,0,namelist[-1],text_font)
plt.xlabel("Composition "+els[1] + "(%)",label_font)
plt.ylabel('Chemical Potential ' +"(eV/atom)",label_font)
ax.legend()
fig.savefig("-".join(els)+"_mucomp",dpi = 200)
plt.show()