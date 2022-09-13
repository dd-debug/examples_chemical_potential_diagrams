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
# from adjustText import adjust_text
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
els = ['O',"Mn"]
elsE = [Element(i) for i in els]
entries = getOrigStableEntriesList(els)
# pd = PhaseDiagram(entries)
# PDPlotter(pd).show()
CPentries = trans_PD_to_ChemPot_entries(entries,els)
label_font = {'fontname':'Calibri', 'size':'20', 'color':'black', 'weight':'normal'}

text_font = {'fontname':'Calibri', 'size':'15', 'weight':'normal'}
# need to sort based on composition
CPentries = sorted(CPentries,key = lambda e: e.entry.composition.get_atomic_fraction(els[1]))

fig = plt.figure(figsize=(9.2, 7))
ax = plt.gca()
ax.tick_params(direction='out',labelsize= 15, length=2, width=2, colors='k',
       grid_color='k', grid_alpha=0.5)
compdict = {}
formEdict = {}
# ax.set_yticks(np.arange(-3.0, 0.25, 0.5))
for e in CPentries:
    comp = e.entry.composition.get_atomic_fraction(els[1])
    compdict[e.name]=comp
    formEdict[e.name]=e.form_E
#     texts = ax.text(comp,e.form_E,e.name)
    print(comp,e.form_E)

# texts = [ax.text(e.entry.composition.get_atomic_fraction(els[1]),e.form_E,e.name,**text_font) for e in CPentries]
# adjust_text(texts)
ax.plot(list(compdict.values()),list(formEdict.values()),linewidth = 3,c= 'k')

ax.set_ylim([-3.25,0.2])
colordict = {'O2': [0.5, 0.0, 1.0, 1.0], 'MnO2': [0.09999999999999998, 0.5877852522924731, 0.9510565162951535, 1.0], 'Mn2O3': [0.30000000000000004, 0.9510565162951535, 0.8090169943749475, 1.0], 'Mn3O4': [0.7, 0.9510565162951536, 0.5877852522924731, 1.0], 'MnO': [1.0, 0.5877852522924732, 0.30901699437494745, 1.0], 'Mn': [1.0, 1.2246467991473532e-16, 6.123233995736766e-17, 1.0]}

for e in compdict:
    ax.plot(compdict[e],formEdict[e],c = colordict[e],marker = 'o',markersize = 15,alpha = 0.9)
plt.xlabel("Composition "+els[1] + "(%)",**label_font)
plt.ylabel('Formation energy' +" (eV/atom)",**label_font)
# fig.savefig("-".join(els)+"_convex_hull",dpi = 200)
plt.show()