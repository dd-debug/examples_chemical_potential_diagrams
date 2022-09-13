'''
Created on 2020.06.03

@author: dd
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
import matplotlib.pyplot as plt
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
els = ["O","Mn"]
elsE = [Element(i) for i in els]
entries = getOrigStableEntriesList(els)
entries = sorted(entries,key = lambda entry: entry.composition.get_atomic_fraction(els[1]))
for e in entries:
    print(e.name)
pd = PhaseDiagram(entries,elsE)
D = {}
for e in entries:
    D[e] = pd.get_form_energy_per_atom(e)
ax = plt.gca()
xx = []
yy = []
name = []
for e in entries:
    formE = D[e]
    ind = entries.index(e)
    if ind != len(entries)-1:
        formEnext = D[entries[ind+1]]
        slope = (formEnext-formE)/(entries[ind+1].composition.get_atomic_fraction(els[1])-e.composition.get_atomic_fraction(els[1]))
        intercept = D[e]-slope*e.composition.get_atomic_fraction(els[1])
        name.append(e.name)
        xx.append(intercept)
        yy.append(intercept+slope)
        ax.text(intercept, intercept+slope,e.name)
ax.plot(xx, yy)

ax.set_xlim([-10,0])
ax.set_ylim([-10,0])
plt.xlabel("mu"+els[0])
plt.ylabel("mu"+els[1])
plt.show()
# PDPlotter(pd).show()