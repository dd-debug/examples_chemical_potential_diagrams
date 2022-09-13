from chemicalDiagram.ChemicalPotentialDiagram import ChemPotEntry,ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter
from myResearch.getOrigStableEntriesList import getEntriesList,getOrigStableEntriesList
from pymatgen.ext.matproj import MPRester
import os
import json
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

newformeN = 2
newformeO = 2
elsList = ["Cu","N"]
PDentries = getEntriesList(elsList)

stableentries = []
for e in PhaseDiagram(PDentries).stable_entries:
    if e.name == "N2":
        entry = PhaseDiagram(PDentries).make_entry_from_formEperatom(e.composition, newformeN)
    stableentries.append(e)

for e in PDentries:
    if e.name == "Cu3N":
        stableentries.append(e)
limits = [[-2,1],[-2,newformeN]]
# limits = [[-2,0],[-1,1.25]]
CPentries = trans_PD_to_ChemPot_entries(stableentries,elsList)
newentries = []
for e in CPentries:
    if e.name != "N2" and e.name != "O2":
        newentries.append(e)

newentries.append(ChemPotEntry(entry,newformeN,elsList))


cp = ChemPotDiagram(newentries,elsList,limits = limits)
# ChemPotPlotter(cp).get_phimumu_plot(alpha = 0.3,limits = limits,
#                          show_label=True,label_equiline=True)
ChemPotPlotter(cp).get_equil_line_on_CP_Halfspace(limits = limits,
    alpha = 0.3,show_polytope = True,label_domain = False, show_label = False)

'''plot the changed convex hull'''
# newpdes = []
# for e in PDentries:
#     if e.name != "N2" and e.name != "O2":
#         newpdes.append(e)
#         
# newpdes.append(entry)
# newpdes.append(entryO)
# new_PDPlotter(PhaseDiagram(newpdes)).show()