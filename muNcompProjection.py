# -*- coding: UTF-8 -*-

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotEntry,ChemPotPlotter,trans_PD_to_ChemPot_entries
from pymatgen.core.composition import Composition
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.analysis.phase_diagram import PhaseDiagram
compound = "Ba4NaAl2B8(ClO6)3"
compound = "Ba2CaV2CoF14"
compound = "NbNO"
elsList = [str(i) for i in Composition(compound).elements]
elsList = ['Ag',"Ge",'O','N']
elsList = ["Ta","O","N"]
# choose a element and its chempot you want to project
projEle = ['O',"N"]
newformeN = 2
newformeO = 2
# elsList = ["Nb","O","N"]
# projEle = ["O","N"]
PDentries = getOrigStableEntriesList(elsList)
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
for e in PhaseDiagram(PDentries).stable_entries:
    if e.name == "N2":
        entry = PhaseDiagram(PDentries).make_entry_from_formEperatom(e.composition, newformeN)
    if e.name == "O2":
        entryO = PhaseDiagram(PDentries).make_entry_from_formEperatom(e.composition, newformeO)
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
limits = [[-10,0],[-10,newformeO],[-10,newformeN]]
newentries = []
for e in CPentries:
    if e.name != "N2" and e.name != "O2" and e.name != "N2O" and e.name != "NO2" and e.name != "N2O5":
        newentries.append(e)
newentries.append(ChemPotEntry(entry,newformeN,elsList))
newentries.append(ChemPotEntry(entryO,newformeO,elsList))
cp = ChemPotDiagram(newentries,elsList,limits = limits)
# limits = []
# for i in range(len(elsList)):
#     limits += [[-8,0]]
# limits = [[-10,0],[-10,0],[-10,0]]
# cp = ChemPotDiagram(CPentries,elsList,limits = limits)

plotter = ChemPotPlotter(cp)
plotter.get_projection_equilLine_plot(projEle, \
    limits = limits,show_polytope = False,
    alpha = 0.2,label_domains = False,label_equilLine=False)
# Cl2 {'Al': 0, 'Si': 0, 'O': 0, 'Cl': 1.0, 'Li': 0}
# hyperplane [0, 0, 0, 1.0, 0, -0.0]
# LiSiNO
# KGeNO
# Sr6Sn2NO