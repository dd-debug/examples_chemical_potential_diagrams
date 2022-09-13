# -*- coding: UTF-8 -*-

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from pymatgen.core.composition import Composition
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList


compound = "Ba4NaAl2B8(ClO6)3"
compound = "Ba2CaV2CoF14"
compound = "Ba2GaAsSe5"
elsList = [str(i) for i in Composition(compound).elements]
# elsList = ['Ba','Mn','N']
# choose a element and its chempot you want to project
# projEle = ['Ba','Mn']
projEle = ['Ba',"Ga"]
PDentries = getOrigStableEntriesList(elsList)
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)

# limits = []
# for i in range(len(elsList)):
#     limits += [[-20,0]]
limits = [[-14,0],[-13,0],[-12,0],[-15,0]]
cp = ChemPotDiagram(CPentries,elsList,limits = limits)
plotter = ChemPotPlotter(cp)
plotter.get_projection_equilLine_plot(projEle, \
    limits = limits,show_polytope = False,
    alpha = 0.1,label_domains = True)
