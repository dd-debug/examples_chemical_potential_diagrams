# -*- coding: UTF-8 -*-

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from pymatgen.core.composition import Composition
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList

compound = "BaMnN2"
elsList = [str(i) for i in Composition(compound).elements]

projEle = ['Ba','N']

PDentries = getOrigStableEntriesList(elsList)
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)

limits = []
for i in range(len(elsList)):
    limits += [[-20,0]]
limits = [[-4,0],[-3,0],[-3,0]]
cp = ChemPotDiagram(CPentries,elsList,limits = limits)

plotter = ChemPotPlotter(cp)
plotter.get_projection_equilLine_mu_x_plot(projEle, \
    limits = limits,show_polytope = False,
    alpha = 0.3,label_equilLine = True)
# [array([[-4.00000000e+00, -2.54679250e-03],
#        [-4.00000000e+00, -2.77200725e-01]]), array([[-4.00000000e+00, -2.54679250e-03],
#        [-9.81853757e-01, -2.54679250e-03]]), array([[-1.80581556, -0.27720073],
#        [-4.        , -0.27720073]]), array([[-1.80581556, -0.27720073],
#        [-0.98185376, -0.00254679]])]