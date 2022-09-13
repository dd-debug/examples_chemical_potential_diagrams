'''
Created on 09.03.2021 laptop

@author: jiadongc
'''
from chemicalDiagram.ChemicalPotentialDiagram import trans_PD_to_ChemPot_entries, ChemPotDiagram, ChemPotPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
import plotly.graph_objects as go
import json

els = ["Li","Co","O"]
el_mu="Li"
entries = getOrigStableEntriesList(els)

limits = [[-6,0],[-10,0],[-10,0]]
# limits = [[-20,0],[-10,0],[-10,0]]
CPentries = trans_PD_to_ChemPot_entries(entries,els)
cp = ChemPotDiagram(CPentries,els,limits = limits)
projEle = [el_mu]
plotter = ChemPotPlotter(cp)
plotter.get_mu_xx_plot(projEle, 
    limits = limits,show_polytope = False,
    alpha = 0.3,label_domains = True)

