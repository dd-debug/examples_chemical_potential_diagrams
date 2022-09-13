from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine
from myResearch.getOrigStableEntriesList import getEntriesList,getOrigStableEntriesList
from pymatgen.ext.matproj import MPRester
import os
import json
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter

elsList = ['Rb',"Ge",'N'] 
PDentries = getOrigStableEntriesList(elsList)
# PDPlotter(PhaseDiagram(PDentries)).show()
limits = [[-3,0],[-2,0],[-2,0]]
# limits = [[-10,0],[-10,0],[-10,0]]
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
cp = ChemPotDiagram(CPentries,elsList,limits = limits)
ChemPotPlotter(cp).get_equil_line_on_CP_Halfspace(limits = limits,
    alpha = 0.3,show_polytope = True, show_label=False,
    label_equilLine = False)
# ['N2', 'Ba2N', 'BaN2', 'BaN6', 'Ba3N', 'Ba']
# [[[[-10.0], [-3.569602414279995]], [[0.0], [0.0]]], [[[-0.28486814500000257], [-0.010507995000000214]], [[-1.0047829123799974], [-1.5535032123800003]]], [[[-0.28486814500000257], [-1.6568497474999968]], [[-1.0047829123799974], [-0.3187921111299996]]], [[[-3.569602414279995], [-1.6568497474999968]], [[0.0], [-0.3187921111299996]]], [[[0.0], [-0.010507995000000214]], [[-1.585027197380001], [-1.5535032123800003]]], [[[0.0], [0.0]], [[-10.0], [-1.585027197380001]]]]
