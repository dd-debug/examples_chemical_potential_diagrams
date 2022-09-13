'''
Created on 2020.06.03

@author: dd
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
els = ["Mn","O"]
elsE = [Element(i) for i in els]
entries = getOrigStableEntriesList(els)

CPentries = trans_PD_to_ChemPot_entries(entries,els)

# need to sort based on composition
CPentries = sorted(CPentries,key = lambda e: e.entry.composition.get_atomic_fraction(els[1]))
limits = [[-7,0],[-6,0]]
cp = ChemPotDiagram(CPentries,els,limits = limits) #,limits= [[-5,0],[-5,0]])
cpPlotter = ChemPotPlotter(cp)

# show cp diagram and equi line
cpPlotter.get_equil_line_on_CP_Intercept(filename = "-".join(els),limits = limits) #(limits= [[-5,0],[-5,0]])

