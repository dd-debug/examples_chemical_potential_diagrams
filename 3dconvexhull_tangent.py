'''
Created on Jun 5, 2021

@author: jiadongc
'''


import os 
import json
from pymatgen.ext.matproj import MPRester
# MPR = MPRester("SZXJWLvi8njBGvA4sT")
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.periodic_table import Element
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
'''import new_PDPlotter'''
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter


elsList = ["Ge","Zr","N"]
elementlist = [Element(el) for el in elsList]
PDentries = getOrigStableEntriesList(elsList)
pd = PhaseDiagram(PDentries,elementlist)

new_PDPlotter(pd).show()
# new_PDPlotter(pd).add_line()

