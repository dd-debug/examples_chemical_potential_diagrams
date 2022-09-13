'''
Created on Jun 5, 2021

@author: jiadongc
'''


import os 
import json
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core.periodic_table import Element


elsList = ["Li","Ge","S"]
elementlist = [Element(el) for el in elsList]
PDentries = getOrigStableEntriesList(elsList)
pd = PhaseDiagram(PDentries,elementlist)
new_PDPlotter(pd).show()


