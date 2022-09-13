from chemicalDiagram.ChemicalPotentialDiagram import ChemPotEntry,ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter
from myResearch.getOrigStableEntriesList import getEntriesList,getOrigStableEntriesList
from pymatgen.ext.matproj import MPRester
import os
import json
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

els = ["Ba","Mn","N"]
entries = getOrigStableEntriesList(els)
pd = PhaseDiagram(entries)
for entry in entries:
    if entry.name =="BaN2":
        print(entry.name, pd.get_form_energy_per_atom(entry))
    if entry.name == "MnN":
        print(entry.name, pd.get_form_energy_per_atom(entry))
new_PDPlotter(pd).show()

# BaN2 (0.3333333333333333, 0.5773502691896257)
# MnN (0.75, 0.4330127018922193)