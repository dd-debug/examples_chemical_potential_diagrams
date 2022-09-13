
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
#from EquilibriumLine import EquilLine
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.ext.matproj import MPRester
import os
import json
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram
from myResearch.A200526_duality import plane
# def getOrigStableEntriesList(els,filename = None):
#     directory = os.path.join('/home/anagha/Documents/UMich Project/Week 1')
#     s_els = els.copy()
#     s_els.sort()
#     if filename == None:
#         filename = '-'.join(s_els)
#     print(filename)
#     cache = os.path.join(directory, filename)
#     if os.path.exists(cache):
#         print('loading from cache.','-'.join(s_els))
#         with open(cache, 'r') as f:
#             dict_entries = json.load(f)
#         list_entries = []
#         for e in dict_entries:
#             list_entries.append(ComputedEntry.from_dict(e))
#         return list_entries
#     else:
#         print('Reading from database.','-'.join(s_els))
#         with MPRester("IWpSi6GcnAqLmtVo") as MPR:
#             entries = MPR.get_entries_in_chemsys(s_els)
#         pd = PhaseDiagram(entries)
#         newentries=[]
#         for e in pd.stable_entries:
# #             print(e.entry_id,':',e.name)
#             newentries.append(e)
#         dict_entries = []
#         for e in newentries:
#             dict_entries.append(e.as_dict())
#         with open(cache,'w') as f:
#             json.dump(dict_entries,f)
#         return newentries
elsList = ["Ba","Mn","N"]
PDentries = getOrigStableEntriesList(elsList)
limits = [[-4,0],[-3,0],[-3,0]]
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
cp = ChemPotDiagram(CPentries,elsList,limits = limits)
ChemPotPlotter(cp).get_equil_line_on_CP_Halfspace(limits = limits,
    alpha = 0.4,show_polytope = False)
