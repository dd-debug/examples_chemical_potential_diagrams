
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
# elsList = ["Ba","Mn","N"]
# # elsList = ["Mn","O"]
# PDentries = getOrigStableEntriesList(elsList)
# limits = [[-4,0],[-3,0],[-3,0]]
# # limits = [[-7,0],[-6,0]]
# CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
# cp = ChemPotDiagram(CPentries,elsList,limits = limits)
#  
# ChemPotPlotter(cp).get_equil_line_on_CP_Halfspace(limits = limits,
#     alpha = 0.4,show_polytope = False)
# ChemPotPlotter(cp).get_chempot_plot(limits = limits,alpha = 0.2).show()
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram,GrandPotPDEntry,\
    PDPlotter,PhaseDiagram
from pymatgen.ext.matproj import MPRester
from pymatgen.core.periodic_table import Element
MPR = MPRester("2d5wyVmhDCpPMAkq")
elsList = ["Ba","Mn","N"]
# entries = MPR.get_entries_in_chemsys(elsList)
entries = getOrigStableEntriesList(elsList)
# pd = PhaseDiagram(entries)
gentries = []
for e in entries:
    gentries.append(GrandPotPDEntry(e,{Element("N"):-8}))
# for e in gentries:
#     print(e.energy)
#     print(e.original_entry)
for e in gentries:
    if e.name == "MnN":
        entry = e
gpd = GrandPotentialPhaseDiagram(gentries,{Element("N"):-8})
print(gpd.el_refs)
for e in gpd.el_refs:
    print("buxingle",gpd.el_refs[e].name)
print("MnN",gpd.get_form_energy_per_atom(entry))
for e in gentries:
    if e.name == "Mn":
        print("Mn",e.energy_per_atom)
        print(e.original_entry.energy_per_atom)
    if e.name == "MnN":
        print("e.energy",e.energy)
        print("e.originenergy",e.original_entry.energy)

for e in gentries:
    print(e.name,e.composition)
    print(e.energy,e.original_entry.energy)
# PDPlotter(gpd).show()
# # Dehn-Sommerville