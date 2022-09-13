
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
elsList = ["Ba","Mn","N"]
from pymatgen.ext.matproj import MPRester
MPR = MPRester("2d5wyVmhDCpPMAkq")
PDentries=MPR.get_entries_in_chemsys(elsList)
# PDentries = getOrigStableEntriesList(elsList)
limits = [[-4,0],[-3,0],[-3,0]]
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
cp = ChemPotDiagram(CPentries,elsList,limits = limits)
ChemPotPlotter(cp).get_equil_line_on_CP_Halfspace(limits = limits,
    alpha = 0.4,show_polytope = True)
