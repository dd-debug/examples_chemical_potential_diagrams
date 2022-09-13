# -*- coding: UTF-8 -*-
import mpl_toolkits.mplot3d as a3
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from scipy.spatial import ConvexHull
from chemicalDiagram.EquilibriumLine import EquilLine
from pymatgen.core.composition import Composition
compound = "Ba4NaAl2B8(ClO6)3"
# compound = "Ba2CaV2CoF14"
# compound = "Ba2GaAsSe5"
# compound = "BaMnN2"
elsList = [str(i) for i in Composition(compound).elements]
# elsList = ["Ba","Mn","N"]
# elsList = ["Al","Si","O","Cl","Li"]
# choose a element and its chempot you want to fix
fixedEle = [('Ba',-6.6),("V",-4)]
fixedEle = [('Ba',-5),("Na",-5),("Al",-5)]
PDentries = getOrigStableEntriesList(elsList)

CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
limits = []
for e in range(6):
    limits += [[-10,0]]
print(limits)

cp = ChemPotDiagram(CPentries,elsList,fixed = fixedEle,limits = limits)

ChemPotPlotter(cp,fixed = fixedEle).get_slice_equilLine_plot(
    limits=limits,show_polytope=True)

