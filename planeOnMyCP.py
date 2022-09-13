# -*- coding: UTF-8 -*-
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram,PDPlotter
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
MPR = MPRester("2d5wyVmhDCpPMAkq")
import numpy as np
import mpl_toolkits.mplot3d as a3
elsList = ['Ba','Mn','N'] 
PDentries = getOrigStableEntriesList(elsList)
# PDentries = MPR.get_entries_in_chemsys(elsList)
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
# limits = []
# for i in range(len(elsList)):
#     limits += [[-10,0]]
cp = ChemPotDiagram(CPentries,elsList)
ChemPotPlotter(cp).show()
# ChemPotPlotter(cp).show(limits = limits)