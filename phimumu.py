# -*- coding: UTF-8 -*-
from pymatgen.ext.matproj import MPRester
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList

MPR = MPRester("2d5wyVmhDCpPMAkq")
# elsList = ['Ba','Fe','Mn'] #no stable materials except single matters
elsList = ['Mn','O'] 
# PDentries = MPR.get_entries_in_chemsys(elsList)
PDentries = getOrigStableEntriesList(elsList)
limits = [[-8,0],[-8,0]]
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
cp = ChemPotDiagram(CPentries,elsList,limits = limits)
plotter = ChemPotPlotter(cp)
limits = [[-8,0],[-8,0],[-10,8]]
# plotter.show()
plotter.get_phimumu_plot(alpha = 0.3,limits = limits,
                         show_label=False,label_equiline=False)
# [-4.05734142 -3.94265858  1.92447197]
