# -*- coding: UTF-8 -*-
from pymatgen.ext.matproj import MPRester
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
import pandas as pd
MPR = MPRester("2d5wyVmhDCpPMAkq")
fomula = 'BaMnN2'
elsList = [str(e) for e in Composition(fomula).elements] 
# elsList = ["La","H","Ni"]
elsList = ["Fe","O","Al"]
# elsList = ['Ni','O'] 
PDentries = getOrigStableEntriesList(elsList)
unshowList = ["AlFe","Al12Fe7","Al6Fe","AlFe3"]
PDentries.remove(PDentries[7])
for e in PDentries:
    print(e.composition)
    if e.name in unshowList:
        print("hahahahaha",e.name)
        PDentries.remove(e)
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
limits = []
for i in range(len(elsList)):
    limits += [[-10,0]]
cp = ChemPotDiagram(CPentries,elsList,limits = limits)
# ChemPotPlotter(cp).show()
ChemPotPlotter(cp).get_equil_line_on_CP_Intercept(limits = limits,PDentries=PDentries)
# ChemPotPlotter(cp).show(limits = limits)

