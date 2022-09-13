'''
Created on Sep 29, 2021

@author: jiadongc
'''
'''I am going to change chemical potential diagram from mu1,mu2 axis 
to mu1-mu2,-mu2 axis'''

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine
from myResearch.getOrigStableEntriesList import getEntriesList,getOrigStableEntriesList
from pymatgen.ext.matproj import MPRester
import os
import numpy as np
import json
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.A200526_duality.convexhullpdplotter import new_PDPlotter

elsList = ["Ba","N"]
PDentries = getOrigStableEntriesList(elsList)
# PDPlotter(PhaseDiagram(PDentries)).show()
limits = [[-4,0],[-4,0],[-4,0]]
limits = [[-4,0],[-4,0]]
CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
cp = ChemPotDiagram(CPentries,elsList,limits = limits)

ww=[]
for i in range(len(elsList)):
    aa=[0 for j in range(len(elsList))]
    aa[i]=1
    aa[-1]=-1
    ww.append(aa)
print(ww)

for e in cp._stable_domain_vertices:
    X_scaled = cp._stable_domain_vertices[e]
    vec_list = []
    for i in range(len(X_scaled)):
        vec_new = (np.array(X_scaled[i])).dot(np.linalg.inv(np.array(ww)))
        vec_list.append(vec_new)
    cp._stable_domain_vertices[e] = np.array(vec_list)

ChemPotPlotter(cp).get_scaled_equil_line_on_CP_Halfspace(ww=ww,limits = limits,
    alpha = 0.2,label_equilLine=False,show_polytope = False, show_label=False)