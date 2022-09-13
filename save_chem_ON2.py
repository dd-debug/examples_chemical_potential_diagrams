# -*- coding: UTF-8 -*-

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotEntry, ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine
from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PhaseDiagram
from myResearch.getOrigStableEntriesList import getEntriesList,getOrigStableEntriesList
from pymatgen.ext.matproj import MPRester
import json
import os
import imageio
import numpy as np
from adjustText import adjust_text
import matplotlib.pyplot as plt
MPR = MPRester("SZXJWLvi8njBGvA4sT")
current_dir = os.path.join(os.path.dirname(__file__))
def get_ON_entries():
    ## Many queries are very large, so this python 
    # method either queries the MP and saves it in the 'cache' file, 
    # or if the cache file exists, it loads it directly from the cache. 
    
    cache = os.path.join(current_dir, 'ternary_N_O_without_nitrates_stable')
    if os.path.exists(cache):
        print("Loading from cache.")
        with open(cache, 'r') as f:
            return json.load(f)
    else:
        print("Reading from db.")
#         from pymatgen.ext.matproj import MPRester
        
        criteria = {'elements':{'$all': ["O","N"]}, 'nelements':3, 'e_above_hull':{'$lte':0.00}}
        # The criteria uses mongodb query language. See here for more details: https://docs.mongodb.com/manual/reference/operator/query/
                
        props = ["material_id",'pretty_formula','e_above_hull','structure',"warnings","formation_energy_per_atom"]
        #The properties and the criteria use MaterialsProject features 
        #You can see what is queryable from the MP API documentation: https://github.com/materialsproject/mapidoc/tree/master/materials 
        
        entries = MPR.query(criteria=criteria, properties=props)
        print(len(entries))
        #Save files are prepared in a 'JSON' file. 
        #Some MP objects are not JSONable, and so they must be turned into a dictionary before they can be saved. 
        new_entries=[]
        for e in entries:
            X=e
            X['structure']=X['structure'].as_dict()
            new_entries.append(X)
        
            
        with open(cache, 'w') as f:
            json.dump(new_entries, f)
        return entries
    

ternary_O_N = get_ON_entries()



projEle = ["O","N"]
newformeN = 2
newformeO = 2
fig = None

print(len(ternary_O_N))

formulas = [ee["pretty_formula"] for ee in ternary_O_N]
formulas = list(dict.fromkeys(formulas))
'''set colordict'''
jet= plt.get_cmap('gist_rainbow')
n=len(formulas)
color=iter(jet(np.linspace(0,1,n)))
colordict = {}
for e in formulas:
    c = next(color)
    colordict[e] = c
print(len(formulas))

import random

# random.shuffle(formulas)
m=0
texts = []
filenames = []
aalist=[]

for ee in formulas:
    phi_dict = {}
    m+=1
    if ee != "TaNO":
        continue

    elsList = [str(el) for el in Composition(ee).elements]
    for el in elsList:
        if el != "O" and el != "N":
            tael = el
    elsList = [tael, "O","N"]
    PDentries = getEntriesList(elsList)

    for e in PhaseDiagram(PDentries).stable_entries:
        if e.name == "N2":
            entry = PhaseDiagram(PDentries).make_entry_from_formEperatom(e.composition, newformeN)
        if e.name == "O2":
            entryO = PhaseDiagram(PDentries).make_entry_from_formEperatom(e.composition, newformeO)
    limits = [[-10,0],[-10,newformeO],[-10,newformeN]]
    CPentries = trans_PD_to_ChemPot_entries(PDentries,elsList)
    
    newentries = []
    for e in CPentries:
        if e.name != "N2" and e.name != "O2":
            newentries.append(e)
    newentries.append(ChemPotEntry(entry,newformeN,elsList))
    newentries.append(ChemPotEntry(entryO,newformeO,elsList))
    eql = EquilLine(newentries,elsList,limits = limits)

    for entry in eql._stable_domain_vertices:
        if entry.name == ee:
            center = np.average(eql._stable_domain_vertices[entry], axis=0)
    for entry in eql._stable_domain_vertices:
        phi = entry.form_E
#         print(entry.name,phi)
        for ind in range(len(elsList)):
            phi -= center[ind]*entry.ncomp[elsList[ind]]
#             print(center[ind],entry.ncomp[elsList[ind]])
#         print(phi)
        phi_dict[entry]=phi
    print()
    print(m, ee)
    phi_dict = {k: v for k, v in sorted(phi_dict.items(), key=lambda item: item[1])}
    for e in phi_dict:
        print(e.name, phi_dict[e])



























