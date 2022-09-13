# -*- coding: UTF-8 -*-

from chemicalDiagram.ChemicalPotentialDiagram import ChemPotEntry,ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
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
def get_ON_entries(n=3,name = "oxynitrides_ternary"):
    ## Many queries are very large, so this python 
    # method either queries the MP and saves it in the 'cache' file, 
    # or if the cache file exists, it loads it directly from the cache. 
    
    cache = os.path.join(current_dir, name)
    if os.path.exists(cache):
        print("Loading from cache.")
        with open(cache, 'r') as f:
            return json.load(f)
    else:
        print("Reading from db.")
#         from pymatgen.ext.matproj import MPRester
        
        criteria = {'elements':{'$all': ["O","N"]}, 'nelements':n, 'e_above_hull':{'$lte':0.00}}
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
    

ternaryon = get_ON_entries()
quanternaryon =get_ON_entries(n=4,name="oxynitrides_quanternary")
t_on = [i['pretty_formula'] for i in ternaryon]
q_on = [i['pretty_formula'] for i in quanternaryon]
for i in ternaryon:
    print(ternaryon.index(i),i['pretty_formula'])
print()
for i in quanternaryon:
    print(quanternaryon.index(i),i['pretty_formula']) 
# print(t_on)
# print(q_on)
 
# with open("text_mined_recipes.json","r") as f:
#     txt_recipes = json.load(f)
# #     for i in txt_recipes:
# #         try:
# #             els = [str(j) for j in Composition(i["targets_string"][0]).elements]
# #             if "O" in els and "N" in els:
# #                 if len(Composition(i["targets_string"][0]).elements)==4:
# #                     print(i["targets_string"][0])
# #                     print(i["reaction_string"])
# #                     print()
# #         except:
# #             a=1
# # #             print(i["targets_string"][0])
#              
#     for i in txt_recipes:
#         if i["targets_string"][0] in t_on or i["targets_string"][0] in q_on:
#             print(i["targets_string"][0])
#             print(i["reaction_string"])
#             print(i["doi"])

























