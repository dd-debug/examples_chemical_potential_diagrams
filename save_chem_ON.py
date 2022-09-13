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
def get_ON_entries():
    ## Many queries are very large, so this python 
    # method either queries the MP and saves it in the 'cache' file, 
    # or if the cache file exists, it loads it directly from the cache. 
    
    cache = os.path.join(current_dir, 'ternary_N_O_without_nitrates')
    if os.path.exists(cache):
        print("Loading from cache.")
        with open(cache, 'r') as f:
            return json.load(f)
    else:
        print("Reading from db.")
#         from pymatgen.ext.matproj import MPRester
        
        criteria = {'elements':{'$all': ["O","N"]}, 'nelements':3, 'e_above_hull':{'$lte':0.05}}
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
projEle = ['O',"N"]
getridlist = ["N2","O2","N2O","NO2","N2O5","NO","N2O3","NO3","NO4","NO6","N4O9"]
for ee in formulas:
    m+=1
#     if ee != "Ge2N2O" and ee != "Zr2N2O" and ee != "TaNO" and ee != "VN3O10" and ee != "Bi5NO10":
#         continue

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
        if e.name not in getridlist:
            newentries.append(e)
    newentries.append(ChemPotEntry(entry,newformeN,elsList))
    newentries.append(ChemPotEntry(entryO,newformeO,elsList))
    cp = ChemPotDiagram(newentries,elsList,limits = limits)
    plotter = ChemPotPlotter(cp)
    try:
        fig = plotter.get_equil_line_on_CP_Halfspace(limits = limits,
            alpha = 0.2,show_polytope = False)
        fn = str(m) + "_" + ee+".png"
        filenames.append(fn)
        print(fn)
        fig.savefig(current_dir+"/figs_chempot_nogas/" + fn )
    except:
        print(m,ee,"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa")
        aalist.append(ee)
print(aalist)
#     fig = plotter.get_ON_ternary_projection_equilLine_plot(projEle, texts, \
#         fig = fig, colordict = colordict, compound = ee, limits = limits,show_polytope = False,
#         alpha = 0.1,label_domains = True)
#     fn = str(m) + "_" + ee+".png"
#     filenames.append(fn)
# plt.show()
#     fig.savefig(current_dir+"/figs_no_nitrates/" + fn )
# # plt.scatter(0,0)
# with open("filenames","w") as f:
#     json.dump(filenames,f)
# writer = []
# for filename in filenames:
#     image = imageio.imread(current_dir+"/figs_no_nitrates/" +filename)
#     writer.append(image)
# exportname = "figs_no_nitrates.gif"
# kargs = { 'duration': 0.1 }
# imageio.mimsave(exportname, writer, 'GIF', **kargs)


























