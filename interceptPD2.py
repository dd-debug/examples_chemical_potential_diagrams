'''
Created on 2020.06.03

@author: dd
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
from pymatgen.core.composition import Composition
import mpl_toolkits.mplot3d as a3
import matplotlib.pyplot as plt
import numpy as np
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
# els = ["Mg","Ge","N"]
els = ["Na","Fe","N"]
els = ["Fe","O","Al"]
elsE = [Element(i) for i in els]
entries = getOrigStableEntriesList(els)
unshowList = ["AlFe","Al12Fe7","Al6Fe","AlFe3"]
entries.remove(entries[7])

for e in entries:
    print(e.composition)
    if e.name in unshowList:
        print("hahahahaha",e.name)
        entries.remove(e)
#     if e.composition.get_atomic_fraction("Fe") == 0.75:
#         print("gaga",e.name)
#         entries.remove(e)
for e in entries:
    print(e.name)
pd = PhaseDiagram(entries,elsE)
pdp = PDPlotter(pd)
pdp.show()
# haha
# cppdata = pdp.tern_facets_chem_pot_data()
# print(cppdata)
# # CPentries = trans_PD_to_ChemPot_entries(entries,els)
# # cp = ChemPotDiagram(CPentries,els)
# # plt = ChemPotPlotter(cp).get_chempot_plot()
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# cppdata=np.array(cppdata)
# # fig = plt.figure(figsize=(9.2, 7))
# ax = a3.Axes3D(fig) 
# for cpdata in cppdata:
#     ax.scatter(cpdata[0],cpdata[1],cpdata[2])
#     
# 
# 
# limits = []
# for i in range(len(els)):
#     limits += [[-4,0]]
# ax.dist=10
# ax.azim=30
# ax.elev=10
# ax.set_xlim(limits[0])
# ax.set_ylim(limits[1])
# ax.set_zlim(limits[2])
# ax.set_xlabel('Chem Pot '+els[0],fontname='Consolas',fontsize = 12)
# ax.set_ylabel('Chem Pot '+els[1],fontname='Consolas',fontsize = 12)
# ax.set_zlabel('Chem Pot '+els[2],fontname='Consolas',fontsize = 12)
# plt.show()
# 
# 
# 

