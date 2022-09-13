'''
Created on Jul 29, 2021

@author: jiadongc
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram,GrandPotentialPhaseDiagram,PDPlotter
from myResearch.A200526_duality.grandpotpdplotter import GrandPotentialPhaseDiagram_E,GPDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.periodic_table import Element
import plotly.graph_objects as go
import json
with open("C:/Users/jiadongc/anaconda3/envs/my_pymatgen/lib/site-packages/pymatgen/util/plotly_pd_layouts.json", "r") as f:
    plotly_layouts = json.load(f)
els = ["Ba","Mn","N"]
el_mu="N"
entries = getOrigStableEntriesList(els)
# PDPlotter(PhaseDiagram(entries)).show()
jj=-5
data = []
annotations_lists = []
for i in range(0,100): 
    gpd = GrandPotentialPhaseDiagram_E(entries,{Element(el_mu):jj})
    jj=jj-0.05
    data += GPDPlotter(gpd).get_plot()[0]
#     annotations_lists += GPDPlotter(gpd).get_plot()[1]

fig = go.Figure(data=data)
# layouts = plotly_layouts["advanced_ternary_layout"].copy()
# layouts["scene"].update({"annotations": annotations_lists})
# fig.layout = layouts
els_copy=[el for el in els if el != el_mu]
fig.update_layout(
    scene = dict(
        xaxis = dict(title=els_copy[1]+"x"+els_copy[0]+"1-x composition"),
                     yaxis = dict(title = "mu_"+el_mu),
                     zaxis = dict(title = "phi=G-mu"+el_mu+"x"+el_mu)),
    margin=dict(r=20, l=10, b=10, t=10))
fig.show()
#do not forget to change get_form_energy back-