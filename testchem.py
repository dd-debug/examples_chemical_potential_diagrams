from pymatgen.ext.matproj import MPRester
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from chemicalDiagram.EquilibriumLine import EquilLine,EquilLinePlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
from pymatgen.core.composition import Composition
els = ["Fe","O","Al"]
# compound = "Ba2GaAsSe5"
# els = [str(i) for i in Composition(compound).elements]
PDentries = getOrigStableEntriesList(els)
fixedEle = [('Ba',-5)]
unshowList = ["AlFe","Al12Fe7","Al6Fe","AlFe3"]
PDentries.remove(PDentries[7])
for e in PDentries:
    print(e.composition)
    if e.name in unshowList:
        print("hahahahaha",e.name)
        PDentries.remove(e)
# fixedEle = [('Fe',-1)]
CPentries = trans_PD_to_ChemPot_entries(PDentries,els)
cp = EquilLine(CPentries,els)#, fixed = fixedEle)
EquilLinePlotter(cp).get_equilLine_plot()

