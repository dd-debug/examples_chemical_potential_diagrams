from pymatgen.analysis.phase_diagram import PhaseDiagram,PDPlotter
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList
entries = getOrigStableEntriesList(["Ba","Mn","N"])
PDPlotter(PhaseDiagram(entries)).show()