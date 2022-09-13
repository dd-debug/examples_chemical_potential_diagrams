'''
Created on Jun 5, 2021

@author: jiadongc
'''
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter
import plotly.graph_objects as go
import json
import math
import itertools
import numpy as np
from pymatgen.analysis.reaction_calculator import ComputedReaction
from chemicalDiagram.EquilibriumLine import EquilLine
from chemicalDiagram.ChemicalPotentialDiagram import ChemPotDiagram,ChemPotPlotter,trans_PD_to_ChemPot_entries
from pymatgen.core.composition import Composition
from itertools import combinations
from myResearch.getOrigStableEntriesList import getOrigStableEntriesList

with open("C:/Users/jdche/anaconda3/envs/my_pymatgen/lib/site-packages/pymatgen/util/plotly_pd_layouts.json", "r") as f:
    plotly_layouts = json.load(f)

def triangular_coord(coord):
    """
    Convert a 2D coordinate into a triangle-based coordinate system for a
    prettier phase diagram.

    Args:
        coord: coordinate used in the convex hull computation.

    Returns:
        coordinates in a triangular-based coordinate system.
    """
    unitvec = np.array([[1, 0], [0.5, math.sqrt(3) / 2]])

    result = np.dot(np.array(coord), unitvec)
    return result.transpose()


class new_PDPlotter(PDPlotter):
    def get_plot(
        self,
        label_stable=True,
        label_unstable=True,
        ordering=None,
        energy_colormap=None,
        process_attributes=False,
        plt=None,
        label_uncertainties=False,
    ):
        """
        :param label_stable: Whether to label stable compounds.
        :param label_unstable: Whether to label unstable compounds.
        :param ordering: Ordering of vertices (matplotlib backend only).
        :param energy_colormap: Colormap for coloring energy (matplotlib backend only).
        :param process_attributes: Whether to process the attributes (matplotlib
            backend only).
        :param plt: Existing plt object if plotting multiple phase diagrams (
            matplotlib backend only).
        :param label_uncertainties: Whether to add error bars to the hull (plotly
            backend only). For binaries, this also shades the hull with the
            uncertainty window.
        :return: go.Figure (plotly) or matplotlib.pyplot (matplotlib)
        """
        fig = None

        if self.backend == "plotly":
            data = [self._create_plotly_lines()]


            if self._dim == 3:
                data.append(self._create_plotly_ternary_support_lines())
                data.append(self._create_plotly_ternary_hull())

            stable_labels_plot = self._create_plotly_stable_labels(label_stable)
            stable_marker_plot, unstable_marker_plot = self._create_plotly_markers(
                label_uncertainties,marksize=10
            )

            if self._dim == 2 and label_uncertainties:
                data.append(self._create_plotly_uncertainty_shading(stable_marker_plot))
            
#             data.append(stable_labels_plot)
            data.append(unstable_marker_plot)
            data.append(stable_marker_plot)
            linelist = self.add_line()
            data += linelist
#             data += self.create_tangent_plane("Ta3NO6")
# #             data += self.create_reaction_plane()
#             comp1 = Composition("Ta9N")
#             comp2 = Composition("NO9")
#             data += self.create_reaction_plane3([comp1, comp2])
#             data += self.get_irpd(comp1, comp2)
            
            
            fig = go.Figure(data=data)
            fig.layout = self._create_plotly_figure_layout()

        elif self.backend == "matplotlib":
            if self._dim <= 3:
                fig = self._get_2d_plot(
                    label_stable,
                    label_unstable,
                    ordering,
                    energy_colormap,
                    plt=plt,
                    process_attributes=process_attributes,
                )
            elif self._dim == 4:
                fig = self._get_3d_plot(label_stable)
        fig.write_html("Ge-Zr-N.html")
        return fig
#     def add_line(self,marksize=10,linewidth=8):
#         els1str = [str(i) for i in self._pd.elements]
#         els1 = combinations(els1str,2)
# #         for i in els1:
# #             print(i)
# #         for eee in els1:
# #             print(ads)
# #         gagaS
#         
#         data1 = []
#         vertices = []
#         for els in els1:
#             gaga
    def add_line(self,marksize=10,linewidth=8):
        els1str = [str(i) for i in self._pd.elements]
        els1 = combinations(els1str,2)
        
        data1 = []
        vertices = []
        for els in els1:
            
            entries = getOrigStableEntriesList(els)
            pd = PhaseDiagram(entries)
            deepest = 0
            for e in entries:
                formE = pd.get_form_energy_per_atom(e)
                if deepest >= formE:
                    deepest = formE
                    deepestentry = e
            
            ver=np.append(triangular_coord([deepestentry.composition.get_atomic_fraction(self._pd.elements[1]),
                              deepestentry.composition.get_atomic_fraction(self._pd.elements[2])]),deepest)

            vertices.append(ver.tolist())
            print(vertices)
            
        vertices = np.array(vertices)
        print(len(vertices))

        for i in range(0,len(vertices)):
            if abs(vertices[i][2] + 0.8471310825000005) < 1e-6:
                c="red"
            else:
                c="blue"
            data1.append(go.Scatter3d(
                x=[vertices[i][1],vertices[i][1]], y=[vertices[i][0],vertices[i][0]], z=[vertices[i][2],0],
#                 marker=dict(
#                     size=0,
#                     opacity=0.5,
#                     color="#4d94ff",
#                     colorscale='Viridis',
#                 ),
                line=dict(
                    color=c,
                    width=linewidth
                )))
        return data1
    
    def _create_plotly_figure_layout(self, label_stable=True):
        """
        Creates layout for plotly phase diagram figure and updates with
        figure annotations.

        :return: Dictionary with Plotly figure layout settings.
        """
        annotations_list = None
        layout = dict()

        if label_stable:
            annotations_list = self._create_plotly_element_annotations()

        if self._dim == 2:
            layout = plotly_layouts["default_binary_layout"].copy()
            layout["annotations"] = annotations_list
        elif self._dim == 3:
            layout = plotly_layouts["default_ternary_layout"].copy()
            layout["scene"].update({"annotations": annotations_list})
        elif self._dim == 4:
            layout = plotly_layouts["default_quaternary_layout"].copy()
            layout["scene"].update({"annotations": annotations_list})

        return layout
    
    def create_reaction_plane(self):
        '''xx yy zz are elemental mu at the intercepts (the order is the elements order)'''
        # BaN2 
        aa = (0.3333333333333333, 0.5773502691896257,-0.7648113232533321)
        # MnN 
        bb = (0.75, 0.4330127018922193,-0.4743314211899978)
        para = 1.5
        coords = np.array([[0.3333333333333333, 0.5773502691896257,0],[0.3333333333333333, 0.5773502691896257,self._min_energy*para],
                           [0.75, 0.4330127018922193,0],[0.75, 0.4330127018922193,self._min_energy*para]])
        print(coords[:,0])
        data1 = [go.Scatter3d(
            x=[aa[1],bb[1]], y=[aa[0],bb[0]], z=[aa[2],bb[2]],
                marker=dict(
                size=8,
                color="darkblue",
                colorscale='Viridis',
            ),
            line=dict(
                color='red',
                width=5
            ))]
        
        facets = [0,1,2,3]
        # i, j and k give the vertices of triangles
        # here we represent the 2 triangles of the rectangle
        data1.append(go.Mesh3d(
            x=list(coords[:, 1]),
            y=list(coords[:, 0]),
            z=list(coords[:, 2]),
            i=[facets[1],facets[1]],
            j=[facets[0],facets[2]],
            k=[facets[2],facets[3]],
            opacity=0.8,
            color = "lightblue",
            hoverinfo="none",
            lighting=dict(diffuse=0.0, ambient=1.0),
            name="Convex Hull (shading)",
            flatshading=True,
            showlegend=True,
        ))
        return data1
    def create_reaction_plane2(self):
        '''xx yy zz are elemental mu at the intercepts (the order is the elements order)'''
        # BaN2 
        aa = (0.3333333333333333, 0.5773502691896257,-0.7648113232533321)
        '''plot a series planes with AxB1-x'''
        para = 2
        data1 = []
        w = 5
        for ii in [0.2,0.4,0.6,0.8]:
            if ii == 0.8:
                c = ["#cc33ff","#e699ff"]
#                 w = 10
            else:
                c = ["#4d94ff","#b3d1ff"]
            coords = np.array([np.append(triangular_coord((ii,0)),0),np.append(triangular_coord((ii,0)),self._min_energy*para),
                               np.append(triangular_coord((0,1)),0),np.append(triangular_coord((0,1)),self._min_energy*para)])
            facets = [0,1,2,3]
            print(coords)
            print(np.append(triangular_coord((0,1)),0))
            # i, j and k give the vertices of triangles
            # here we represent the 2 triangles of the rectangle
            for aa,bb in zip([0,0,2,1],[1,2,3,3]):
                data1.append(go.Scatter3d(
                x=[coords[aa][1],coords[bb][1]], y=[coords[aa][0],coords[bb][0]], z=[coords[aa][2],coords[bb][2]],
                    marker=dict(
                    size=8,
                    opacity=0.5,
                    color=c[0],
                    colorscale='Viridis',
                ),
                line=dict(
                    color=c[0],
                    width=w
                )))
            if ii == 0.8:
                data1.append(go.Mesh3d(
                    x=list(coords[:, 1]),
                    y=list(coords[:, 0]),
                    z=list(coords[:, 2]),
                    i=[facets[1],facets[1]],
                    j=[facets[0],facets[2]],
                    k=[facets[2],facets[3]],
                    opacity=0.2,
                    color = c[1],
                    hoverinfo="none",
                    lighting=dict(diffuse=0.0, ambient=1.0),
                    name="Convex Hull (shading)",
                    flatshading=True,
                    showlegend=True,
                ))
        return data1

    def create_reaction_plane3(self,compositions,marksize=10,linewidth=8):
        '''input compositions of two entries you want to link in a rectangle'''

        '''plot a series planes with AxB1-x'''
        para = 2
        data1 = []
        comp_vers = [[comp.get_atomic_fraction(self._pd.elements[1]),
                      comp.get_atomic_fraction(self._pd.elements[2])] for comp in compositions]
        ii = comp_vers[0]
        jj = comp_vers[1]
        c = ["#cc33ff","#e699ff"]
#                 w = 10
#         else:
#             c = ["#4d94ff","#b3d1ff"]
        coords = np.array([np.append(triangular_coord((ii[0],ii[1])),0),np.append(triangular_coord((ii[0],ii[1])),self._min_energy*para),
                           np.append(triangular_coord((jj[0],jj[1])),0),np.append(triangular_coord((jj[0],jj[1])),self._min_energy*para)])
        facets = [0,1,2,3]
        print(coords)
        print(np.append(triangular_coord((0,1)),0))
        # i, j and k give the vertices of triangles
        # here we represent the 2 triangles of the rectangle
        for aa,bb in zip([0,0,2,1],[1,2,3,3]):
            data1.append(go.Scatter3d(
            x=[coords[aa][1],coords[bb][1]], y=[coords[aa][0],coords[bb][0]], z=[coords[aa][2],coords[bb][2]],
                marker=dict(
                size=marksize,
                opacity=0.5,
                color=c[0],
                colorscale='Viridis',
            ),
            line=dict(
                color=c[0],
                width=linewidth
            )))

        data1.append(go.Mesh3d(
            x=list(coords[:, 1]),
            y=list(coords[:, 0]),
            z=list(coords[:, 2]),
            i=[facets[1],facets[1]],
            j=[facets[0],facets[2]],
            k=[facets[2],facets[3]],
            opacity=0.2,
            color = c[1],
            hoverinfo="none",
            lighting=dict(diffuse=0.0, ambient=1.0),
            name="Convex Hull (shading)",
            flatshading=True,
            showlegend=True,
        ))
        return data1
    
    def get_irpd(self,comp1=Composition("N2"),comp2=Composition("BaMn4"),marksize=10,linewidth=8):
        c = ["#cc33ff","#e699ff"]
        cricomps = self._pd.get_critical_compositions(comp1,comp2)
        print()
        vertices = []
        for comp in cricomps:
            formE = (self._pd.get_hull_energy(comp)- sum([comp[el] * self._pd.el_refs[el].energy_per_atom for el in comp.elements]))/comp.num_atoms
            ver = np.append(triangular_coord([comp.get_atomic_fraction(self._pd.elements[1]),
                              comp.get_atomic_fraction(self._pd.elements[2])]),formE)
            print(comp)
            print(ver)
            vertices.append(ver.tolist())
        vertices = np.array(vertices)
        print(vertices)
        print(vertices[:,0])
        data1 = []
        for i in range(0,len(vertices)-1):
            cc = c[0]
#             if i==2:
#                 cc="#4d94ff"
            print(vertices[:,0][i:i+2])
            data1.append(go.Scatter3d(
                x=vertices[:,1][i:i+2], y=vertices[:,0][i:i+2], z=vertices[:,2][i:i+2],
                    marker=dict(
                    size=marksize,
                    color=cc,
                    colorscale='Viridis',
                ),
                line=dict(
                    color=c[0],
                    width=linewidth
                )))
        return data1
    

    def create_tangent_plane(self,entryname = "Ge2N2O",marksize=8,linewidth=8):

        '''xx yy zz are elemental mu at the intercepts (the order is the elements order)'''
        CPentries = trans_PD_to_ChemPot_entries(self._pd.stable_entries, 
                                                [str(el) for el in self._pd.elements])
        cp=EquilLine(CPentries, [str(el) for el in self._pd.elements])
        for entry, vertices in cp._stable_domain_vertices.items():
            print(entry.name)
            if entry.name == entryname:
                tarentry = entry
                xx,yy,zz = np.average(vertices, axis=0)
        print(xx,yy,zz)
        
        data1 = []
        coords = np.array([[0,0,xx],[1,0,yy],[0.5,0.8660254037844386,zz]])
        plotline = True
        if plotline:
            Ba = coords[0]
            Mn = coords[1]
            N = coords[2]
            point1 = (Ba-N)*(1-tarentry.composition.get_atomic_fraction(self._pd.elements[-1]))+N
            point2 = (Mn-N)*(1-tarentry.composition.get_atomic_fraction(self._pd.elements[-1]))+N
            
            data1.append(go.Scatter3d(
                x=[point1[1],point2[1]], y=[point1[0],point2[0]], z=[point1[2],point2[2]],
                    marker=dict(
                    size=marksize,
                    opacity=0.5,
                    color="#4d94ff",
                    colorscale='Viridis',
                ),
                line=dict(
                    dash='dash',
                    color="#4d94ff",
                    width=linewidth
                )))
        print(coords[:,0])
        facets = [0,1,2]
        for aa,bb in zip([0,1,2],[1,2,0]):
            data1.append(go.Scatter3d(
            x=[coords[aa][1],coords[bb][1]], y=[coords[aa][0],coords[bb][0]], z=[coords[aa][2],coords[bb][2]],
                marker=dict(
                size=marksize,
                opacity=0.5,
                color="#4d94ff",
                colorscale='Viridis',
            ),
            line=dict(
                color="#4d94ff",
                width=linewidth
            )))
        data1.append(go.Mesh3d(
            x=list(coords[:, 1]),
            y=list(coords[:, 0]),
            z=list(coords[:, 2]),
            i=[facets[1]],
            j=[facets[0]],
            k=[facets[2]],
            opacity=0.2,
            color = "#4d94ff",
            hoverinfo="none",
            lighting=dict(diffuse=0.0, ambient=1.0),
            name="Convex Hull (shading)",
            flatshading=True,
            showlegend=True,
        ))
        return data1
        

    def _create_plotly_ternary_hull(self):
        """
        Creates shaded mesh plot for coloring the ternary hull by formation energy.

        :return: go.Mesh3d plot
        """
        facets = np.array(self._pd.facets)
        coords = np.array(
            [
                triangular_coord(c)
                for c in zip(self._pd.qhull_data[:-1, 0], self._pd.qhull_data[:-1, 1])
            ]
        )
        for c in zip(self._pd.qhull_data[:-1, 0], self._pd.qhull_data[:-1, 1]):
            print("haha",c, triangular_coord(c))
            
        energies = np.array(
            [self._pd.get_form_energy_per_atom(e) for e in self._pd.qhull_entries]
        )
#         for a,b in zip(facets, coords):
#             print(b, a)

        return go.Mesh3d(
            x=list(coords[:, 1]),
            y=list(coords[:, 0]),
            z=list(energies),
            i=list(facets[:, 1]),
            j=list(facets[:, 0]),
            k=list(facets[:, 2]),
            opacity=0.8,
            intensity=list(energies),
            colorscale=plotly_layouts["stable_colorscale"],
            colorbar=dict(title="Formation energy<br>(eV/atom)", x=0.9, len=0.75),
            hoverinfo="none",
            lighting=dict(diffuse=0.0, ambient=1.0),
            name="Convex Hull (shading)",
            flatshading=True,
            showlegend=True,
        )



    def _create_plotly_lines(self):
        """
        Creates Plotly scatter (line) plots for all phase diagram facets.

        :return: go.Scatter (or go.Scatter3d) plot
        """
        line_plot = None
        x, y, z, energies = [], [], [], []

        for line in self.pd_plot_data[0]:
            x.extend(list(line[0]) + [None])
            y.extend(list(line[1]) + [None])

            if self._dim == 3:
                z.extend(
                    [
                        self._pd.get_form_energy_per_atom(self.pd_plot_data[1][coord])
                        for coord in zip(line[0], line[1])
                    ]
                    + [None]
                )

            elif self._dim == 4:
                energies.extend(
                    [
                        self._pd.get_form_energy_per_atom(self.pd_plot_data[1][coord])
                        for coord in zip(line[0], line[1], line[2])
                    ]
                    + [None]
                )
                z.extend(list(line[2]) + [None])

        plot_args = dict(
            mode="lines",
            hoverinfo="none",
            line={"color": "rgba(0,0,0,1.0)", "width": 3.0},
            showlegend=False,
        )

        if self._dim == 2:
            line_plot = go.Scatter(x=x, y=y, **plot_args)
        elif self._dim == 3:
            line_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            line_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)

        return line_plot


    def _create_plotly_ternary_support_lines(self):
        """
        Creates support lines which aid in seeing the ternary hull in three
        dimensions.

        :return: go.Scatter3d plot of support lines for ternary phase diagram.
        """
        stable_entry_coords = dict(map(reversed, self.pd_plot_data[1].items()))
        for e in stable_entry_coords:
            print(e.name,stable_entry_coords[e])
        print()
        elem_coords = [stable_entry_coords[e] for e in self._pd.el_refs.values()]
        print(elem_coords)

        # add top and bottom triangle guidelines
        para = 2
        x, y, z = [], [], []
        for line in itertools.combinations(elem_coords, 2):
            x.extend([line[0][0], line[1][0], None] * 2)
            y.extend([line[0][1], line[1][1], None] * 2)
            z.extend([0, 0, None, self._min_energy*para, self._min_energy*para, None])

        # add vertical guidelines
        for elem in elem_coords:
            x.extend([elem[0], elem[0], None])
            y.extend([elem[1], elem[1], None])
            z.extend([0, self._min_energy*para, None])

        return go.Scatter3d(
            x=list(y),
            y=list(x),
            z=list(z),
            mode="lines",
            hoverinfo="none",
            line=dict(color="rgba (0, 0, 0, 0.4)", dash="solid", width=5.0),
            showlegend=False,
        )

    def _create_plotly_markers(self, label_uncertainties=False,marksize=5):
        """
        Creates stable and unstable marker plots for overlaying on the phase diagram.

        :return: Tuple of Plotly go.Scatter (or go.Scatter3d) objects in order: (
            stable markers, unstable markers)
        """

        def get_marker_props(coords, entries, stable=True):
            """ Method for getting marker locations, hovertext, and error bars
            from pd_plot_data"""
            x, y, z, texts, energies, uncertainties = [], [], [], [], [], []
            colors = []
            colordict = {'Mn2N': [1.0, 0.0, 0.16, 1.0], 'N2': [1.0, 0.325384207737149, 0.0, 1.0], 'MnN': [1.0, 0.81293057763646, 0.0, 1.0], 'Mn4N': [0.6995230524642289, 1.0, 0.0, 1.0], 'Mn': [0.19077901430842603, 1.0, 0.0, 1.0], 'BaN6': [0.0, 1.0, 0.29517183217372983, 1.0], 'BaMnN2': [0.0, 1.0, 0.7800969850305717, 1.0], 'BaN2': [0.0, 0.7320971867007668, 1.0, 1.0], 'Ba3MnN3': [0.0, 0.22058823529411742, 1.0, 1.0], 'Ba2N': [0.2696078431372551, 0.0, 1.0, 1.0], 'Ba3N': [0.7598039215686277, 0.0, 1.0, 1.0], 'Ba': [1.0, 0.0, 0.75, 1.0]}

            for coord, entry in zip(coords, entries):
                energy = round(self._pd.get_form_energy_per_atom(entry), 3)
                colors.append("k")
                entry_id = getattr(entry, "entry_id", "no ID")
                comp = entry.composition

                if hasattr(entry, "original_entry"):
                    comp = entry.original_entry.composition

                formula = comp.reduced_formula
                clean_formula = self._htmlize_formula(formula)
                label = f"{clean_formula} ({entry_id}) <br> " f"{energy} eV/atom"

                if not stable:
                    e_above_hull = round(self._pd.get_e_above_hull(entry), 3)
                    if e_above_hull > self.show_unstable:
                        continue
                    label += f" (+{e_above_hull} eV/atom)"
                    energies.append(e_above_hull)
                else:
                    uncertainty = 0
                    if (
                        hasattr(entry, "correction_uncertainty_per_atom")
                        and label_uncertainties
                    ):
                        uncertainty = round(entry.correction_uncertainty_per_atom, 4)
                        label += f"<br> (Error: +/- {uncertainty} eV/atom)"

                    uncertainties.append(uncertainty)
                    energies.append(energy)

                texts.append(label)

                x.append(coord[0])
                y.append(coord[1])

                if self._dim == 3:
                    z.append(energy)
                elif self._dim == 4:
                    z.append(coord[2])

            return {
                "x": x,
                "y": y,
                "z": z,
                "texts": texts,
                "energies": energies,
                "uncertainties": uncertainties,
                "colors":colors
            }

        stable_coords, stable_entries = (
            self.pd_plot_data[1].keys(),
            self.pd_plot_data[1].values(),
        )
        unstable_entries, unstable_coords = (
            self.pd_plot_data[2].keys(),
            self.pd_plot_data[2].values(),
        )

        stable_props = get_marker_props(stable_coords, stable_entries)

        unstable_props = get_marker_props(
            unstable_coords, unstable_entries, stable=False
        )

        stable_markers, unstable_markers = dict(), dict()

        if self._dim == 2:
            stable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["x"]),
                    y=list(stable_props["y"]),
                    name="Stable",
                    marker=dict(
                        color="darkgreen", size=11, line=dict(color="black", width=2)
                    ),
                    opacity=0.9,
                    hovertext=stable_props["texts"],
                    error_y=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="gray",
                        thickness=2.5,
                        width=5,
                    ),
                )
            )

            unstable_markers = plotly_layouts["default_binary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=list(unstable_props["x"]),
                    y=list(unstable_props["y"]),
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=6,
                        symbol="diamond",
                    ),
                    hovertext=unstable_props["texts"],
                )
            )

        elif self._dim == 3:
            stable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["y"]),
                    y=list(stable_props["x"]),
                    z=list(stable_props["z"]),
#                     z=[0 for ii in range(len(list(stable_props["z"])))],
                    name="Stable",
                    marker=dict(
                        color="#005ce6",
#                         color = stable_props["colors"],
                        size=marksize,
                        opacity=0.6,
                        line=dict(width=3,color="black"),
                    ),
#                     marker_symbol="square",
                    hovertext=stable_props["texts"],
                    error_z=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="darkgray",
                        width=10,
                        thickness=5,
                    ),
                )
            )

            unstable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            unstable_markers.update(
                dict(
                    x=unstable_props["y"],
                    y=unstable_props["x"],
                    z=unstable_props["z"],
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=6,
                        symbol="diamond",
                        colorbar=dict(
                            title="Energy Above Hull<br>(eV/atom)", x=0.05, len=0.75
                        ),
                    ),
                    hovertext=unstable_props["texts"],
                )
            )

        elif self._dim == 4:
            stable_markers = plotly_layouts["default_quaternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=stable_props["x"],
                    y=stable_props["y"],
                    z=stable_props["z"],
                    name="Stable",
                    marker=dict(
                        color=stable_props["energies"],
                        colorscale=plotly_layouts["stable_markers_colorscale"],
                        size=8,
                        opacity=0.9,
                    ),
                    hovertext=stable_props["texts"],
                )
            )

            unstable_markers = plotly_layouts[
                "default_quaternary_marker_settings"
            ].copy()
            unstable_markers.update(
                dict(
                    x=unstable_props["x"],
                    y=unstable_props["y"],
                    z=unstable_props["z"],
                    name="Above Hull",
                    marker=dict(
                        color=unstable_props["energies"],
                        colorscale=plotly_layouts["unstable_colorscale"],
                        size=5,
                        symbol="diamond",
                        colorbar=dict(
                            title="Energy Above Hull<br>(eV/atom)", x=0.05, len=0.75
                        ),
                    ),
                    hovertext=unstable_props["texts"],
                    visible="legendonly",
                )
            )

        stable_marker_plot = (
            go.Scatter(**stable_markers)
            if self._dim == 2
            else go.Scatter3d(**stable_markers)
        )
        unstable_marker_plot = (
            go.Scatter(**unstable_markers)
            if self._dim == 2
            else go.Scatter3d(**unstable_markers)
        )

        return stable_marker_plot, unstable_marker_plot
        









