'''
Created on Jun 5, 2021

@author: jiadongc
'''
from pymatgen.analysis.phase_diagram import GrandPotentialPhaseDiagram,PhaseDiagram, PDPlotter
import plotly.graph_objects as go
import json
import math
from pymatgen.core.periodic_table import Element
import itertools
import numpy as np
from pymatgen.analysis.reaction_calculator import ComputedReaction
with open("C:/Users/jiadongc/anaconda3/envs/my_pymatgen/lib/site-packages/pymatgen/util/plotly_pd_layouts.json", "r") as f:
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

class GrandPotentialPhaseDiagram_E(GrandPotentialPhaseDiagram):
    def get_form_energy(self, entry):
        """
        Returns the formation energy for an entry (NOT normalized) from the
        elemental references.

        Args:
            entry (PDEntry): A PDEntry-like object.

        Returns:
            Formation energy from the elemental references.
        """
        c = entry.composition
#         for el in c.elements:
#             print()
#             print(el, c[el], self.el_refs[el].energy_per_atom)
        return entry.energy
    
class GPDPlotter(PDPlotter):
    def show(self, *args, **kwargs):
        r"""
        Draw the phase diagram using Plotly (or Matplotlib) and show it.

        Args:
            *args: Passed to get_plot.
            **kwargs: Passed to get_plot.
        """
        fig = self.get_plot(*args, **kwargs)
        fig.update_yaxes(range = [-23,0])
        fig.show()
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

#             stable_labels_plot = self._create_plotly_stable_labels(label_stable)
            stable_marker_plot, unstable_marker_plot = self._create_plotly_markers(
                label_uncertainties
            )

            if self._dim == 2 and label_uncertainties:
                data.append(self._create_plotly_uncertainty_shading(stable_marker_plot))

#             data.append(stable_labels_plot)
#             data.append(unstable_marker_plot)
            data.append(stable_marker_plot)
            
#             fig = go.Figure(data=data)
#             fig.layout = self._create_plotly_figure_layout()
            layout = self._create_plotly_figure_layout()
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

        return data,layout

    def _create_plotly_element_annotations(self):
        """
        Creates terminal element annotations for Plotly phase diagrams.

        :return: list of annotation dicts.
        """
        annotations_list = []
        x, y, z = None, None, None
        mu = None
        for coords, entry in self.pd_plot_data[1].items():
            if not entry.composition.is_element:
                continue

            x, y = coords[0], coords[1]
            mu = list(entry.chempots.values())[0]
            if self._dim == 3:
                z = self._pd.get_form_energy_per_atom(entry)
            elif self._dim == 4:
                z = coords[2]

            if entry.composition.is_element:
                clean_formula = str(entry.composition.elements[0])
                if hasattr(entry, "original_entry"):
                    orig_comp = entry.original_entry.composition
                    clean_formula = self._htmlize_formula(orig_comp.reduced_formula)

                font_dict = {"color": "#000000", "size": 24.0}
                opacity = 1.0

            annotation = plotly_layouts["default_annotation_layout"].copy()
            for d in ["xref", "yref"]:
                annotation.pop(d)  # Scatter3d cannot contain xref, yref
            annotation.update(
                {
                    "x": x,
                    "y": mu,
                    "z": y,
                    "font": font_dict,
                    "text": clean_formula,
                    "opacity": opacity,
                }
            )

            if self._dim == 3 or self._dim == 4:
                for d in ["xref", "yref"]:
                    annotation.pop(d)  # Scatter3d cannot contain xref, yref
                    if self._dim == 3:
                        annotation.update({"x": y, "y": x})
                        if entry.composition.is_element:
                            z = 0.9 * self._min_energy  # place label 10% above base

                annotation.update({"z": z})

            annotations_list.append(annotation)

        # extra point ensures equilateral triangular scaling is displayed
        if self._dim == 3:
            annotations_list.append(dict(x=1, y=1, z=0, opacity=0, text=""))

        return annotations_list
    
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
        return annotations_list
#         if self._dim == 2:
# 
#             layout = plotly_layouts["advanced_ternary_layout"].copy()
#             layout["scene"].update({"annotations": annotations_list})
# 
#         elif self._dim == 3:
#             layout = plotly_layouts["default_ternary_layout"].copy()
#             layout["scene"].update({"annotations": annotations_list})
#         elif self._dim == 4:
#             layout = plotly_layouts["default_quaternary_layout"].copy()
#             layout["scene"].update({"annotations": annotations_list})

#         return layout
    def _create_plotly_markers(self, label_uncertainties=False):
        """
        Creates stable and unstable marker plots for overlaying on the phase diagram.
 
        :return: Tuple of Plotly go.Scatter (or go.Scatter3d) objects in order: (
            stable markers, unstable markers)
        """
 
        def get_marker_props(coords, entries, stable=True):
            """ Method for getting marker locations, hovertext, and error bars
            from pd_plot_data"""
            x, y, z, texts, energies, uncertainties = [], [], [], [], [], []
            mu = []
            for coord, entry in zip(coords, entries):
                energy = round(self._pd.get_form_energy_per_atom(entry), 3)
 
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
                mu.append(list(entry.chempots.values())[0])
                if self._dim == 3:
                    z.append(energy)
                elif self._dim == 4:
                    z.append(coord[2])
 
            return {
                "x": x,
                "y": mu,
                "z": y,
                "texts": texts,
                "energies": energies,
                "uncertainties": uncertainties,
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
            stable_markers = plotly_layouts["default_ternary_marker_settings"].copy()
            stable_markers.update(
                dict(
                    x=list(stable_props["x"]),
                    y=list(stable_props["y"]),
                    z=list(stable_props["z"]),
                    name="Stable",
                    marker=dict(
                        color="darkgreen", size=6, line=dict(color="black", width=2)
                    ),
                    opacity=0.8,
                    hovertext=stable_props["texts"],
                    error_z=dict(
                        array=list(stable_props["uncertainties"]),
                        type="data",
                        color="gray",
                        thickness=2.5,
                        width=5,
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
                    name="Stable",
                    marker=dict(
                        color="black",
                        size=12,
                        opacity=0.8,
                        line=dict(color="black", width=3),
                    ),
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
            go.Scatter3d(**stable_markers)
            if self._dim == 2
            else go.Scatter3d(**stable_markers)
        )
        unstable_marker_plot = (
            go.Scatter3d(**unstable_markers)
            if self._dim == 2
            else go.Scatter3d(**unstable_markers)
        )
 
        return stable_marker_plot, unstable_marker_plot
     
    def _create_plotly_lines(self):
        """
        Creates Plotly scatter (line) plots for all phase diagram facets.
  
        :return: go.Scatter (or go.Scatter3d) plot
        """
        line_plot = None
        x, y, z, energies = [], [], [], []
        mu = []
        (lines, stable_entries, unstables) = self.pd_plot_data
          
        for e in stable_entries:
            print(stable_entries[e].name,stable_entries[e].chempots)
            muN = list(stable_entries[e].chempots.values())[0]
            break
        for line in lines:
            x.extend(list(line[0]) + [None])
            y.extend(list(line[1]) + [None])
            mu.extend([muN,muN,None])
  
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
            line={"color": "rgba(0,0,0,1.0)", "width": 7.0},
            showlegend=False,
        )
  
        if self._dim == 2:
            line_plot = go.Scatter3d(x=x, y=mu, z=y, **plot_args)
        elif self._dim == 3:
            line_plot = go.Scatter3d(x=y, y=x, z=z, **plot_args)
        elif self._dim == 4:
            line_plot = go.Scatter3d(x=x, y=y, z=z, **plot_args)
  
        return line_plot
 
 
 
 
 

