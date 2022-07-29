"""
Load crossovers and flightlines from CSVs and plot them for comparison.

This script was made to verify crossover calculation following the implementation of
segment geom simplification in the database.

# Dependencies

A GDAL whl can be obtained from https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal.
    Get the whl that matches your version of python and install it with
    `pip install ./GDAL-3.4.2-cp37-cp37m-win32.whl`
    replacing the file with the path to the one you downloaded.
Likewise, a Fiona whl can obtained from https://www.lfd.uci.edu/~gohlke/pythonlibs/#fiona.
    Install it in the same manner.
Install the remaining dependencies with `pip install -r requirements.txt`


Author: Reece Mathews
"""
import csv
import os
import ctypes

from pathlib import Path
from typing import List

import geopandas
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import matplotlib as mpl

from shapely import wkb
from shapely.ops import transform
import pyproj

load_ruler = True

try:
    from matplotlib_tools import Ruler
except ModuleNotFoundError:
    load_ruler = False
    print("Not loading ruler as the matplotlib_tools.py module was not found")
    print("This file can be obtained from https://github.com/terranjp/matplotlib-tools/blob/master/matplotlib_tools/tools.py")
    print("I have opted not to include this module in the repo as no license is provided by the author.")

# Default csv field size limit too small for segment geoms
# Set to max size supported by system
# https://stackoverflow.com/a/54517228
csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))

# Set working directory to crossover-test-suite
if Path(os.getcwd()).name == "cresis-toolbox":
    os.chdir("python/crossover-test-suite")


DATA_DIR = Path("data")
WIDGETS = []

COLORS = {
        "1m segments":   "C1",
        "1m crossovers": "C2",
        "0m segments":   "C3",
        "0m crossovers": "C4",
    }

DEFAULT_LIMITS = [-3.5e6, -3.5e6, 1e6, -0.5e6]

wgs84 = pyproj.CRS('EPSG:4326')
greenland = pyproj.CRS('EPSG:3413')
projection = pyproj.Transformer.from_crs(wgs84, greenland, always_xy=True).transform


def reproject(geom_obj):
    """Convert a geometry to the 3413 CRS."""
    return transform(projection, geom_obj)


def convert_geom(geom: str):
    """Convert geom strs from the Postgis HEXEWKB format to WKB."""
    geom_obj = wkb.loads(geom, hex=True)
    geom_proj = reproject(geom_obj)
    return geom_proj


def load_data():
    """Load each CSV from the DATA_DIR."""
    data = {}

    for file in os.listdir(DATA_DIR):
        name = Path(file).stem
        with open(DATA_DIR / file, newline='') as f:
            data[name] = list(csv.DictReader(f))

        for row in data[name]:
            for field in row:
                if field.endswith("geom"):
                    geom_str = row[field]
                    row[field] = convert_geom(geom_str)

    return data


def plot_map():
    """Make a base plot of the Greenland map."""
    map_df = geopandas.read_file("maps/GRL_adm0.shp")
    map_df = map_df.set_crs("EPSG:4326")
    map_df = map_df.to_crs("EPSG:3413")
    return map_df.boundary.plot(color='black')


def plot_geoms(geoms: List[str], base, color=None, zorder=1):
    """Plot a list of geoms string on the map base with the given color."""

    # Pandas plot method returns the base axes rather than the object plotted like Matplotlib does
    # Keep note of children of base axes before plotting so that we can find the new child after plotting
    old_children = set(base.get_children())
    df = geopandas.GeoDataFrame(geometry=geoms)
    ax = df.plot(ax=base, color=color, zorder=zorder)

    # Find the child we just plotted to return it
    plotted_child = list(set(ax.get_children()) - old_children)
    if plotted_child:
        return plotted_child[0]

    return None


class VisibilityState():
    showing_0m = True
    showing_1m = True

    def __init__(self, plot_seg_0m, plot_seg_1m, plot_cx_0m, plot_cx_1m):
        self.plot_seg_0m = plot_seg_0m
        self.plot_seg_1m = plot_seg_1m
        self.plot_cx_0m = plot_cx_0m
        self.plot_cx_1m = plot_cx_1m

    def toggle_0m(self, clk=None):
        self.showing_0m = not self.showing_0m
        if self.plot_seg_0m is not None:
            self.plot_seg_0m.set_visible(self.showing_0m)
        if self.plot_cx_0m is not None:
            self.plot_cx_0m.set_visible(self.showing_0m)
        # plt.draw()

    def toggle_1m(self, clk=None):
        self.showing_1m = not self.showing_1m
        if self.plot_seg_1m is not None:
            self.plot_seg_1m.set_visible(self.showing_1m)
        if self.plot_cx_1m is not None:
            self.plot_cx_1m.set_visible(self.showing_1m)
        # plt.draw()

    def set_state(self, state):
        self.plot_seg_1m.set_visible(state)
        self.plot_seg_0m.set_visible(state)
        self.plot_cx_0m.set_visible(state)
        self.plot_cx_1m.set_visible(state)


def plot_from_data(map_base, data):
    mpl.style.use("seaborn")

    if load_ruler:
        # Monkey-patch annotate method to handle erroneous use of "s" param instead of "text" by the matplotlib_tools module
        original_annotate = map_base.annotate
        map_base.annotate = lambda *args, **kwargs: original_annotate(text=kwargs["s"], *args, **{k: v for k, v in kwargs.items() if k != "s"}) if "s" in kwargs else original_annotate(*args, **kwargs)
        # Ruler is immediately garbage collected if not assigned to anything
        # Matplotlib widgets are maintained with only weakrefs internally. The user is expected to maintain their own strong reference
        WIDGETS.append(Ruler(ax=map_base, useblit=True, lineprops={"color": "red", "linewidth": "2"}))

    map_base.set_xlabel('Meters East of North pole')
    map_base.set_ylabel('Meters North of North pole')

    plot_seg_1m = plot_geoms([row["geom"] for row in data["1m segments"]], map_base, COLORS["1m segments"])
    plot_cx_1m = plot_geoms([row["cx_geom"] for row in data["1m crossovers"]], map_base, COLORS["1m crossovers"], 2)
    plot_seg_0m = plot_geoms([row["geom"] for row in data["0m segments"]], map_base, COLORS["0m segments"])
    plot_cx_0m = plot_geoms([row["cx_geom"] for row in data["0m crossovers"]], map_base, COLORS["0m crossovers"], 2)

    # Create visibility toggle buttons
    state = VisibilityState(plot_seg_0m, plot_seg_1m, plot_cx_0m, plot_cx_1m)

    ax_toggle_0 = plt.axes((0.9, 0.8, 0.1, 0.075))
    ax_toggle_1 = plt.axes((0.9, 0.7, 0.1, 0.075))

    b_toggle_0 = Button(ax_toggle_0, 'Toggle Full Res')
    b_toggle_0.on_clicked(state.toggle_0m)
    WIDGETS.append(b_toggle_0)  # Necessary to avoid garbage collection
    b_toggle_1 = Button(ax_toggle_1, 'Toggle 1m')
    b_toggle_1.on_clicked(state.toggle_1m)
    WIDGETS.append(b_toggle_1)

    map_base.legend([plot_seg_1m, plot_cx_1m, plot_seg_0m, plot_cx_0m], 
                    ['1m Segments', '1m Crossovers', 'Full Res Segments', 'Full Res Crossovers',])

    return state

    # plt.draw()


def get_segments(name, data):
    """Find the rows in the data files for a given segment by name."""
    segments = {}
    for row in data["0m segments"]:
        if row["name"] == name:
            segments["0m segments"] = row
            break
    for row in data["1m segments"]:
        if row["name"] == name:
            segments["1m segments"] = row
            break

    return segments


def plot_pair(pair, map_base, data):
    """Plot a pair of crossovers."""
    elements = []
    
    segments = pair.segment_pair.split(" ")
    for segment_name in segments:
        segment_objs = get_segments(segment_name, data)
        for segment_file in segment_objs:
            elements.append(plot_geoms([segment_objs[segment_file]["geom"]], map_base, COLORS[segment_file]))
    
    if pair.cx_pair[0] is not None:
        elements.append(plot_geoms([pair.cx_pair[0]["cx_geom"]], map_base, COLORS["0m crossovers"], 2))
    if pair.cx_pair[1] is not None:
        elements.append(plot_geoms([pair.cx_pair[1]["cx_geom"]], map_base, COLORS["1m crossovers"], 2))


    return elements
    


class DistanceState():

    def __init__(self, map_base, data, cx_distances, visibility_state):
        self.map_base = map_base
        self.data = data
        self.cx_distances = cx_distances
        self.selected_cx_idx = None
        self.visibility_state = visibility_state
        self.cx_plot_elements = None

    def view_prev(self, clk=None):
        if self.selected_cx_idx is None:
            self.selected_cx_idx = 1
        if self.selected_cx_idx == 0:
            self.selected_cx_idx = len(self.cx_distances)
        
        self.selected_cx_idx -= 1

        self.show_cx()
    
    def view_next(self, clk=None):
        if self.selected_cx_idx is None:
            self.selected_cx_idx = len(self.cx_distances) - 2
        if self.selected_cx_idx == len(self.cx_distances) - 1:
            self.selected_cx_idx = -1
        
        self.selected_cx_idx += 1

        self.show_cx()

    def show_cx(self):
        self.visibility_state.set_state(False)
        self.hide_cx()

        pair = self.cx_distances[self.selected_cx_idx]

        self.cx_plot_elements = plot_pair(pair, self.map_base, data)

        cx_1_bounds = []
        cx_2_bounds = []
        if pair.cx_pair[0] is not None and pair.cx_pair[1] is not None:
            cx_1_bounds = pair.cx_pair[0]["cx_geom"].bounds
            cx_2_bounds = pair.cx_pair[1]["cx_geom"].bounds
    
        elif pair.cx_pair[0] is not None:
            cx_1_bounds = pair.cx_pair[0]["cx_geom"].bounds
            cx_2_bounds = cx_1_bounds
        elif pair.cx_pair[1] is not None:
            cx_2_bounds = pair.cx_pair[1]["cx_geom"].bounds
            cx_1_bounds = cx_2_bounds

        minx = min(cx_1_bounds[0], cx_2_bounds[0]) - 1
        maxx = max(cx_1_bounds[2], cx_2_bounds[2]) + 1
        miny = min(cx_1_bounds[1], cx_2_bounds[1]) - 1
        maxy = max(cx_1_bounds[3], cx_2_bounds[3]) + 1

        width = maxx - minx
        height = maxy - miny

        # Keep the viewport square
        xmargin = 0
        ymargin = 0
        if width < height:
            xmargin = (height - width) // 2
        elif height < width:
            ymargin = (width - height) // 2

        self.map_base.set_xlim(
            minx - xmargin,
            maxx + xmargin
            )
        self.map_base.set_ylim(
            miny - ymargin,
            maxy + ymargin
            )

        # plt.draw()

    def hide_cx(self):
        if self.cx_plot_elements is not None:
            for element in self.cx_plot_elements:
                element.remove()
        self.cx_plot_elements = None

    def view_home(self, clk=None):
        self.hide_cx()
        self.visibility_state.set_state(True)
        self.map_base.set_xlim(DEFAULT_LIMITS[0], DEFAULT_LIMITS[2])
        self.map_base.set_ylim(DEFAULT_LIMITS[1], DEFAULT_LIMITS[3])


def plot_dist_analyzer(map_base, data, cx_distances, visibility_state):

    ax_next_cx = plt.axes((0.9, 0.6, 0.1, 0.075))
    ax_prev_cx = plt.axes((0.9, 0.5, 0.1, 0.075))
    ax_home_cx = plt.axes((0.9, 0.4, 0.1, 0.075))

    state = DistanceState(map_base, data, cx_distances, visibility_state)

    b_next_cx = Button(ax_next_cx, 'Next Large Distance')
    b_next_cx.on_clicked(state.view_next)
    WIDGETS.append(b_next_cx)
    b_prev_cx = Button(ax_prev_cx, 'Prev Large Distance')
    b_prev_cx.on_clicked(state.view_prev)
    WIDGETS.append(b_prev_cx)
    b_home = Button(ax_home_cx, 'Home')
    b_home.on_clicked(state.view_home)
    WIDGETS.append(b_home)

    # plt.draw()


if __name__ == "__main__":
    from compare_crossovers import find_differences

    data = load_data()
    map_base = plot_map()

    visibility_state = plot_from_data(map_base, data)
    cx_distances = find_differences(data["0m crossovers"], data["1m crossovers"])
    plot_dist_analyzer(map_base, data, cx_distances, visibility_state)

    plt.show()
