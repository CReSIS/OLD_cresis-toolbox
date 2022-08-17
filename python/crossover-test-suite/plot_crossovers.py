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

Map shape files can be obtained from https://maps.princeton.edu/catalog/stanford-sd368wz2435
GRL_adm0.shp and GRL_adm0.shx should be placed in a maps folder.

The data files csvs can be obtained by running the data_queries.py file on virtual boxes 
which are already simplified to the desired resolutions. Use the postgres-conn.sample.json
to produce a postgres-conn.json file which points to the vbox DB.

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
from matplotlib.offsetbox import AnchoredText
import fiona.errors

from shapely import wkb
from shapely.ops import transform
import pyproj


os.chdir(Path(__file__).parent)

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
TARGETA = "ops0/0m"
TARGETB = "ops0/0.1m"
# TARGETA = "vbox/15m"
# TARGETB = "vbox/1m"
IGNORED_SEGMENTS = ['20190512_02', '20190515_01', '20190516_02', '20190508_01']

COLORS = {
        f"{TARGETB} segments":   "C1",
        f"{TARGETB} crossovers": "C2",
        f"{TARGETA} segments":   "C3",
        f"{TARGETA} crossovers": "C5",
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


def load_data(data_types=("segments", "crossovers")):
    """Load each CSV from the DATA_DIR."""
    data = {}

    for target in (TARGETA, TARGETB):
        for file_type in data_types:
            name = f"{target} {file_type}"
            file = DATA_DIR / f"{name}.csv"
            with open(file, newline='') as f:
                data[name] = list(csv.DictReader(f))

            for row in data[name]:
                for field in row:
                    if field.endswith("geom"):
                        geom_str = row[field]
                        row[field] = convert_geom(geom_str)

    # Remove ignored segments
    for file in data:
        if file.endswith("segments"):
            data[file] = [row for row in data[file] if row["name"] not in IGNORED_SEGMENTS]
        if file.endswith("crossovers"):
            data[file] = [row for row in data[file] if row["seg1_name"] not in IGNORED_SEGMENTS 
                                                    and row["seg2_name"] not in IGNORED_SEGMENTS]

    return data


def plot_map():
    """Make a base plot of the Greenland map."""
    try:
        map_df = geopandas.read_file("maps/GRL_adm0.shp")
    except fiona.errors.DriverError:
        raise FileNotFoundError("No map files found")
    map_df = map_df.set_crs("EPSG:4326")
    map_df = map_df.to_crs("EPSG:3413")

    map_base = map_df.boundary.plot(color='black')
    return map_base


def plot_geoms(geoms: List[str], base, color=None, zorder=1):
    """Plot a list of geoms string on the map base with the given color."""

    # Pandas plot method returns the base axes rather than the object plotted like Matplotlib does
    # Keep note of children of base axes before plotting so that we can find the new child after plotting
    old_children = set(base.get_children())
    df = geopandas.GeoDataFrame(geometry=geoms)
    ax = df.plot(ax=base, color=color, zorder=zorder, markersize=30)

    # Find the child we just plotted to return it
    plotted_child = list(set(ax.get_children()) - old_children)
    if plotted_child:
        return plotted_child[0]

    return None


class VisibilityState():
    """Control display of all segments."""
    showing_A = True
    showing_B = True

    def __init__(self, data, map_base):
        self.plot_seg_A = None
        self.plot_seg_B = None
        self.plot_cx_A = None
        self.plot_cx_B = None
        self.data = data
        self.map_base = map_base

        self.toggle_button_A = None
        self.toggle_button_B = None

    def first_plot(self):
        self.plot_seg_A = plot_geoms([row["geom"] for row in data[f"{TARGETA} segments"]], map_base, COLORS[f"{TARGETA} segments"])
        self.plot_cx_A = plot_geoms([row["cx_geom"] for row in data[f"{TARGETA} crossovers"]], map_base, COLORS[f"{TARGETA} crossovers"], 2)
        self.plot_seg_B = plot_geoms([row["geom"] for row in data[f"{TARGETB} segments"]], map_base, COLORS[f"{TARGETB} segments"])
        self.plot_cx_B = plot_geoms([row["cx_geom"] for row in data[f"{TARGETB} crossovers"]], map_base, COLORS[f"{TARGETB} crossovers"], 2)
        self.map_base.legend([self.plot_seg_B, self.plot_cx_B, self.plot_seg_A, self.plot_cx_A], 
                             [f"{TARGETB} Segments", f'{TARGETB} Crossovers', f'{TARGETA} Segments', f'{TARGETA} Crossovers',])


    def toggle_A(self, clk=None):
        self.showing_A = not self.showing_A
        if self.plot_seg_A is not None:
            self.plot_seg_A.set_visible(self.showing_A)
        if self.plot_cx_A is not None:
            self.plot_cx_A.set_visible(self.showing_A)

    def toggle_B(self, clk=None):
        self.showing_B = not self.showing_B
        if self.plot_seg_B is not None:
            self.plot_seg_B.set_visible(self.showing_B)
        if self.plot_cx_B is not None:
            self.plot_cx_B.set_visible(self.showing_B)

    def set_state(self, state):
        if self.plot_seg_A is None:
            if not state:
                return 
            self.first_plot()
        self.plot_seg_A.set_visible(state)
        self.plot_seg_B.set_visible(state)
        self.plot_cx_A.set_visible(state)
        self.plot_cx_B.set_visible(state)

    def set_button_state(self, state):
        self.toggle_button_A.set_visible(state)
        self.toggle_button_B.set_visible(state)


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

    # Create visibility toggle buttons
    state = VisibilityState(data, map_base)

    ax_toggle_A = plt.axes((0.9, 0.8, 0.1, 0.075))
    ax_toggle_B = plt.axes((0.9, 0.7, 0.1, 0.075))

    b_toggle_A = Button(ax_toggle_A, f'Toggle {TARGETA}')
    b_toggle_A.on_clicked(state.toggle_A)
    WIDGETS.append(b_toggle_A)  # Necessary to avoid garbage collection
    b_toggle_B = Button(ax_toggle_B, f'Toggle {TARGETB}')
    b_toggle_B.on_clicked(state.toggle_B)
    WIDGETS.append(b_toggle_B)

    state.toggle_button_A = ax_toggle_A
    state.toggle_button_B = ax_toggle_B

    return state


def get_segments(name, data):
    """Find the rows in the data files for a given segment by name."""
    segments = {}
    for row in data[f"{TARGETA} segments"]:
        if row["name"] == name:
            segments[f"{TARGETA} segments"] = row
            break
    for row in data[f"{TARGETB} segments"]:
        if row["name"] == name:
            segments[f"{TARGETB} segments"] = row
            break

    return segments


def plot_pair(pair, map_base, data):
    """Plot a pair of crossovers."""
    elements = []
    
    segments = pair.segment_pair.split(" ")
    for i, segment_name in enumerate(segments):
        segment_objs = get_segments(segment_name, data)
        for segment_file in segment_objs:
            elements.append(plot_geoms([segment_objs[segment_file]["geom"]], map_base, COLORS[segment_file]))
            if i == 0:
                elements[-1].set_label(segment_file)

    if pair.cx_pair[0] is not None:
        elements.append(plot_geoms([pair.cx_pair[0]["cx_geom"]], map_base, COLORS[f"{TARGETA} crossovers"], 2))
        elements[-1].set_label(f"{TARGETA} crossovers")
    if pair.cx_pair[1] is not None:
        elements.append(plot_geoms([pair.cx_pair[1]["cx_geom"]], map_base, COLORS[f"{TARGETB} crossovers"], 2))
        elements[-1].set_label(f"{TARGETB} crossovers")

    map_base.legend()
    return elements
    

class DistanceState():
    """Control display of crossover comparisons"""

    def __init__(self, map_base, data, cx_distances, visibility_state):
        self.map_base = map_base
        self.data = data
        self.cx_distances = cx_distances
        self.selected_cx_idx = None
        self.visibility_state = visibility_state
        self.cx_plot_elements = None
        self.home_button = None
        self.zoom_toggle_button = None
        self.cx_text = None
        self.bounds = None
        self.zoomed = False

    def view_prev(self, clk=None):
        if self.selected_cx_idx is None:
            self.zoomed = True

        if self.selected_cx_idx is None or self.selected_cx_idx == 0:
            self.selected_cx_idx = len(self.cx_distances)
        
        self.selected_cx_idx -= 1

        self.show_cx()
    
    def view_next(self, clk=None):
        if self.selected_cx_idx is None:
            self.zoomed = True

        if self.selected_cx_idx is None or self.selected_cx_idx == len(self.cx_distances) - 1:
            self.selected_cx_idx = -1
        
        self.selected_cx_idx += 1

        self.show_cx()

    def show_cx(self):
        self.visibility_state.set_state(False)
        self.home_button.set_visible(True)
        self.zoom_toggle_button.set_visible(True)
        self.visibility_state.set_button_state(False)
        self.hide_cx()

        pair = self.cx_distances[self.selected_cx_idx]

        self.cx_plot_elements = plot_pair(pair, self.map_base, self.data)

        # Get info for annotation box
        cx_1_pp1 = None
        cx_1_pp2 = None
        cx_2_pp1 = None
        cx_2_pp2 = None
        if pair.cx_pair[0] is not None:
            cx_1_pp1 = pair.cx_pair[0]["pp1_id"]
            cx_1_pp2 = pair.cx_pair[0]["pp2_id"]
        if pair.cx_pair[1] is not None:
            cx_2_pp1 = pair.cx_pair[1]["pp1_id"]
            cx_2_pp2 = pair.cx_pair[1]["pp2_id"]
        
        self.cx_text.txt.set_text(f"Crossover {self.selected_cx_idx}\nDistance {pair.distance:.3} m\nSegments {pair.segment_pair}\nFull Res point path ids {cx_1_pp1} {cx_1_pp2}\n1m point path ids {cx_2_pp1} {cx_2_pp2}")

        # Set window viewport
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
            xmargin = (height - width) / 2
        elif height < width:
            ymargin = (width - height) / 2

        self.bounds = [minx - xmargin, miny - ymargin, maxx + xmargin, maxy + ymargin]

        if self.zoomed:
            self.zoom_in()

    def zoom_in(self):
        if self.bounds is not None:
            self.map_base.set_xlim(
                self.bounds[0],
                self.bounds[2]
                )
            self.map_base.set_ylim(
                self.bounds[1],
                self.bounds[3]
                )
            self.zoomed = True

    def zoom_out(self):
        self.map_base.set_xlim(DEFAULT_LIMITS[0], DEFAULT_LIMITS[2])
        self.map_base.set_ylim(DEFAULT_LIMITS[1], DEFAULT_LIMITS[3])
        self.zoomed = False

    def hide_cx(self):
        if self.cx_plot_elements is not None:
            for element in self.cx_plot_elements:
                element.remove()
        self.cx_plot_elements = None

    def view_home(self, clk=None):
        self.hide_cx()
        self.visibility_state.set_state(True)
        self.visibility_state.set_button_state(True)
        self.home_button.set_visible(False)
        self.zoom_toggle_button.set_visible(False)
        self.zoom_out()
        self.cx_text.txt.set_text("")

    def zoom_toggle(self, clk=None):
        self.zoomed = not self.zoomed
        if self.zoomed:
            self.zoom_in()
        else:
            self.zoom_out()


def plot_dist_analyzer(map_base, data, cx_distances, visibility_state):

    ax_next_cx = plt.axes((0.9, 0.6, 0.1, 0.075))
    ax_prev_cx = plt.axes((0.9, 0.5, 0.1, 0.075))
    ax_zoom_toggle_cx = plt.axes((0.9, 0.4, 0.1, 0.075))
    ax_home_cx = plt.axes((0.9, 0.3, 0.1, 0.075))

    state = DistanceState(map_base, data, cx_distances, visibility_state)

    b_next_cx = Button(ax_next_cx, 'Show Next CX')
    b_next_cx.on_clicked(state.view_next)
    WIDGETS.append(b_next_cx)
    b_prev_cx = Button(ax_prev_cx, 'Show Prev CX')
    b_prev_cx.on_clicked(state.view_prev)
    WIDGETS.append(b_prev_cx)
    b_zoom_toggle = Button(ax_zoom_toggle_cx, 'Toggle Zoom')
    b_zoom_toggle.on_clicked(state.zoom_toggle)
    WIDGETS.append(b_zoom_toggle)
    b_home = Button(ax_home_cx, 'Show all')
    b_home.on_clicked(state.view_home)
    WIDGETS.append(b_home)

    cx_text = AnchoredText(
        "", prop={"size": 10}, frameon=False, loc='lower left')
    map_base.add_artist(cx_text)

    state.home_button = ax_home_cx
    state.home_button.set_visible(False)
    state.zoom_toggle_button = ax_zoom_toggle_cx
    state.zoom_toggle_button.set_visible(False)

    state.cx_text = cx_text
    state.view_next()


if __name__ == "__main__":
    from compare_crossovers import find_differences

    print("Loading Data")
    data = load_data()
    print("Plotting map")
    map_base = plot_map()

    visibility_state = plot_from_data(map_base, data)
    print("Comparing crossovers")
    cx_distances = find_differences(data[f"{TARGETA} crossovers"], data[f"{TARGETB} crossovers"])
    print("Plotting")
    plot_dist_analyzer(map_base, data, cx_distances, visibility_state)

    plt.show()


# TODO[Reece]: Do not assume the name or order of the files
