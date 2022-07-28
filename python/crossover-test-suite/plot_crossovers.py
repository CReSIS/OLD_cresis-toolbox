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


def plot_geoms(geoms: List[str], base, color):
    """Plot a list of geoms string on the map base with the given color."""

    # Pandas plot method returns the base axes rather than the object plotted like Matplotlib does
    # Keep note of children of base axes before plotting so that we can find the new child after plotting
    old_children = set(base.get_children())
    df = geopandas.GeoDataFrame(geometry=geoms)
    ax = df.plot(ax=base, color=color)

    # Find the child we just plotted to return it
    plotted_child = list(set(ax.get_children()) - old_children)[0]

    return plotted_child


class State():
    showing_0m = True
    showing_1m = True

    def __init__(self, plot_seg_0m, plot_seg_1m, plot_cx_0m, plot_cx_1m):
        self.plot_seg_0m = plot_seg_0m
        self.plot_seg_1m = plot_seg_1m
        self.plot_cx_0m = plot_cx_0m
        self.plot_cx_1m = plot_cx_1m

    def toggle_0m(self, clk):
        self.showing_0m = not self.showing_0m
        self.plot_seg_0m.set_visible(self.showing_0m)
        self.plot_cx_0m.set_visible(self.showing_0m)
        plt.draw()

    def toggle_1m(self, clk):
        self.showing_1m = not self.showing_1m
        self.plot_seg_1m.set_visible(self.showing_1m)
        self.plot_cx_1m.set_visible(self.showing_1m)
        plt.draw()


def plot_from_data():
    mpl.style.use("seaborn")

    data = load_data()
    map_base = plot_map()

    if load_ruler:
        # Monkey-patch annotate method to handle erroneous use of "s" param instead of "text" by the matplotlib_tools module
        original_annotate = map_base.annotate
        map_base.annotate = lambda *args, **kwargs: original_annotate(text=kwargs["s"], *args, **{k: v for k, v in kwargs.items() if k != "s"}) if "s" in kwargs else original_annotate(*args, **kwargs)
        # Ruler is immediately garbage collected if not assigned to anything
        # Matplotlib widgets are maintained with only weakrefs internally. The user is expected to maintain their own strong reference
        ruler = Ruler(ax=map_base, useblit=True, lineprops={"color": "red", "linewidth": "2"})

    map_base.set_xlabel('Meters East of North pole')
    map_base.set_ylabel('Meters North of North pole')

    plot_seg_1m = plot_geoms([row["geom"] for row in data["1m segments"]], map_base, "C0")
    plot_cx_1m = plot_geoms([row["cx_geom"] for row in data["1m crossovers"]], map_base, "C1")
    plot_seg_0m = plot_geoms([row["geom"] for row in data["0m segments"]], map_base, "C2")
    plot_cx_0m = plot_geoms([row["cx_geom"] for row in data["0m crossovers"]], map_base, "C3")

    # Create visibility toggle buttons
    state = State(plot_seg_0m, plot_seg_1m, plot_cx_0m, plot_cx_1m)

    ax_toggle_0 = plt.axes((0.9, 0.8, 0.1, 0.075))
    ax_toggle_1 = plt.axes((0.9, 0.7, 0.1, 0.075))

    b_toggle_0 = Button(ax_toggle_0, 'Toggle Full Res')
    b_toggle_0.on_clicked(state.toggle_0m)
    b_toggle_1 = Button(ax_toggle_1, 'Toggle 1m')
    b_toggle_1.on_clicked(state.toggle_1m)

    map_base.legend([plot_seg_1m, plot_cx_1m, plot_seg_0m, plot_cx_0m], 
                    ['1m Segments', '1m Crossovers', 'Full Res Segments', 'Full Res Crossovers',])

    plt.show()


# TODO[Reece]: Find large differences in crossovers and circle them on map
# Or only plot them maybe


if __name__ == "__main__":
    plot_from_data()
