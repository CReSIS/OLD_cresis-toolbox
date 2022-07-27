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

import geopandas
import matplotlib.pyplot as plt

from shapely import wkb


# Default csv field size limit too small for segment geoms
# Set to max size supported by system
# https://stackoverflow.com/a/54517228
csv.field_size_limit(int(ctypes.c_ulong(-1).value // 2))

# Set working directory to crossover-test-suite
if Path(os.getcwd()).name == "cresis-toolbox":
    os.chdir("python/crossover-test-suite")


DATA_DIR = Path("data")


def load_data():
    """Load each CSV from the DATA_DIR."""
    data = {}

    for file in os.listdir(DATA_DIR):
        name = Path(file).stem
        with open(DATA_DIR / file, newline='') as f:
            data[name] = list(csv.DictReader(f))

    return data


def plot_map():
    """Make a base plot of the Greenland map."""
    map_df = geopandas.read_file("maps/GRL_adm0.shp")
    return map_df.boundary.plot(color='black')


def plot_geom(geom: str, base, color):
    
    geom_obj = wkb.loads(geom, hex=True)
    df = geopandas.GeoDataFrame(geometry=[geom_obj])
    df.plot(ax=base, color=color)


def plot_geoms():
    data = load_data()
    map_base = plot_map()
    for row in data["0m segments"]:
        plot_geom(row["geom"], map_base, "b")
    for row in data["1m segments"]:
        plot_geom(row["geom"], map_base, "r")

    plt.show()



if __name__ == "__main__":
    plot_geoms()
