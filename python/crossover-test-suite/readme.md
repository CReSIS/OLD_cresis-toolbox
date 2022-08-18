# Crossover Test Suite

Author: Reece Mathews

This directory holds a few scripts for comparing crossovers. These scripts were made to verify crossover calculation following the implementation of
segment geom simplification in the database.

## Dependencies

A GDAL whl can be obtained from [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/#gdal).
    Get the whl that matches your version of python and install it with
    `pip install ./GDAL-3.4.2-cp37-cp37m-win32.whl`
    replacing the file with the path to the one you downloaded.
Likewise, a Fiona whl can obtained from [here](https://www.lfd.uci.edu/~gohlke/pythonlibs/#fiona).
    Install it in the same manner.
Install the remaining dependencies with `pip install -r requirements.txt`

Map shape files can be obtained from [here](https://maps.princeton.edu/catalog/stanford-sd368wz2435)
GRL_adm0.shp and GRL_adm0.shx should be placed in the maps folder (create a maps folder under crossover-test-suite if not present).

## Setup

First you'll have to run the `data_queries.py` file against two separate databases (or the same DB before and after simplification). This process is described at the top of the `data_queries.py` file.

Once you have the crossover and segments CSVs for two separate databases you wish to compare, follow the directions in the `plot_crossovers.py` file to match crossovers together and produce a plot comparing shifts in location between the databases.

Loading large amounts of data can take several minutes but the UI seems to be fairly response once it's loaded.

## Usage

The `plot_crossovers.py` script attempts to match crossovers between the two database download targets and plot them one pair at a time.

You can click "Show Next CX" (or use the arrow keys) to iterate through each pair.

The "jump to crossover" box lets you enter a crossover number (ordered by greatest distance) to jump to. Press enter to submit the number.

"Toggle zoom" removes the axis limits and zooms out to a full view of the continent. This makes it easier to then zoom back in to the desired scale using Matplotlib's zoom tool. Press this button again to jump back in to the close-up view of the crossovers.

The "Toggle `TARGETA`" button hides the first database's segments and crossovers. This can be done for the other database's as well.

"Show All" zooms out and plots all segments and crossovers simultaneously. This can be slow.

If you have installed the `matplotlib_tools` module (*the install process is explained when running `plot_crossovers.py` without having first installed it*), you can drag on the plot to measure a distance in meters (distance is the `L` measurement in the top left of the plot.)

Matplotlib has a zoom tool accesible by pressing `O` and a pan tool from `P`.

## Notes

These tools were made with Greenland missions. To use Anatactica, a shapefile would have to be downloaded to replace the Greenland map and, in `plot_crossovers.py`, the `projection` variable would have to be swapped to Antartica's EPSG and the `plot_map` function would have to be altered to reflect these changes and plot the right map.
