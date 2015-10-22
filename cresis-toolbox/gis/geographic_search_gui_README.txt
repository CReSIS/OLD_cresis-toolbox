
CURRENT UPDATE DATE: 05/31/2012

- If you find any bugs please notify us at cresis_data@cresis.ku.edu (Subject: GeoSearch)
- Please see the "Planned Updates" under "Current Version" to see if we are aware of your issue.

run_geographic_search_gui.m

  Interactive Tool for downloading level-II CReSIS data.

YOU MUST HAVE geographic_search_gui.m in the same directory.

DEPENDS ON THE FOLLOWING FUNCTIONS:
geographic_search_gui.m
geodetic_to_along_track.m
get_equal_alongtrack_spacing_idxs.m

YOU MUST DOWNLOAD AND STORE ALL OF THE ABOVE FUNCTIONS.

Known Bugs:

  MAY NOT BE AN ISSUE WITH VERSION v3+
  You may need to set the Java Heap size property in Matlab to a larger value
  (e.g. 256 MB or more) to download large files. The error message from
  urlread when Java runs out of memory indicates that the URL is bad or
  the connection is bad (i.e. the error message for too small of a heap
  size is misleading).

TOOLBOXES USED:
-mapping
-stats
-images

How to Use:

(1) Download BOTH files (.m files)
(2) Open the run_**** .m file
(2) READ the comments in the user input section of the script and fill out the values


Note: If you have issues contact us at cresis_data@cresis.ku.edu (Subject: GeoSearch).

_______________________________

Datasets Included for FTP Data

 ** Not Included in this version 
_________________________________________

Greenland and NE Canada

1993_Greenland_P3
1995_Greenland_P3
1996_Greenland_P3
1997_Greenland_P3
1998_Greenland_P3
1999_Greenland_P3
2001_Greenland_P3
2002_Greenland_P3
2003_Greenland_P3
2005_Greenland_TO
2006_Greenland_TO
2007_Greenland_P3
2008_Greenland_TO
** 2008_Greenland_Ground
2009_Greenland_TO
2010_Greenland_DC8
2010_Greenland_P3
2011_Greenland_P3
2011_Greenland_TO


Antarctica

2002_Antarctica_P3chile
2004_Antarctica_P3chile
2009_Antarctica_DC8
2009_Antarctica_TO
2010_Antarctica_DC8

_________________________________________

Datasets Included for CReSIS Local Data

 ** Not Included in this version 
_________________________________________

Greenland and NE Canada

1993_Greenland_P3
1995_Greenland_P3
1996_Greenland_P3
1997_Greenland_P3
1998_Greenland_P3
1999_Greenland_P3
2001_Greenland_P3
2002_Greenland_P3
2003_Greenland_P3
2005_Greenland_TO
2006_Greenland_TO
2007_Greenland_P3
2008_Greenland_TO
2008_Greenland_Ground
2009_Greenland_TO
2010_Greenland_DC8
2010_Greenland_P3
2011_Greenland_P3
2011_Greenland_TO


Antarctica

2002_Antarctica_P3chile
2004_Antarctica_P3chile
2009_Antarctica_DC8
2009_Antarctica_TO
2010_Antarctica_DC8