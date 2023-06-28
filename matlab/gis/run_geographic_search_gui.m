% Script Information/Details
% =============================================
% Interactive Geographic Search/Download
% Author: Kyle Purdon
% Contributors: John Paden, Aric Beaver
% Version: V3.7 01/05/2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "Geographic Search Tool"
% Toolboxes Used: Standard, Mapping, Stats, Images
% Known Bugs: None
% Planned Updates: None
% Additional Information: See M-File
% see also geographic_search_gui.m
% ==================================================

% Datasets Included for FTP Data
%  ** Not Included in this version
% Greenland and NE Canada
% 1993_Greenland_P3
% 1995_Greenland_P3
% 1996_Greenland_P3
% 1997_Greenland_P3
% 1998_Greenland_P3
% 1999_Greenland_P3
% 2001_Greenland_P3
% 2002_Greenland_P3
% 2003_Greenland_P3
% 2005_Greenland_TO
% 2006_Greenland_TO
% 2007_Greenland_P3
% ** 2008_Greenland_TO
% ** 2008_Greenland_Ground
% ** 2009_Greenland_TO
% 2010_Greenland_DC8
% 2010_Greenland_P3
% 2011_Greenland_P3
% Antarctica
% 2002_Antarctica_P3chile
% 2004_Antarctica_P3chile
% 2009_Antarctica_DC8
% 2009_Antarctica_TO
% 2010_Antarctica_DC8
% ----------------------------------------------------
% Datasets Included for CReSIS Local Data
%  ** Not Included in this version
% Greenland and NE Canada
% 1993_Greenland_P3
% 1995_Greenland_P3
% 1996_Greenland_P3
% 1997_Greenland_P3
% 1998_Greenland_P3
% 1999_Greenland_P3
% 2001_Greenland_P3
% 2002_Greenland_P3
% 2003_Greenland_P3
% 2005_Greenland_TO
% 2006_Greenland_TO
% 2007_Greenland_P3
% 2008_Greenland_TO
% 2008_Greenland_Ground
% 2009_Greenland_TO
% 2010_Greenland_DC8
% 2010_Greenland_P3
% 2011_Greenland_P3
% Antarctica
% 2002_Antarctica_P3chile
% 2004_Antarctica_P3chile
% 2009_Antarctica_DC8
% 2009_Antarctica_TO
% 2010_Antarctica_DC8

% Accum/Snow/KUBand Data Included
% 2011_Greenland_P3

% see also geographic_search_gui.m
% =============================================

%% GeoSearch Directions (How to use the script)

% ----- STEP 1 -----
% You have 2 options when running this script
%    (1) Use Local CReSIS Data
%    (2) Use FTP CReSIS Data
% To use local CReSIS Data you must be on a CReSIS machine.
% If you are not within the CReSIS network you must use the FTP Data.


% ----- STEP 2 -----
% You have 2 options when running this script
%    (1) Enter all the variables below in "User Input"
%    (2) Respond to "prompts" in the command window
% If you set runType to 2 press F5 now and follow the
% prompts in the command window. Otherwise go to STEP 3.

% Read the instructions here and fill out
% your choices in "User Input" below.

% ----- STEP 3 -----
% Set your location of interest. This should be:
%    (1) Greenland
%    (2) Antarctica
%    (3) Canada

% ----- STEP 4 -----
% Set your radar of interest. This should be:
%    (1) RDS
%    (2) Snow
%    (3) Accum
%    (4) KUBand
%
% If you set this value to "1" RDS, go to STEP 5
% IF you set this value to "2-4" go to STEP 11a
% **No Download Support for Non-RDS Radars

% ----- STEP 5 -----
% Do you want the files in your StudyArea to be downloaded?
% Answer can be "True" OR "False"
% If you set this value to "True" go to STEP 6.
% If you set this value to "False" go to STEP 11.

% ----- STEP 6 -----
% Set the folder your files will be saved to. Can be a relative or absolute
% path. Can be LINUX or WINDOWS path. See Examples below.
% Current directory (download_dir = '.';)
% Other Path (download_dir = 'C:\Users\SomeFolder\Output';)

% ----- STEP 7 -----
% Specify the output filename (No Spaces)
% param.file_name = 'Someglacier_Output';

% ----- STEP 8 -----
% Set the type of file to output
% CSV File (param.extID = 1;)
% MAT File (param.extID = 2;)
% BOTH Files (param.extID = 3;)

% ----- STEP 9 -----
% Set the type of files to download.
%    (1) csv_good
%    (2) csv
% csv: all data cresis has (Including bad data)
% csv_good: all data cresis has with a surface and bottom value.
% Default: csv_good (See Examples)
% CSV_GOOD (download_tid = 1;)
% CSV (download_tid = 2;)

% ----- STEP 10a -----
% Do you want the download files to be downsampled?
% Makes the file easier to use in ArcGIS and other programs
% Answer can be "True" OR "False"
% If you set this value to "True" go to STEP 10b.
% If you set this value to "False" go to STEP 11.

% ----- STEP 10b -----
% Set the factor to downsample to. (meters)
% (Should be between 300-700m)
% 500m (dsamp_meters = 500;)

% ----- STEP 11a -----
% Set option to use a default or custom GeoTIFF
% Answer can be:
%    (1) Use the Default GeoTiff (See list below)
%    (2) Use a custom GeoTIff
% If you answer (1)STEP 12  or (2)STEP 11b

% DEFAULT IMAGES:
% Greenland: Landsat-7 Greenland_natural_150m.tif
% Canada: Landsat-7 Canada_250m.tif
% Antarctica: Landsat-7 Antarcitca_LIMA_peninsula.tif

% ----- STEP 11b -----
% Specify the path to a custom GeoTIFF Image
% Example (geotiff_fn = 'C:\Users\SomeFolder\SomeImage.tif';)
% MUST be a GeoTiff Image

% ----- STEP 12 -----
% Specify the base path for TEMP files (Must be able to hold ~200mb - 2gb)
% The files will be removed after the script runs.
% Example (param.tPath = 'C:\Users\AppData\Temp\';)

% ----- STEP 13 -----
% Make sure you have entered all the parameters correctly.
% Press F5 to save and run the script.
% Follow any on-screen prompts.

%% User Input

% Please read the above instructions before filling out the user input.

%Type of Data to Get (1)=Local (2)=FTP
param.dataGetType = 1;
%Interactive or Script Input Mode (1)=Script (2)=Interactive
param.runType = 2;  
%Location of interest (1)=Greenland (2)=Antarctica (3)=Canada
param.location_en = 1;
%Radar of Interest (1)=RDS, (2)=Snow or (3)=Accum (4)=KUBand
param.radarID = 1;
%Download files?
param.download_en = true;
%Set download output path (No Spaces)
param.download_dir = 'C:\Users\SomeFolder';
%Set download output file name (With NO extension) (No Spaces)
param.file_name = 'FileOutputName';
%Set the output file type(s) (1)=CSV (2)=MAT (3)=BOTH
param.extID = 1;
%Set the type of files to download (1)=CSV_GOOD (2)=CSV
param.download_tid = 1;
% Downsample files?
param.down_sample = true;
%Downsampling factor (units in meters)
param.dsamp_meters = 500;
%Set GeoTiff Option (1)=Default (2)=Custom
param.image_en = 1;
%Custom GeoTiff Path (No Spaces) (ONLY IF param.image_en = 2)
param.geotiff_fn = 'C:\Users\SomeFolder\SomeImage.tif';
%Temporary Base Path (No Spaces)
param.tPath = 'C:\Users\AppData\Temp\';

%% Automated Section

geographic_search_gui(param);

