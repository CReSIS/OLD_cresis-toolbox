% Script Information/Details
% =============================================
% Prepares CReSIS Data for ArcGIS Gridding
% Author: Kyle Purdon
% Version: 03/06/2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "gis data prep tool"
% Toolboxes Used: 
% Known Bugs: None
% Planned Updates: None
% Additional Information:
% =============================================

%% User Input

% Input File
% Must have fields:
%   'LAT','LON','UTCTime','THICK','ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY','SEASON'
% Created using geographic_search_gui.m
param.input = 'SomePath\SomeInputFile.csv';

% Data Value Removal Parameters
param.remove_high_elev = true;         % Remove Elevation > 10,000ft
param.remove_negative_thick = true;    % Remove Thickness < 0m
param.remove_derived = true;           % Remove Quality < 3 (Derived)

% Years of Data to Keep (Empty keeps all years)
% Ex. param.keep_years = {'2010','2011'};  - Keeps only 2010 or 2011 data
% Ex. param.keep_years = {};               - Keeps all years of data
param.keep_years = {'2006','2008','2009','2010','2011','2012'};

% Sync with NASA ATM Lidar Data
param.sync_atm = true;

% Split the FRAME variable.
% FRAME splits into YYYYMMDD,SEGMENT,FRAME
param.split_frame = true;

% Output file (Absolute path with filename and .csv/.txt extention)
param.output = 'SomePath\SomeOutputFile.txt';

%% Automated Section
gis_dataprep(param);
clear param;
