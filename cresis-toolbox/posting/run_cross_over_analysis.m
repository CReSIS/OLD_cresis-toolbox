% Script Information/Details
% =============================================================
% Title: run_cross_over_analysis.m
% About: Calculates the crossover error for a CReSIS CSV file
%        and plots the resulting errors.
% Author: Brady Maasen, John Paden, Kyle Purdon
% Version Date: 02/21/2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: crossover analysis tool"
% Toolboxes Used: Standard, Mapping, Simulink
% Known Bugs: None
% Planned Updates: None
% Additional Information:
% see also plot_cross_overs.m cross_over_analysis.m
% =============================================================

%% User Input

% Filename that you want output .mat and .csv to be created under
% The .txt, .csv and .mat will be appended to this.
% DO NOT include an extention in this variable.
param.fout = 'C:\Users\SomeFolder\SomeOutputFile';

% Name of input file or path
param.fin = 'C:\Users\SomeFolder\SomeInputFile.csv';

% Specify a geotiff for plotting cross-overs
% These can be downloaded at:
%   ftp://data.cresis.ku.edu/data/picker/GIS_data/

% GREENLAND
%geotiff_fn = '/cresis/projects/GIS_data/greenland/Landsat-7/Greenland_natural_150m.tif';
geotiff_fn = 'P:\GIS_data\greenland\Landsat-7\Greenland_natural_150m.tif';

% CANADA
% geotiff_fn = '/cresis/projects/GIS_data/canada/Landsat-7/Canada_150m.tif';
% geotiff_fn = 'P:\GIS_data\canada\Landsat-7\Canada_150m.tif';

% ANTARCTICA
% geotiff_fn = '/cresis/projects/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
% geotiff_fn = 'P:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif';

%Param Structure for DataType (CReSIS CSV Default)
param.data_type.type = 4; 
%1=CSV 
%2=MATLAB 
%3=GISCSV(Split Frame, ATM Sync)
%4=GISTXT(GISCSV + Interpolation Processing + IceFree)

if param.data_type.type == 1
  param.data_type.format = '%f%f%f%f%f%s%f%f%f%s'; %Input String
  param.data_type.delim = ',';
  param.data_type.headerlines = 1;
  param.data_type.col_lat = 1;
  param.data_type.col_lon = 2;
  param.data_type.col_time =3;
  param.data_type.col_thickness = 4;
  param.data_type.col_frame = 6;

elseif param.data_type.type == 2
  param.data_type.segs = {}; % Leave empty to do all segments

elseif param.data_type.type == 3
  param.data_type.format = '%f%f%f%f%f%s%s%s%f%f%d%s%f%f%s'; %Input String
  param.data_type.delim = ',';
  param.data_type.headerlines = 1;
  param.data_type.col_lat = 1;
  param.data_type.col_lon = 2;
  param.data_type.col_time =3;
  param.data_type.col_thickness = 4;
  param.data_type.col_YYYYMMDD = 6;
  param.data_type.col_segment = 7;
  param.data_type.col_frame = 8;
  
  elseif param.data_type.type == 4
  param.data_type.format = '%f%f%f%s%s%f%f%f%f%s%s%s%f%f%d'; %Input String
  param.data_type.delim = ',';
  param.data_type.headerlines = 1;
  param.data_type.col_lat = 6;
  param.data_type.col_lon = 7;
  param.data_type.col_time =8;
  param.data_type.col_thickness = 3;
  param.data_type.col_YYYYMMDD = 10;
  param.data_type.col_segment = 11;
  param.data_type.col_frame = 12;
end

% --------------------------------------------------------
% If using CReSIS Standardized Data (With ~15m Sampling):
% The following parameters can be left to default values.
% --------------------------------------------------------

% Run points through a sorter which sorts them by flight line and
% in the order that the data was collected. This is useful if the
% data points you are loading are randomly sorted (e.g. an export
% from ARC will typically produce unsorted points). This is NOT
% needed when loading .mat layer files.
param.run_sort = false;

% The resolution of data (meters) after cross_over_analysis is run
% 20 or 25 meters is a good "standard". Increase this number for faster proc.
param.min_sample_spacing = 20;

% Return all cross overs that include data points with in dist meters
% where dist = (param.min_sample_spacing * param.scale_check)
param.scale_check = 1.5;

% Do not allow two flight lines two have two cross overs right next to
% each other. The minimum number of points in between cross overs
% is set by param.min_idx_sep
param.min_idx_sep = 20;

% Force cross overs to be at least this many seconds apart in time
param.UTC_check = 100;

% Flightline plot decimation (10 = Keep every 10th Point)
% End goal is ~ 250-300m sampling. (15m original, plot_samp = 20, 
% becomes 300m) Must be a postive integer. (Ex. 1 2 3 4 5 6 ...)
param.plot_samp = 10;

% Debug level (1 is default)
param.debug_level = 1;


%% Automated Section

% Call Run Function
fprintf('Running cross_over_analysis.m\n');
cross_over_analysis(param);

% Call Plot Function
fprintf('Running plot_cross_overs.m\n');
plot_cross_overs(param.fout,geotiff_fn,param);

