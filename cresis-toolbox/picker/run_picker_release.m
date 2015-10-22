% MATLAB script run_picker_release
%
% This script shows how to use the image browser or picker.
%
% Steps to run:
%  1. Please set the base_path, gis_path, and image_browser_path variables
%     at the start of the script.
%  2. Copy CSARP_$PROC_TYPES and CSARP_layerData folders into the base_path
%     directory where $PROC_TYPES is all the processing types you want
%     to browse (e.g. csarp-combined, standard, qlook, and mvdr).
%  3. Copy any GIS files you would like to use into the gis_path variable.
%     If you use something other than the GIS files provided, then you
%     will need to add them to the geotiff_fns variable in this file.
%  4. If you plan to use the landmarks tool, make sure param.landmarks
%     points to a valid file.
%
% This tool does not require any toolboxes to run, but some of the functionality
% will be missing if you do no have the following toolboxes:
%  1. Mapping toolbox: you will not be able to load the geotiffs
%  2. Signal Processing toolbox: you will not be able to use the averaged
%     and detrended viewing options (only raw works)
%
% This is version 1.4 of the software. We have tested the software with Linux
% and Windows.
%
% Please report any problems to cresis_data@cresis.ku.edu
%
% Author: John Paden

% ==========================================================
% User Settings
% ==========================================================

if ispc
  base_path = 'Z:\mdce\mcords\2010_Greenland_DC8\';
  gis_path = 'V:\GIS_data\';
  image_browser_path = 'C:\Users\radar\Documents\scripts\image_browser_v1_4\';
else
  base_path = '/mnt/scratch1/output/mcords/2010_Greenland_P3/CSARP_post';
  gis_path = '/mnt/scratch1/GIS_data/';
  image_browser_path = '/mnt/scratch1/scripts/image_browser_v1_4/';
end

% Add the path (comment out this line if image_browser already added to the path)
addpath(genpath(image_browser_path));

% source_data = cell array of the location of the folders with the .mat files
%   which contain the echograms -- listing directories that do not exist
%   does not effect the operation
source_data      = {fullfile(base_path,'CSARP_standard'),
  fullfile(base_path,'CSARP_csarp-combined'),
  fullfile(base_path,'CSARP_mvdr'),
  fullfile(base_path,'CSARP_qlook')};

% layer_data = location of the .mat files which contain the layer information
layer_data       = fullfile(base_path,'CSARP_layerData');

% geotiff_fns = cell array of the location of each .tif files you want to be able to load
%   into the Map view
geotiff_fns         = {fullfile(gis_path,'greenland','Landsat-7','Greenland_natural_150m.tif'),
  fullfile(gis_path,'greenland','Landsat-7','Greenland_natural_90m.tif'),
  fullfile(gis_path,'antarctica','Landsat-7','Antarctica_LIMA.tif'),
  fullfile(gis_path,'antarctica','Landsat-7','Antarctica_LIMA_peninsula.tif')};

% param.source_data_mask = cell vector of strings which each frame ID
%   will be matched against to see if it should be included. A match occurs
%   when the frame ID contains the string. Leaving the
%   field undefined or empty causes all frames to be loaded.
% param.source_data_mask = {'20080702','20080706'};

% param.fast_load.en = enable fast loading
param.fast_load.en = false;

% param.fast_load.tmp_file = location of the temporary file for fast loading
param.fast_load.tmp_file = 'C:\tmp\picker_fast_load.mat';

% param.fast_load.recreate = true forces the file to be recreated every time
%   false uses the previously created file (software will create the file
%   if it does not exist)
%   If the source/layer files change, you MUST set this to true or delete
%   the previous tmp_file, otherwise the old (and incorrect) fast load
%   information will be loaded.
param.fast_load.recreate = true;

% param.landmarks = file path for location to store landmarks file (used with
%   the lankmarks tool)
param.landmarks = 'C:\tmp\picker_landmarks.mat';

% Run the picker
picker(source_data,layer_data,geotiff_fns,param);

return;

