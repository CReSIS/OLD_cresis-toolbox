% script run_picker
%
% Sets up parameters for running the picker. Make a copy in your local
% directory.
%
% Author: John Paden
%
% See also picker.m

clear param;
param.radar_name = 'kuband';
param.season_name = '2011_Greenland_P3';

% Set to empty if not in the post directory, otherwise set to 'CSARP_post'
post_dir = '';

source_data = {};
source_data{end+1} = fullfile(ct_filename_out(param,'','',1),post_dir,'CSARP_qlook');

param.source_data_mask = {'20110406'};

layer_data = fullfile(ct_filename_out(param,'','',1),post_dir,'CSARP_layerData');

geotiff_fns = {};
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','mzl7geo_90m_lzw.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','mzl7geo_45m_lzw.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','Greenland_natural_250m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','Greenland_natural_150m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','Greenland_natural_90m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'canada'),'Landsat-7','Canada_250m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'canada'),'Landsat-7','Canada_150m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'canada'),'Landsat-7','Canada_90m.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'antarctica'),'Landsat-7','Antarctica_LIMA.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'antarctica'),'Landsat-7','Antarctica_LIMA_480.tif');
geotiff_fns{end+1} = fullfile(ct_filename_gis(param,'antarctica'),'Landsat-7','Antarctica_LIMA_peninsula.tif');

param.fast_load.en = false;
param.fast_load.recreate = true;
param.fast_load.tmp_file = '';

param.landmarks = ct_filename_tmp(param,'','picker','landmarks.mat');

picker(source_data,layer_data,geotiff_fns,param);

return;
