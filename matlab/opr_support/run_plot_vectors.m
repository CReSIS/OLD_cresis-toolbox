% script run_plot_vectors
%
% Shows how to run plot_vectors. Copy and paste a snippet of code
% from below (select and hit F9) to plot the particular season
% of interest.
%
% Author: John Paden
%
% See also: plot_vectors.m

%% All files for a season
clear param;
param.radar_name = 'snow';
param.season_name = '2014_Greenland_P3';
vectors_path = ct_filename_support(param, '', 'vectors');
filenames = get_filenames(vectors_path,'vectors','','.mat');
geoTiff = fullfile(gRadar.gis_path,'arctic','NaturalEarth_Data','Arctic_NaturalEarth.tif');
% geoTiff = fullfile(gRadar.gis_path,'greenland','Landsat-7','mzl7geo_90m_lzw.tif');
% geoTiff = fullfile(gRadar.gis_path,'antarctica','Landsat-7','Antarctica_LIMA_480m.tif');
plot_vectors(filenames,[],geoTiff);

%% Single day for a season
clear param;
param.radar_name = 'snow';
param.season_name = '2014_Greenland_P3';
vectors_path = ct_filename_support(param, '', 'vectors');
filenames = get_filenames(vectors_path,'vectors_20140514','','.mat');
geoTiff = fullfile(gRadar.gis_path,'arctic','NaturalEarth_Data','Arctic_NaturalEarth.tif');
% geoTiff = fullfile(gRadar.gis_path,'greenland','Landsat-7','mzl7geo_90m_lzw.tif');
% geoTiff = fullfile(gRadar.gis_path,'antarctica','Landsat-7','Antarctica_LIMA_480m.tif');
plot_vectors(filenames,[],geoTiff);

%% Single segment for a season
clear param;
param.radar_name = 'mcords3';
param.season_name = '2017_Greenland_P3';
vectors_path = ct_filename_support(param, '', 'vectors');
filenames = get_filenames(vectors_path,'vectors_20170403_01','','.mat');
% geoTiff = fullfile(gRadar.gis_path,'arctic','NaturalEarth_Data','Arctic_NaturalEarth.tif');
geoTiff = fullfile(gRadar.gis_path,'greenland','Landsat-7','mzl7geo_90m_lzw.tif');
% geoTiff = fullfile(gRadar.gis_path,'antarctica','Landsat-7','Antarctica_LIMA_480m.tif');
plot_vectors(filenames,[],geoTiff);

