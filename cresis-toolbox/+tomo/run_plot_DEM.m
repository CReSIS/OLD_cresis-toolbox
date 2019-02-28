% script run_plot_DEM
%
% Run script for tomo.plot_surf_geotiff
%
% Authors: Victor Berger, John Paden
%
% See also: plot_DEM.m

param.radar_name      = 'rds';
param.season_name     = '2014_Greenland_P3';
param.DEM_source      = 'DEM';
param.location        = 'arctic';
param.day_seg         = '20140401_03';
param.frm             = 33;
param.geotiff_fn      = ct_filename_gis(param,fullfile('canada','Landsat-7','Canada_90m.tif'));

global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Automated section
[fig_h, ax, cax] = tomo.plot_DEM(param,param_override, 'width', 800, 'height', 640,...
  'crossover', 'on', 'title', 'on')

