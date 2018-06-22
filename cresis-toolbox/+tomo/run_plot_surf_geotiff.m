% script run_plot_surf_geotiff
%
% Run script for tomo.plot_surf_geotiff
%
% Authors: Victor Berger, John Paden
%
% See also: plot_surf_geotiff.m

param.radar_name      = 'rds';
param.season_name     = '2014_Greenland_P3';
param.DEM_source      = 'DEM';
param.location        = 'arctic';
param.day_seg         = '20140325_06';
param.frm             = 1;
param.geotiff_fn      = ct_filename_gis(param,fullfile('canada','Landsat-7','Canada_90m.tif'));

global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

%% Automated section
[fig_h, ax, cax] = plot_surf_geotiff(param,param_override);

