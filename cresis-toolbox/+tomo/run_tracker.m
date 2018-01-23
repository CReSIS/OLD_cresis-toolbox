%% Script run_tracker
%
% Script for running tracker
%
% Authors: Victor Berger, John Paden
%
% See also: tracker.m

fprintf('\n\n========================================================\n');
fprintf('run tracker\n');
fprintf('========================================================\n');

%% User Settings

params             = read_param_xls(ct_filename_param('rds_param_2017_Greenland_P3.xls'),'','post');
% params             = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_DC8.xls'),'20091102_01','post');
% params.cmd.generic = 1;
% params.cmd.frms    = [];
options.name       = 'CSARP_post/mvdr';
options.debug      = false;
% options.name       = 'standard1';
geotiff_fn         = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
% geotiff_fn         = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif';
% geotiff2_fn        = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_rockmask.tif';

geotiff_fn  = ct_filename_gis([],geotiff_fn);
geotiff2_fn = '';
% geotiff2_fn = ct_filename_gis([],geotiff2_fn);

%% Automated section
global gRadar;
tomo.tracker(params, options, geotiff_fn, geotiff2_fn);
