% script run_fmcw_set_nyquist_zone_Tsys_GPSoffset
%
% Runs fmcw_set_nyquist_zone_Tsys_GPSoffset.m
%
% Author: John Paden

%% User Settings
% param_fn = ct_filename_param('kuband_param_2009_Antarctica_DC8.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2009_Antarctica_DC8.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2010_Greenland_P3.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2010_Antarctica_DC8.xls'); max_nz = 2;
param_fn = ct_filename_param('snow_param_2011_Antarctica_DC8.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls'); max_nz = 2;
% param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls'); max_nz = 4;
% param_fn = ct_filename_param('snow_param_2015_Greenland_Polar6.xls'); max_nz = 4;

use_lidar_data = true; % Usually true
combine_elev_lidar_en = true; % Usually true
combine_elev_lidar_max_diff = 400e-9;
combine_surface_land_dems = true; % Usually true

debug_level = 0; % Set to zero to remove stops

save_records_en = true;

refine_Tsys_en = false;

%lidar_interp_gaps_dist = [1000 250];
lidar_interp_gaps_dist = [150 75];
lidar_source = 'atm';

radar_twtt_ratio = 1; % Usually one: this value will be multiplied onto the radar twtt
radar_twtt_offset = 0e-9; % Usually zero: this value will be added to the radar twtt

% =========================================================================
%% Automated Section
% =========================================================================

params = read_param_xls(param_fn,'','post');

fmcw_set_nyquist_zone_Tsys_GPSoffset;

return;
