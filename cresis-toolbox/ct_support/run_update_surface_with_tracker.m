% script run_update_surface_with_tracker
%
% Runs update_surface_with_tracker.m
%
% Author: John Paden

%% User Settings
% ----------------------------------------------------------------------
% params = read_param_xls(ct_filename_param('snow_param_2009_Greenland_P3.xls'),'20090421_04','post');
% params = read_param_xls(ct_filename_param('snow_param_2010_Greenland_DC8.xls'),'20100326_01','post');
% params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),'20110325_02','post');
% params = read_param_xls(ct_filename_param('snow_param_2012_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2015_Greenland_Polar6.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'),'20150413_01','post');
params = read_param_xls(ct_filename_param('rds_param_2015_Greenland_C130.xls'),'20150413_01','post');

params.cmd.generic = 1;
params.cmd.frms = [1]; % Specify specific frames (or leave empty/undefined to do all frames)
debug_level = 1;

echogram_source = 'qlook'; % location of echogram data used for tracking (and for saving)
% save_sources: cell array of strings indicating which layer sources
%   should be updated (options are ops, layerData, and echogram)
save_sources = {'ops'};
save_ops_layer = 'surface'; % OPS save_source only (layer name)
layerdata_source = 'layerData'; % layerData save_source only (ct_filename_location parameter)

surf_override = [];
surf_override.en = true;
surf_override.min_bin = 0;
surf_override.manual = false;

%% Enable one set of parameters
if 1
  % RDS threshold method
  surf_override.method = 'threshold';
  surf_override.noise_rng = [0 -50 -10];
  surf_override.min_bin = 1.9e-6;
  surf_override.threshold = 9;
  surf_override.sidelobe	= 15;
  surf_override.max_diff	= 100e-9;
  surf_override.filter_len	= 1;
  surf_override.medfilt = 3;
  surf_override.search_rng	= 0:2;
%   surf_override.init.method	= 'snake';
%   surf_override.init.search_rng	= [-10:10];
  surf_override.feedthru.time = 1e-9*[-75390765878 1115131099362];
  surf_override.feedthru.power_dB = [-74 -139];
elseif 0
  % RDS max method
  surf_override.method = 'max';
  surf_override.min_bin = 0.75e-6;
  surf_override.search_rng	= [-11:0];
  surf_override.threshold = 13;
  surf_override.medfilt = 3;
  % FEEDTHRU IS SEASON DEPENDENT (this example is for RDS 2009 Greenland TO)
  surf_override.feedthru.time = 1e-9*[102;1069;1803;2436;3003;3603];
  surf_override.feedthru.power_dB = [-50.9;-82.0;-106.3;-128.9;-131.7;-152.1];
  
elseif 0
  % Accum
  surf_override.method = 'threshold';
  surf_override.noise_rng = [200 -300 -100];
  surf_override.min_bin = 0.1e-6;
  surf_override.threshold = 9;
  surf_override.sidelobe	= 12;
  surf_override.max_diff	= inf;
  surf_override.filter_len	= 5;
  surf_override.search_rng	= 0:1;
  
elseif 0
  % FMCW Land Ice
  surf_override.method = 'threshold';
  surf_override.noise_rng = [100 -700 -400];
  surf_override.noise_override = 1;
  surf_override.min_bin = 2e-6;
  surf_override.threshold = 5;
  surf_override.sidelobe	= 13;
  surf_override.max_diff	= 75e-9;
  surf_override.filter_len	= [3 13];
  surf_override.search_rng	= 0:2;
  surf_override.detrend = 0;
  surf_override.init.method	= 'medfilt';
  surf_override.init.medfilt	= 11;
  surf_override.init.method	= 'dem';
  surf_override.init.dem_offset = 0e-9;
  surf_override.init.method	= 'snake';
  surf_override.init.search_rng	= [-240:240];
  surf_override.method = 'snake';
  surf_override.search_rng	= -90:90;
  
elseif 1
  % FMCW Sea Ice
  surf_override.method = 'threshold';
  surf_override.noise_rng = [100 -400 -100];
  surf_override.threshold = 9;
  surf_override.sidelobe	= 13;
  surf_override.max_diff	= 45e-9;
  surf_override.filter_len	= 7;
  surf_override.search_rng	= 0:10;
  surf_override.detrend = 2;
%   surf_override.init.method	= 'dem';
%   surf_override.init.dem_offset = 0e-9;
  surf_override.init.method	= 'medfilt';
  surf_override.init.medfilt	= 51;
end

%% Automated Section
% ----------------------------------------------------------------------
  
update_surface_with_tracker;

return
