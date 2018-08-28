% script run_update_surface_with_tracker
%
% Runs update_surface_with_tracker.m
%
% Author: John Paden

%% User Settings
% ----------------------------------------------------------------------
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'','post');
% params = read_param_xls(ct_filename_param('rds_param_2017_Greenland_P3.xls'),'','post');
% params = read_param_xls(ct_filename_param('snow_param_2017_Greenland_P3.xls'),'','post');

% params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.frms',[]); % Specify specific frames (or leave empty/undefined to do all frames)


param_override.update_surface.debug_level = 1;
param_override.update_surface.echogram_img = 0; % To choose an image besides the base (0) image
% echogram_source: location of echogram data used for tracking
param_override.update_surface.echogram_source = 'qlook';

% layer_params: structure of layer references of where to store the output
param_override.update_surface.layer_params = []; idx = 0;
% idx = idx + 1;
% param_override.update_surface.layer_params(idx).name = 'surface';
% param_override.update_surface.layer_params(idx).source = 'echogram';
% param_override.update_surface.layer_params(idx).echogram_source = 'qlook';
idx = idx + 1;
param_override.update_surface.layer_params(idx).name = 'surface';
param_override.update_surface.layer_params(idx).source = 'layerdata';
param_override.update_surface.layer_params(idx).layerdata_source = 'layerData';
param_override.update_surface.layer_params(idx).echogram_source = 'qlook';
% idx = idx + 1;
% param_override.update_surface.layer_params(idx).name = 'surface';
% param_override.update_surface.layer_params(idx).source = 'ops';


surf_override = [];
surf_override.en = true;
surf_override.min_bin = 0;
surf_override.manual = false;

%% Enable one set of parameters
if 0
  % RDS threshold method
  debug_time_guard = 2e-6;
  surf_override.method = 'threshold';
  surf_override.noise_rng = [0 -50 -10];
  surf_override.min_bin = 2e-6;
  surf_override.threshold = 15;
  surf_override.sidelobe	= 15;
  surf_override.max_diff	= inf;
  surf_override.filter_len	= 5;
  surf_override.medfilt = 3;
  surf_override.search_rng	= 0:2;
%   surf_override.init.method	= 'snake';
%   surf_override.init.search_rng	= [-10:10];
%   surf_override.feedthru.time = 1e-5*[-0.0084    0.0471    0.0768    0.1342    0.1749    0.2039];
%   surf_override.feedthru.power_dB = [-95.7021  -69.3129  -39.6250  -40.4855  -58.8433  -96.5627];
elseif 0
  % RDS max method
  debug_time_guard = 2e-6;
  surf_override.method = 'max';
  surf_override.min_bin = 0.75e-6;
  surf_override.search_rng	= [-11:0];
  surf_override.threshold = 13;
  surf_override.medfilt = 3;
  % FEEDTHRU IS SEASON DEPENDENT (this example is for RDS 2009 Greenland TO)
  surf_override.feedthru.time = 1e-9*[102;1069;1803;2436;3003;3603];
  surf_override.feedthru.power_dB = [-50.9;-82.0;-106.3;-128.9;-131.7;-152.1];
  
elseif 1
  % Accum
  debug_time_guard = 2e-6;
  surf_override.method = 'threshold';
  surf_override.noise_rng = [200 -300 -100];
  surf_override.min_bin = 0e-6;
  surf_override.threshold = 9;
  surf_override.sidelobe	= 12;
  surf_override.max_diff	= inf;
  surf_override.filter_len	= 5;
  surf_override.search_rng	= 0:1;
  
elseif 0
  % FMCW Land Ice
  debug_time_guard = 50e-9;
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
%   surf_override.init.method	= 'dem';
%   surf_override.init.dem_offset = 0e-9;
%   surf_override.init.method	= 'snake';
%   surf_override.init.search_rng	= [-240:240];
  surf_override.method = 'snake';
  surf_override.search_rng	= -90:90;
  
elseif 1
  % FMCW Sea Ice
  debug_time_guard = 50e-9;
  surf_override.method = 'threshold';
  surf_override.method = '';
  surf_override.noise_rng = [100 -400 -100];
  surf_override.threshold = 13;
  surf_override.sidelobe	= 25;
  surf_override.max_diff	= 45e-9;
  surf_override.filter_len	= 7;
  surf_override.search_rng	= 0:10;
  surf_override.detrend = 0;
  surf_override.init.method	= 'dem';
  surf_override.init.dem_offset = 67e-9;
%   surf_override.init.lidar_source = 'atm';
%   surf_override.init.method	= 'medfilt';
%   surf_override.init.medfilt	= 51;
end
param_override.update_surface.surf = surf_override;

%% Automated Section
% ----------------------------------------------------------------------
  

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    %update_surface_with_tracker(param,param_override);
    update_surface_with_tracker;
  end
end
