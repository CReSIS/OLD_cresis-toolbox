function track = layer_tracker_profile(param,profile_str)
% track = layer_tracker_profile(param,profile_str)
%
% Returns a default set of tracking parameters for layer tracking with
% layer_tracker.m.
%
% param: parameter spreadsheet structure
% profile_str: optional string containing the profile name to load (default
%   is "default"
%
% track: Parameters used to control the tracker. See layer_tracker.m.
%
% Example:
%   Should only be used from layer_tracker.m
%
% Author: John Paden

%% Check input arguments
if ~exist('profile_str','var') || isempty(profile_str)
  profile_str = 'default';
end
track.profile = profile_str;

%% Default track settings
track.data_noise_en = false;
track.detrend = 0;
track.en = true;
track.init = [];
track.init.method = 'max';
track.init.snake_rng = [-2e-7 2e-7];
track.init.dem_layer = '';
track.init.max_diff = inf;
track.init.max_diff_method = 'interp_finite';
track.filter_mocomp = false;
track.filter = [1 1];
track.filter_trim = [0 0];
track.fixed_value = 0;
track.min_bin = 0;
track.max_bin = inf;
track.max_rng = [0 0];
track.max_rng_units = 'time';
track.medfilt = 0;
track.medfilt_threshold = 0;
track.method = 'threshold';
track.prefilter_trim = [0 0];
track.sidelobe_dB = [];
track.sidelobe_rows = [];
track.snake_rng = [-2e-7 2e-7];
track.threshold = 15;
track.threshold_noise_rng = [0 -inf -1];

if strcmpi(profile_str,'default')
  %% Default profile
  
elseif strcmpi(profile_str,'RDS_OIB')
  %% RDS_OIB profile
  track.debug_time_guard = 2e-6;
  track.filter	= [1 5];
  track.init.method	= 'medfilt';
  track.init.medfilt	= 11;
  track.init.max_diff = 0.5e-6;
  track.max_rng	= [0 1];
  track.max_rng_units = 'bins';
  track.medfilt = 11;
  track.medfilt_threshold = 30;
  track.method = 'threshold';
  track.min_bin = 1.6e-6;
  track.threshold = 10;
  track.threshold_noise_rng = [0 -2e-6 -0.2e-6];
  track.threshold_rel_max = -9;
  track.threshold_rng = 5;
  
elseif strcmpi(profile_str,'ACCUM')
  %% ACCUM profile
  track.debug_time_guard = 2e-6;
  track.filter	= [3 3];
  track.filter_trim = [3 3];
  track.init.method	= 'medfilt';
  track.init.medfilt	= 11;
  track.init.max_diff = 0.5e-6;
  track.max_rng	= [0 2];
  track.max_rng_units = 'bins';
  track.medfilt = 11;
  track.medfilt_threshold = 30;
  track.method = 'threshold';
  track.min_bin = 1.6e-6;
  track.threshold = 10;
  track.threshold_noise_rng = [0 -1e-6 -0.3e-6];
  track.threshold_rel_max = -9;
  track.threshold_rng = 5;
  
elseif strcmpi(profile_str,'SNOW_AWI')
  %% SNOW_AWI profile
  track.debug_time_guard = 50e-9;
  track.min_bin = 0.1e-6;
  track.prefilter_trim = [0 0];
  track.filter = [5 3];
  track.filter_trim = [10 10];
  track.init.method	= 'medfilt';
  track.init.medfilt	= 51;
  track.init.max_diff = 0.3e-6;
  track.max_rng	= [0 9];
  track.max_rng_units = 'bins';
  track.medfilt = 11;
  track.medfilt_threshold = 100;
  track.method = 'threshold';
  track.threshold = 8;
  track.threshold_noise_rng = [15e-9 -75e-9 -30e-9];
  track.threshold_rel_max = -9;
  
elseif strcmpi(profile_str,'DEM_LIDAR')
  %% DEM profile
  track.debug_time_guard = 100e-9;
  track.min_bin = 0e-6;
  track.init.method  = 'dem';
  track.init.dem_offset = 6.5e-9;
  track.init.dem_layer.name = 'surface';
  track.init.dem_layer.source = 'lidar';
  track.init.dem_layer.lidar_source = 'atm';
  track.init.max_diff_method = 'merge_vectors';
  track.init.max_diff = 0e-6; % Force output layer to equal DEM
  track.max_rng	= [0 0];
  track.method = ''; % No data dependent tracking
  
else
  error('Invalid profile selected: %s\n', profile_str);
end
