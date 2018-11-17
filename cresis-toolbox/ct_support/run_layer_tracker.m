% script run_layer_tracker
%
% Runs layer_tracker.m
%
% Author: John Paden

%% User Settings
% ----------------------------------------------------------------------
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2011_Greenland_P3.xls'),'','post');

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20110426');
params = ct_set_params(params,'cmd.frms',[]); % Specify specific frames (or leave empty/undefined to do all frames)

param_override.layer_tracker.debug_level = 1;

param_override.layer_tracker.echogram_img = 0; % To choose an image besides the base (0) image
% echogram_source: location of echogram data used for tracking
param_override.layer_tracker.echogram_source = 'CSARP_post/qlook';

% layer_params: structure of layer references of where to store the output
param_override.layer_tracker.layer_params = []; idx = 0;
% idx = idx + 1;
% param_override.layer_tracker.layer_params(idx).name = 'surface';
% param_override.layer_tracker.layer_params(idx).source = 'echogram';
% param_override.layer_tracker.layer_params(idx).echogram_source = 'deconv';
idx = idx + 1;
param_override.layer_tracker.layer_params(idx).name = 'surface';
param_override.layer_tracker.layer_params(idx).source = 'layerdata';
param_override.layer_tracker.layer_params(idx).layerdata_source = 'layerData';
param_override.layer_tracker.layer_params(idx).echogram_source = 'qlook';
% idx = idx + 1;
% param_override.layer_tracker.layer_params(idx).name = 'surface';
% param_override.layer_tracker.layer_params(idx).source = 'ops';

%% Enable one set of parameters
track_override = [];
track_override.en = true;
track_override.min_bin = 0;
track_override.manual = false;
switch ct_output_dir(params(1).radar_name)
  case 'rds'
    debug_time_guard = 2e-6;
    track_override.method = 'threshold';
    track_override.noise_rng = [0 -50 -10];
    track_override.min_bin = 2e-6;
    track_override.threshold = 15;
    track_override.sidelobe	= 15;
    track_override.max_diff	= inf;
    track_override.filter_len	= 5;
    track_override.medfilt = 3;
    track_override.search_rng	= 0:2;
    
  case 'accum'
    debug_time_guard = 2e-6;
    track_override.method = 'threshold';
    track_override.noise_rng = [200 -300 -100];
    track_override.min_bin = 0.1e-6;
    track_override.threshold = 9;
    track_override.sidelobe	= 12;
    track_override.max_diff	= inf;
    track_override.filter_len	= 5;
    track_override.search_rng	= 0:1;
    
  case {'snow','kuband','kaband'}
    % FMCW
    param_override.layer_tracker.debug_time_guard = 50e-9;
    track_override.method = 'threshold';
    track_override.threshold_noise_rng = [15e-9 -75e-9 -30e-9];
    track_override.threshold = 8;
    sidelobe = load('/cresis/snfs1/dataproducts/ct_data/ct_tmp/echogram_stats/snow/2011_Greenland_P3/stats_20110329_01.mat','sidelobe_rows','sidelobe_dB','sidelobe_vals');
    track_override.sidelobe_rows = [sidelobe.sidelobe_rows(1:194)];
    track_override.sidelobe_dB = [-sidelobe.sidelobe_dB(1:194,end)]-sidelobe.sidelobe_vals(end)+4.5;
    track_override.min_bin = 0.4e-6;
    track_override.prefilter_trim = [0 0];
    track_override.filter = [5 3];
    track_override.filter_trim = [10 10];
    track_override.max_rng	= [0 0.5e-9];
    track_override.detrend = '';
    if 0
      track_override.init.method	= '';
    elseif 0
      track_override.init.method	= 'dem';
      track_override.init.dem_offset = 39e-9;
      track_override.init.dem_layer.name = 'surface';
      track_override.init.dem_layer.source = 'lidar';
      track_override.init.dem_layer.lidar_source = 'atm';
    elseif 1
      track_override.init.method	= 'medfilt';
      track_override.init.medfilt	= 51;
    end
    track_override.init.max_diff = 0.5e-6;
    track_override.medfilt = 11;
    track_override.medfilt_threshold = 100;
    
end
param_override.layer_tracker.track = track_override;

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
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    layer_tracker(param,param_override);
  end
end
