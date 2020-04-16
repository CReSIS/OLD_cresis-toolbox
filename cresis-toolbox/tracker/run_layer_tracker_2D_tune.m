% script run_layer_tracker_2D
%
% Runs layer_tracker_2D.m

%% User Settings
% ----------------------------------------------------------------------
param_override = [];

% params = read_param_xls(ct_filename_param('rds_param_2011_Greenland_P3.xls'));
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));
% params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'));

params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20140313_08');
params = ct_set_params(params,'cmd.frms',[1 2]); % Specify specific frames (or leave empty/undefined to do all frames)
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20110331_02');
% params = ct_set_params(params,'cmd.frms',19); % Specify specific frames (or leave empty/undefined to do all frames)

param_override.layer_tracker.debug_plots = {'tracked_images'};
% param_override.layer_tracker.debug_plots = {'tracked_images','visible'}; % Uncomment for debugging

param_override.layer_tracker.echogram_img = 0; % To choose an image besides the base (0) image
% echogram_source: location of echogram data used for tracking
% param_override.layer_tracker.echogram_source = 'CSARP_post/qlook';
% param_override.layer_tracker.echogram_source = 'CSARP_post/mvdr';
param_override.layer_tracker.echogram_source = 'CSARP_post/standard';

% layer_params: layerparams structure of where to store the output using
% opsCopyLayers.m
param_override.layer_tracker.layer_params = [];
% Uncomment to enable layerdata storage
param_override.layer_tracker.layer_params.layerdata_source = 'layer_paden';
% Uncomment to enable OPS storage
% param_override.layer_tracker.layer_params.source = 'ops';

% block_size_frms: Number of frames to be loaded at a time
param_override.layer_tracker.block_size_frms = 1;

% track_per_task: Number of tracks per task
param_override.layer_tracker.track_per_task = 1;

%% param.layer_tracker.track options
track_idx = 0;
for y = 160:20:300
  for dy = [5 10 20 40]
    track = [];
    
    %% Enable one set of parameters
    track.en = true;
    % RDS
    track.profile = 'rds_OIB';
    
    %% LSM User Settings
    track.method           = 'lsm';
    track.lsm.lyrtop       = 'lsm_top'; %layername, layer_dest.name
    track.lsm.lyrbot       = 'lsm_bot';
    track.lsm.y            = y; % = '' for y = mean(SURF)
    track.lsm.dy           = dy;
    track.lsm.storeIter    = [25:25:400];
    track.init.max_diff    = inf;
    track.detrend          = [];
    track.norm.scale       = [-40 90];
    
    track_idx = track_idx + 1;
    param_override.layer_tracker.track{track_idx} = track;
  end
end

% param_override.layer_tracker.surf_layer = struct('name','surface','source','layerdata','layerdata_source','layer');

% param_override.layer_tracker.crossover_layer = struct('name','bottom','source','ops');

% dbstop if error;
% param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
param_override.cluster.type = 'debug';
% param_override.cluster.type = 'slurm';
% param_override.cluster.rerun_only = true;
% param_override.cluster.desired_time_per_job  = 240*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

%% Automated Section
% ----------------------------------------------------------------------

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

ctrl_chain = {};
% Process each of the segments
for param_idx = 1:length(params)
  param = params(param_idx);
  if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
    [ctrl_chain{end+1},param] = layer_tracker_2D(param,param_override);
    % Since we are tuning, save the parameters
    out_param_fn = [ct_filename_ct_tmp(param,'','layer_tracker','') ...
      sprintf('_%s_t%03d_%s.mat',datestr(now,'yyyymmdd_HHMMSS'), ...
        length(param.layer_tracker.track), param.layer_tracker.track{1}.method)];
    fprintf('Saving tuning parameters %s\n',out_param_fn);
    ct_save(out_param_fn,'param');
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
