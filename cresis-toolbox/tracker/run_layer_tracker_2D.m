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
param_override.layer_tracker.echogram_source = 'CSARP_post/mvdr';

% layer_params: layerparams structure of where to store the output using
% opsCopyLayers.m
param_override.layer_tracker.layer_params = [];
% Uncomment to enable layerdata storage
param_override.layer_tracker.layer_params.layerdata_source = 'layer_test';
% Uncomment to enable OPS storage
% param_override.layer_tracker.layer_params.source = 'ops';

% block_size_frms: Number of frames to be loaded at a time
param_override.layer_tracker.block_size_frms = 2;

%% param.layer_tracker.track options
track_override = [];

% if 1 % If using GeoTIFF file for ice mask
%   track_override.binary_icemask = false;
%   track_override.icemask_fn = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
%   track_override.icemask_fn = ct_filename_gis([], track_override.icemask_fn);
%                        
%   % Useful for Antarctica seasons:
%   % track_override.icemask_fn  = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_icemask_grounded_and_shelves.tif';
%   % track_override.icemask2_fn = 'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_rockmask.tif';
%   %   if isfield(track_override, 'icemask2_fn') && ~isempty(track_override.icemask2_fn)
%   %     track_override.icemask2_fn = ct_filename_gis([],track_override.icemask2_fn);
%   %   end
% else
%   track_override.binary_icemask = true;
%   track_override.icemask_fn = '/cresis/snfs1/dataproducts/GIS_data/canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.bin';
%   [track_override.ice_mask_fn_dir,track_override.ice_mask_fn_name] = fileparts(track_override.icemask_fn);
%   track_override.ice_mask_mat_fn = fullfile(track_override.ice_mask_fn_dir,[track_override.ice_mask_fn_name '.mat']);
% end

if 1 % If using GeoTIFF file for ice mask
  if strcmpi(params(1).post.ops.location,'arctic')
    if 1
      % Greenland
      track_override.binary_icemask = false;
      track_override.icemask_fn = 'greenland/IceMask/GimpIceMask_90m_v1.1.tif';
      track_override.icemask_fn = ct_filename_gis([], track_override.icemask_fn);
    else
      % Canada
      track_override.binary_icemask = true;
      track_override.icemask_fn = '/cresis/snfs1/dataproducts/GIS_data/canada/ice_mask/03_rgi50_ArcticCanadaNorth/03_rgi50_ArcticCanadaNorth.bin';
      [track_override.ice_mask_fn_dir,track_override.ice_mask_fn_name] = fileparts(track_override.icemask_fn);
      track_override.ice_mask_mat_fn = fullfile(track_override.ice_mask_fn_dir,[track_override.ice_mask_fn_name '.mat']);
    end
  else
    % Useful for Antarctica seasons:
    track_override.binary_icemask = false;
    track_override.icemask_fn = ct_filename_gis([], 'greenland/IceMask/GimpIceMask_90m_v1.1.tif');
  end
end

%% Enable one set of parameters
track_override.en = true;
switch ct_output_dir(params(1).radar_name)
  case 'rds'
    % RDS
    track_override.profile = 'rds_OIB';

    % Override default filter settings
    if 0
      track_override.filter	= [3 3];
      track_override.filter_trim	= [3 3];
      track_override.threshold = 10;
      track_override.max_rng	= [0 2];
    end
    
    % Use sidelobe rejection
    if 0
      % run_get_echogram_stats output
      sidelobe = load('/N/dcwan/projects/cresis/output/ct_tmp/echogram_stats/rds/2018_Greenland_P3/stats_20180421_01.mat','sidelobe_rows','sidelobe_dB','sidelobe_vals');
      track_override.sidelobe_rows = [sidelobe.sidelobe_rows(75:98)];
      track_override.sidelobe_dB = -(sidelobe.sidelobe_dB(75:98,1)-max(sidelobe.sidelobe_dB(:,1))+21);
      track_override.sidelobe_dB(track_override.sidelobe_dB<9) = 9;
      track_override.threshold_rel_max = -max(track_override.sidelobe_dB);
      track_override.data_noise_en = true;
    end
    
    % Use feedthrough rejection
    if 0
      % run_get_echogram_stats output
      feedthru = load('/N/dcwan/projects/cresis/output/ct_tmp/echogram_stats/rds/2018_Greenland_P3/stats_20180421_01.mat');
      track_override.feedthru.time = feedthru.dt*feedthru.bins;
      track_override.feedthru.power_dB = feedthru.min_means+20;
      bin_mask = track_override.feedthru.time<2e-6;
      track_override.feedthru.time = track_override.feedthru.time(bin_mask);
      track_override.feedthru.power_dB = track_override.feedthru.power_dB(bin_mask);
      track_override.feedthru.power_dB(end) = -inf;
      track_override.min_bin = 0.5e-6;
      track_override.data_noise_en = true;
    end
    
    % Override default init method
    if 0
      track_override.init.method	= 'dem';
      track_override.init.dem_offset = 0;
      track_override.init.dem_layer.name = 'surface';
      track_override.init.dem_layer.source = 'lidar';
      track_override.init.dem_layer.lidar_source = 'atm';
      track_override.init.max_diff = 1e-6;
    elseif 0
      track_override.init.method	= 'snake';
      track_override.init.snake_rng	= [-0.5e-6 0.5e-6];
      track_override.init.max_diff	= 0.5e-6;
    end
    
    % Override default method
    if 0
      track_override.method = 'snake';
      track_override.snake_rng	= [-0.15e-6 0.15e-6];
    end
    
    %% Viterbi
    if 1
      %% Viterbi User Settings
      track_override.method                 = 'viterbi';
      track_override.viterbi.crossoverload  = true;
      track_override.viterbi.layername      = 'viterbi_bot'; %surface or bottom
      track_override.viterbi.detrending     = true;
      track_override.viterbi.top_sup        = true;
      track_override.viterbi.mult_sup       = true;
      track_override.viterbi.use_surf_for_slope = true;
      track_override.viterbi.custom_combine = false;
      track_override.viterbi.DIM_matrix     = fullfile('+tomo', 'Layer_tracking_2D_parameters_Matrix.mat');

      track_override.viterbi.surf_weight    = -1;
      track_override.viterbi.mult_weight    = -1;
      track_override.viterbi.mult_weight_decay       = -1;
      track_override.viterbi.mult_weight_local_decay = -1;
      track_override.viterbi.manual_slope   = 0;
      track_override.viterbi.max_slope      = -1;
      track_override.viterbi.transition_weight = 1;
      track_override.viterbi.image_mag_weight = 1;
      track_override.viterbi.gt_weight = 1;
      track_override.init.max_diff    = inf;
      track_override.viterbi.gt_cutoff = -1;
    end
    
    %% MCMC
    if 0
      %% MCMC User Settings
      track_override.method      = 'mcmc';
      track_override.mcmc.lyrtop = 'mcmc_top'; %layername, layer_dest.name
      track_override.mcmc.lyrbot = 'mcmc_bot';
      track_override.mcmc.alg    = 'MCMC';
      track_override.init.max_diff    = inf;
    end
    
    %% LSM
    if 0
      %% LSM User Settings
      track_override.method           = 'lsm';
      track_override.lsm.lyrtop       = 'lsm_top'; %layername, layer_dest.name
      track_override.lsm.lyrbot       = 'lsm_bot';
      track_override.lsm.y            = 220; % = '' for y = mean(SURF)
      track_override.lsm.dy           = 10;
      track_override.lsm.storeIter    = [200 400];
      track_override.init.max_diff    = inf;
      track_override.detrend          = [];
      track_override.norm.scale       = [-40 90];

    end
    
    %% Stereo
    if 0
      %% Stereo User Settings
      track_override.method               = 'stereo';
      track_override.stereo.lyrtop        = 'stereo_top'; %layername, layer_dest.name
      track_override.stereo.lyrbot        = 'stereo_bot';
      track_override.stereo.surfaceload   = true;
      track_override.stereo.crossoverload = true;
      track_override.stereo.top_smooth    = 1000;
      track_override.stereo.bottom_smooth = 1000;
      track_override.stereo.top_peak      = 0.5;
      track_override.stereo.bottom_peak   = 0.5;
      track_override.stereo.repulsion     = 10;
      track_override.stereo.alg           = 'HMM';
      track_override.init.max_diff    = inf;
    end
    
    %% Threshold
    if 0
      track_override.method = 'threshold';
    end
    
    %% Fixed
    if 0
      track_override.method = 'fixed';
    end
  case 'accum'
    % ACCUM
    track_override.profile = 'ACCUM';
    
    % Override default init method
    if 0
      track_override.init.method	= 'dem';
      track_override.init.dem_offset = 0;
      track_override.init.dem_layer.name = 'surface';
      track_override.init.dem_layer.source = 'lidar';
      track_override.init.dem_layer.lidar_source = 'atm';
      track_override.init.max_diff = 0.3e-6;
    end
    
    % Override default method
    if 0
      track_override.method = 'snake';
      track_override.snake_rng	= [-0.15e-6 0.15e-6];
    end
    
  case {'snow','kuband','kaband'}
    % FMCW
    track_override.profile = 'snow_AWI';
    
    % Use sidelobe rejection
    if 0
      % run_get_echogram_stats output
      sidelobe = load('/cresis/snfs1/dataproducts/ct_data/ct_tmp/echogram_stats/snow/2011_Greenland_P3/stats_20110329_01.mat','sidelobe_rows','sidelobe_dB','sidelobe_vals');
      track_override.sidelobe_rows = [sidelobe.sidelobe_rows(1:194)];
      track_override.sidelobe_dB = [-sidelobe.sidelobe_dB(1:194,end)]-sidelobe.sidelobe_vals(end)+4.5;
      track_override.threshold_rel_max = -max(track_override.sidelobe_dB);
    end
    
    % Override default init method
    if 0
      track_override.init.method  = 'dem';
      track_override.init.dem_offset = 0;
      track_override.init.dem_layer.name = 'surface';
      track_override.init.dem_layer.source = 'lidar';
      track_override.init.dem_layer.lidar_source = 'atm';
      track_override.init.max_diff = 0.3e-6;
    elseif 0
      track_override.init.method  = 'snake';
      track_override.init.snake_rng = [-15e-9 15e-9];
      track_override.init.max_diff  = 0.3e-6;
    end
    
end

param_override.layer_tracker.track = track_override;

% param_override.layer_tracker.surf_layer = struct('name','surface','source','layerdata','layerdata_source','layer');

% param_override.layer_tracker.crossover_layer = struct('name','bottom','source','ops');


% dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
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
    
    ctrl_chain{end+1} = layer_tracker_2D(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);
