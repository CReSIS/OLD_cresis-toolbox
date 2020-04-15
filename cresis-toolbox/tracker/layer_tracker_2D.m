function ctrl_chain = layer_tracker_2D(param,param_override)
% ctrl_chain = layer_tracker_2D(param,param_override)
%
% Check input parameters and create tasks for layer_tracker.
% layer_tracker_task does the actual tracking.

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks: cmd
% =====================================================================

% Remove frames that do not exist from param.cmd.frms list
frames = frames_load(param);
param.cmd.frms = frames_param_cmd_frms(param,frames);

%% Input Checks: layer_tracker
% =====================================================================

physical_constants;

%  .layer_tracker: parameter structure controlling layer tracker
if ~isfield(param,'layer_tracker') || isempty(param.layer_tracker)
  param.layer_tracker = [];
end

%  .copy_param: Final output opsCopyLayer parameter struct
if ~isfield(param.layer_tracker,'copy_param') || isempty(param.layer_tracker.copy_param)
  param.layer_tracker.copy_param = [];
end
%  .copy_param.copy_method
if ~isfield(param.layer_tracker.copy_param,'copy_method') || isempty(param.layer_tracker.copy_param.copy_method)
  param.layer_tracker.copy_param.copy_method = 'overwrite';
end
%  .copy_param.gaps_fill
if ~isfield(param.layer_tracker.copy_param,'gaps_fill') || isempty(param.layer_tracker.copy_param.gaps_fill)
  param.layer_tracker.copy_param.gaps_fill = [];
end
%  .copy_param.gaps_fill.method
if ~isfield(param.layer_tracker.copy_param.gaps_fill,'method') || isempty(param.layer_tracker.copy_param.gaps_fill.method)
  param.layer_tracker.copy_param.gaps_fill.method = 'preserve_gaps';
end

%  .debug_plots: cell array of strings
if ~isfield(param.layer_tracker,'debug_plots') || isempty(param.layer_tracker.debug_plots)
  param.layer_tracker.debug_plots = {};
end

if ~isfield(param.layer_tracker,'debug_out_dir') || isempty(param.layer_tracker.debug_out_dir)
  param.layer_tracker.debug_out_dir = 'layer_tracker';
end
debug_out_dir = param.layer_tracker.debug_out_dir;

%  .echogram_img: To choose an image besides the base (0) image
if ~isfield(param.layer_tracker,'echogram_img') || isempty(param.layer_tracker.echogram_img)
  param.layer_tracker.echogram_img = 0;
end

%  .echogram_source: string containing location of echogram data used for
%  tracking
if ~isfield(param.layer_tracker,'echogram_source') || isempty(param.layer_tracker.echogram_source)
  param.layer_tracker.echogram_source = 'qlook';
end

%  .layer_params: layerparams structure of where to store the output using
%  opsCopyLayers.m
if ~isfield(param.layer_tracker,'layer_params') || isempty(param.layer_tracker.layer_params)
  param.layer_tracker.layer_params = [];
end
if ~isfield(param.layer_tracker.layer_params,'source') || isempty(param.layer_tracker.layer_params.source)
  param.layer_tracker.layer_params.source = 'layerdata';
end
if ~isfield(param.layer_tracker.layer_params,'layerdata_source') || isempty(param.layer_tracker.layer_params.layerdata_source)
  param.layer_tracker.layer_params.layerdata_source = 'layer';
end

%  .surf_layer: layer parameter structure for loading a layer with
%  opsLoadLayers.m
if ~isfield(param.layer_tracker,'surf_layer') || isempty(param.layer_tracker.surf_layer)
  param.layer_tracker.surf_layer = [];
end

%  .crossover_layer: layer parameter structure for loading crossovers with
%  opsGetCrossovers.m
if ~isfield(param.layer_tracker,'crossover_layer') || isempty(param.layer_tracker.crossover_layer)
  param.layer_tracker.crossover_layer = [];
end
if ~isempty(param.layer_tracker.crossover_layer)
  if ~isfield(param.layer_tracker.crossover_layer,'name') || isempty(param.layer_tracker.crossover_layer.name)
    param.layer_tracker.crossover_layer.name = 'bottom';
  end
  if ~isfield(param.layer_tracker.crossover_layer,'source') || isempty(param.layer_tracker.crossover_layer.source)
    param.layer_tracker.crossover_layer.source = 'ops';
  end
end

%% Input Checks: layer_tracker.track field
% ======================================================================

if ~isfield(param.layer_tracker,'track') || isempty(param.layer_tracker.track)
  param.layer_tracker.track = [];
end

method_idx = 1;

track = merge_structs(param.qlook.surf,param.layer_tracker.track);

% profile: default is no profile, otherwise loads preset configuration
if ~isfield(param.layer_tracker.track,'profile') || isempty(param.layer_tracker.track.profile)
  param.layer_tracker.track.profile = '';
end

track = merge_structs(layer_tracker_profile(param,track.profile), track);

if ~isfield(track,'en') || isempty(track.en)
  % If true, tracking will be done on this segment. If false, then no
  % tracking is done. Default is true.
  track.en = true;
end
if ~track.en
  return;
end

if ~isfield(track,'data_noise_en') || isempty(track.data_noise_en)
  % If true, then a matrix data_noise will be created which will go through
  % all the same operations as data except no "sidelobe" and no "feedthru"
  % masking will be done. This is sometimes necessary to enable when
  % sidelobe or feedthru remove too much data so that the noise estimatation
  % process in the tracker does not function properly. Currently only
  % tracker_threshold makes use of data_noise.
  track.data_noise_en = false;
end

% debug_time_guard: Vertical band around layer in debug plot
if ~isfield(track,'debug_time_guard') || isempty(track.debug_time_guard)
  track.debug_time_guard = 50e-9;
end
param.layer_tracker.debug_time_guard = track.debug_time_guard;

if ~isfield(track,'detrend') || isempty(track.detrend)
  track.detrend = [];
end

if ~isfield(track,'feedthru') || isempty(track.feedthru)
  track.feedthru = [];
end

if ~isfield(track,'filter') || isempty(track.filter)
  track.filter = [1 1];
end
if length(track.filter) == 1
  warning('Deprecated surf.filter format. Should specify 2 element vector that specifies the multilooks in [cross-track along-track].');
  track.filter = [1 track.filter(1)];
end
if any(mod(track.filter,2) == 0)
  error('Surface filter lengths must be odd. layer_tracker.track.filter = [%d %d].', layer_tracker.track.filter);
end

if ~isfield(track,'filter_trim') || isempty(track.filter_trim)
  track.filter_trim = [0 0];
end

if ~isfield(track,'fixed_value') || isempty(track.fixed_value)
  track.fixed_value = 0;
end

if ~isfield(track,'init') || isempty(track.init)
  track.init = [];
end
if ~isfield(track.init,'method') || isempty(track.init.method)
  track.init.method = 'max';
end
if ~isfield(track.init,'snake_rng') || isempty(track.init.snake_rng)
  track.init.snake_rng = [-2e-7 2e-7];
end
if ~isfield(track.init,'dem_layer') || isempty(track.init.dem_layer)
  track.init.dem_layer = '';
end
if ~any(strcmpi(track.init.method,{'max','nan','snake','medfilt','dem'}))
  error('Unsupported surface init method %s. Options are max, nan, snake, medfilt, or dem. max is default.', track.init.method);
end
if ~isfield(track.init,'max_diff') || isempty(track.init.max_diff)
  track.init.max_diff = inf;
end
if ~isfield(track.init,'max_diff_method') || isempty(track.init.max_diff_method)
  if strcmpi(track.init.method,{'dem'}) || ~isempty(track.init.dem_layer)
    % If the initial surface is from a dem or reference layer, the default
    % method for outliers is to use merge_vectors to fill the outliers in.
    track.init.max_diff_method = 'merge_vectors';
  else
    % Otherwise the default method is just to interpolate between the good
    % points that exist to fill the outliers in.
    track.init.max_diff_method = 'interp_finite';
  end
end
if ~any(strcmpi(track.init.max_diff_method,{'merge_vectors','interp_finite'}))
  error('Unsupported max diff method %s. Options are merge_vectors, interp_finite. The default is interp_finite unless a dem or reference layer is provided.', track.init.max_diff_method);
end

if ~isfield(track,'max_bin') || isempty(track.max_bin)
  track.max_bin = inf;
end

if ~isfield(track,'max_rng') || isempty(track.max_rng)
  track.max_rng = [0 0];
end

if ~isfield(track,'max_rng_units') || isempty(track.max_rng_units)
  track.max_rng_units = 'time';
end

if ~isfield(track,'medfilt') || isempty(track.medfilt)
  % Like medfilt1 except it handles the edges of the frame better
  track.medfilt = 0;
end
if ~isfield(track,'medfilt_threshold') || isempty(track.medfilt_threshold)
  % medfilt_threshold: the point is compared to the medfilt1 result and if
  % the point is > medfilt_threshold from medfilt1, then the point is
  % replaced with the medfilt1 result. The two extremes are:
  %   0 makes medfilt act like medfilt1
  %   inf effectively disables the medfilt operation
  track.medfilt_threshold = 0;
end

if ~isfield(track,'method') || isempty(track.method)
  track.method = '';
end

if ~isfield(track,'min_bin') || isempty(track.min_bin)
  track.min_bin = 0;
end

if ~isfield(track,'name') || isempty(track.name)
  track.name = sprintf('m%03d', method_idx);
end

if ~isfield(track,'prefilter_trim') || isempty(track.prefilter_trim)
  track.prefilter_trim = [0 0];
end

if ~isfield(track,'sidelobe_rows') || isempty(track.sidelobe_rows) || ~isfield(track,'sidelobe_dB') || isempty(track.sidelobe_dB)
  track.sidelobe_dB = [];
  track.sidelobe_rows = [];
end

if ~isfield(track,'snake_rng') || isempty(track.snake_rng)
  track.snake_rng = [-2e-7 2e-7];
end

if ~isfield(track,'threshold') || isempty(track.threshold)
  track.threshold = 15;
end

if ~isfield(track,'threshold_noise_rng') || isempty(track.threshold_noise_rng)
  track.threshold_noise_rng = [0 -inf -1];
end

param.layer_tracker.track = track;

%% Set up Cluster
% ===================================================================

ctrl = cluster_new_batch(param);
%cluster_compile({'layer_tracker_task'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);
cluster_compile({'layer_tracker_task','layer_tracker_combine_task'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);
ctrl_chain = {};

%% layer_tracker
% =========================================================================
% =========================================================================

% Cluster setup
% -------------------------------------------------------------------------

sparam.argsin{1} = param;
sparam.task_function = 'layer_tracker_task';
sparam.num_args_out = 1;
sparam.argsin{1}.load.echogram_img = param.layer_tracker.echogram_img;
sparam.cpu_time = 60;
sparam.mem = 500e6;
sparam.notes = '';

cpu_time_mult = 0;
mem_mult = 0;
for method_idx = 1:length(param.layer_tracker.track)
  switch param.layer_tracker.track.method(method_idx)
    case 'viterbi'
      cpu_time_mult = cpu_time_mult + 11e-6;
      mem_mult = max(mem_mult,64);
      
    case 'lsm'
      cpu_time_mult = cpu_time_mult + 5.5e-5;
      mem_mult = max(mem_mult,64);
      
    otherwise
      cpu_time_mult = cpu_time_mult + 11e-6;
      mem_mult = max(mem_mult,64);
  end
end

%% layer_tracker: Loop to create tasks
% -------------------------------------------------------------------------
in_fn_dir = ct_filename_out(param,param.layer_tracker.echogram_source,'');
if strcmp(param.layer_tracker.layer_params.source,'ops')
  tmp_out_fn_dir_dir = ct_filename_out(param,'ops','layer_tracker_tmp');
  param.layer_tracker.layer_params.layerdata_source = 'ops'; % Only used for stdout
else
  tmp_out_fn_dir_dir = ct_filename_out(param,param.layer_tracker.layer_params.layerdata_source,'layer_tracker_tmp');
end
mem_combine = 0;
cputime_combine = 0;
frm_idx = 1;
while frm_idx <= length(param.cmd.frms)
  dparam = [];
  dparam.file_success = {};
  dparam.argsin{1}.layer_tracker.frms = [];
  Nx = 0;
  Nt = 0;
  
  for subblock_idx = 1:param.layer_tracker.block_size_frms
    frm = param.cmd.frms(frm_idx);
    if ~any(frm == param.cmd.frms)
      break;
    end
    % Add frame to this block
    frm_idx = frm_idx + 1;
    dparam.argsin{1}.layer_tracker.frms(end+1) = frm;
    
    % File Success
    % ---------------------------------------------------------------------
    % Outputs stored in:
    % /cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_layer_tracker_tmp/CSARP_layer_test/20140313_08/
    %
    % Comparing four different methods:
    %   layer_tracker_001/m001_lsm.mat, ..., layer_tracker_00N/m001_lsm.mat
    %   layer_tracker_001/m002_mcmc.mat, ..., layer_tracker_00N/m002_mcmc.mat
    %   layer_tracker_001/m003_stereo.mat, ..., layer_tracker_00N/m003_stereo.mat
    %   layer_tracker_001/m004_viterbi.mat, ..., layer_tracker_00N/m004_viterbi.mat
    % Layers in the files (all combined into one file during combine):
    %   m001_lsm_surface_001, ..., m001_lsm_surface_016, m001_lsm_bottom_001, ..., m001_lsm_bottom_016
    %   m002_mcmc_surface, m002_mcmc_bottom
    %   m003_stereo_surface, m003_stereo_bottom
    %   m004_viterbi_bottom
    %
    % Comparing the same method with four different sets of parameters:
    %   layer_tracker_001/m001_lsm.mat, ..., layer_tracker_00N/m001_lsm.mat
    %   layer_tracker_001/m002_lsm.mat, ..., layer_tracker_00N/m001_lsm.mat
    %   layer_tracker_001/m003_lsm.mat, ..., layer_tracker_00N/m001_lsm.mat
    %   layer_tracker_001/m004_lsm.mat, ..., layer_tracker_00N/m001_lsm.mat
    % Layers in the files (all combined into one file during combine):
    %   m001_lsm_surface_001, ..., m001_lsm_surface_016, m001_lsm_bottom_001, ..., m001_lsm_bottom_016
    %   m002_lsm_surface_001, ..., m002_lsm_surface_016, m002_lsm_bottom_001, ..., m002_lsm_bottom_016
    %   m003_lsm_surface_001, ..., m003_lsm_surface_016, m003_lsm_bottom_001, ..., m003_lsm_bottom_016
    %   m004_lsm_surface_001, ..., m004_lsm_surface_016, m004_lsm_bottom_001, ..., m004_lsm_bottom_016
    for method_idx = 1:length(param.layer_tracker.track)
      tmp_out_fn_name = sprintf('%s_%s.mat', param.layer_tracker.track(method_idx).name, param.layer_tracker.track.method);
      tmp_out_fn = fullfile(tmp_out_fn_dir_dir,sprintf('layer_tracker_%03d', frm),tmp_out_fn_name);
      dparam.file_success{end+1} = tmp_out_fn;
      if ~ctrl.cluster.rerun_only && exist(tmp_out_fn,'file')
        delete(tmp_out_fn);
      end
    end
    
    % CPU time and memory
    % ---------------------------------------------------------------------
    data_fn = fullfile(in_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    mdata = load(data_fn, 'GPS_time','Time');
    if (subblock_idx==1)
      max_time = mdata.Time(end);
      min_time = mdata.Time(1);
    else
      if(max_time <= mdata.Time(end))
        max_time = mdata.Time(end);
      end
      if(min_time >= mdata.Time(1))
        min_time = mdata.Time(1);
      end
    end
    Nx = Nx + length(mdata.GPS_time);
  end
  dt = mdata.Time(2) - mdata.Time(1);
  Nt = 1 + (max_time-min_time)/dt;
    
  % Rerun only check
  % ---------------------------------------------------------------------
  if ~ctrl.cluster.rerun_only
    if ~cluster_file_success(dparam.file_success)
      fprintf('  Already exists [rerun_only skipping]: %s (%s)\n', ...
        dparam.notes, datestr(now));
      continue;
    end
  end
  
  % CPU time and memory
  % ---------------------------------------------------------------------
  dparam.cpu_time = cpu_time_mult * Nx * Nt;
  dparam.mem = 800e6 + mem_mult * Nx * Nt;
  mem_combine = mem_combine + 256*Nx;
  cputime_combine = cputime_combine + 1e-1*Nx;
  
  % Notes
  % ---------------------------------------------------------------------
  dparam.notes = sprintf('%s %s:%s:%s %s %s %s %d-%d (%d of %d)', ...
    sparam.task_function, param.radar_name, param.season_name, ...
    param.layer_tracker.echogram_source, param.layer_tracker.layer_params.layerdata_source, ...
    param.layer_tracker.track(1).method, param.day_seg, ...
    dparam.argsin{1}.layer_tracker.frms([1 end]), frm_idx, length(param.cmd.frms));
  
  % Create task
  % ---------------------------------------------------------------------
  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
end

ctrl = cluster_save_dparam(ctrl);
ctrl_chain{end+1} = ctrl;
fprintf('Done %s\n',datestr(now));

%% layer_tracker_combine
% =========================================================================
% =========================================================================
ctrl = cluster_new_batch(param);

sparam = [];
sparam.argsin{1} = param;
sparam.task_function = 'layer_tracker_combine_task';
sparam.num_args_out = 1;

sparam.cpu_time = 30 + cputime_combine;
sparam.mem = 300e6 + mem_combine;
sparam.notes = '';

if strcmp(param.layer_tracker.layer_params.source,'ops')
  sparam.file_success = {};
else
  sparam.file_success = {};
  out_fn_dir = ct_filename_out(param,'',param.layer_tracker.layer_params.layerdata_source);
  for frm = param.cmd.frms
    out_fn = fullfile(out_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    sparam.file_success{end+1} = out_fn;
  end
end

ctrl = cluster_new_task(ctrl,sparam,[]);
ctrl_chain{end+1} = ctrl;
fprintf('Done %s\n',datestr(now));

