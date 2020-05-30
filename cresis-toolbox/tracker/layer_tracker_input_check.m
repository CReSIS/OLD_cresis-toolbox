% script layer_tracker_input_check
%
% Support function for layer_tracker.m

%% Input Checks: layer_tracker
% =====================================================================

physical_constants;

%  .layer_tracker: parameter structure controlling layer tracker
if ~isfield(param,'layer_tracker') || isempty(param.layer_tracker)
  param.layer_tracker = [];
end

%  .block_size_frms: Number of data frames to load and process at a time.
%  Default is 1 frame at a time.
if ~isfield(param.layer_tracker,'block_size_frms') || isempty(param.layer_tracker.block_size_frms)
  param.layer_tracker.block_size_frms = 1;
end
param.layer_tracker.block_size_frms = min(length(param.cmd.frms), param.layer_tracker.block_size_frms);

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

%% Input Checks: layer_tracker.track field
% ======================================================================

if ~isfield(param.layer_tracker,'track') || isempty(param.layer_tracker.track)
  param.layer_tracker.track = {};
end

for track_idx = 1:length(param.layer_tracker.track)
  
  track = merge_structs(param.qlook.surf,param.layer_tracker.track{track_idx});
  
  % profile: default is the profile named after the system (e.g. accum,
  % kaband, kuband, rds, or snow), otherwise loads preset configuration
  if ~isfield(track,'profile') || isempty(track.profile)
    track.profile = ct_output_dir(param.radar_name);
  end
  
  track = merge_structs(layer_tracker_profile(param,track.profile), track);
  
  if ~isfield(track,'en') || isempty(track.en)
    % If true, tracking will be done on this segment. If false, then no
    % tracking is done. Default is true.
    track.en = true;
  end
  if ~track.en
    continue;
  end
  
  %  .crossover: struct controlling crossover loading
  if ~isfield(track,'crossover') || isempty(track.crossover)
    track.crossover = [];
  end
  %  .crossover.en: enable loading of crossovers
  if ~isfield(track.crossover,'en') || isempty(track.crossover.en)
    track.crossover.en = false;
  end
  %  .crossover.name: layer name to load crossovers for
  if ~isfield(track.crossover,'name') || isempty(track.crossover.name)
    track.crossover.name = 'bottom';
  end
  %  .crossover.season_names_bad: cell array of strings with bad seasons to
  %  not include in crossovers
  if ~isfield(track.crossover,'season_names_bad') || isempty(track.crossover.season_names_bad)
    track.crossover.season_names_bad = {};
  end
  %  .crossover.gps_time_good_eval: function which returns good/bad based
  %  on crossover gps time.
  if ~isfield(track.crossover,'gps_time_good_eval') || isempty(track.crossover.gps_time_good_eval)
    track.crossover.gps_time_good_eval = @(x) true;
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
  
  %  .ice_mask: struct controlling ice mask loading
  if ~isfield(track,'ice_mask') || isempty(track.ice_mask)
    track.ice_mask = [];
  end
  %  .ice_mask.en: enable loading of ice mask
  if ~isfield(track.ice_mask,'en') || isempty(track.ice_mask.en)
    track.ice_mask.en = false;
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
  elseif isstruct(track.max_bin)
    if ~isfield(track.max_bin,'existence_check')
      track.max_bin.existence_check = false;
    end
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
  elseif isstruct(track.min_bin)
    if ~isfield(track.min_bin,'existence_check')
      track.min_bin.existence_check = false;
    end
  end
  
  %  .mult_suppress: struct controlling surface multiple suppression
  if ~isfield(track,'mult_suppress') || isempty(track.mult_suppress)
    track.mult_suppress = [];
  end
  %  .mult_suppress.en: enable surface multiple suppression
  if ~isfield(track.mult_suppress,'en') || isempty(track.mult_suppress.en)
    track.mult_suppress.en = false;
  end
  %  .mult_suppress.param: param field to pass to mult_suppress.m
  if ~isfield(track.mult_suppress,'param') || isempty(track.mult_suppress.param)
    track.mult_suppress.param = [];
  end
  
  if ~isfield(track,'name') || isempty(track.name)
    track.name = sprintf('t%03d', track_idx);
  end
  
  if ~isfield(track,'track_per_task') || isempty(track.track_per_task)
    track.track_per_task = inf;
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
  
  param.layer_tracker.track{track_idx} = track;
end
