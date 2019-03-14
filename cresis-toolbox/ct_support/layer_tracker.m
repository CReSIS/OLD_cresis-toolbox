function Surface = layer_tracker(param,param_override)
% Surface = layer_tracker(param,param_override)
%
% Tracks a layer using one of the layer tracking methods.
%
% track: Parameters used to control the tracker. Only general parameters
%   used within update_surface_with_tracker are included here. For more
%   details about tracker specific parameters, look at the help for:
%   tracker_snake_simple, tracker_threshold, tracker_max
%
% track structure fields:
% -------------------------------------------------------------------------
% .en: must be true for update_surface_with_tracker to run on a segment
% .medfilt: scalar which specifies the length of a median filter to apply
%   to the tracked data. Similar to medfilt1, but it handles edge
%   conditions and allows for a threshold parameter to be set.
% .medfilt_threshold: scalar used with medfilt. Median filter will only
%   update in the difference of the median filter is more than this
%   threshold. Default is 0 which means every pixel is updated.
% .feedthru: all values in the radar image which exceed the power mask set
%   by the time and power_dB vectors will be set to zero power. To get
%   these values a typical procedure is:
%     plot(mdata.Time, lp(mean(mdata.Data,2)))
%     [track.feedthru.time,track.feedthru.power_dB] = ginput; % Press enter to end the input
%   .time: N length vector of two way travel times
%   .power_dB: N length vector of power in dB
% .method: string containing the method to use for tracking. Options:
%    'threshold': runs tracker_threshold (generally the best for surface altimetry)
%    'snake': runs tracker_snake_simple
%    'max': tracker_max
% .max_diff: used by tracker routines (optional, default is inf)
% .min_bin: used by tracker routines
% .max_bin: used by tracker routines (optional, default is not used)
% .init.method: used by tracker routines (optional, default is not used)
% .init.dem_layer: dem layer struct to be used with opsLoadLayers
%   (usually lidar surface)
%
% Surface: If there is an output argument, then param.cmd.frms must be set
% to a single frame and the surface for this frame will be returned and no
% outputs will be saved.
%
% For more details about tracker specific parameters:
% tracker_snake_simple, tracker_threshold, tracker_max
%
% Example:
%   See run_layer_tracker.m to run in regular mode.
%   See qlook_combine_task.m to run in single frame mode.
%
% Author: John Paden


%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

%% Input Checks
% =====================================================================

physical_constants;

if ~isfield(param,'layer_tracker') || isempty(param.layer_tracker)
  param.layer_tracker = [];
end

if ~isfield(param.layer_tracker,'debug_plots') || isempty(param.layer_tracker.debug_plots)
  param.layer_tracker.debug_plots = {};
end
if ~isempty(param.layer_tracker.debug_plots)
  h_fig = get_figures(1,true);
end
enable_debug_plot = any(strcmp('debug',param.layer_tracker.debug_plots));

% echogram_img: To choose an image besides the base (0) image
if ~isfield(param.layer_tracker,'echogram_img') || isempty(param.layer_tracker.echogram_img)
  param.layer_tracker.echogram_img = 0;
end
echogram_img = param.layer_tracker.echogram_img;

% echogram_source: location of echogram data used for tracking
if ~isfield(param.layer_tracker,'echogram_source') || isempty(param.layer_tracker.echogram_source)
  error('An echogram_source must be specified.');
end
echogram_source = param.layer_tracker.echogram_source;

if ~isfield(param.layer_tracker,'layer_params') || isempty(param.layer_tracker.layer_params)
  param.layer_tracker.layer_params = [];
end
layer_params = param.layer_tracker.layer_params;

%% General Setup Continued
% =====================================================================

if isstruct(echogram_source)
  fprintf('  %s: %s (%s)\n', mfilename, 'param.layer_tracker.echogram_source', datestr(now));
else
  fprintf('=====================================================================\n');
  fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
  fprintf('=====================================================================\n');
end

%% Input Checks: track field
% =====================================================================

if ~isfield(param.layer_tracker,'track') || isempty(param.layer_tracker.track)
  param.layer_tracker.track = [];
end
track = merge_structs(param.qlook.surf,param.layer_tracker.track);

% profile: default is no profile, otherwise loads preset configuration
if ~isfield(track,'profile') || isempty(track.profile)
  track.profile = '';
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
debug_time_guard = track.debug_time_guard;

if ~isfield(track,'detrend') || (isempty(track.detrend) && ~ischar(track.detrend))
  track.detrend = 0;
end
if ischar(track.detrend)
  % Load detrend file generated by run_get_echogram_stats
  % e.g. ct_tmp/echogram_stats/snow/2011_Greenland_P3/stats_20110329_02.mat
  if isempty(track.detrend)
    track.detrend = [ct_filename_ct_tmp(param,'','echogram_stats','stats') '.mat'];
  end
  detrend = load(track.detrend,'dt','bins','min_means');
  detrend.time = detrend.dt*detrend.bins;
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

%% Load in ocean mask, land DEM, and sea surface DEM
if isfield(track,'init') && strcmpi(track.init.method,'dem')
  global gdem;
  if isempty(gdem) || ~isa(gdem,'dem_class') || ~isvalid(gdem)
    gdem = dem_class(param,500);
  end
  gdem.set_res(500);
  gdem.ocean_mask_mode = 'l';

  gdem_str = sprintf('%s:%s:%s',param.radar_name,param.season_name,param.day_seg);
  if ~strcmpi(gdem_str,gdem.name)
    % Load records file
    records_fn = ct_filename_support(param,'','records');
    records = load(records_fn);
    gdem.set_vector(records.lat,records.lon,gdem_str);
  end
end

%% Determine valid frames to process
if nargout ~= 1
  frames_fn = ct_filename_support(param,'','frames');
  if ~exist(frames_fn,'file')
    warning('Cannot verify param.cmd.frms because frames file does not exist: %s.', frames_fn);
  else
    load(frames_fn);
    if isempty(param.cmd.frms)
      param.cmd.frms = 1:length(frames.frame_idxs);
    end
    % Remove frames that do not exist from param.cmd.frms list
    [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
    if length(valid_frms) ~= length(param.cmd.frms)
      bad_mask = ones(size(param.cmd.frms));
      bad_mask(keep_idxs) = 0;
      warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
        param.cmd.frms(find(bad_mask,1)));
      param.cmd.frms = valid_frms;
    end
  end
end

%% Load reference surface
if isfield(track,'init') && isfield(track.init,'dem_layer') ...
    && ~isempty(track.init.dem_layer)
  layers = opsLoadLayers(param,track.init.dem_layer);
  
  % Ensure that layer gps times are monotonically increasing
  lay_idx = 1;
  layers_fieldnames = fieldnames(layers(lay_idx));
  [~,unique_idxs] = unique(layers(lay_idx).gps_time);
  for field_idx = 1:length(layers_fieldnames)-1
    if ~isempty(layers(lay_idx).(layers_fieldnames{field_idx}))
      layers(lay_idx).(layers_fieldnames{field_idx}) = layers(lay_idx).(layers_fieldnames{field_idx})(unique_idxs);
    end
  end
end

%% Track
orig_track = track;
for frm_idx = 1:length(param.cmd.frms)
  %% Track: Load echogram data
  frm = param.cmd.frms(frm_idx);
  track = orig_track;
  if isstruct(echogram_source)
    % echogram_source is the data structure
    mdata = echogram_source;
    data_fn_name = 'param.layer_tracker.echogram_source';
    
  else
    % echogram_source is a file path indicator
    if echogram_img == 0
      data_fn = fullfile(ct_filename_out(param,echogram_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    else
      data_fn = fullfile(ct_filename_out(param,echogram_source,''), ...
        sprintf('Data_img_%02d_%s_%03d.mat', echogram_img, param.day_seg, frm));
    end
    [~,data_fn_name] = fileparts(data_fn);
    fprintf('%d of %d %s (%s)\n', frm_idx, length(param.cmd.frms), data_fn, datestr(now));
    
    if ~exist(data_fn,'file')
      warning('  Missing file\n');
      continue;
    end
    
    if strcmpi(track.method,'') && ~enable_debug_plot
      mdata = load_L1B(data_fn,'GPS_time','Latitude','Longitude','Elevation','Time');
    else
      mdata = load_L1B(data_fn);
    end
  end
  
  data = lp(mdata.Data);
  Nx = size(mdata.Data,2);
  if size(mdata.Data,1) < 2
    Surface = nan(1,Nx);
    new_quality = ones(1,Nx);
    
  else
    
    %% Track: Interpolate GIMP and Geoid
    if strcmpi(track.init.method,'dem')
      gdem.set_vector(mdata.Latitude,mdata.Longitude);
      [land_dem,msl,ocean_mask] = gdem.get_vector_dem();
      
      % Merge land surface and sea surface DEMs
      track.dem = land_dem;
      track.dem(ocean_mask) = msl(ocean_mask);
      track.dem = (mdata.Elevation - track.dem) / (c/2);
      track.dem = interp1(mdata.Time,1:length(mdata.Time),track.dem + track.init.dem_offset,'linear','extrap');
      track.dem = interp_finite(track.dem,1);
    end
    
    %% Track: Prefilter trim
    % Also set leading/following zeros or ~isfinite to NaN
    for rline = 1:Nx
      start_bin = find(isfinite(data(:,rline)),1);
      if ~isempty(start_bin)
        stop_bin = min(size(data,1), start_bin+track.prefilter_trim(1)-1);
        data(1:stop_bin,rline) = NaN;
      end
      stop_bin = find(isfinite(data(:,rline)),1,'last');
      if ~isempty(stop_bin)
        start_bin = max(1, stop_bin-track.prefilter_trim(2)+1);
        data(start_bin:end,rline) = NaN;
      end
    end
    
    if track.data_noise_en
      data_noise = data;
    end
    
    %% Track: Feed through removal
    if ~isempty(track.feedthru)
      % Interpolate feed through power levels on to data time axis
      feedthru_threshold = interp1(track.feedthru.time,track.feedthru.power_dB,mdata.Time);
      feedthru_threshold = interp_finite(feedthru_threshold,-inf,[],@(x) ~isnan(x));
      
      % Set all data to zero that does not exceed the feed through
      % threshold power
      for rline=1:Nx
        data(data(:,rline) <= feedthru_threshold,rline) = NaN;
      end
    end
    
    %% Track: Detrend
    if ischar(track.detrend)
      % A detrend file was passed in
      detrend_curve = interp_finite(interp1(detrend.time,interp_finite(detrend.min_means),mdata.Time),NaN);
      if all(isnan(detrend_curve))
        error('Detrend curve is all NaN.');
      end
      if 0
        % Debug
        rline = 200;
        figure(1); clf;
        plot(data(:,rline))
        hold on
        mean_power = nanmean(data,2);
        plot(mean_power)
        plot(detrend_curve);
        keyboard
      end
      data = bsxfun(@minus,data,detrend_curve);
      if track.data_noise_en
        data_noise = bsxfun(@minus,data_noise,detrend_curve);
      end
      
    elseif track.detrend > 0
      poly_x = (-size(data,1)/2+(1:size(data,1))).';
      mean_power = nanmean(data,2);
      good_mask = isfinite(mean_power);
      p = polyfit(poly_x(good_mask),mean_power(good_mask),track.detrend);
      detrend_curve = polyval(p,poly_x);
      detrend_curve(~good_mask) = NaN;
      detrend_curve = interp_finite(detrend_curve,0);
      if 0
        % Debug
        rline = 200;
        figure(1); clf;
        plot(data(:,rline))
        hold on
        plot(mean_power)
        plot(detrend_curve);
        keyboard
      end
      data = bsxfun(@minus,data,detrend_curve);
      if track.data_noise_en
        data_noise = bsxfun(@minus,data_noise,detrend_curve);
      end
    end
    
    %% Track: Sidelobe
    if ~isempty(track.sidelobe_rows)
      mask = sidelobe_mask_mex(single(data),int32(track.sidelobe_rows),single(track.sidelobe_dB));
      data(mask) = NaN;
    end
    
    %% Track: Filter
    if track.filter(1) ~= 1
      % Multilooking in cross-track/fast-time
      data = lp(nan_fir_dec(10.^(data/10).',ones(1,track.filter(1))/track.filter(1),1,[],[],[],[],2.0).');
      if track.data_noise_en
        data_noise = lp(nan_fir_dec(10.^(data_noise/10).',ones(1,track.filter(1))/track.filter(1),1,[],[],[],[],2.0).');
      end
    end
    if track.filter(2) ~= 1
      % Multilooking in along-track
      data = lp(nan_fir_dec(10.^(data/10),ones(1,track.filter(2))/track.filter(2),1,[],[],[],[],2.0));
      if track.data_noise_en
        data_noise = lp(nan_fir_dec(10.^(data_noise/10),ones(1,track.filter(2))/track.filter(2),1,[],[],[],[],2.0));
      end
    end
    
    %% Track: Post-filter trim
    for rline = 1:Nx
      start_bin = find(isfinite(data(:,rline)),1);
      if ~isempty(start_bin)
        stop_bin = min(size(data,1), start_bin+track.filter_trim(1)-1);
        data(start_bin:stop_bin,rline) = NaN;
        if track.data_noise_en
          data_noise(start_bin:stop_bin,rline) = NaN;
        end
      end
      stop_bin = find(isfinite(data(:,rline)),1,'last');
      if ~isempty(stop_bin)
        start_bin = max(1, stop_bin-track.filter_trim(2)+1);
        data(start_bin:stop_bin,rline) = NaN;
        if track.data_noise_en
          data_noise(start_bin:stop_bin,rline) = NaN;
        end
      end
    end
    
    %% Track: min_bin/max_bin + time to bins conversions
    % Convert from two way travel time to bins
    track.min_bin = find(mdata.Time >= orig_track.min_bin, 1);
    track.max_bin = find(mdata.Time <= orig_track.max_bin, 1, 'last');
    dt = mdata.Time(2) - mdata.Time(1);
    track.init.max_diff = orig_track.init.max_diff/dt;
    if strcmpi(track.max_rng_units,'bins')
      track.max_rng = orig_track.max_rng(1) : orig_track.max_rng(end);
    else
      track.max_rng = round(orig_track.max_rng(1)/dt) : round(orig_track.max_rng(end)/dt);
    end
    data = data(track.min_bin:track.max_bin,:);
    if track.data_noise_en
      data_noise = data_noise(track.min_bin:track.max_bin,:);
    end
    
    %% Track: Create Initial Surface
    if strcmpi(track.init.method,'dem')
      % Correct for min_bin removal
      track.dem = track.dem - track.min_bin + 1;
    elseif strcmp(track.init.method,'snake')
      track.init.search_rng = round(orig_track.init.snake_rng(1)/dt) : round(orig_track.init.snake_rng(2)/dt);
      track.dem = tracker_snake_simple(data,track.init);
    elseif strcmp(track.init.method,'nan')
      track.dem = nan(1,Nx);
    else
      % max or medfilt
      [~,track.dem] = max(data,[],1);
      if strcmp(track.init.method,'medfilt')
        track.dem = medfilt1(track.dem,track.init.medfilt);
      end
    end
    
    %% Track: Merge DEM and reference layer
    if ~isempty(track.init.dem_layer)
      % Interpolate reference layer onto radar GPS time
      ref_interp_gaps_dist = [150 75];
      ops_layer = [];
      ops_layer{1}.gps_time = layers(lay_idx).gps_time;
      ops_layer{1}.type = layers(lay_idx).type;
      ops_layer{1}.quality = layers(lay_idx).quality;
      ops_layer{1}.twtt = layers(lay_idx).twtt;
      ops_layer{1}.type(isnan(ops_layer{1}.type)) = 2;
      ops_layer{1}.quality(isnan(ops_layer{1}.quality)) = 1;
      lay = opsInterpLayersToMasterGPSTime(mdata,ops_layer,ref_interp_gaps_dist);
      dem_layer = lay.layerData{1}.value{2}.data;
      dem_layer = interp1(mdata.Time,1:length(mdata.Time),dem_layer);
      track.dem = merge_vectors(dem_layer, track.dem);
    end
    
    %% Track: Tracking
    if strcmpi(track.method,'threshold')
      track.threshold_noise_rng = round(orig_track.threshold_noise_rng/dt);
      if track.data_noise_en
        track.data_noise = data_noise;
      else
        track.data_noise = data;
      end
      [new_layer,new_quality] = tracker_threshold(data,track);
    elseif strcmpi(track.method,'viterbi')
      new_layer = tracker_viterbi(data,track);
      new_quality = ones(1,Nx);
    elseif strcmpi(track.method,'max')
      new_layer = tracker_max(data,track);
      new_quality = ones(1,Nx);
    elseif strcmpi(track.method,'snake')
      track.search_rng = round(orig_track.snake_rng(1)/dt) : round(orig_track.snake_rng(2)/dt);
      new_layer = tracker_snake_simple(data,track);
      new_quality = ones(1,Nx);
    elseif strcmpi(track.method,'fixed')
      new_layer = ones(size(mdata.GPS_time)) * track.fixed_value;
      new_quality = ones(1,Nx);
    elseif isempty(track.method)
      new_layer = track.dem;
      new_quality = ones(1,Nx);
    else
      error('Not a supported layer tracking method.');
    end
    
    %% Track: max_diff
    new_layer(abs(new_layer  - track.dem) > track.init.max_diff) = NaN;
    switch (track.init.max_diff_method)
      case 'merge_vectors'
        new_layer = merge_vectors(new_layer, track.dem);
      case 'interp_finite'
        new_layer = interp_finite(new_layer,NaN);
    end
    
    %% Track: Median filtering
    if isfield(track,'medfilt') && ~isempty(track.medfilt)
      % This median filter is designed to only operate when the point in
      % question is an outlier exceeding track.medfilt_threshold
      for rline=1:Nx
        rlines = rline + (-track.medfilt:track.medfilt);
        rlines = rlines(rlines>=1 & rlines<=length(new_layer));
        if abs(new_layer(rline) - nanmedian(new_layer(rlines))) > track.medfilt_threshold
          new_layer(rline) = nanmedian(new_layer(rlines));
        end
      end
    end
    
    %% Track: Max search
    if (length(track.max_rng) > 1 || track.max_rng ~= 0)
      % Find the next peak after the threshold
      for rline = 1:Nx
        search_bins = round(new_layer(rline)) + track.max_rng;
        search_bins = search_bins(find(search_bins >= 1 & search_bins <= size(data,1)));
        [~,offset] = max(data(search_bins,rline));
        if ~isempty(offset)
          new_layer(rline) = search_bins(offset);
        end
      end
    end
    
    %% Track: Convert bins to twtt
    Surface = interp1(1:length(mdata.Time), mdata.Time, new_layer + track.min_bin - 1,'linear','extrap');
    
  end
  
  % Some layer sources may not be "double", but we require that Surface be
  % double:
  Surface = double(Surface);
  
  %% Track: Debug plot
  if enable_debug_plot
    clf(h_fig(1));
    set(h_fig(1),'name',sprintf('layer_tracker %s',data_fn_name));
    h_axes(1) = axes('parent',h_fig(1));
    imagesc([],mdata.Time,lp(mdata.Data), 'parent', h_axes(1));
    colormap(h_axes(1), 1-gray(256));
    hold(h_axes(1),'on');
    plot(h_axes(1),find(new_quality==1),Surface(new_quality==1),'g.');
    plot(h_axes(1),find(new_quality==3),Surface(new_quality==3),'r.');
    if strcmpi(track.init.method,{'dem'}) || ~isempty(track.init.dem_layer)
      plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time,track.dem+track.min_bin-1),'m--')
      plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time, ...
        track.dem+track.min_bin-1-track.init.max_diff),'r--')
      plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time, ...
        track.dem+track.min_bin-1+track.init.max_diff),'b--')
    end
    hold(h_axes(1),'off');
    if ~isempty(mdata.Time)
      ylims = [max(mdata.Time(1),min(Surface)-debug_time_guard) min(mdata.Time(end),max(Surface)+debug_time_guard)];
      if ylims(end)>ylims(1)
        ylim(h_axes(1),ylims);
      end
    end
    title(h_axes(1),sprintf('%s',regexprep(data_fn_name,'_','\\_')));
    keyboard
  end
  
  %% Track: Save
  if nargout == 1
    return;
  end
  for layer_idx = 1:length(layer_params)
    layer_param = layer_params(layer_idx);
    
    if strcmpi(layer_param.source,'echogram')
      if isempty(layer_param.echogram_source)
        layer_param.echogram_source = echogram_source;
      end
      
      data_fn = fullfile(ct_filename_out(param,layer_param.echogram_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      fprintf('  Saving %s (%s)\n', data_fn, datestr(now));
      save(data_fn,'-append','Surface');
    end
    
    if strcmpi(layer_param.source,'layerdata')
      layer_fn = fullfile(ct_filename_out(param,layer_param.layerdata_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      if ~exist(layer_fn,'file')
        fprintf('  Create  %s (%s)\n', layer_fn, datestr(now));
        
        lay.GPS_time = mdata.GPS_time;
        lay.Latitude = mdata.Latitude;
        lay.Longitude = mdata.Longitude;
        lay.Elevation = mdata.Elevation;
        
        lay.layerData{1}.quality = interp1(mdata.GPS_time,new_quality,lay.GPS_time,'nearest');
        lay.layerData{1}.value{1}.data = nan(size(lay.GPS_time));
        lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
        lay.layerData{1}.value{2}.data = interp_finite(lay.layerData{1}.value{2}.data,NaN);
        
        lay.layerData{2}.quality = ones(size(lay.GPS_time));
        lay.layerData{2}.value{1}.data = nan(size(lay.GPS_time));
        lay.layerData{2}.value{2}.data = nan(size(lay.GPS_time));
        
        layer_fn_dir = fileparts(layer_fn);
        if ~exist(layer_fn_dir,'dir')
          mkdir(layer_fn_dir);
        end
        save(layer_fn,'-struct','lay');
        
      else
        % Load the layerData file
        lay = load(layer_fn);
        % Update the surface auto picks
        lay.layerData{1}.quality = interp1(mdata.GPS_time,new_quality,lay.GPS_time,'nearest');
        lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
        lay.layerData{1}.value{2}.data = interp_finite(lay.layerData{1}.value{2}.data,NaN);
        % Append the new results back to the layerData file
        fprintf('  Saving %s (%s)\n', layer_fn, datestr(now));
        save(layer_fn,'-append','-struct','lay','layerData');
      end
    end
    
    if strcmpi(layer_param.source,'ops')
      % Get all the frames for this segment
      if any(strcmpi({layer_params.source},'ops'))
        opsAuthenticate(param,false);
        sys = ct_output_dir(param.radar_name);
        ops_param = struct('properties',[]);
        ops_param.properties.season = param.season_name;
        ops_param.properties.segment = param.day_seg;
        [status,ops_seg_data] = opsGetSegmentInfo(sys,ops_param);
      end
      
      % OPS query to get the point path ID's
      ops_param = struct('properties',[]);
      ops_param.properties.location = param.post.ops.location;
      ops_param.properties.season = param.season_name;
      ops_param.properties.start_gps_time = ops_seg_data.properties.start_gps_time(frm);
      ops_param.properties.stop_gps_time = ops_seg_data.properties.stop_gps_time(frm);
      
      sys = ct_output_dir(param.radar_name);
      [status,data] = opsGetPath(sys,ops_param);
      
      % Write the new layer information to these point path ID's
      ops_param = struct('properties',[]);
      ops_param.properties.point_path_id = data.properties.id;
      ops_param.properties.twtt = interp_finite(interp1(mdata.GPS_time,Surface,data.properties.gps_time));
      ops_param.properties.type = 2*ones(size(ops_param.properties.twtt));
      ops_param.properties.quality = interp1(mdata.GPS_time,new_quality,data.properties.gps_time,'nearest');
      ops_param.properties.lyr_name = layer_param.name;
      
      opsCreateLayerPoints(sys,ops_param);
    end
    
  end
  
end
