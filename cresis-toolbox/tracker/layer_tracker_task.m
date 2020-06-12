function success = layer_tracker_task(param)
% success = layer_tracker_task(param)
%
% layer_tracker_task is used to load and prepare the data and metadata and
% run the tracker. See layer_tracker.m.
%
% param: parameter structure from parameter spreadsheet
%   param.layer_tracker.echogram_source: The normal cluster mode is for
%   this to be a ct_filename_out compatible string. The normal
%   qlook_combine_task mode is for this to be a structure and this function
%   returns the tracked surface and does not save the result.
%
% success:
%   In cluster mode: logical scalar, true when task completes successfully
%   In qlook_combine_task mode: surface twtt corresponding to the
%   param.layer_tracker.echogram_source.GPS_time.

%% Input Checks: track field
% =====================================================================

physical_constants('c');

if isstruct(param.layer_tracker.echogram_source)
  % echogram_source is the data structure
  param.layer_tracker.track = {param.qlook.surf};
  layer_tracker_input_check;
end

tracked_images_en = any(strcmp('tracked_images',param.layer_tracker.debug_plots));
visible_en = any(strcmp('visible',param.layer_tracker.debug_plots));
if tracked_images_en
  h_fig = get_figures(1,visible_en);
end

%% Create output directory paths
if strcmp(param.layer_tracker.layer_params.source,'ops')
  tmp_out_fn_dir_dir = ct_filename_out(param,'ops','layer_tracker_tmp');
  param.layer_tracker.layer_params.layerdata_source = 'ops'; % Only used for stdout
else
  tmp_out_fn_dir_dir = ct_filename_out(param,param.layer_tracker.layer_params.layerdata_source,'layer_tracker_tmp');
end

%% Load echogram data
if isstruct(param.layer_tracker.echogram_source)
  % echogram_source is the data structure
  mdata = param.layer_tracker.echogram_source;
  frm_str = {sprintf('%s_%03d',param.day_seg,param.cmd.frms)};
  frm_start = 1;
  frm_stop = length(mdata.GPS_time);
  param.layer_tracker.tracks_in_task = 1;
  
else
  % Create input directory paths
  in_fn_dir = ct_filename_out(param,param.layer_tracker.echogram_source,'');
  
  mdata = [];
  mdata.GPS_time = [];
  mdata.Data = [];
  mdata.Time = [];
  mdata.Elevation=[];
  mdata.Latitude = [];
  mdata.Longitude = [];
  mdata.Roll = [];
  % Keep track of start/stop index for each frame
  frm_start = zeros(size(param.layer_tracker.frms));
  frm_stop = zeros(size(param.layer_tracker.frms));
  frm_str = cell(size(param.layer_tracker.frms));
  for frm_idx = 1:length(param.layer_tracker.frms)
    frm = param.layer_tracker.frms(frm_idx);
    
    % Load the current frame
    frm_str{frm_idx} = sprintf('%s_%03d',param.day_seg,frm);
    data_fn = fullfile(in_fn_dir, sprintf('Data_%s.mat',frm_str{frm_idx}));
    if frm_idx == 1
      mdata = load_L1B(data_fn);
      frm_start(frm_idx) = 1;
      frm_stop(frm_idx) = length(mdata.GPS_time);
      
    else
      tmp_data = load_L1B(data_fn);
      
      Nx = length(tmp_data.GPS_time);
      frm_start(frm_idx) = length(mdata.GPS_time)+1;
      frm_stop(frm_idx) = length(mdata.GPS_time)+Nx;
      
      mdata.GPS_time(1,end+(1:Nx)) = tmp_data.GPS_time;
      mdata.Elevation(1,end+(1:Nx)) = tmp_data.Elevation;
      mdata.Latitude(1,end+(1:Nx)) = tmp_data.Latitude;
      mdata.Longitude(1,end+(1:Nx)) = tmp_data.Longitude;
      mdata.Roll(1,end+(1:Nx)) = tmp_data.Roll;
      mdata.Surface(1,end+(1:Nx)) = tmp_data.Surface;
      
      % Handle time axis that changes from one frame to the next
      mdata.Time = tmp_data.Time;
      mdata.Data(:,end+(1:Nx)) = tmp_data.Data;
    end
  end
end
Nx = size(mdata.Data,2);
data = lp(mdata.Data);

%% Track
% =========================================================================
for track_idx = param.layer_tracker.tracks_in_task
  track = param.layer_tracker.track{track_idx};
  orig_track = track;
  
  %% Load reference surface
  if isfield(track,'init') && isfield(track.init,'dem_layer') ...
      && ~isempty(track.init.dem_layer)
    layers = opsLoadLayers(param,track.init.dem_layer);
    
    % Ensure that layer gps times are monotonically increasing
    lay_idx = 1;
    layers_fieldnames = {'gps_time','twtt','elev','lat','lon','type','quality','frm'};
    [~,unique_idxs] = unique(layers(lay_idx).gps_time);
    for field_idx = 1:length(layers_fieldnames)-1
      if ~isempty(layers(lay_idx).(layers_fieldnames{field_idx}))
        layers(lay_idx).(layers_fieldnames{field_idx}) = layers(lay_idx).(layers_fieldnames{field_idx})(unique_idxs);
      end
    end
    mdata.Surface = interp_finite(interp1(layers(lay_idx).gps_time,layers(lay_idx).twtt,mdata.GPS_time));
  end
  
  %% Load in ocean mask, land DEM, and sea surface DEM
  if isfield(track,'init') && strcmpi(track.init.method,'dem')
    global gdem;
    if isempty(gdem) || ~isa(gdem,'dem_class') || ~isvalid(gdem)
      gdem = dem_class(param,500);
    end
    gdem.set_res(500);
    gdem.ocean_mask_mode = 'l';
    
    gdem_str = sprintf('%s:%s:%s_%03d_%03d',param.radar_name,param.season_name,param.day_seg,param.layer_tracker.frms([1 end]));
    if ~strcmpi(gdem_str,gdem.name)
      gdem.set_vector(mdata.Latitude,mdata.Longitude,gdem_str);
    end
  end
  
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
  
  %% Ice mask calculation
  if track.ice_mask.en
    if strcmp(track.ice_mask.type,'bin')
      mask = load(track.ice_mask.mat_fn,'R','X','Y','proj');
      [fid,msg] = fopen(track.ice_mask.fn,'r');
      if fid < 1
        fprintf('Could not open file %s\n', track.ice_mask.fn);
        error(msg);
      end
      mask.maskmask = logical(fread(fid,[length(mask.Y),length(mask.X)],'uint8'));
      fclose(fid);
    else
      [mask.maskmask,mask.R,~] = geotiffread(track.ice_mask.fn);
      mask.proj = geotiffinfo(track.icemask_fn);
    end
    [mask.x, mask.y] = projfwd(mask.proj, mdata.Lat, mdata.Lon);
    mask.X = mask.R(3,1) + mask.R(2,1)*(1:size(mask.maskmask,2));
    mask.Y = mask.R(3,2) + mask.R(1,2)*(1:size(mask.maskmask,1));
    [mask.X,mask.Y] = meshgrid(mask.X,mask.Y);
    ice_mask.mask = round(interp2(mask.X, mask.Y, double(mask.maskmask), mask.x, mask.y));
    ice_mask.mask(isnan(ice_mask.mask)) = 1;
    track.ice_mask = ice_mask.mask;
  else
    track.ice_mask = ones(1,Nx);
  end

  %% Crossover loading
  if track.crossover.en
    sys = ct_output_dir(param.radar_name);
    ops_param = [];
    ops_param.properties.search_str = frm_str{frm_idx};
    ops_param.properties.location = param.post.ops.location;
    ops_param.properties.season = param.season_name;
    [status,ops_data_segment] = opsGetFrameSearch(sys,ops_param);
    if status == 1
      ops_param = [];
      ops_param.properties.location = param.post.ops.location;
      ops_param.properties.lyr_name = track.crossover.name;
      ops_param.properties.frame = frm_str;
      ops_param.properties.segment_id = ops_data_segment.properties.segment_id;
      [status,ops_data] = opsGetCrossovers(sys,ops_param);
      if status == 1
        ops_data.properties.twtt = ops_data.properties.twtt(:).';
        ops_data.properties.source_point_path_id = double(ops_data.properties.source_point_path_id(:).');
        ops_data.properties.cross_point_path_id = double(ops_data.properties.cross_point_path_id(:).');
        ops_data.properties.source_elev = ops_data.properties.source_elev.';
        ops_data.properties.cross_elev = ops_data.properties.cross_elev.';
        
        ops_param = [];
        ops_param.properties.location = param.post.ops.location;
        ops_param.properties.season = param.season_name;
        ops_param.properties.segment_id = ops_data_segment.properties.segment_id;
        ops_param.properties.point_path_id = ops_data.properties.source_point_path_id;
        ops_param.properties.return_geom = 'geog';
        ops_param.properties.lyr_name = track.crossover.name;
        [status,ops_data_gps_time] = opsGetLayerPoints(sys,ops_param);
        
        ops_param = [];
        ops_param.properties.location = param.post.ops.location;
        ops_param.properties.season = param.season_name; % This field is actually ignored
        ops_param.properties.segment_id = ops_data_segment.properties.segment_id;
        ops_param.properties.point_path_id = ops_data.properties.cross_point_path_id;
        ops_param.properties.return_geom = 'geog';
        ops_param.properties.lyr_name = track.crossover.name;
        [status,ops_data_cross_gps_time] = opsGetLayerPoints(sys,ops_param);
        
        [~,crossover_idx,gps_time_idx] = intersect(ops_data.properties.source_point_path_id.',ops_data_gps_time.properties.point_path_id);
        
        source_gps_time = zeros(size(ops_data.properties.source_point_path_id));
        crossover_gps_time = zeros(size(ops_data.properties.source_point_path_id));
        good_mask = true(size(ops_data.properties.source_point_path_id));
        for idx = 1:length(ops_data.properties.source_point_path_id)
          match_idx = find(ops_data.properties.source_point_path_id(idx) == ops_data_gps_time.properties.point_path_id);
          source_gps_time(idx) = ops_data_gps_time.properties.gps_time(match_idx);
          match_idx = find(ops_data.properties.cross_point_path_id(idx) == ops_data_cross_gps_time.properties.point_path_id);
          crossover_gps_time(idx) = ops_data_cross_gps_time.properties.gps_time(match_idx);
          if any(strcmp(ops_data.properties.season_name{idx}, track.crossover.season_names_bad))
            good_mask(idx) = false;
          end
        end
        % Ensure 
        good_mask = good_mask & track.crossover.gps_time_good_eval(crossover_gps_time);
        
        % Convert crossover twtt to source twtt and then convert from twtt
        % to rows/bins
        cols = round(interp1(mdata.GPS_time,1:length(mdata.GPS_time), source_gps_time(good_mask)));
        rows = round(interp1(mdata.Time,1:length(mdata.Time), ops_data.properties.twtt(good_mask) ...
          + (ops_data.properties.source_elev(good_mask) - ops_data.properties.cross_elev(good_mask))*2/c ));
        track.crossovers = [cols; rows];
      end
    end
  else
    track.crossovers = zeros(2,0);
  end
  
  %% Track: Multiple Suppression
  if track.mult_suppress.en
    data = 10*log10(echo_mult_suppress(mdata));
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
    if ~isempty(stop_bin)%
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
  if isfield(track,'detrend') && ~isempty(track.detrend)
    data = echo_detrend(data,track.detrend);
  end
  
  %% Track: Fasttime Cross Correlation (xcorr)
  if isfield(track,'xcorr') && ~isempty(track.xcorr)
    data = echo_xcorr(data,track.xcorr);
  end
  
  %% Track: Sidelobe
  if ~isempty(track.sidelobe_rows)
    mask = sidelobe_mask_mex(single(data),int32(track.sidelobe_rows),single(track.sidelobe_dB));
    data(mask) = NaN;
  end
  
  %% Track: Max Range Filter
  if ~isequal(track.max_rng_filter,track.filter)
    % Multilooking in cross-track/fast-time
    data_max_rng = lp(nan_fir_dec(10.^(data/10).',ones(1,track.max_rng_filter(1))/track.max_rng_filter(1),1,[],[],[],[],2.0).');
    % Multilooking in along-track
    data_max_rng = lp(nan_fir_dec(10.^(data_max_rng/10),ones(1,track.max_rng_filter(2))/track.max_rng_filter(2),1,[],[],[],[],2.0));
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
      if ~isequal(track.max_rng_filter,track.filter)
        data_max_rng(start_bin:stop_bin,rline) = NaN;
      end
    end
    stop_bin = find(isfinite(data(:,rline)),1,'last');
    if ~isempty(stop_bin)
      start_bin = max(1, stop_bin-track.filter_trim(2)+1);
      data(start_bin:stop_bin,rline) = NaN;
      if track.data_noise_en
        data_noise(start_bin:stop_bin,rline) = NaN;
      end
      if ~isequal(track.max_rng_filter,track.filter)
        data_max_rng(start_bin:stop_bin,rline) = NaN;
      end
    end
  end
  
  %% Track: min_bin/max_bin + time to bins conversions
  if isnumeric(track.min_bin)
    % Convert from two way travel time to bins
    track.min_bin = ones(1,Nx) * find(mdata.Time >= orig_track.min_bin, 1);
  elseif isstruct(track.min_bin)
    % Load layer
    track.min_bin = opsLoadLayers(param, orig_track.min_bin);
    % Interpolate min_bin_layer onto echogram GPS times
    track.min_bin = layerdata.interp(mdata,track.min_bin);
    track.min_bin = interp1(mdata.Time,1:length(mdata.Time),track.min_bin.twtt_ref);
  end
  
  if isnumeric(track.max_bin)
    % Convert from two way travel time to bins
    track.max_bin = find(mdata.Time >= orig_track.max_bin, 1);
    if isempty(track.max_bin)
      track.max_bin = ones(1,Nx) * length(mdata.Time);
    else
      track.max_bin = ones(1,Nx) * track.max_bin;
    end
    
  elseif isstruct(track.max_bin)
    % Load layer
    track.max_bin = opsLoadLayers(param, orig_track.max_bin);
    % Interpolate max_bin_layer onto echogram GPS times
    track.max_bin = layerdata.interp(mdata,track.max_bin);
    track.max_bin = interp1(mdata.Time,1:length(mdata.Time),track.max_bin.twtt_ref);
  end
  
  % Fix invalid min_bin and max_bin by setting to the midpoint
  invalid_min_max = track.max_bin < track.min_bin;
  track.min_bin(invalid_min_max) = round((track.max_bin(invalid_min_max) + track.min_bin(invalid_min_max))/2);
  track.max_bin(invalid_min_max) = track.min_bin(invalid_min_max);
  
  min_min_bin = round(max(1,min(track.min_bin)));
  max_max_bin = round(min(size(mdata.Data,1),max(track.max_bin)));
  
  dt = mdata.Time(2) - mdata.Time(1);
  track.init.max_diff = orig_track.init.max_diff/dt;
  if strcmpi(track.max_rng_units,'bins')
    track.max_rng = orig_track.max_rng(1) : orig_track.max_rng(end);
  else
    track.max_rng = round(orig_track.max_rng(1)/dt) : round(orig_track.max_rng(end)/dt);
  end
  data = data(min_min_bin:max_max_bin,:);
  if track.data_noise_en
    data_noise = data_noise(min_min_bin:max_max_bin,:);
  end
  if ~isequal(track.max_rng_filter,track.filter)
    data_max_rng = data_max_rng(min_min_bin:max_max_bin,:);
  end
  track.zero_bin = floor(-mdata.Time(min_min_bin)/dt + 1);
  
  %% Track: Create Initial Surface
  if strcmpi(track.init.method,'dem')
    % Correct for min_bin removal
    track.dem = track.dem - min_min_bin + 1;
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
  
  %% Track: Normalization
  if isfield(track,'norm') && ~isempty(track.norm)
    data = echo_norm(data,track.norm);
  end
  
  %% Track: Tracking
  if strcmpi(track.method,'threshold')
    track.threshold_noise_rng = round(orig_track.threshold_noise_rng/dt);
    if track.data_noise_en
      track.data_noise = data_noise;
    else
      track.data_noise = data;
    end
    [new_layers,new_quality] = tracker_threshold(data,track);
  elseif strcmpi(track.method,'max')
    new_layers = tracker_max(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'viterbi')
    new_layers = tracker_viterbi(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'lsm')
    new_layers = tracker_lsm(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'stereo')
    [new_layers,big_matrix] = tracker_stereo(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'mcmc')
    new_layers = tracker_mcmc(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'snake')
    track.search_rng = round(orig_track.snake_rng(1)/dt) : round(orig_track.snake_rng(2)/dt);
    new_layers = tracker_snake_simple(data,track);
    new_quality = ones(1,Nx);
  elseif strcmpi(track.method,'fixed')
    new_layers = ones(size(mdata.GPS_time)) * track.fixed_value;
    new_quality = ones(1,Nx);
  elseif isempty(track.method)
    new_layers = track.dem;
    new_quality = ones(1,Nx);
  else
    error('Not a supported layer tracking method.');
  end
  
  % Trackers may output multiple layers (different rows in new_layers). For
  % loop to post-process each of these layers.
  for layer_idx = 1:size(new_layers,1)
    new_layer = new_layers(layer_idx,:);
    
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
        if ~isequal(track.max_rng_filter,track.filter)
          [~,offset] = max(data_max_rng(search_bins,rline));
        else
          [~,offset] = max(data(search_bins,rline));
        end
        if ~isempty(offset)
          new_layer(rline) = search_bins(offset);
        end
      end
    end
    
    %% Track: Convert bins to twtt
    new_layer = interp1(1:length(mdata.Time), mdata.Time, new_layer + min_min_bin - 1,'linear','extrap');
    % Some layer sources may not be "double", but we require that layers be double type:
    new_layer = double(new_layer);
    
    %% Track: Debug plot
    if tracked_images_en
      clf(h_fig(1));
      figure_name = sprintf('layer_tracker %s %d-%d layer %d',param.day_seg, param.layer_tracker.frms([1 end]), layer_idx);
      set(h_fig(1),'name',figure_name);
      h_axes(1) = axes('parent',h_fig(1));
      imagesc([],mdata.Time,lp(mdata.Data), 'parent', h_axes(1));
      colormap(h_axes(1), 1-gray(256));
      hold(h_axes(1),'on');
      plot(h_axes(1),find(new_quality==1),new_layer(new_quality==1),'g');
      plot(h_axes(1),find(new_quality==3),new_layer(new_quality==3),'r');
      if strcmpi(track.init.method,{'dem'}) || ~isempty(track.init.dem_layer)
        plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time,track.dem+min_min_bin-1),'m--')
        plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time, ...
          track.dem+min_min_bin-1-track.init.max_diff),'r--')
        plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time, ...
          track.dem+min_min_bin-1+track.init.max_diff),'b--')
      end
      hold(h_axes(1),'off');
      if ~isempty(mdata.Time)
        ylims = [max(mdata.Time(1),min(new_layer)-track.debug_time_guard) min(mdata.Time(end),max(new_layer)+track.debug_time_guard)];
        if ylims(end)>ylims(1)
          ylim(h_axes(1),ylims);
        end
      end
      title(h_axes(1),regexprep(figure_name,'_','\\_'));
      if visible_en
        fprintf('Debug plots: review "%s" and then run "dbcont" to continue.\n', figure_name);
        keyboard
      end
      
      fig_fn = fullfile(ct_filename_ct_tmp(param,'',param.layer_tracker.debug_out_dir,''), ...
        param.layer_tracker.layer_params.layerdata_source, ...
        sprintf('%s_%s_l%03d_frm_%03d.jpg',track.name,track.method,layer_idx,param.layer_tracker.frms(1)));
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      ct_saveas(h_fig(1),fig_fn);
    end
    
    new_layers(layer_idx,:) = new_layer;
  end
  
  %% Track: Save
  if isstruct(param.layer_tracker.echogram_source)
    % Qlook mode
    
    % Return twtt in "success" variable
    success = new_layers(:, frm_start(1):frm_stop(1));
    
  else
    % Cluster Mode
    for frm_idx = 1:length(param.layer_tracker.frms)
      frm = param.layer_tracker.frms(frm_idx);
      
      % Get just the current frames layer data
      gps_time = mdata.GPS_time(frm_start(frm_idx):frm_stop(frm_idx));
      twtt = new_layers(:, frm_start(frm_idx):frm_stop(frm_idx));
      param_layer_tracker = param;
      file_version = '1';
      file_type = 'layer_tracker';
      
      tmp_out_fn_name = sprintf('t%03d_%s.mat', track_idx, track.method);
      tmp_out_fn = fullfile(tmp_out_fn_dir_dir,sprintf('layer_tracker_%03d', frm),tmp_out_fn_name);
      fprintf('  Saving %s (%s)\n', tmp_out_fn, datestr(now));
      tmp_out_fn_dir = fileparts(tmp_out_fn);
      if ~exist(tmp_out_fn_dir,'dir')
        mkdir(tmp_out_fn_dir);
      end
      ct_save(tmp_out_fn,'twtt','gps_time','param_layer_tracker','file_type','file_version');
    end
    success = true;
  end
end

fprintf('Done %s\n',datestr(now));
