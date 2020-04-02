function success = layer_temp_task(param)
% layer_tracker_2D_task is used to set the parameters (options) for the
% different tracking algorithms. Look into calling the algorithms and
% running it on the cluster
success = {};
%% Load reference surface
if any(strcmp(param.layer_tracker.track.method,'snake')) || any(strcmp(param.layer_tracker.track.method,'threshold'))|| any(strcmp(param.layer_tracker.track.method,'max')) || any(strcmp(param.layer_tracker.track.method,'fixed')) || any(strcmp(param.layer_tracker.track.method,'viterbi')) || any(strcmp(param.layer_tracker.track.method,'lsm')) || any(strcmp(param.layer_tracker.track.method,'stereo')) || any(strcmp(param.layer_tracker.track.method,'mcmc'))
  
  %intialising twtt, gps_time to save layer information
  layer_params = param.layer_tracker.cmds.layer_params;
  param.layer_tracker.fname = [];
  if param.layer_tracker.track.track_data
    for iter = 1:16
      for layers_idx = 1:length(layer_params)
        twtt{iter,layers_idx} = [];
      end
    end
  else
    for layer_idx = 1:length(layer_params)
      twtt{layer_idx} = [];
    end
  end
  %
  
  
  %
  gps_time = [];

  
  if isfield(param.layer_tracker.tracker,'init') && isfield(param.layer_tracker.tracker.init,'dem_layer') ...
      && ~isempty(param.layer_tracker.tracker.init.dem_layer)
    layers = opsLoadLayers(param,param.layer_tracker.tracker.init.dem_layer);
    
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
  
  orig_track = param.layer_tracker.tracker;
  %%
  
  save_name = '/cresis/snfs1/scratch/anjali/CSARP_layerData_Anjali';
  dir_name1 = 'CSARP_layer_tracker_tmp';
  dir_name2 = 'CSARP_layer_tracker';
  frm_dir = sprintf('layer_tracker_%03d_%03d', param.frm_nums(1),param.frm_nums(end));
  if ~exist(fullfile(save_name,param.season_name,dir_name1,dir_name2,param.day_seg,frm_dir))
    mkdir(fullfile(save_name,param.season_name,dir_name1,dir_name2,param.day_seg,frm_dir));
  end
  tmp_out_fn_dir = fullfile(save_name,param.season_name,dir_name1,dir_name2,param.day_seg,frm_dir);
  
  
  mdata = [];
  mdata.GPS_time = [];
  mdata.Data = [];
  mdata.Time = [];
  mdata.Elevation=[];
  mdata.Latitude = [];
  mdata.Longitude = [];
  mdata.data_length = [];
  mdata.Bottom = [];
  
  for frm_idx = 1:length(param.frm_nums)
    frm = param.frm_nums(frm_idx);
    data_fn = fullfile(ct_filename_out(param,param.layer_tracker.echogram_source,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    lay=load(data_fn);
    param.layer_tracker.fname{end+1} = data_fn;
    mdata.GPS_time = cat(2,mdata.GPS_time,lay.GPS_time);
    mdata.Data = cat(2,mdata.Data,lay.Data);
    mdata.Elevation = cat(2,mdata.Elevation,lay.Elevation);
    mdata.Latitude = cat(2,mdata.Latitude,lay.Latitude);
    mdata.Longitude = cat(2, mdata.Longitude,lay.Longitude);
    param.layer_tracker.tracker.frm = frm;
    mdata.data_length = cat(2,mdata.data_length,length(lay.Bottom));
    mdata.Bottom = cat(2,mdata.Bottom,lay.Bottom);
  end
  mdata.Time = lay.Time;
  orig_track = param.layer_tracker.tracker;
  track = orig_track;
  Nx = size(mdata.Data,2);
  data = lp(mdata.Data);
  
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
  
  %% Track: Prefilter trim%
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
  if ischar(track.detrend)
    % A detrend file was passed in
    track.detrend_struct = detrend;
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
      elseif strcmpi(track.method,'max')
        new_layer = tracker_max(data,track);
        new_quality = ones(1,Nx);
      elseif strcmpi(track.method,'viterbi')
        new_layer = tracker_viterbi(mdata,param);
        new_quality = ones(1,Nx);
      elseif strcmpi(track.method,'lsm')
        if param.layer_tracker.track.track_data
          %new_layer = tracker_lsm_track(mdata,param);
          %new_quality = ones(1,Nx);
          new_layer = load('/cresis/snfs1/projects/LSM_anjali/LSM_result_part9.mat');
          new_quality = ones(1,Nx);
        else
          new_layer = tracker_lsm(mdata,param);
          new_quality = ones(1,Nx);
        end
      elseif strcmpi(track.method,'stereo')
        [new_layer,big_matrix] = tracker_stereo(mdata,param);
        new_quality = ones(1,Nx);
      elseif strcmpi(track.method,'mcmc')
        new_layer = tracker_mcmc(mdata,param);
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
      
  new_temp= new_layer;
  
  for layer_idx = 1:length(layer_params)
    
    layer_param = layer_params(layer_idx);
    
    new_layer = new_temp(layer_idx,:);
    track.init.max_diff = inf;
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
    %Bottom = interp1(1:length(mdata.Time), mdata.Time, new_layer + track.min_bin - 1,'linear','extrap');
    %end x check if needed or not
    % Some layer sources may not be "double", but we require that Surface be double:
    Surface = double(Surface);
    %Bottom = double(Bottom);
    
    %% Track: Debug plot
    if ~isempty(param.layer_tracker.debug_plots)
      h_fig = get_figures(1,true);
    end
    if param.layer_tracker.enable_debug_plot
      clf(h_fig(1));
      %set(h_fig(1),'name',sprintf('layer_tracker %s',data_fn_name));
      h_axes(1) = axes('parent',h_fig(1));
      imagesc([],mdata.Time,lp(mdata.Data), 'parent', h_axes(1));
      colormap(h_axes(1), 1-gray(256));
      hold(h_axes(1),'on');
      plot(h_axes(1),find(new_quality==1),Surface(new_quality==1),'g');
      plot(h_axes(1),find(new_quality==3),Surface(new_quality==3),'r');
      if strcmpi(track.init.method,{'dem'}) || ~isempty(track.init.dem_layer)
        plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time,track.dem+track.min_bin-1),'m--')
        plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time, ...
          track.dem+track.min_bin-1-track.init.max_diff),'r--')
        plot(h_axes(1),interp1(1:length(mdata.Time),mdata.Time, ...
          track.dem+track.min_bin-1+track.init.max_diff),'b--')
      end
      hold(h_axes(1),'off');
      if ~isempty(mdata.Time)
        ylims = [max(mdata.Time(1),min(Surface)-track.debug_time_guard) min(mdata.Time(end),max(Surface)+track.debug_time_guard)];
        %         if ylims(end)>ylims(1)
        %           ylim(h_axes(1),ylims);
        %         end
      end
      %title(h_axes(1),sprintf('%s',regexprep(data_fn_name,'_','\\_')));
      keyboard
    end
    
    %% Track: Save
    if param.layer_tracker.track.save_layerData
      %     if nargout == 1
      %       return;
      %     end
      %     layer_params = param.layer_tracker.layer_params;
      %     for layer_idx = 1:length(layer_params)
      %       layer_param = layer_params(layer_idx);
      try
        lay_param.source = layer_param{1}.source;
        lay_param.echogram_source = layer_param{1}.echogram_source;
        lay_param.layerdata_source = layer_param{1}.layerdata_source;
      catch ME
        lay_param.source = layer_param.source;
        lay_param.echogram_source = layer_param.echogram_source;
        lay_param.layerdata_source = layer_param.layerdata_source;
      end
      
      if strcmpi(lay_param.source,'echogram')
        if isempty(lay_param.echogram_source)
          lay_param.echogram_source = echogram_source;
        end
        
        data_fn = fullfile(ct_filename_out(param,lay_param.echogram_source,''), ...
          sprintf('Data_%s_%03d.mat', param.day_seg, frm));
        fprintf('  Saving %s (%s)\n', data_fn, datestr(now));
        %save(data_fn,'-append','Surface');
      end
      
      if strcmpi(lay_param.source,'layerdata')
        layer_fn = fullfile(ct_filename_out(param,lay_param.layerdata_source,''), ...
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
          %save(layer_fn,'-struct','lay');
          
        else
          % Load the layerData file
          %           lay = load(layer_fn);
          %           % Update the surface auto picks
          %           lay.layerData{layer_idx}.quality = interp1(mdata.GPS_time,new_quality,lay.GPS_time,'nearest');
          %           lay.layerData{layer_idx}.value{2}.data = interp1(mdata.GPS_time,Surface,lay.GPS_time);
          %           %lay.layerData{1}.value{2}.data = interp1(mdata.GPS_time,new_temp.(sprintf('layer_%s_%03d',param.day_seg,frm)).top,lay.GPS_time);
          %           lay.layerData{layer_idx}.value{2}.data = interp_finite(lay.layerData{layer_idx}.value{2}.data,NaN);
          % Append the new results back to the layerData file
          fprintf('  Saving %s (%s)\n', layer_fn, datestr(now));
          %filename = fullfile(save_name,param.season_name,dir_name,param.day_seg,frm_dir, sprintf('layers_%s_%03d.mat',param.layer_tracker.track.method,layer_idx));
          filename = fullfile(save_name,param.season_name,dir_name1,dir_name2,param.day_seg,frm_dir,sprintf('layer_%s_%s',param.layer_tracker.track.method,param.layer_tracker.name));
          tmp_name = sprintf('%s_%s_%03d',param.layer_tracker.track.method,param.layer_tracker.name,layer_idx);
          %quality{layer_idx} = interp1(mdata.GPS_time,new_quality,lay.GPS_time,'nearest');
          
          %           tmp.layerData{1} = cat(2,tmp.layerData{1},lay.layerData);
          %           tmp.GPS_time{1} = cat(2,tmp.GPS_time{1},lay.GPS_time);
          %           tmp.param_layer_tracker{1} = cat(2,tmp.tracker{1},param);
          twtt{layer_idx} = cat(2,twtt{layer_idx},Surface);
          
          %           tmp.layerData = lay.layerData;
          %           tmp.GPS_time  = cat(2,tmp.GPS_time,lay.GPS_time);
          %           tmp.param_layer_tracker = param;
          %           if param.ct_file_lock
          %             tmp.file_version = '1L';
          %           else
          %             tmp.file_version = '1';
          %           end
          %           tracker_data{layer_idx} = struct(tmp_name,tmp);
          %save(filename, '-struct','lay','layerData');
          %save(layer_fn,'-append','-struct','lay','layerData');
        end
      end
      
      if strcmpi(lay_param.source,'ops')
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
      %end
    end
  end
  
  gps_time = cat(2,gps_time,mdata.GPS_time);
  param_layer_tracker = param;
  if param.ct_file_lock
    file_version = '1L';
  else
    file_version = '1';
  end
    save(filename,'twtt','gps_time','param_layer_tracker','file_version');%,'indexes','twtt_test');
  success = true;
end