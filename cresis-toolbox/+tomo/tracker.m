function tracker (params, options, geotiff_fn, geotiff2_fn)
% function tracker (params, options)
%
% Builds a large matrix containing multiple segments of 2D data and applies
%  the Viterbi layer-tracking program on it.
% Loads crossovers from OPS to use as ground truth for the tracking.
% Other functionality includes multiple supression and
%  a simple detrending technique.
%
% See also: run_tracker.m
%
% Authors: Victor Berger, John Paden


%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar;

physical_constants;
clear('param_override');

% Input checking
if ~exist('params','var')
  error('Use run_tracker: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar, param_override);
else
  param_override = gRadar;
end

warning('off','all');

bottom_bin    = -1;
egt_weight    = -1;
mu_size       = 31;
mu            = log10(exp(-(-(mu_size-1)/2 : (mu_size-1)/2).^4/1));
mu(mu<-30)    = -30;
mu            = mu - mean(mu);
sigma         = sum(abs(mu))/10*ones(1,mu_size);
smooth_var    = inf;
repulsion     = 150000;
smooth_weight = 100;
ice_bin_thr   = 10;

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  clear big_matrix Surface;
  
  param = merge_structs(param,param_override);
  
  dbstack_info = dbstack;
  fprintf('=====================================================================\n');
  fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
  fprintf('=====================================================================\n');
  
  big_matrix          = {};
  big_matrix.Data     = [];
  big_matrix.GPS_time = [];
  big_matrix.Lat      = [];
  big_matrix.Lon      = [];

  % Load frames file
  load(ct_filename_support(param,'','frames'));
  
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
  
  data_fn_dir = ct_filename_out(param, options.name, '');
  for frm = param.cmd.frms
    data_fn_name        = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
    data_fn             = fullfile(data_fn_dir, data_fn_name);
    cur_matrix          = load(data_fn);
    big_matrix.GPS_time = horzcat(big_matrix.GPS_time, cur_matrix.GPS_time);
    big_matrix.Lat      = horzcat(big_matrix.Lat, cur_matrix.Latitude);
    big_matrix.Lon      = horzcat(big_matrix.Lon, cur_matrix.Longitude);
    big_matrix.Time     = cur_matrix.Time;
    %big_matrix.Data     = horzcat(big_matrix.Data, cur_matrix.Data);
    
    %% Image combine
    load(data_fn, 'param_combine');
    param.combine.imgs                   = param_combine.combine.imgs;
    combine                              = param.combine;
    combine.img_comb_mult                = 1.15;
    combine.img_comb_bins                = 50;
    update_img_combine_param.mode        = 'combine';
    update_img_combine_param.update_surf = false;
    combine.img_comb_weights             = [];
    try
      layer_params.name   = 'surface';
      layer_params.source = 'ops';
      layers              = opsLoadLayers(param,layer_params);
    catch ME
      warning(ME.getReport);
      continue;
    end
    combine.out_path              = data_fn_dir;
    combine.frm                   = frm;
    combine.imb_comb_surf         = 0;
    combine.img_comb_weights_mode = 'get_heights';  
    [Data, big_matrix.Time]       = img_combine(param, combine, layers);    
    big_matrix.Data               = horzcat(big_matrix.Data, Data);

    big_matrix.Data = horzcat(big_matrix.Data, cur_matrix.Data);
    
  end
  
  %%
  clear Surface;
  layer_params.name   = 'surface';
  layer_params.source = 'ops';
  try
    Surface = opsLoadLayers(param,layer_params);
  catch ME
    warning(ME.getReport);
    continue;
  end
  
  % Catch error with layer not having correctly loaded from OPS
  if(isempty(Surface.gps_time) || any(isnan(Surface.gps_time)) || any(isnan(Surface.twtt)))
    warning('on', 'all');
    warning('PROBLEM PROCESSING FRAME %s', param.day_seg);
    warning('off', 'all');
    continue;
  end
  
  %%
  opsAuthenticate(param,false);
  layer_name                   = 'bottom';
  sys                          = ct_output_dir(param.radar_name);
  ops_param                    = struct('properties',[]);
  ops_param.properties.season  = param.season_name;
  ops_param.properties.segment = param.day_seg;
  [status,ops_frames]          = opsGetSegmentInfo(sys,ops_param);
  
  ops_param = struct('properties',[]);
  ops_param.properties.location = param.post.ops.location;
  ops_param.properties.season = param.season_name;
  ops_param.properties.start_gps_time = ops_frames.properties.start_gps_time(1);
  ops_param.properties.stop_gps_time = ops_frames.properties.stop_gps_time(end);
  ops_param.properties.nativeGeom = true;
  [~,ops_data] = opsGetPath(sys,ops_param);
  
  big_matrix.Data = interp1(big_matrix.GPS_time, big_matrix.Data.', ops_data.properties.gps_time);
  for rbin=1:size(big_matrix.Data,2)
    big_matrix.Data(:,rbin) = interp_finite(big_matrix.Data(:,rbin));
  end
  big_matrix.Data = big_matrix.Data.';
  big_matrix.Lat = interp_finite(interp1(big_matrix.GPS_time, big_matrix.Lat.', ops_data.properties.gps_time));
  big_matrix.Lon = interp_finite(interp1(big_matrix.GPS_time, big_matrix.Lon.', ops_data.properties.gps_time));
  big_matrix.GPS_time = ops_data.properties.gps_time;
  big_matrix.Surface = interp_finite(interp1(Surface.gps_time,Surface.twtt,big_matrix.GPS_time));
  
  %% Multiple suppression
  big_matrix.Data = lp(big_matrix.Data);
  Tmultiple_offset = -200e-9;
  threshold = -60;
  surface_guard = 10;
  surf_bins = round(interp1(big_matrix.Time, 1:length(big_matrix.Time), big_matrix.Surface));
  for rline = 1:size(big_matrix.Data,2)
    % Threshold and shift the large values
    mult = zeros(size(big_matrix.Data,1),1);
    mult(big_matrix.Data(:,rline)>threshold) = big_matrix.Data(big_matrix.Data(:,rline)>threshold,rline) - threshold;
    % Ignore time before the surface bin since no multiples can occur at these
    % times and the feed through does occur and could create problems if we do not suppress it.
    mult(1:surf_bins(rline)-surface_guard) = 0;
    % Interpolate to the surface multiple time
    mult = interp1((big_matrix.Time+Tmultiple_offset)*2,mult,(big_matrix.Time+Tmultiple_offset));
    big_matrix.Data(:,rline) = big_matrix.Data(:,rline) - mult;
  end
  
  %%
  Nx             = size(big_matrix.Data, 2);
  slope          = round(diff(big_matrix.Surface));
  viterbi_weight = ones([1 Nx]);
  mc             = -1 * ones(1, Nx);
  mc_weight      = 0;
  
  if isempty(geotiff_fn)
    mask = inf * ones([1 Nx]);
  else
    [mask,R,~] = geotiffread(geotiff_fn);
    proj = geotiffinfo(geotiff_fn);
    [points.x, points.y] = projfwd(proj, big_matrix.Lat, big_matrix.Lon);
    X = R(3,1) + R(2,1)*(1:size(mask,2));
    Y = R(3,2) + R(1,2)*(1:size(mask,1));
    if 0
      figure;
      imagesc(X,Y,mask);
      hold on;
      plot(points.x,points.y)
    end
    [X,Y] = meshgrid(X,Y);
    fl_mask = round(interp2(X, Y, double(mask), points.x, points.y));
    mask = fl_mask;
    mask(isnan(mask)) = 1;
    
    % TEMPORARILY HARDCODED FOR ANTARTICA:
    if 0
      [mask,R,~] = geotiffread(geotiff_fn);
      [mask2, R2, ~] = geotiffread(geotiff2_fn);
      proj = geotiffinfo(geotiff_fn);
      [points.x, points.y] = projfwd(proj, big_matrix.Lat, big_matrix.Lon);
      X = R(3,1) + R(2,1)*(1:size(mask,2));
      Y = R(3,2) + R(1,2)*(1:size(mask,1));
      if 0
        figure;
        imagesc(X,Y,mask);
        hold on;
        plot(points.x,points.y)
      end
      [X,Y] = meshgrid(X,Y);
      fl_mask = interp2(X, Y, double(mask), points.x, points.y);
      fl_mask = round(interp2(X, Y, double(mask), points.x, points.y));
      mask = fl_mask;
      mask(isnan(mask)) = 1;
      fl_mask2 = interp2(X, Y, double(mask2), points.x, points.y);
      fl_mask2 = round(interp2(X, Y, double(mask2), points.x, points.y));
      mask2 = fl_mask2;
      mask2(isnan(mask2)) = 1;
      f_mask = zeros(size(mask));
      for i = 1:length(f_mask)
        if ((mask(1, i) == 0 || mask(1, i) == 1) && (mask2(1, i) == 127))
          f_mask(1, i) = 1;
        end
      end
      mask = f_mask;
      if max(max(mask)) > 100
        tmp_mask = zeros(size(mask));
        tmp_mask(mask == 0) = 1;
        mask = tmp_mask;
      end
    end
  end
  
  query = sprintf('SELECT rds_segments.id FROM rds_seasons,rds_segments where rds_seasons.name=''%s'' and rds_seasons.id=rds_segments.season_id and rds_segments.name=''%s''',param.season_name,param.day_seg);
  [status,tables] = opsQuery(query);
  segment_id = tables{1};
  
  ops_param                       = struct('properties',[]);
  ops_param.properties.location   = param.post.ops.location;
  ops_param.properties.lyr_name   = layer_name;
  ops_param.properties.frame      = ops_frames.properties.frame;
  ops_param.properties.segment_id = ones(size(ops_param.properties.frame)) ...
    *double(segment_id);
  
  [status,data] = opsGetCrossovers(sys,ops_param);
  
  %%
  rline = [];
  rows  = [];
  cols  = [];
  gps_time = [];
  season_name = {};
  for i = 1 : length(data.properties.source_point_path_id)
    if ~isnan(data.properties.twtt(i))
      new_rline = find(ops_data.properties.id == data.properties.source_point_path_id(i));
      new_gps_time = ops_data.properties.gps_time(new_rline);
      new_season_name = data.properties.season_name{i};
      if big_matrix.GPS_time(1) <= new_gps_time ...
          && big_matrix.GPS_time(end) >= new_gps_time ...
          && str2double(new_season_name(1:4)) >= 2006
        rline(end+1) = new_rline;
        gps_time(end+1) = new_gps_time;
        season_name{end+1} = new_season_name;
        [~, cols(end+1)] = min(abs(big_matrix.GPS_time - ops_data.properties.gps_time(rline(end))));
        twtt = data.properties.twtt(i);
        twtt = twtt + (data.properties.source_elev(i)-data.properties.cross_elev(i))/(c/2);
        rows(end+1) = round(interp1(big_matrix.Time, 1:length(big_matrix.Time), twtt));
      end
    end
  end
  
  %% Set variable echogram tracking parameters
  gt     = [cols(:).'; rows(:).'];
  bounds = [0 Nx];
  mask   = 90*fir_dec(double(mask), ones(1,5)/3.7);
  mask(mask>=90) = inf;
  
  %% Detrending routine
  if 1
    detect_data = big_matrix.Data;
    % Along track filtering
    detect_data = fir_dec(detect_data,ones(1,5)/5,1);
    % Estimate noise level
    noise_value = mean(mean(detect_data(end-80:end-60,:)));
    % Estimate trend
    trend = mean(detect_data,2);
    trend(trend<noise_value) = noise_value;
    % Subtract trend
    detect_data = bsxfun(@minus,detect_data,trend);
    % Remove bad circular convolution wrap around at end of record
    detect_data(end-70:end,:) = 0;
    big_matrix.Data = detect_data;
  end
  
  %% Call viterbi.cpp
  tic
  big_matrix.Labels = tomo.viterbi(double(detect_data), double(surf_bins), ...
    double(bottom_bin), double(gt), double(mask), double(mu), ...
    double(sigma), double(egt_weight), double(smooth_weight), ...
    double(smooth_var), double(slope), int64(bounds), ...
    double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
    double(mc), double(mc_weight));
  toc
  
  % Used to debug, stop and check tracker output before saving result
  if options.debug
    figure; imagesc(big_matrix.Data); hold on; plot(surf_bins, 'g'); plot(big_matrix.Labels, 'r');
    keyboard
  end
  
  % Interpolate from row number to TWTT
  big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, big_matrix.Labels);
  big_matrix.TWTT(~mask) = big_matrix.Surface(~mask);
  big_matrix.TWTT(isnan(mask)) = NaN;
  
  %% Load labels into OPS using opsCopyLayers
  copy_param = [];
  copy_param.layer_source.existence_check = false;
  copy_param.layer_dest.existence_check = false;
  
  % Set the source (choose one)
  copy_param.layer_source.source = 'custom';
  copy_param.layer_source.gps_time = big_matrix.GPS_time;
  copy_param.layer_source.twtt = big_matrix.TWTT;
  
  copy_param.layer_dest.name = 'mcmc_bottom';
  copy_param.layer_dest.source = 'ops';
  copy_param.layer_dest.existence_check = false;
  
  copy_param.copy_method = 'overwrite';
  
  if strcmpi(copy_param.layer_dest.name,'surface')
    copy_param.gaps_fill.method = 'interp_finite';
  else
    copy_param.gaps_fill.method = 'preserve_gaps';
    copy_param.gaps_fill.method_args = [40 20];
  end
  
  param = merge_structs(param,gRadar);
  fprintf('opsCopyLayers %s (%s)\n', param.day_seg, datestr(now));
  opsCopyLayers(param,copy_param);
  fprintf('  Complete (%s)\n', datestr(now));
  
end

warning('on','all');