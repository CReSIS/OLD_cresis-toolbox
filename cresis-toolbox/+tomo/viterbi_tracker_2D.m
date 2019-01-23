function [labels, surf_bins] = viterbi_tracker_2D (params, options, data_struct)
% function viterbi_tracker_2D (params, options, data_struct)
%
% Builds a large matrix containing multiple segments of 2D data and applies
%  the Viterbi layer-tracking program on it.
% Loads crossovers from OPS to use as ground truth for the tracking.
% Other functionality includes multiple supression and
%  a simple detrending technique.
%
% See also: run_tracker_2D.m
%
% Authors: Victor Berger

%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar;

physical_constants;
clear('param_override');

% Input checking
if ~exist('params','var')
  error('Use run_tracker_2D: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar, param_override);
else
  param_override = gRadar;
end

%% Check parameters and set to default if needed
if isfield(options.viterbi, 'bottom_bin')
  bottom_bin = options.viterbi.bottom_bin;
else
  bottom_bin = -1;
end

if isfield(options.viterbi, 'egt_weight')
  egt_weight = options.viterbi.egt_weight;
else
  egt_weight = -1;
end

if isfield(options.viterbi, 'mu_size')
  mu_size = options.viterbi.mu_size;
else
  mu_size = 31;
end

if isfield(options.viterbi, 'mu_thr')
  mu_thr = options.viterbi.mu_thr;
else
  mu_thr = -30;
end

if isfield(options.viterbi, 'mu')
  mu = options.viterbi.mu;
else
  mu = log10(exp(-(-(mu_size-1)/2 : (mu_size-1)/2).^4/1));
  mu(mu < mu_thr) = mu_thr;
  mu              = mu - mean(mu);
end

if isfield(options.viterbi, 'sigma')
  sigma = options.viterbi.sigma;
else
  sigma = sum(abs(mu))/10*ones(1, mu_size);
end

if isfield(options.viterbi, 'smooth_var')
  smooth_var = options.viterbi.smooth_var;
else
  smooth_var = Inf;
end

if isfield(options.viterbi, 'repulsion')
  repulsion = options.viterbi.repulsion;
else
  repulsion = 150000;
end

if isfield(options.viterbi, 'smooth_weight')
  smooth_weight = options.viterbi.smooth_weight;
else
  smooth_weight = 5;
end

if isfield(options.viterbi, 'ice_bin_thr')
  ice_bin_thr = options.viterbi.ice_bin_thr;
else
  ice_bin_thr = 10;
end

%% Distance-to-Ice-Margin model
% Obtained from geostatistical analysis of 2014 Greenland P3
prob.DIM_means = [    16.0582   27.3381   34.0492   40.5714...
  46.9463   52.3826   58.4664   63.7750   70.5143   76.3140...
  81.3519   86.9523   92.2285   98.1430   102.8310  107.0558...
  112.4538  116.3923  121.0109  125.1486  128.7611  133.3286...
  136.3500  139.5058  142.3812  146.1565  148.3781  151.1398...
  153.5642  155.2606  157.4824  159.8529  161.0239  163.1799...
  164.2849  166.1013  166.0728  167.3132  168.1448  169.1323...
  169.7559  170.3869  171.4264  171.8926  171.5201  171.8870...
  171.6407  172.3505  171.3533  171.6161];

prob.DIM_vars  = [     20.5619   24.9082   28.0037   32.0840...
  35.6021   38.9544    42.4034   45.3588   48.9714   52.1360...
  55.0826   57.5144    60.3847   63.1485   65.1199   67.3734...
  69.3662   71.2849    72.9471   74.3759   75.5521   76.9737...
  77.9961   79.3596    79.9999   81.0342   81.6340   82.2424...
  82.9658   83.4794    84.1632   84.4168   85.0014   85.3065...
  85.7757   86.1880    86.3563   86.5577   86.7289   87.0748...
  87.1360   87.2473    87.2828   86.6350   86.5453   86.2989...
  86.4736   86.7318   87.1606   87.6966];

costm = ones(1201, 101);
fd = 18 * fir_dec(1:101, ones(1,5)/3.7.');

for t = 1:1200
  for i = 1:101
    costm(t,i) = 200 * exp(-0.075 .* t) - 200 * exp(-0.075 * 50);
    if t > fd(i)
      costm(t,i) = 200;
    end
  end
end

prob.DIM_means = [prob.DIM_means prob.DIM_means(end)* ones(1,50)];
prob.DIM_vars  = [prob.DIM_vars  prob.DIM_vars(end) * ones(1,50)];

DIM_costmatrix = ones(1200, 100);
for DIM = 1 : 100
  for T = 1 : 1200
    var = DIM * 0.05 * prob.DIM_vars(end);
    cost = 0.1 ./ normpdf(T, prob.DIM_means(DIM), var);
    cost(cost > 800) = 800;
    DIM_costmatrix(T, DIM) = cost;
  end
end

for k = 1:100
  DIM_costmatrix(1:50, k) = DIM_costmatrix(1:50, k) + 0.06 * k * costm(1:50, k);
end

DIM_costmatrix(DIM_costmatrix < 0) = 0;
DIM_costmatrix(DIM_costmatrix > 800) = 800;

DIM_costmatrix = DIM_costmatrix(3:end, :);
DIM_costmatrix(end+1, :) = DIM_costmatrix(end,:);
DIM_costmatrix(end+1, :) = DIM_costmatrix(end,:);
DIM_costmatrix = DIM_costmatrix ./ 4;

% Visualization of DIM_costmatrix
if 0
  figure; (imagesc(DIM_costmatrix)); colorbar; hold on;
  xlabel('Distance to nearest ice-margin [m]');
  ylabel('Ice thickness distribution [range-bin]')
  title('Added cost');
end

labels = {};

%% Process each of the segments
% =====================================================================
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
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
  
  viterbi_tic = tic;
  
  for frm = param.cmd.frms
    if options.debug
      fprintf('Running frame %s_%03d.\n',param.day_seg,frm);
    end
    
    try
      cur_matrix = data_struct.(sprintf('data_%s_%03d',param.day_seg,frm));
    catch ME
      continue;
    end
    
    big_matrix.GPS_time = horzcat(big_matrix.GPS_time, cur_matrix.GPS_time);
    big_matrix.Lat      = horzcat(big_matrix.Lat, cur_matrix.Latitude);
    big_matrix.Lon      = horzcat(big_matrix.Lon, cur_matrix.Longitude);
    big_matrix.Time     = cur_matrix.Time;
    
    %% Image combine
    if 1
      big_matrix.Data   = horzcat(big_matrix.Data, cur_matrix.Data);
    else
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
    end
  end
  
  if options.debug
    fprintf('\nDone: image combine (%s)', datestr(now,'HH:MM:SS'));
  end
  
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
    fprintf('\nProblem processing frame %s', param.day_seg);
    continue;
  end
  
  if options.debug
    fprintf('\nDone: opsLoadLayers (%s)', datestr(now,'HH:MM:SS'));
  end
  
  %%
  opsAuthenticate(param,false);
  layer_name                   = 'bottom';
  sys                          = ct_output_dir(param.radar_name);
  ops_param                    = struct('properties',[]);
  ops_param.properties.season  = param.season_name;
  ops_param.properties.segment = param.day_seg;
  [~,ops_frames]               = opsGetSegmentInfo(sys,ops_param);
  
  ops_param = struct('properties',[]);
  ops_param.properties.location = param.post.ops.location;
  ops_param.properties.season = param.season_name;
  ops_param.properties.start_gps_time = ops_frames.properties.start_gps_time(1);
  ops_param.properties.stop_gps_time = ops_frames.properties.stop_gps_time(end);
  ops_param.properties.nativeGeom = true;
  [~,ops_data] = opsGetPath(sys,ops_param);
  
  big_matrix.Surface = interp_finite(interp1(Surface.gps_time,Surface.twtt,big_matrix.GPS_time));
  
  if options.debug
    fprintf('\nDone: opsAuthenticate, opsGetSegmentInfo, time-range bin interpolation (%s)', datestr(now,'HH:MM:SS'));
  end
  
  big_matrix.Data = lp(big_matrix.Data);
  surf_bins = round(interp1(big_matrix.Time, 1:length(big_matrix.Time), big_matrix.Surface));
  
  %% Top suppression
  if options.viterbi.top_sup
    topbuffer = 10;
    botbuffer = 30;
    filtvalue = 50;
    for rline = 1 : size(big_matrix.Data, 2)
      column_chunk = big_matrix.Data(round(surf_bins(rline) - topbuffer) : ...
        round(surf_bins(rline) + botbuffer), rline);
      big_matrix.Data(round(surf_bins(rline) - topbuffer) : ...
        round(surf_bins(rline) + botbuffer), rline) = imgaussfilt(column_chunk, filtvalue);
    end
    
    if options.debug
      fprintf('\nDone: top suppression (%s)', datestr(now,'HH:MM:SS'));
    end
  end
  
  %% Multiple suppression
  if options.viterbi.mult_sup
    topbuffer = 10;
    botbuffer = 5;
    filtvalue = 100;
    for rline = 1 : size(big_matrix.Data, 2)
      try
        column_chunk = big_matrix.Data(round(2*surf_bins(rline) - topbuffer) : ...
          round(2*surf_bins(rline) + botbuffer), rline);
        big_matrix.Data(round(2*surf_bins(rline) - topbuffer) : ...
          round(2*surf_bins(rline) + botbuffer), rline) = imgaussfilt(column_chunk, filtvalue);
      catch ME
        continue;
      end
    end
    
    if options.debug
      fprintf('\nDone: multiple suppression (%s)', datestr(now,'HH:MM:SS'));
    end
  end
  
  Nx             = size(big_matrix.Data, 2);
  slope          = round(diff(big_matrix.Surface));
  viterbi_weight = ones([1 Nx]);
  
  %% Ice mask calculation
  if isempty(options.icemask_fn)
    ice_mask.mask = inf * ones([1 Nx]);
  else
    if ~options.binary_icemask
      try
        [ice_mask.mask,ice_mask.R,~] = geotiffread(options.icemask_fn);
        ice_mask.proj = geotiffinfo(options.icemask_fn);
      catch ME
        fprintf('\nProblem loading ice mask file, check.\n');
        keyboard
      end
      [points.x, points.y] = projfwd(ice_mask.proj, big_matrix.Lat, big_matrix.Lon);
      ice_mask.X = ice_mask.R(3,1) + ice_mask.R(2,1)*(1:size(ice_mask.mask,2));
      ice_mask.Y = ice_mask.R(3,2) + ice_mask.R(1,2)*(1:size(ice_mask.mask,1));
      if 0
        figure;
        imagesc(ice_mask.X,ice_mask.Y,ice_mask.mask);
        hold on;
        plot(points.x,points.y)
      end
      [ice_mask.X,ice_mask.Y] = meshgrid(ice_mask.X,ice_mask.Y);
      ice_mask.fl_mask = round(interp2(ice_mask.X, ice_mask.Y, double(ice_mask.mask), points.x, points.y));
      ice_mask.mask = ice_mask.fl_mask;
      ice_mask.mask(isnan(ice_mask.mask)) = 1;
      
      % Useful for Antarctica seasons:
      if 0
        [ice_mask.mask,ice_mask.R,~] = geotiffread(options.icemask_fn);
        [ice_mask.mask2, ~, ~] = geotiffread(options.icemask2_fn);
        ice_mask.proj = geotiffinfo(options.icemask_fn);
        [points.x, points.y] = projfwd(ice_mask.proj, big_matrix.Lat, big_matrix.Lon);
        ice_mask.X = ice_mask.R(3,1) + ice_mask.R(2,1)*(1:size(ice_mask.mask,2));
        ice_mask.Y = ice_mask.R(3,2) + ice_mask.R(1,2)*(1:size(ice_mask.mask,1));
        if 0
          figure;
          imagesc(ice_mask.X,ice_mask.Y,ice_mask.mask);
          hold on;
          plot(points.x,points.y)
        end
        [ice_mask.X,ice_mask.Y] = meshgrid(ice_mask.X,ice_mask.Y);
        ice_mask.fl_mask = interp2(ice_mask.X, ice_mask.Y, double(ice_mask.mask), points.x, points.y);
        ice_mask.fl_mask = round(interp2(ice_mask.X, ice_mask.Y, double(ice_mask.mask), points.x, points.y));
        ice_mask.mask = ice_mask.fl_mask;
        ice_mask.mask(isnan(ice_mask.mask)) = 1;
        ice_mask.fl_mask2 = interp2(ice_mask.X, ice_mask.Y, double(ice_mask.mask2), points.x, points.y);
        ice_mask.fl_mask2 = round(interp2(ice_mask.X, ice_mask.Y, double(ice_mask.mask2), points.x, points.y));
        ice_mask.mask2 = ice_mask.fl_mask2;
        ice_mask.mask2(isnan(ice_mask.mask2)) = 1;
        ice_mask.f_mask = zeros(size(ice_mask.mask));
        for i = 1:length(ice_mask.f_mask)
          if ((ice_mask.mask(1, i) == 0 || ice_mask.mask(1, i) == 1) && (ice_mask.mask2(1, i) == 127))
            ice_mask.f_mask(1, i) = 1;
          end
        end
        ice_mask.mask = ice_mask.f_mask;
        if max(max(ice_mask.mask)) > 100
          ice_mask.tmp_mask = zeros(size(ice_mask.mask));
          ice_mask.tmp_mask(ice_mask.mask == 0) = 1;
          ice_mask.mask = ice_mask.tmp_mask;
        end
      end
    else
      [ice_mask_fn_dir, ice_mask_fn_name] = fileparts(options.icemask_fn);
      ice_mask_mat_fn = fullfile(ice_mask_fn_dir,[ice_mask_fn_name '.mat']);
      ice_mask = load(ice_mask_mat_fn,'R','X','Y','proj');
      
      [fid,msg] = fopen(options.icemask_fn,'r');
      if fid < 1
        fprintf('Could not open file %s\n', ice_mask_bin_fn);
        error(msg);
      end
      ice_mask.mask = logical(fread(fid,[length(ice_mask.Y),length(ice_mask.X)],'uint8'));
      fclose(fid);
      [points.x, points.y] = projfwd(ice_mask.proj, big_matrix.Lat, big_matrix.Lon);
      
      if 0
        figure;
        imagesc(ice_mask.X,ice_mask.Y,ice_mask.mask);
        hold on;
        plot(points.x,points.y)
      end
      
      ice_mask.mask = round(interp2(ice_mask.X, ice_mask.Y, double(ice_mask.mask), points.x, points.y));
      ice_mask.mask(isnan(ice_mask.mask)) = 1;
      ice_mask.mask_dist = round(bwdist(ice_mask.mask == 0));
    end
  end
  
  if options.debug
    fprintf('\nDone: ice mask calculation (%s)', datestr(now,'HH:MM:SS'));
  end
  
  %% Crossover loading
  if options.viterbi.crossoverload
    query = sprintf('SELECT rds_segments.id FROM rds_seasons,rds_segments where rds_seasons.name=''%s'' and rds_seasons.id=rds_segments.season_id and rds_segments.name=''%s''',param.season_name,param.day_seg);
    [~,tables] = opsQuery(query);
    segment_id = tables{1};
    
    ops_param                       = struct('properties',[]);
    ops_param.properties.location   = param.post.ops.location;
    ops_param.properties.lyr_name   = layer_name;
    ops_param.properties.frame      = ops_frames.properties.frame;
    ops_param.properties.segment_id = ones(size(ops_param.properties.frame)) ...
      *double(segment_id);
    
    [~,data] = opsGetCrossovers(sys,ops_param);
    
    %% Load and align crossovers
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
    gt     = [cols(:).'; rows(:).'];
    if options.debug
      fprintf('\nDone: opsGetCrossovers (%s)', datestr(now,'HH:MM:SS'));
    end
  else
    gt = '';
  end
  
  %% Set variable echogram tracking parameters
  bounds = [0 Nx];
  ice_mask.mask_dist = round(bwdist(ice_mask.mask == 0));
  ice_mask.mask   = 90*fir_dec(double(ice_mask.mask), ones(1,5)/3.7);
  ice_mask.mask(ice_mask.mask>=90) = Inf;
  transition_mu = -1 * ones(1, size(big_matrix.Data, 2));
  transition_sigma = -1 * ones(1, size(big_matrix.Data, 2));
  
  %% Detrending routine
  if options.viterbi.detrending
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
  
  if options.debug
    fprintf('\nDone: detrending (%s)', datestr(now,'HH:MM:SS'));
  end
  
  %% Call viterbi.cpp
  if options.debug
    fprintf('\nProceeding with Viterbi call... ');
  end
 
  labels_wholeseg = tomo.viterbi(double(big_matrix.Data), double(surf_bins), ...
    double(bottom_bin), double(gt), double(ice_mask.mask), double(mu), ...
    double(sigma), double(egt_weight), double(smooth_weight), ...
    double(smooth_var), double(slope), int64(bounds), ...
    double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
    double(ice_mask.mask_dist), double(DIM_costmatrix), ...
    double(transition_mu), double(transition_sigma));
  
  viterbi_toc = toc(viterbi_tic);
  if options.debug
    fprintf(' done. (%s)\n', datestr(now));
  end
  
  % Stop and check tracker output before saving result
  if options.debug
    figure; imagesc(big_matrix.Data); hold on; plot(surf_bins, 'g'); plot(labels_wholeseg, 'r');
    legend('Ice-surface', 'Ice-bottom');
    colormap(1-gray(256));
    keyboard
  end
  
  if options.ops_write
    warning('off');
    % Interpolate from row number to TWTT
    big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, labels_wholeseg);
    big_matrix.TWTT(~ice_mask.mask) = big_matrix.Surface(~ice_mask.mask);
    big_matrix.TWTT(isnan(ice_mask.mask)) = NaN;
    
    %% Load labels into OPS using opsCopyLayers
    copy_param = [];
    copy_param.layer_source.existence_check = false;
    copy_param.layer_dest.existence_check = false;
    
    % Set the source
    copy_param.layer_source.source = 'custom';
    copy_param.layer_source.gps_time = big_matrix.GPS_time;
    copy_param.layer_source.twtt = big_matrix.TWTT;
    
    % Set the destination
    copy_param.layer_dest.name = options.viterbi.layername;
    copy_param.layer_dest.source = options.layer_dest_source;
    
    if strcmp(copy_param.layer_dest.source, 'layerdata')
      copy_param.layer_dest.layerdata_source = options.layer_dest_layerdata_source;
      copy_param.layer_dest.echogram_source = options.layer_dest_echogram_source;
    end
    
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
    warning('on');
  end
  
  %% Generate labels struct with separate frame results
  idx_ctr = 0;
  for frm = param.cmd.frms
    try
      labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).bot = ...
        labels_wholeseg(idx_ctr + 1 : idx_ctr + length(data_struct.(sprintf('data_%s_%03d',param.day_seg,frm)).Bottom));
      labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).toc = ...
        viterbi_toc/length(param.cmd.frms);
      idx_ctr = idx_ctr + length(data_struct.(sprintf('data_%s_%03d',param.day_seg,frm)).Bottom);
    catch ME
      continue;
    end
  end
end