function labels_wholeseg = tracker_viterbi(data_struct,param)
% function viterbi_tracker_2D (params, options, data_struct)
% ICE MASK OPTION, layername(444), image_path(420)
% Builds a large matrix containing multiple segments of 2D data and applies
%  the Viterbi layer-tracking program on it.
% Loads crossovers from OPS to use as ground truth for the tracking.
% Other functionality includes multiple supression and
%  a simple detrending technique.
%

%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar;
param = merge_structs(param,gRadar);
physical_constants;

%% Check parameters and set to default if needed
for layer_idx = 1:length(param.layer_tracker.cmds.layer_params)
  
  if isfield(param.layer_tracker.track.viterbi, 'surf_layer_params') && ~isempty(param.layer_tracker.track.viterbi.surf_layer_params)
    surf_layer_params = param.layer_tracker.track.viterbi.surf_layer_params;
  else
    surf_layer_params = [];
  end
  
  if ~isfield(surf_layer_params, 'name') || ~isempty(surf_layer_params.name)
    surf_layer_params.name = 'surface';
  end
%  keyboard;
  %% Distance-to-Ice-Margin model
  if ~isfield(param.layer_tracker.track.viterbi, 'DIM_matrix') || isempty(param.layer_tracker.track.viterbi.DIM_matrix)
    prob.DIM_means = [6.908407 14.603709 22.077745 29.333980 36.375879 43.206908 49.830531 56.250214 62.469423 68.491622 74.320277 79.958853 85.410815 90.679629 95.768760 100.681673 105.421833 109.992706 114.397756 118.640450 122.724252 126.652627 130.429042 134.056960 137.539848 140.881171 144.084393 147.152981 150.090399 152.900113 155.585588 158.150289 160.597681 162.931230 165.154401 167.270659 169.283470 171.196299 173.012610 174.735870 176.369543 177.917095 179.381991 180.767696 182.077676 183.315395 184.484320 185.587915 186.629645 187.612977 188.541374 189.418303 190.247228 191.031615 191.774929 192.480636 193.152200 193.793087 194.406763 194.996691 195.566338 196.119170 196.658650 197.188245 197.711419 198.231639 198.752368 199.277073 199.809219 200.352271 200.909693 201.484953 202.081514 202.702842 203.352402 204.033660 204.750081 205.505130 206.302272 207.144972 208.036696 208.980910 209.981077 211.040664 212.163136 213.351958 214.610596 215.942514 217.351178 218.840053 220.412604 222.072297 223.822597 225.666969 227.608878 229.651790 231.799170 234.054483 236.421194 238.902770];
    prob.DIM_vars  = [19.032807 25.055615 28.640839 32.753438 36.144273 39.329838 42.688930 45.511082 49.036818 52.177650 55.028380 57.516627 60.519048 63.217206 65.211203 67.459337 69.609678 71.543557 73.182822 74.615772 75.628159 77.127086 78.155483 79.447090 80.011376 81.108576 81.618789 82.287856 82.979740 83.561585 84.281769 84.648076 85.290095 85.566969 86.052342 86.487424 86.675812 86.959733 87.181337 87.641261 87.674246 87.947628 87.895269 87.286380 87.202972 86.878606 87.151259 87.477659 88.049960 88.587946 88.515276 89.070799 88.756636 88.345201 87.754785 87.689382 87.240118 86.800999 86.164340 86.085916 85.803664 85.356194 85.831974 85.264038 85.222428 84.898093 84.652262 84.332790 84.249144 83.871931 83.552786 83.233334 82.842279 82.658637 82.008042 81.694151 81.421515 80.901673 80.885452 81.070003 80.524210 80.776716 80.320438 80.445820 80.085639 79.751146 79.557559 78.923447 78.522063 77.525973 77.426494 76.624448 76.855826 77.277564 76.777165 76.716292 75.970217 77.149291 76.900846 76.890210];
    gauss = @(x, mean, var)((1 / (sqrt(2 * pi * (var.^2)))) * (exp(-((x - mean).^2)/(2 * var.^ 2))));
    costm = ones(50, 100);
    for t = 1:50
      for i = 1:100
        costm(t,i) = 10 * exp(-0.01 .* t) - 10 * exp(-0.01 * 50);
      end
    end
    alpha = 0.025;
    beta = 0.05;
    
    DIM_costmatrix = ones(1000, 100);
    for DIM = 1 : 100;
      for T = 1 : 1000
        if alpha / gauss(T, prob.DIM_means(DIM), DIM * beta * prob.DIM_vars(DIM)) >= 200
          DIM_costmatrix(T, DIM) = 200;
        else
          DIM_costmatrix(T, DIM) = alpha/ gauss(T, prob.DIM_means(DIM), DIM * beta * prob.DIM_vars(DIM));
        end
      end
    end
    
    DIM_costmatrix = DIM_costmatrix(15:end,:);
    
    for k = 1:100
      DIM_costmatrix(1:50, k) = DIM_costmatrix(1:50, k) + k * costm(1:50, k);
    end
  else
    clear DIM DIM_costmatrix;
    %
    % DIM = load(fullfile(gRadar.path, param.layer_tracker.tracker.viterbi.DIM_matrix));
    DIM = load(fullfile(gRadar.path, param.layer_tracker.track.viterbi.DIM_matrix));
    DIM_costmatrix = DIM.Layer_tracking_2D_parameters;
    DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));
  end
  
  labels = {};
  big_matrix          = {};
  big_matrix.Data     = [];
  big_matrix.GPS_time = [];
  big_matrix.Lat      = [];
  big_matrix.Lon      = [];
  good_frms           = param.cmd.frms;
  
  try
    Surface = opsLoadLayers(param,surf_layer_params);
  catch ME
    warning(ME.getReport);
    keyboard
  end
  
  for frm = good_frms
%   frm = param.layer_tracker.track.frm;
    if param.layer_tracker.track.debug
      fprintf('\nRunning frame %s_%03d.\n',param.day_seg,frm);
    end
  %  
%     try
%       cur_matrix = mdata(frm);
%     catch ME
%       cur_matrix = mdata;
%     end

    cur_matrix = data_struct;
    
    big_matrix.GPS_time = horzcat(big_matrix.GPS_time, cur_matrix.GPS_time);
    big_matrix.Lat      = horzcat(big_matrix.Lat, cur_matrix.Latitude);
    big_matrix.Lon      = horzcat(big_matrix.Lon, cur_matrix.Longitude);
    big_matrix.Time     = cur_matrix.Time;
    
    %% Image combine
    if ~param.layer_tracker.track.viterbi.custom_combine
      big_matrix.Data = horzcat(big_matrix.Data, cur_matrix.Data);%horzcat(big_matrix.Data, data_struct); % horzcat(big_matrix.Data, cur_matrix.Data); 
    else
      try
        mode = 'get_heights';
        % param.(mode).out_path = param.layer_tracker.track.name;
        param.(mode).out_path = param.layer_tracker.cmds.layer_params(layer_idx).name;
        param.load.frm = frm;
        param.(mode).img_comb = param.combine.img_comb;
        param.(mode).img_comb_mult = 1.15;
        param.(mode).img_comb_bins = 50;
        [Data, big_matrix.Time]    = img_combine(param, mode, Surface);
        big_matrix.Data            = horzcat(big_matrix.Data, Data);
      catch ME
        if param.layer_tracker.track.debug
          fprintf('\nProblem during custom image combining');
        end
        big_matrix.Data = horzcat(big_matrix.Data, cur_matrix.Data);%horzcat(big_matrix.Data, data_struct); % horzcat(big_matrix.Data, cur_matrix.Data);
      end
    end
  end%frm for loop
  
  if param.layer_tracker.track.debug && param.layer_tracker.track.viterbi.custom_combine
    fprintf('\nDone: image combine (%s)', datestr(now,'HH:MM:SS'));
  end
  
  % Catch error with layer not having correctly loaded from OPS
  if(isempty(Surface.gps_time) || any(isnan(Surface.gps_time)) || any(isnan(Surface.twtt)))
    fprintf('\nProblem processing frame %s', param.day_seg);
    return;
  end
  
  if param.layer_tracker.track.debug
    fprintf('\nDone: opsLoadLayers (%s)', datestr(now,'HH:MM:SS'));
  end
  
  %%
  
  big_matrix.Surface = interp_finite(interp1(Surface.gps_time,Surface.twtt,big_matrix.GPS_time));
  
  if param.layer_tracker.track.debug
    fprintf('\nDone: opsAuthenticate, opsGetSegmentInfo, time-range bin interpolation (%s)', datestr(now,'HH:MM:SS'));
  end
  
  big_matrix.Data = lp(big_matrix.Data); %(needed when not using data_struct but instead using mdata)
  surf_bins = round(interp1(big_matrix.Time, 1:length(big_matrix.Time), big_matrix.Surface));
  
  Nx = size(big_matrix.Data, 2);
  slope = diff(surf_bins);
    
  try
    surf_weight = param.layer_tracker.track.viterbi.surf_weight;
  catch ME
    surf_weight = 1000;
  end
  try
    mult_weight = param.layer_tracker.track.viterbi.mult_weight;
  catch ME
    mult_weight = 100;
  end
  try
    mult_weight_decay = param.layer_tracker.track.viterbi.mult_weight_decay;
  catch ME
    mult_weight_decay = 0;
  end
  try
    mult_weight_local_decay = param.layer_tracker.track.viterbi.mult_weight_local_decay;
  catch ME
    mult_weight_local_decay = .8;
  end
  try
    manual_slope = param.layer_tracker.track.viterbi.manual_slope;
  catch ME
    manual_slope = 0;
  end
  try
    max_slope = param.layer_tracker.track.viterbi.max_slope;
  catch ME
    max_slope = -1;
  end
  try
    transition_weight = param.layer_tracker.track.viterbi.transition_weight;
  catch ME
    transition_weight = 1;
  end
  try
    image_mag_weight = param.layer_tracker.track.viterbi.image_mag_weight;
  catch ME
    image_mag_weight = 1;
  end
  try
    gt_weight = -param.layer_tracker.track.viterbi.gt_weight;
  catch ME
    gt_weight = -1;
  end
  try
    gt_cutoff = param.layer_tracker.track.viterbi.gt_cutoff;
  catch ME
    gt_cutoff = 5;
  end
  
  transition_weights = ones(1, Nx-1) * transition_weight;
  manual_slope = ones(1, Nx-1) * manual_slope;
  if ~param.layer_tracker.track.viterbi.use_surf_for_slope
    slope = manual_slope;
  end

  dt = big_matrix.Time(2) - big_matrix.Time(1);
  zero_bin = floor(-big_matrix.Time(1)/dt + 1);

  %% Ice mask calculation
  if param.layer_tracker.track.binary_icemask
    mask = load(param.layer_tracker.track.ice_mask_mat_fn,'R','X','Y','proj');
    [fid,msg] = fopen(param.layer_tracker.track.icemask_fn,'r');
    if fid < 1
      fprintf('Could not open file %s\n', param.layer_tracker.track.ice_mask_bin_fn);
      error(msg);
    end
    mask.maskmask = logical(fread(fid,[length(mask.Y),length(mask.X)],'uint8'));
    fclose(fid);
  else
    [mask.maskmask,mask.R,~] = geotiffread(param.layer_tracker.track.icemask_fn);
    mask.proj = geotiffinfo(param.layer_tracker.track.icemask_fn);
  end
  [mask.x, mask.y] = projfwd(mask.proj, big_matrix.Lat, big_matrix.Lon);
  mask.X = mask.R(3,1) + mask.R(2,1)*(1:size(mask.maskmask,2));
  mask.Y = mask.R(3,2) + mask.R(1,2)*(1:size(mask.maskmask,1));
  [mask.X,mask.Y] = meshgrid(mask.X,mask.Y);
  ice_mask.mask = round(interp2(mask.X, mask.Y, double(mask.maskmask), mask.x, mask.y));
  ice_mask.mask(isnan(ice_mask.mask)) = 1;
  
  if param.layer_tracker.track.debug
    fprintf('\nDone: ice mask calculation (%s)', datestr(now,'HH:MM:SS'));
  end
  
  %% Crossover loading
  if param.layer_tracker.track.viterbi.crossoverload
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
    
    query = sprintf('SELECT rds_segments.id FROM rds_seasons,rds_segments where rds_seasons.name=''%s'' and rds_seasons.id=rds_segments.season_id and rds_segments.name=''%s''',param.season_name,param.day_seg);
    query
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
            && big_matrix.GPS_time(end) >= new_gps_time %...
          %&& str2double(new_season_name(1:4)) >= 2006
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
    if param.layer_tracker.track.debug
      fprintf('\nDone: opsGetCrossovers (%s)', datestr(now,'HH:MM:SS'));
    end
  else
    gt = '';
  end
  
  %% Set variable echogram tracking parameters
  bounds = [0 Nx];
  ice_mask.mask_dist = round(bwdist(ice_mask.mask == 0));
  
  %% Detrending routine
   if 1
%       big_matrix.Data = echo_norm(big_matrix.Data,struct('scale',param.layer_tracker.track.norm.scale));
         % Along track filtering
         big_matrix.Data = fir_dec(big_matrix.Data,ones(1,5)/5,1);
         % Estimate noise level
         noise_value = mean(mean(big_matrix.Data(end-80:end-60,:)));
         % Estimate trend
         trend = mean(big_matrix.Data,2);
         trend(trend<noise_value) = noise_value;
         % Subtract trend
         big_matrix.Data = bsxfun(@minus,big_matrix.Data,trend);
         % Remove bad circular convolution wrap around at end of record
         big_matrix.Data(end-70:end,:) = 0;
      if param.layer_tracker.track.debug
        fprintf('\nDone: detrending (%s)', datestr(now,'HH:MM:SS'));
      end
   end
  
  ice_mask.mask_dist = round(bwdist(ice_mask.mask == 0));
  ice_mask.mask_dist = round(ice_mask.mask_dist ./ 45);
  
  DIM_costmatrix = 18 .* DIM_costmatrix;
  
  %% Call viterbi.cpp
  if param.layer_tracker.track.debug
    fprintf('\nProceeding with Viterbi call... ');
  end
  
  gt_layer = ones(1, size(surf_bins, 2)) * NaN;
  gt_layer(gt(1, :)) = gt(2, :);
  layers = [surf_bins; gt_layer];
  
  surf_costs = ones(1, size(surf_bins, 2)) * surf_weight;
  gt_costs = ones(1, size(surf_bins, 2)) * gt_weight;
  layer_costs = [surf_costs; gt_costs];
  
  surf_cutoffs = ones(1, size(surf_bins, 2)) * -1;
  gt_cutoffs = ones(1, size(surf_bins, 2)) * gt_cutoff;
  layer_cutoffs = [surf_cutoffs; gt_cutoffs];
    
  viterbi_tic = tic;

  labels_wholeseg = tomo.viterbi(double(big_matrix.Data), double(layers), ...
    double(layer_costs), double(layer_cutoffs), double(ice_mask.mask), ...
    double(image_mag_weight), double(slope), double(max_slope), ...
    int64(bounds), double(ice_mask.mask_dist), double(DIM_costmatrix), ...
    double(transition_weights), double(mult_weight), ...
    double(mult_weight_decay), double(mult_weight_local_decay), int64(zero_bin));

  viterbi_toc = toc(viterbi_tic);
  
  fprintf('\nAverage time elapsed per frame: %.2f seconds', ...
    viterbi_toc/length(param.cmd.frms));
  
  if param.layer_tracker.track.debug
    figure; imagesc(big_matrix.Data); hold on; plot(surf_bins, 'g'); plot(labels_wholeseg, 'r');
    legend('Ice-surface', 'Ice-bottom');
    colormap(1-gray(256));
    
    ice_plot = surf_bins;
    ice_plot(~ice_mask.mask) = NaN;
    ice_plot = 2 + (ice_plot ./ ice_plot);
    no_ice_plot = surf_bins;
    no_ice_plot(~isnan(ice_plot)) = NaN;
    no_ice_plot = 2 + (no_ice_plot ./ no_ice_plot);
    
    plot(ice_plot, 'w', 'LineWidth', 3); plot(no_ice_plot, 'k', 'LineWidth', 3);
  end
 
  %% Save image
  if param.layer_tracker.track.save_img
    rbin_ctr = 1;
    for frm = good_frms
      close all;
      frm_size = size(data_struct.(sprintf('data_%s_%03d',param.day_seg,frm)).Data, 2);
      startpt = rbin_ctr; endpt = startpt + frm_size - 1; rbin_ctr = endpt+1;
      f = figure; imagesc(lp(data_struct.(sprintf('data_%s_%03d',param.day_seg,frm)).Data)); hold on;
      colormap(1-gray(256)); plot(labels_wholeseg.data(startpt:endpt), 'r');
      title(sprintf('Viterbi -- %s_%03d -- Time elapsed: %.2f seconds', param.day_seg, frm, ...
        (viterbi_toc/length(param.cmd.frms))), 'Interpreter', 'none');
      path = [param.layer_tracker.track.save_img_path, sprintf('%s_%03d',param.day_seg,frm)];
      print(f, path, param.layer_tracker.track.save_img_format);
    end
  end
  
  param.layer_tracker.track.ops_write = false;
  if param.layer_tracker.track.ops_write
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
    %   copy_param.layer_dest.name = param.layer_tracker.track.viterbi.layername;
    %   copy_param.layer_dest.source = param.layer_tracker.track.layer_dest_source;
    %
    %   if strcmp(copy_param.layer_dest.source, 'layerdata')
    %     copy_param.layer_dest.layerdata_source = param.layer_tracker.track.layer_dest_layerdata_source;
    %     copy_param.layer_dest.echogram_source = param.layer_tracker.track.layer_dest_echogram_source;
    %   end
    
    copy_param.layer_dest = param.layer_tracker.cmds.layer_params(layer_idx);
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
  
  %% Write additional file
  if param.layer_tracker.track.save_add_f
    rbin_ctr = 1;
    for frm = good_frms
      frm_size = size(data_struct.(sprintf('data_%s_%03d',param.day_seg,frm)).Data, 2);
      startpt = rbin_ctr; endpt = startpt + frm_size - 1; rbin_ctr = endpt+1;
      labels.bot = labels_wholeseg(startpt:endpt);
      labels.toc = viterbi_toc/length(param.cmd.frms);
      path = [param.layer_tracker.track.save_add_f_path, sprintf('%s_%03d.mat',param.day_seg,frm)];
      save(path, 'labels');
      fprintf('\nSaved additional file to %s', path);
    end
  end
end

