function labels = mcmc_tracker_2D (params, options, data_struct)
% function mcmc_tracker_2D (params, options, data_struct)
%
% See also: run_tracker_2D.m
%
% Authors: Victor Berger

%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar;
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

if isfield(options.mcmc, 'top_smooth')
  top_smooth = options.mcmc.top_smooth;
else
  top_smooth = 1000;
end

if isfield(options.mcmc, 'bottom_smooth')
  bottom_smooth = options.mcmc.bottom_smooth;
else
  bottom_smooth = 1000;
end

if isfield(options.mcmc, 'top_peak')
  top_peak = options.mcmc.top_peak;
else
  top_peak = 0.5;
end

if isfield(options.mcmc, 'bottom_peak')
  bottom_peak = options.mcmc.bottom_peak;
else
  bottom_peak = 0.5;
end

if isfield(options.mcmc, 'repulsion')
  repulsion = options.mcmc.repulsion;
else
  repulsion = 10;
end

labels = {};

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
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
  
  for frm = param.cmd.frms
    
    mcmc_tic = tic;
    
    if options.debug
      fprintf('\nMCMC: Running frame %s_%03d (%s)',param.day_seg, frm, datestr(now,'HH:MM:SS'));
    end
    
    big_matrix      = data_struct.(sprintf('data_%s_%03d',param.day_seg,frm));
    big_matrix.Data = lp(big_matrix.Data);
    big_matrix.Data = 255  * ((big_matrix.Data - min(big_matrix.Data(:))) ...
      ./ (max(big_matrix.Data(:) - min(big_matrix.Data(:)))));
    big_matrix.Data = uint8(repmat(big_matrix.Data, [1 1 3]));
    pts1            = [];
    pts2            = [];
    
    %% Ice mask calculation from geotiff
    if isempty(options.geotiff_fn)
      mask = inf * ones([1 Nx]);
    else
      [mask,R,~] = geotiffread(options.geotiff_fn);
      proj = geotiffinfo(options.geotiff_fn);
      [points.x, points.y] = projfwd(proj, big_matrix.Latitude, big_matrix.Longitude);
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
      
      % Useful for Antarctica seasons:
      if 0
        [mask,R,~] = geotiffread(options.geotiff_fn);
        [mask2, ~, ~] = geotiffread(options.geotiff2_fn);
        proj = geotiffinfo(options.geotiff_fn);
        [points.x, points.y] = projfwd(proj, big_matrix.Latitude, big_matrix.Longitude);
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
    
    if strcmp(options.mcmc.alg, 'MCMC')
      if(numel(pts1) + numel(pts2) > 0 )
        [big_matrix.Labels, big_matrix.lower, big_matrix.upper] = ...
          tomo.RJ_MCMC(double(big_matrix.Data(:,:,1)), pts1, pts2);
      else
        tic
        [big_matrix.Labels, big_matrix.lower, big_matrix.upper] = ...
          tomo.RJ_MCMC(double(big_matrix.Data(:,:,1)));
        toc
      end
    elseif strcmp(options.mcmc.alg, 'HMM')
      opts = [top_smooth bottom_smooth top_peak bottom_peak repulsion];
      if(numel(pts1) + numel(pts2) > 0 )
        [~, big_matrix.Labels] = tomo.stereo(1, double(big_matrix.Data(:,:,1)), opts, pts1, pts2);
      else
        [~, big_matrix.Labels] = tomo.stereo(1, double(big_matrix.Data(:,:,1)), opts);
      end
    else
      fprintf('\nUnrecognized algorithm (MCMC or HMM)');
      keyboard
    end
    
    mcmc_toc = toc(mcmc_tic);
    
    if options.debug
      figure; imshow(big_matrix.Data); hold on;
      plot(big_matrix.Labels(1, :), 'g'); plot(big_matrix.Labels(2, :), 'r');
      legend('Ice-surface', 'Ice-bottom');
      keyboard
    end
    
    if options.ops_write
      %% %% Write surface layer
      % Interpolate from row number to TWTT
      big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, big_matrix.Labels(1,:));
      big_matrix.TWTT(~mask) = big_matrix.Surface(~mask);
      big_matrix.TWTT(isnan(mask)) = NaN;
      
      %% Load labels into OPS using opsCopyLayers
      copy_param = [];
      copy_param.layer_source.existence_check = false;
      copy_param.layer_dest.existence_check = false;
      
      % Set the source
      copy_param.layer_source.source = 'custom';
      copy_param.layer_source.gps_time = big_matrix.GPS_time;
      copy_param.layer_source.twtt = big_matrix.TWTT;
      
      copy_param.layer_dest.name = options.mcmc.layername;
      copy_param.layer_dest.source = 'ops';
      
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
      
      %% Write bottom layer
      % Interpolate from row number to TWTT
      big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, big_matrix.Labels(2,:));
      big_matrix.TWTT(~mask) = big_matrix.Surface(~mask);
      big_matrix.TWTT(isnan(mask)) = NaN;
      
      %% Load labels into OPS using opsCopyLayers
      copy_param = [];
      copy_param.layer_source.existence_check = false;
      copy_param.layer_dest.existence_check = false;
      
      % Set the source
      copy_param.layer_source.source = 'custom';
      copy_param.layer_source.gps_time = big_matrix.GPS_time;
      copy_param.layer_source.twtt = big_matrix.TWTT;
      
      copy_param.layer_dest.name = options.mcmc.layername;
      copy_param.layer_dest.source = 'ops';
      
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
    
    labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).top = ...
      big_matrix.Labels(1,:);
    labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).bot = ...
      big_matrix.Labels(2,:);
    labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).toc = ...
      mcmc_toc;
    
  end
end
end