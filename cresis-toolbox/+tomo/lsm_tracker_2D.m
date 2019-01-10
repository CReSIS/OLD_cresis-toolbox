function labels = lsm_tracker_2D (params, options, data_struct)
% function lsm_tracker_2D (params, options)
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
  
  data_fn_dir = ct_filename_out(param, options.name, '');
  for frm = param.cmd.frms
    
    lsm_tic = tic;
    
    data_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
    data_fn      = fullfile(data_fn_dir, data_fn_name);
    data         = data_struct.(sprintf('data_%s_%03d',param.day_seg,frm));
    imds         = imageDatastore(data_fn, 'FileExtensions', '.mat');
    obj          = tomo.LSMObject(imds.Files);
    obj.setLSMOptions('y', 240, 'dy', 10, 'outerIter', options.lsm.numOuterIter);
    
    %% Ice mask calculation from geotiff
    if isempty(options.geotiff_fn)
      mask = inf * ones([1 Nx]);
    else
      [mask,R,~] = geotiffread(options.geotiff_fn);
      proj = geotiffinfo(options.geotiff_fn);
      [points.x, points.y] = projfwd(proj, data.Latitude, data.Longitude);
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
        [points.x, points.y] = projfwd(proj, data.Latitude, data.Longitude);
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
    
    [flag, Labels.top, Labels.bot] = obj.runLSM();
    fprintf('\nDone: runLSM() call. ');

    if flag == 0
      warning('LSM failed to converge to two layers after %d inner iterations. Trying again with %d inner iterations...', 2*options.lsm.numOuterIter, 2*options.lsm.maxOuterIter);
      obj = tomo.LSMObject(imds.Files);
      obj.setLSMOptions('y', 240, 'dy', 10, 'outerIter', options.lsm.maxOuterIter);
      [flag, Labels.top, Labels.bot] = obj.runLSM();
      fprintf('\nDone: second runLSM() call. ');
      if flag == 0
        warning('LSM failed to converge to two layers after 2000 inner iterations. Skipping...');
        continue;
      end
    end

    Labels.top.y = interp1(Labels.top.x, Labels.top.y, 1:length(data.Bottom));
    Labels.bot.y = interp1(Labels.bot.x, Labels.bot.y, 1:length(data.Bottom));
    
    lsm_toc = toc(lsm_tic);
    
    if options.debug
      figure; imagesc(lp(data.Data)); colormap(1 - gray(256)); hold on;
      plot(Labels.top.y, 'g'); plot(Labels.bot.y, 'r');
      legend('Ice-surface', 'Ice-bottom');
      keyboard
    end
    
    if options.ops_write
      %% Write surface layer
      % Interpolate from row number to TWTT
      data.TWTT = interp1(1:length(data.Time), data.Time, Labels.top.y);
      data.TWTT(~mask) = data.Surface(~mask);
      data.TWTT(isnan(mask)) = NaN;
      
      %% Load labels into OPS using opsCopyLayers
      copy_param = [];
      copy_param.layer_source.existence_check = false;
      copy_param.layer_dest.existence_check = false;
      
      % Set the source
      copy_param.layer_source.source = 'custom';
      copy_param.layer_source.gps_time = data.GPS_time;
      copy_param.layer_source.twtt = data.TWTT;
      
      copy_param.layer_dest.name = options.lsm.lyrtop;
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
      data.TWTT = interp1(1:length(data.Time), data.Time, Labels.bot.y);
      data.TWTT(~mask) = data.Surface(~mask);
      data.TWTT(isnan(mask)) = NaN;
      
      %% Load labels into OPS using opsCopyLayers
      copy_param = [];
      copy_param.layer_source.existence_check = false;
      copy_param.layer_dest.existence_check = false;
      
      % Set the source
      copy_param.layer_source.source = 'custom';
      copy_param.layer_source.gps_time = data.GPS_time;
      copy_param.layer_source.twtt = data.TWTT;
      
      copy_param.layer_dest.name = options.lsm.lyrbot;
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
      Labels.top.y;
    labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).bot = ...
      Labels.bot.y;
    labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).toc = ...
      lsm_toc;
    
  end
end
end