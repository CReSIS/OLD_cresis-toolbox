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
  frames = frames_load(param);
  
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
    
    fprintf('\nLSM: Running frame %s_%03d (%s)\n',param.day_seg, frm, datestr(now,'HH:MM:SS'));
    
    lsm_tic = tic;
    
    try
      data = data_struct.(sprintf('data_%s_%03d',param.day_seg,frm));
    catch ME
      fprintf('\nProblem with frame %s_%03d, verify.\n' ,param.day_seg,frm);
      continue;
    end

    data_fn_name = sprintf('Data_%s_%03d.mat',param.day_seg,frm);
    data_fn      = fullfile(data_fn_dir, data_fn_name);
    imds         = imageDatastore(data_fn, 'FileExtensions', '.mat');
    obj          = tomo.LSMObject(imds.Files);
    obj.setLSMOptions('y', 240, 'dy', 10, 'outerIter', options.lsm.numOuterIter);
    
    [flag, Labels.top, Labels.bot] = obj.runLSM();
    
    if flag == 0
      warning('LSM failed to converge to two layers after %d inner iterations. Trying again with %d inner iterations...', 2*options.lsm.numOuterIter, 2*options.lsm.maxOuterIter);
      obj = tomo.LSMObject(imds.Files);
      obj.setLSMOptions('y', 240, 'dy', 10, 'outerIter', options.lsm.maxOuterIter);
      [flag, Labels.top, Labels.bot] = obj.runLSM();
      fprintf('\nDone: second runLSM() call. ');
    end
    
    Labels.top.y = interp1(Labels.top.x, Labels.top.y, 1:length(data.Bottom), 'linear', 'extrap');
    Labels.bot.y = interp1(Labels.bot.x, Labels.bot.y, 1:length(data.Bottom), 'linear', 'extrap')';
    
    lsm_toc = toc(lsm_tic);
    
    if options.debug
      figure; imagesc(lp(data.Data)); colormap(1 - gray(256)); hold on;
      plot(Labels.top.y, 'g'); plot(Labels.bot.y, 'r');
      legend('Ice-surface', 'Ice-bottom');
      keyboard
    end

    if options.ops_write
      warning('off')
      %% Set write options
      % Load labels into OPS using opsCopyLayers
      copy_param = [];
      copy_param.layer_source.existence_check = false;
      copy_param.layer_dest.existence_check = false;
      
      % Set the source
      copy_param.layer_source.source = 'custom';      
      copy_param.layer_dest.source = options.layer_dest_source;
      
      if strcmp(copy_param.layer_dest.source, 'layerdata')
        copy_param.layer_dest.layerdata_source = options.layer_dest_layerdata_source;
        copy_param.layer_dest.echogram_source = options.layer_dest_echogram_source;
      end
      
      copy_param.copy_method = 'overwrite';
      
      copy_param.gaps_fill.method = 'preserve_gaps';
      copy_param.gaps_fill.method_args = [40 20];
      
      param = merge_structs(param,gRadar);
      
      %% Write surface layer
      % Interpolate from row number to TWTT
      data.TWTT = interp1(1:length(data.Time), data.Time, Labels.top.y);      
      copy_param.layer_source.gps_time = data.GPS_time;
      copy_param.layer_source.twtt = data.TWTT;
      
      % Load labels into OPS using opsCopyLayers
      copy_param.layer_dest.name = options.lsm.lyrtop;
      fprintf('\nopsCopyLayers %s (%s) [TOP]', param.day_seg, datestr(now)); 
      opsCopyLayers(param,copy_param);
      
      %% Write bottom layer
      % Interpolate from row number to TWTT
      data.TWTT = interp1(1:length(data.Time), data.Time, Labels.bot.y);
      copy_param.layer_source.gps_time = data.GPS_time;
      copy_param.layer_source.twtt = data.TWTT;
      
      % Load labels into OPS using opsCopyLayers
      copy_param.layer_dest.name = options.lsm.lyrbot;
      fprintf('\nopsCopyLayers %s (%s) [BOTTOM]', param.day_seg, datestr(now)); 
      opsCopyLayers(param,copy_param);
      
      fprintf('\n  Complete (%s)\n', datestr(now));
      warning('on');
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