function labels = tracker_lsm(data_struct,param)
%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar
% params = merge_structs(params,gRadar);
  % Load frames file
  load(ct_filename_support(param,'','frames'));
  
  if isempty(param.cmd.frms)
    param.cmd.frms = 1:length(frames.frame_idxs);
  end
  
  data_fn_dir = ct_filename_out(param, param.layer_tracker.track.name, '');
  
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
    imds         = datastore(data_fn, 'Type','image','FileExtensions', '.mat');
    obj          = tomo.LSMObject(imds.Files);
    obj.setLSMOptions('y', param.layer_tracker.track.lsm.y, 'dy', param.layer_tracker.track.lsm.dy, 'outerIter', param.layer_tracker.track.lsm.numOuterIter);
    
    [flag, Labels.top, Labels.bot] = obj.runLSM();
    
    if flag == 0
      warning('LSM failed to converge to two layers after %d inner iterations. Trying again with %d inner iterations...', 2*param.layer_tracker.track.lsm.numOuterIter, 2*param.layer_tracker.track.lsm.maxOuterIter);
      obj = tomo.LSMObject(imds.Files);
      obj.setLSMOptions('y', 240, 'dy', 10, 'outerIter', param.layer_tracker.track.lsm.maxOuterIter);
      [flag, Labels.top, Labels.bot] = obj.runLSM();
      fprintf('\nDone: second runLSM() call. ');
    end
    
    Labels.top.y = interp1(Labels.top.x, Labels.top.y, 1:length(data.Bottom), 'linear', 'extrap');
    Labels.bot.y = interp1(Labels.bot.x, Labels.bot.y, 1:length(data.Bottom), 'linear', 'extrap')';
    
    lsm_toc = toc(lsm_tic);
    
    if param.layer_tracker.track.debug
      figure; imagesc(lp(data.Data)); colormap(1 - gray(256)); hold on;
      plot(Labels.top.y, 'g'); plot(Labels.bot.y, 'r');
      legend('Ice-surface', 'Ice-bottom');
      keyboard
    end

    if param.layer_tracker.track.ops_write
      warning('off')
      %% Set write options
      % Load labels into OPS using opsCopyLayers
      copy_param = [];
      copy_param.layer_source.existence_check = false;
      copy_param.layer_dest.existence_check = false;
      
      % Set the source
      copy_param.layer_source.source = 'custom';      
      copy_param.layer_dest.source = param.layer_tracker.track.layer_dest_source;
      
      if strcmp(copy_param.layer_dest.source, 'layerdata')
        copy_param.layer_dest.layerdata_source = param.layer_tracker.track.layer_dest_layerdata_source;
        copy_param.layer_dest.echogram_source = param.layer_tracker.track.layer_dest_echogram_source;
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
      copy_param.layer_dest.name = param.layer_tracker.track.lsm.lyrtop;
      fprintf('\nopsCopyLayers %s (%s) [TOP]', param.day_seg, datestr(now)); 
      opsCopyLayers(param,copy_param);
      
      %% Write bottom layer
      % Interpolate from row number to TWTT
      data.TWTT = interp1(1:length(data.Time), data.Time, Labels.bot.y);
      copy_param.layer_source.gps_time = data.GPS_time;
      copy_param.layer_source.twtt = data.TWTT;
      
      % Load labels into OPS using opsCopyLayers
      copy_param.layer_dest.name = param.layer_tracker.track.lsm.lyrbot;
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