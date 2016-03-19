function make_layer_files(layer_param)
% make_layer_files(layer_param)
%
% Makes layer files for the picker.m program. This should be run
% from run_make_layer_files.m (that script sets up all the control
% variables).
%
% Author: John Paden
%
% See also: run_make_layer_files

physical_constants;

if ~isfield(layer_param,'post_dir')
  layer_param.post_dir = '';
end

if ~isfield(layer_param,'frm_types') || isempty(layer_param.frm_types)
  layer_param.frm_types = {-1,0,-1,-1,-1};
end

params = read_param_xls(layer_param.param_fn);

dbstack_info = dbstack;
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isnumeric(params(param_idx).cmd.generic) && ~isnumeric(params(param_idx).cmd.generic)
    continue;
  end
  if ~param.cmd.generic
    continue;
  end
  
  fprintf('=====================================================================\n');
  fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
  fprintf('=====================================================================\n');
  
  % Load frames file
  load(ct_filename_support(param,param.records.frames_fn,'frames'));
  
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
  
  in_fn_dir = ct_filename_out(param,layer_param.echogram_input,'');
  out_fn_dir = ct_filename_out(param,layer_param.layer_output,'');

  % Create each layer file
  for frm_idx = 1:length(param.cmd.frms)
    frm = param.cmd.frms(frm_idx);
    if ~ct_proc_frame(frames.proc_mode(frm),layer_param.frm_types)
      continue;
    end
    
    in_fn = fullfile(in_fn_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    
    out_fn = fullfile(out_fn_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));

    if ~exist(out_fn_dir,'dir')
      mkdir(out_fn_dir);
    end
    
    if ~exist(in_fn,'file')
      warning('File not found %s\n', in_fn)
      continue;
    end
    
    % Load data file
    warning off;
    lyr = load(in_fn,'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom', ...
      'Elevation_Correction','Truncate_Bins','Time');
    warning on;
    lyr = uncompress_echogram(lyr);
    Nx = length(lyr.GPS_time);
    
    % ---------------------------------------------------------------------
    % Create a layer struct
    for layer_idx = 1:2
      % Manually picked points
      %  inf/nan: no pick
      %  finite: propagation time to target
      lyr.layerData{layer_idx}.value{1}.data ...
        = inf*ones(1,Nx);
      % Automatically generated points
      %  inf/nan: no pick
      %  finite: propagation time to target
      if layer_idx == 1 && isfield(lyr,'Surface')
        lyr.layerData{layer_idx}.value{2}.data = lyr.Surface;
      elseif layer_idx == 2 && isfield(lyr,'Bottom')
        lyr.layerData{layer_idx}.value{2}.data = lyr.Bottom;
      else
        lyr.layerData{layer_idx}.value{2}.data ...
          = inf*ones(1,Nx);
      end
      % Quality control level
      %  1: good
      %  2: moderate
      %  3: derived
      lyr.layerData{layer_idx}.quality ...
        = ones(1,Nx);
    end
    
    fprintf('  %s\n', in_fn);
    if exist(out_fn,'file') && layer_param.do_not_overwrite_layer_files
      if layer_param.update_gps
        fprintf('    Updating GPS, re-interpolating layers to new GPS time\n');
        % ---------------------------------------------------------------------
        % Interpolate old layer picks onto the new GPS_time axes
        % and adjust for changes in elevation
        old_lyr = load(out_fn,'GPS_time','layerData','Elevation');
        
        if layer_param.adjust_gps_time
          % Modify the GPS time so that the min/max time match between the
          % old and new layer files. This is useful when timing offsets are
          % present, but the data themselves don't have an offset
%           keyboard
%           GPS_time([1 end]) - lyr.GPS_time([1 end])
%           plot(lyr.GPS_time)
          old_lyr.GPS_time = (old_lyr.GPS_time-min(old_lyr.GPS_time)) ...
            / (max(old_lyr.GPS_time)-min(old_lyr.GPS_time)) ...
            * (max(lyr.GPS_time)-min(lyr.GPS_time)) + min(lyr.GPS_time);
%           hold on;
%           plot(lyr.GPS_time,'r:')
%           hold off;
%           GPS_time([1 end]) - lyr.GPS_time([1 end])
        end
        
        % Create fast time correction vector that compensates for
        % differences in elevation
        old_lyr.Elevation = interp1(old_lyr.GPS_time, ...
          old_lyr.Elevation, lyr.GPS_time, 'linear', 'extrap');
        fast_time_correction = (lyr.Elevation - old_lyr.Elevation)/(c/2);
        
        for layer_idx = 1:2
          % Special mapping for manual points
          lyr.layerData{layer_idx}.value{1}.data = inf*ones(1,Nx);
          for rline = 1:length(old_lyr.layerData{layer_idx}.value{1}.data)
            if old_lyr.layerData{layer_idx}.value{1}.data(rline) ~= inf
              [min_val min_idx] = min(abs(old_lyr.GPS_time(rline) - lyr.GPS_time));
              lyr.layerData{layer_idx}.value{1}.data(min_idx) = old_lyr.layerData{lyr.layer_idx}.value{1}.data(rline) + fast_time_correction(min_idx);
            end
          end
          % Linear interpolation for automated points
          layerData{lyr.layer_idx}.value{2}.data = interp1(old_lyr.GPS_time, ...
            old_lyr.layerData{layer_idx}.value{2}.data, lyr.GPS_time, 'linear', 'extrap')  + fast_time_correction;
          % Nearest neighbor interpolation for quality data
          lyr.layerData{layer_idx}.quality = interp1(old_lyr.GPS_time, ...
            old_lyr.layerData{layer_idx}.quality, lyr.GPS_time, 'nearest', 'extrap');
        end

        if layer_param.save_changes
          fprintf('    Update: %s\n', out_fn);
          save(out_fn,'-v6','-struct','lyr','layerData','Latitude','Longitude','Elevation','GPS_time');
        else
          fprintf('  Not saving information (TEST MODE)\n');
        end
      else
        fprintf('    Already exists (not overwriting)\n');
      end
    else
      if layer_param.save_changes
        fprintf('    New: %s\n', out_fn);
        save(out_fn,'-v6','-struct','lyr','layerData','Latitude','Longitude','Elevation','GPS_time');
      else
        fprintf('  Not saving information (TEST MODE)\n');
      end
    end
  end
end


return;





