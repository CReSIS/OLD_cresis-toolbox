function make_layer_files(param,param_override)
% make_layer_files(param,param_override)
%
% Makes layer files for the picker.m program. This should be run
% from run_make_layer_files.m (that script sets up all the control
% variables).
%
% Author: John Paden
%
% See also: run_make_layer_files

param = merge_structs(param,param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Determine which frames we will operate on
% Load frames file
load(ct_filename_support(param,'','frames'));

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
% valid_frms = ones(1,length(param.cmd.frms));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end
physical_constants;

if ~isfield(param.make_layer_files,'post_dir')
  param.make_layer_files.post_dir = '';
end

if ~isfield(param.make_layer_files,'echogram_img')
  param.make_layer_files.echogram_img = 0;
end

if ~isfield(param.make_layer_files,'frm_types') || isempty(param.make_layer_files.frm_types)
  param.make_layer_files.frm_types = {-1,0,-1,-1,-1};
end

if ~isfield(param.make_layer_files,'skip_phrase') || isempty(param.make_layer_files.skip_phrase)
  param.make_layer_files.skip_phrase = ''; % Often set to 'do not process'
end

if ~isfield(param.make_layer_files,'save_changes') || isempty(param.make_layer_files.save_changes)
  param.make_layer_files.save_changes = true; % Useful for debugging
end

% Flag which prevents overwriting layer files which already exist
if ~isfield(param.make_layer_files,'do_not_overwrite_layer_files') || isempty(param.make_layer_files.do_not_overwrite_layer_files)
  param.make_layer_files.do_not_overwrite_layer_files = false;
end

% Flag for updating GPS values (useful when you do not want to overwrite
% the files, but you do want the GPS info to be updated). Note that
% existing layers will be re-interpolated to the new GPS time so the
% layer data does change some.
if ~isfield(param.make_layer_files,'update_gps') || isempty(param.make_layer_files.update_gps)
  param.make_layer_files.update_gps = true;
end

% Flag for modifying the GPS time so that the min/max time match between the
% old and new layer files. This is useful when timing offsets are
% present, but the data themselves don't have an offset.
if ~isfield(param.make_layer_files,'adjust_gps_time') || isempty(param.make_layer_files.adjust_gps_time)
  param.make_layer_files.adjust_gps_time = false;
end

in_fn_dir = ct_filename_out(param,param.make_layer_files.echogram_input,'');
out_fn_dir = ct_filename_out(param,param.make_layer_files.layer_output,'');

% Create each layer file
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  if ~ct_proc_frame(frames.proc_mode(frm),param.make_layer_files.frm_types)
    continue;
  end
  
  if param.make_layer_files.echogram_img == 0
    in_fn = fullfile(in_fn_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
  else
    in_fn = fullfile(in_fn_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
      param.make_layer_files.echogram_img, param.day_seg, frm));
  end
  
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
  if exist(out_fn,'file') && param.make_layer_files.do_not_overwrite_layer_files
    if param.make_layer_files.update_gps
      fprintf('    Updating GPS, re-interpolating layers to new GPS time\n');
      % ---------------------------------------------------------------------
      % Interpolate old layer picks onto the new GPS_time axes
      % and adjust for changes in elevation
      old_lyr = load(out_fn,'GPS_time','layerData','Elevation');
      
      if param.make_layer_files.adjust_gps_time
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
        layerData{layer_idx}.value{2}.data = interp1(old_lyr.GPS_time, ...
          old_lyr.layerData{layer_idx}.value{2}.data, lyr.GPS_time, 'linear', 'extrap')  + fast_time_correction;
        % Nearest neighbor interpolation for quality data
        lyr.layerData{layer_idx}.quality = interp1(old_lyr.GPS_time, ...
          old_lyr.layerData{layer_idx}.quality, lyr.GPS_time, 'nearest', 'extrap');
      end
      
      if param.make_layer_files.save_changes
        fprintf('    Update: %s\n', out_fn);
        save(out_fn,'-v7.3','-struct','lyr','layerData','Latitude','Longitude','Elevation','GPS_time');
      else
        fprintf('  Not saving information (TEST MODE)\n');
      end
    else
      fprintf('    Already exists (not overwriting)\n');
    end
  else
    if param.make_layer_files.save_changes
      fprintf('    New: %s\n', out_fn);
      save(out_fn,'-v7.3','-struct','lyr','layerData','Latitude','Longitude','Elevation','GPS_time');
    else
      fprintf('  Not saving information (TEST MODE)\n');
    end
  end
end
