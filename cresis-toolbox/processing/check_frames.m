function check_frames(param)
% check_frames(param)
%
% Checks the fields in the frames files for potential errors.
%
% param = parameter structure indicating which radar and season to load.
%   This causes the standard parameter spreadsheet filename to be created
%   and loaded.  Therefore it will not work if you don't have the parameter
%   spreadsheet in the standard location with the standard name.
%   ALTERNATE MODE: pass in a string containing the filename of a specific
%   record to load
%
% Example:
% param = [];
% param.radar_name = 'mcrds';
% param.season_name = '2009_Antarctica_TO';
% check_frames(param)
%
% param = [];
% param.radar_name = 'icards';
% param.season_name = '1993_Greenland_P3';
% check_frames(param)
%
% param = [];
% param.radar_name = 'snow3';
% param.season_name = '2013_Antarctica_Basler';
% check_frames(param)

if ischar(param)
  % param is a string containing the filename
  
  records_fn = ct_filename_support(param,'','records');
  records = load(param,'param_records');
  
  frames_fn = ct_filename_support(records.param_records,'','frames');
  
  if ~isfield(records,'param_records')
    error('This mode only supported for new records files with param_records field.');
  end
  check_records_support_func(param,records.param_records);
  
elseif isstruct(param)
  % param is a struct indicating which radar/season to check (checks all segments)
  [output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

  param_fn = ct_filename_param(sprintf('%s_param_%s.xls', output_dir, param.season_name));
  params = read_param_xls(param_fn);
  
  for param_idx = 1:length(params)
    param = params(param_idx);
    
    records_fn = ct_filename_support(param,'','records');
    
    if ~isempty(regexpi(param.cmd.notes,'do not process'))
      continue;
    end
    
    % DEBUG OPTION TO JUST CHECK RECORDS BASED ON GENERIC COLUMN IN SPREADSHEET
    % Uses the generic column of the parameter spreadsheet to determine which segments
    % to check.
    %if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    %  continue;
    %end

    fprintf('%d of %d: Checking frames %s\n', param_idx, length(params), records_fn);
    
    check_frames_support_func(records_fn,param);
  end
  
else
  error('Invalid argument')
end

end

function check_frames_support_func(records_fn,param)
% check_frames_support_func(records_fn,param)
%
% Support function which does the actual checking of the frames

if ~exist(records_fn,'file')
  warning('Missing records file %s\n', records_fn);
  return;
end

frames_fn = ct_filename_support(param,'','frames');
records = load(records_fn,'gps_time');
load(frames_fn);

if frames.frame_idxs(1) ~= 1
  warning('Frames starts with %d instead of 1', frames.frame_idxs(1));
end

if any(frames.frame_idxs > length(records.gps_time))
  warning('frame_idxs has entries (e.g. %d) extending beyond records %d', ...
    max(frames.frame_idxs), length(records.gps_time));
end

fixable_error = false;

if length(frames.frame_idxs) ~= length(frames.nyquist_zone)
  warning('nyquist_zone field length %d does not match frame_idxs length %d', ...
    length(frames.nyquist_zone), length(frames.frame_idxs));
  fixable_error = true;
  if length(frames.frame_idxs) > length(frames.nyquist_zone)
    frames.nyquist_zone(end+1 : length(frames.frame_idxs)) = NaN;
  else
    frames.nyquist_zone = frames.nyquist_zone(1:length(frames.frame_idxs));
  end
end

if length(frames.frame_idxs) ~= length(frames.proc_mode)
  warning('proc_mode field length %d does not match frame_idxs length %d', ...
    length(frames.proc_mode), length(frames.frame_idxs));
  fixable_error = true;
  if length(frames.frame_idxs) > length(frames.proc_mode)
    frames.proc_mode(end+1 : length(frames.frame_idxs)) = NaN;
  else
    frames.proc_mode = frames.proc_mode(1:length(frames.frame_idxs));
  end
end

if 1 && fixable_error
  fprintf('Save %s\n', frames_fn);
  save(frames_fn,'frames');
end

end
