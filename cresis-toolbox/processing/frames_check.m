function frames_check(param)
% frames_check(param)
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
% frames_check(param)
%
% param = [];
% param.radar_name = 'icards';
% param.season_name = '1993_Greenland_P3';
% frames_check(param)
%
% param = [];
% param.radar_name = 'snow3';
% param.season_name = '2013_Antarctica_Basler';
% frames_check(param)

if ischar(param)
  % param is a string containing the filename
  
  records_fn = ct_filename_support(param,'','records');
  records = load(param,'param_records');
  
  frames_fn = ct_filename_support(records.param_records,'','frames');
  
  if ~isfield(records,'param_records')
    error('This mode only supported for new records files with param_records field.');
  end
  records_check_support_func(param,records.param_records);
  
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

    fprintf('%s\tChecking\t%s\n', param.day_seg, records_fn);
    
    frames_create_support_func(records_fn,param);
  end
  
else
  error('Invalid argument')
end

end

function frames_create_support_func(records_fn,param)
% frames_create_support_func(records_fn,param)
%
% Support function which does the actual checking of the frames

if ~exist(records_fn,'file')
  fprintf(2,'%s\tno_records_file\t%s\n', param.day_seg, records_fn);
  return;
end

frames_fn = ct_filename_support(param,'','frames');
if ~exist(frames_fn,'file')
  fprintf(2,'%s\tno_frames_file\t%s\n', param.day_seg, records_fn);
  return;
end

records = load(records_fn,'lat','lon');
along_track = geodetic_to_along_track(records.lat,records.lon);
load(frames_fn);

if frames.frame_idxs(1) ~= 1
  fprintf(2,'%s\tframe_idxs_1_not_1\t%d instead of 1\n', param.day_seg, frames.frame_idxs(1));
end

if any(frames.frame_idxs > length(records.lat))
  fprintf(2,'%s\tframe_idxs_too_large\tindex %d > number of records %d\n', ...
    param.day_seg, max(frames.frame_idxs), length(records.lat));
end

fixable_error = false;

if length(frames.frame_idxs) ~= length(frames.nyquist_zone)
  fprintf(2,'%s\tnyquist_zone_length_mismatch\tlength %d does not match frame_idxs length %d\n', ...
    param.day_seg, length(frames.nyquist_zone), length(frames.frame_idxs));
  fixable_error = true;
  if length(frames.frame_idxs) > length(frames.nyquist_zone)
    frames.nyquist_zone(end+1 : length(frames.frame_idxs)) = NaN;
  else
    frames.nyquist_zone = frames.nyquist_zone(1:length(frames.frame_idxs));
  end
end

if length(frames.frame_idxs) ~= length(frames.proc_mode)
  fprintf(2,'%s\tproc_mode_length_mismatch\tlength %d does not match frame_idxs length %d\n', ...
    param.day_seg, length(frames.proc_mode), length(frames.frame_idxs));
  fixable_error = true;
  if length(frames.frame_idxs) > length(frames.proc_mode)
    frames.proc_mode(end+1 : length(frames.frame_idxs)) = NaN;
  else
    frames.proc_mode = frames.proc_mode(1:length(frames.frame_idxs));
  end
end

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
if any(strcmpi(output_dir,'rds'))
  default_frame_len = 50000;
  default_min_rec = 1000;
elseif any(strcmpi(output_dir,'accum'))
  default_frame_len = 20000;
  default_min_rec = 1000;
elseif any(strcmpi(output_dir,{'kaband','kuband','snow'}))
  default_frame_len = 5000;
  default_min_rec = 1000;
end

for frm = 1:length(frames.frame_idxs)
  if frm < length(frames.frame_idxs)
    frame_len = along_track(frames.frame_idxs(frm+1)) - along_track(frames.frame_idxs(frm));
    frame_num_recs = frames.frame_idxs(frm+1) - frames.frame_idxs(frm);
  else
    frame_len = along_track(end) - along_track(frames.frame_idxs(frm));
    frame_num_recs = length(records.lat) - frames.frame_idxs(frm);
  end
  % Check length of frame
  if length(frames.frame_idxs) == 1 && frame_len < default_frame_len/10
    fprintf(2,'%s\tsingle_frame_too_short\t%g km < %g km for frame \t%d\n', param.day_seg, frame_len/1e3, 1/5*default_frame_len/1e3, frm);
  elseif length(frames.frame_idxs) > 1 && frame_len < default_frame_len/5
    fprintf(2,'%s\ttoo_short\t%g km < %g km for frame \t%d\n', param.day_seg, frame_len/1e3, 1/5*default_frame_len/1e3, frm);
  end
  if frame_len > 2*default_frame_len
    fprintf(2,'%s\ttoo_long\t%g km > %g km for frame \t%d\n', param.day_seg, frame_len/1e3, 2*default_frame_len/1e3, frm);
  end
  % Check number of records in frame
  if frame_num_recs < default_min_rec
    fprintf(2,'%s\tFrame has very few records\t%d\t%d\n', param.day_seg, frame_num_recs, frm);
  end
end

if 1 && fixable_error
  fprintf('Save %s\n', frames_fn);
  save(frames_fn,'frames');
end

end

