function frames_update(param,param_override)
% frames_update(param,param_override)
%
% This function updates the frames according to the param struct. The
% primary purpose is to update an old frames file to the new file format.
%
% param: parameter spreadsheet structure array
%
% Examples: See run_all_frames_update.m
%
% Author: John Paden
%
% See also: run_all_frames_update, run_frames_update, frames_update

param = merge_structs(param,param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Input checks
frames_fn = ct_filename_support(param,'','frames');

% Check file format
var_list = whos('-file',frames_fn);
if any(strcmpi('frames',{var_list.name}))
  % Old format had a variable "frames" which the new format does not have
  load(frames_fn,'frames');
  
  % Load records file
  records = records_load(param,'gps_time');
  
  frames.gps_time = [records.gps_time(frames.frame_idxs), records.gps_time(end)];
  Nfrms = length(frames.frame_idxs);
  frames.notes = cell(1,Nfrms);
  if ~isfield(frames,'quality')
    frames.quality = zeros(1,Nfrms);
  end
  if ~isfield(frames,'proc_mode')
    frames.proc_mode = zeros(1,Nfrms);
  end
  frames.Nx = length(records.gps_time);
  frames.param.day_seg = param.day_seg;
  frames.param.season_name = param.season_name;
  frames.param.radar_name = param.radar_name;
  frames.param.sw_version = param.sw_version;
  
  if isfield(param,'ct_file_lock') && param.ct_file_lock
    frames.file_version = '1L';
  else
    frames.file_version = '1';
  end
  frames.file_type = 'frames';
  
  fprintf('Saving updated frames file: %s\n', frames_fn);
  ct_save(frames_fn,'-struct','frames');
else
  % New format
  frames = load(frames_fn);
  
  % Load records file
  records = records_load(param,'gps_time');

  frames.gps_time = [records.gps_time(frames.frame_idxs), records.gps_time(end)];
  Nfrms = length(frames.frame_idxs);
  if ~isfield(frames,'notes')
    frames.notes = cell(1,Nfrms);
  end
  if ~isfield(frames,'quality')
    frames.quality = zeros(1,Nfrms);
  end
  if ~isfield(frames,'proc_mode')
    frames.proc_mode = zeros(1,Nfrms);
  end
  frames.Nx = length(records.gps_time);
  frames.param.day_seg = param.day_seg;
  frames.param.season_name = param.season_name;
  frames.param.radar_name = param.radar_name;
  frames.param.sw_version = param.sw_version;
  
  if param.ct_file_lock
    frames.file_version = '1L';
  else
    frames.file_version = '1';
  end
  frames.file_type = 'frames';
  
  fprintf('Saving updated frames file: %s\n', frames_fn);
  ct_save(frames_fn,'-struct','frames');
end
