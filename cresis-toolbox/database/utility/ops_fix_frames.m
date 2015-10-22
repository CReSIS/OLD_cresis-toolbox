% OPS FIX FRAMES
%
% GOES THROUGH GIVEN SEGMENTS (PARAM SHEET) AND INSERTS ALL FRAME
% INFORMATION IF THE FRAME IS MISSING FROM THE DATABASE.
%
% Author: Kyle Purdon

% =========================================
% USER INPUT SECTION
% =========================================

param_fn = 'C:\Users\kpurdon\Documents\scripts\params-cr1\snow_param_2012_Antarctica_DC8.xls';
location = 'antarctic';
sys_name = 'snow';

% =========================================
% AUTOMATED SECTION
% =========================================

% MAKE USER CONFIRM TASK
confirm_button = questdlg('DOUBLE CHECK INPUT. (ESPECIALLY LOCATION AND SYSTEM)','Confirm Settings','YES:FIXIT','NO:QUIT','NO:QUIT');
switch confirm_button
  case 'NO:QUIT'
    error('PROCESS STOPPED BY USER.');
end

missing_segments = {};
existing_frames = {};

% QUERY FOR location_id
query = sprintf('SELECT location_id FROM %s_locations WHERE location_name = ''%s'';',sys_name,location);
[status,data] = ops_query(query);
if status == 1
  location_id = data{1};
else
  error(data);
end

% LOAD PARAM SHEET AND CHECK INPUTS
params = read_param_xls(param_fn);
radar_name = params(1).radar_name;
season_name = params(1).season_name;
fprintf('Checking for file inputs for each segment...\n');

for param_idx = 1:length(params)
  param = params(param_idx);
  if param.cmd.generic == 1
    if ~isempty(regexpi(param.cmd.notes,'do not process'))
      warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
      keyboard
    end
    records_fn = ct_filename_support(param,'','records');
    frames_fn = ct_filename_support(param,'','frames');
    layer_dir = ct_filename_out(param,settings.layer_post_directory,param.day_seg);
    if ~exist(frames_fn,'file')
      error('  %s: missing %s\n', param.day_seg, frames_fn);
    elseif ~exist(records_fn,'file')
      error('  %s: missing %s\n', param.day_seg, records_fn);
    else
      fprintf('  %s checked\n', param.day_seg);
    end
  end
  
  param = params(param_idx);
  if param.cmd.generic ~= 1
    continue;
  end
  
  segment = param.day_seg;
  
  fprintf('Fixing segment %s ...\n',segment);
  
  
  % QUERY FOR SEGMENT ID
  query = sprintf('SELECT segment_id FROM %s_segments WHERE segment_name = ''%s'' and season_id = (SELECT season_id FROM %s_seasons WHERE season_name = ''%s'');',sys_name,segment,sys_name,season_name);
  [status,data] = ops_query(query);
  if status == 1
    segment_id = data{1};
  elseif status == 2
    fprintf('\t Segment %s was not found in the database. Skipping.\n',segment);
    missing_segments{end+1} = segment;
    continue;
  else
    error(data);
  end
  
  % LOAD THE FRAMES AND RECORDS FILE
  load(frames_fn);
  records = load(records_fn);
  if isfield(records,'records')
    % Adding this check to support old record (pre-cr1) file format
    records = records.records;
  end
  
  frame_start_gps_time = zeros(1,length(frames.frame_idxs));
  frame_stop_gps_time = zeros(1,length(frames.frame_idxs));
  
  % Get start/stop gps time for each frame
  for frame_idx = 1:length(frames.frame_idxs)
    % Handle the end segment (do this first to handle single frame segments)
    if frame_idx == length(frames.frame_idxs)
      frame_start_gps_time(1,frame_idx) = records.gps_time(frames.frame_idxs(frame_idx));
      frame_stop_gps_time(1,frame_idx) = records.gps_time(end);
    else
      frame_start_gps_time(1,frame_idx) = records.gps_time(frames.frame_idxs(frame_idx));
      frame_stop_gps_time(1,frame_idx) = records.gps_time(frames.frame_idxs(frame_idx+1));
    end
  end
  
  %  (frame_name,start_gps_time,stop_gps_time + segment_id,location_id)
  
  for frame_idx = 1:length(frames.frame_idxs)
    
    frame_name = strcat(segment,sprintf('_%03d',frame_idx));
    fprintf('\tFixing frame %s ...\n',frame_name);

    % CHECK IF THE FRAME IS IN THE DATABASE
    query = sprintf('SELECT count(frame_id) FROM %s_frames WHERE frame_name = ''%s'';',sys_name,frame_name);
    [status,data] = ops_query(query);
    if status == 1
      if data{1} ~= 0
        fprintf('\t Frame %s already exists in the database. Skipping.\n',frame_name);
        existing_frames{end+1} = frame_name;
        continue;
      end
    else
      error(data{1});
    end
    
    % FOR EACH FRAME INSERT INTO THE TABLE RETURNING FRAME ID
    query = sprintf('INSERT INTO %s_frames (segment_id,frame_name,start_gps_time,stop_gps_time,location_id) VALUES (%d,''%s'',%2.7f,%2.7f,%d) RETURNING frame_id;',...
      sys_name,segment_id,frame_name,frame_start_gps_time(frame_idx),frame_stop_gps_time(frame_idx),location_id);
    [status,data] = ops_query(query);
    if status == 1
      frame_id = data{1};
    else
      error(data);
    end
    
    % FOR THE INSERTED FRAME INSERT THE ECHOGRAM INFORMATION
    query = sprintf('INSERT INTO %s_echograms (frame_id,echogram_url,location_id) VALUES (%d,''%s'',%d) RETURNING echogram_id;',sys_name,frame_id,'NO ECHOGRAM ONLINE',location_id);
    [status,data] = ops_query(query);
    if status == 1
      echogram_id = data{1};
    else
      error(data);
    end
  end
end

fprintf('Missing Segments\n')
for idx = 1:length(missing_segments)
  fprintf('%s\n',missing_segments{idx});
end
fprintf('Existing Frames\n')
for idx = 1:length(existing_frames)
  fprintf('%s\n',existing_frames{idx});
end
