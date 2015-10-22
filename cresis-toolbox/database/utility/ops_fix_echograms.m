% OPS FIX ECHOGRAMS
%
% GOES THROUGH GIVEN SEGMENTS (PARAM SHEET) AND UPDATES THE ECHOGRAM URLS
%
% Author: Kyle Purdon

% =========================================
% USER INPUT SECTION
% =========================================

param_fn = 'C:\Users\kpurdon\Documents\scripts\params-cr1\rds_param_2013_Greenland_P3.xls';
location = 'arctic';
sys_name = 'rds';

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

missing_segments = {};

for param_idx = 1:length(params)
  param = params(param_idx);
  if param.cmd.generic == 1
    if ~isempty(regexpi(param.cmd.notes,'do not process'))
      warning('You have enabled a segment with ''do not process'' in the cmd.notes, dbcont to continue');
      keyboard
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
    
    % GET ALL frame_id,frame_name FOR CURRENT SEGMENT
    query = sprintf('SELECT frame_id,frame_name FROM %s_frames WHERE segment_id = %d;',sys_name,segment_id);
    [status,data] = ops_query(query);
    
    % UPDATE THE ECHOGRAM URL FOR EACH FRAME RETURNED
    for frame_idx = 1:size(data,2)
      
      % GET CURRENT FRAME INFORMATION
      [frame_id,frame_name] = data{:,frame_idx};
      
      fprintf('\tUpdating echogram for frame %s ...\n',frame_name);
      
      % CONSTRUCT ECHOGRAM URL (ex. 20130402_01_001_1echo.jpg)
      echogram_url = strcat('ftp://data.cresis.ku.edu/data/',sys_name,'/',season_name,'/images/',segment,'/',frame_name,'_1echo.jpg');
      fprintf('\t\t%s\n',echogram_url);
      
      % UPDATE echograms (echogram_url) WITH NEW VALUES
      query = sprintf('UPDATE %s_echograms SET echogram_url = ''%s'' WHERE frame_id = %d RETURNING echogram_id;',sys_name,echogram_url,frame_id);
      [~,~] = ops_query(query);
      
    end
  end
end

fprintf('Missing Segments\n')
if isempty(missing_segments)
  fprintf('\tNone.\n');
else
  for idx = 1:length(missing_segments)
    fprintf('\t%s\n',missing_segments{idx});
  end
end