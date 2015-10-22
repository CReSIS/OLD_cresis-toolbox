% script ops_update_gps_time
%
% This script is used to update the GPS values (lat,lon,elev,gps_time,roll,pitch,heading,etc)
% in the database.
%
% This script is for the following situations:
% 1. GPS time correction was applied to the records/data files
% 2. New GPS data has been applied to the records/data files (e.g. motion
%    compensated or better GPS results)

params = read_param_xls('/users/paden/scripts/branch/params-cr1/rds_param_2013_Greenland_P3.xls','','post');
% params = read_param_xls('/users/paden/scripts/branch/params-cr1/rds_param_2013_Greenland_P3.xls','20130426_01','post');
% params = read_param_xls('C:\users/dangermo\Documents/scripts/branch/params-cr1/rds_param_2013_Greenland_P3.xls','20130426_01','post');

% Set these to the values you would use in run_ops_bulk_insert.m
settings.path_decimation_spacing = 100;
location_name = 'arctic';

%% Automated Section
% =========================================================================

% confirm_button = questdlg('Are the get_heights lever arm parameters and img parameters set correctly?','Confirm Get Heights Parameters','YES: Update','NO: Quit','YES: Update');
% if ~strcmp(confirm_button,'YES: Update')
%   return;
% end

delta_time_offset_guard = 0.2;

for param_idx = 1:length(params)
  param = params(param_idx);
  delta_time_offset = []; % Used to keep track of GPS time corrections that have been applied to the records file and not the database
  
  if ~param.cmd.generic
    continue;
  end
  
  fprintf('Updating GPS in OPS database for %s\n', param.day_seg);
  
  
  %% Load records and frames files
  records_fn = ct_filename_support(param,'','records');
  fprintf('Loading %s\n', records_fn);
  records = load(records_fn);
  records.along_track = geodetic_to_along_track(records.lat,records.lon,records.elev);
  
  frames_fn = ct_filename_support(param,'','frames');
  fprintf('Loading %s\n', frames_fn);
  load(frames_fn);
  
  %% Apply lever arm to GPS data
  fprintf('Applying lever arm to records GPS data (%s)\n', datestr(now,'HH:MM:SS'));
  if isempty(param.get_heights.lever_arm_fh)
    error('Mode without lever arm not supported');
  end
  
  % Default values to use
  img_idx = 1;
  wf = abs(param.get_heights.imgs{img_idx}(1,1));
  adc = abs(param.get_heights.imgs{img_idx}(1,2));
  lever_arm_fh = param.get_heights.lever_arm_fh;
  trajectory_param = struct('rx_path', param.radar.wfs(wf).rx_paths(adc), ...
    'tx_weights', param.radar.wfs(wf).tx_weights, 'lever_arm_fh', lever_arm_fh);
  for tmp_wf_adc_idx = 2:size(param.get_heights.imgs{1},1)
    tmp_wf = abs(param.get_heights.imgs{img_idx}(tmp_wf_adc_idx,1));
    tmp_adc = abs(param.get_heights.imgs{img_idx}(tmp_wf_adc_idx,2));
    trajectory_param.rx_path(tmp_wf_adc_idx) = param.radar.wfs(tmp_wf).rx_paths(tmp_adc);
  end
    
  trajectory_param.gps_source = records.gps_source;
  trajectory_param.radar_name = param.radar_name;
  trajectory_param.season_name = param.season_name;
  records = trajectory_with_leverarm(records,trajectory_param);
  
  %% Get season_id, segment_id, and location_id from OPS database
  fprintf('Get database info (%s)\n', datestr(now,'HH:MM:SS'));
  system_name = ct_output_dir(param.radar_name)
  [status,data] = ops_query(sprintf('SELECT season_id FROM %s_seasons where season_name = ''%s''', system_name, param.season_name));
  season_id = data{1};
  [status,data] = ops_query(sprintf('SELECT segment_id FROM %s_segments where segment_name = ''%s''', system_name, param.day_seg));
  segment_id = data{1};
  [status,data] = ops_query(sprintf('SELECT location_id FROM %s_locations where location_name = ''%s''', system_name, location_name));
  location_id = data{1};
  
  %% Get point path from OPS database
  [status,data] = ops_query(sprintf('SELECT gps_time,ST_X(point_path),ST_Y(point_path),ST_Z(point_path),roll,pitch,heading,point_path_id FROM %s_point_paths where segment_id = %d', system_name, segment_id));
  point_path.gps_time = cell2mat(data(1,:));
  point_path.lat = cell2mat(data(3,:));
  point_path.lon = cell2mat(data(2,:));
  point_path.elev = cell2mat(data(4,:));
  point_path.roll = cell2mat(data(5,:));
  point_path.pitch = cell2mat(data(6,:));
  point_path.heading = cell2mat(data(7,:));
  point_path.point_path_id = cell2mat(data(7,:));
  
  [point_path.gps_time sort_idxs] = sort(point_path.gps_time);
  point_path.lat = point_path.lat(sort_idxs);
  point_path.lon = point_path.lon(sort_idxs);
  point_path.elev = point_path.elev(sort_idxs);
  point_path.roll = point_path.roll(sort_idxs);
  point_path.pitch = point_path.pitch(sort_idxs);
  point_path.heading = point_path.heading(sort_idxs);
  point_path.point_path_id = point_path.point_path_id(sort_idxs);
  
  figure(1); clf;
  plot(point_path.lon,point_path.lat)
  
  delta_time_offset(end+1) = records.gps_time(1) - point_path.gps_time(1);
%   delta_time_offset(end+1) = records.gps_time(end) - point_path.gps_time(end);
  
  %% Get start/stop gps time for each fram and verify that the frames have not changed
  [status,data] = ops_query(sprintf('SELECT frame_id,frame_name,start_gps_time,stop_gps_time FROM %s_frames WHERE segment_id = %f;\n', system_name, segment_id));
  frame_start_gps_time = [];
  frame_stop_gps_time = [];
  frame_id = [];
  for frame_idx = 1:length(frames.frame_idxs)
    % Handle the end segment (do this first to handle single frame segments)
    if frame_idx == length(frames.frame_idxs)
      frame_start_gps_time(1,frame_idx) = records.gps_time(frames.frame_idxs(frame_idx));
      frame_stop_gps_time(1,frame_idx) = records.gps_time(end);
    else
      frame_start_gps_time(1,frame_idx) = records.gps_time(frames.frame_idxs(frame_idx));
      frame_stop_gps_time(1,frame_idx) = records.gps_time(frames.frame_idxs(frame_idx+1));
    end
    frm_name = sprintf('%s_%03d', param.day_seg, frame_idx);
    match_idx = strcmp(frm_name,data(2,:));
    if isempty(match_idx)
      warning('Missing frame %s: code does not know how to handle this case', frm_name);
      keyboard
      return;
    end
    frame_id(1,frame_idx) = data{1,match_idx};
    
    delta_time_offset(end+1) = frame_start_gps_time(1,frame_idx) - data{3,match_idx};
    delta_time_offset(end+1) = frame_stop_gps_time(1,frame_idx) - data{4,match_idx};
  end
  
  %% Estimate the delta_time offset (i.e. checks to see if a GPS time correction was applied and reports this)
  delta_time_offset_est = median(delta_time_offset);
  fprintf('Estimating a delta time offset of %f sec\n', delta_time_offset_est);
  if max(delta_time_offset) - min(delta_time_offset) > delta_time_offset_guard
    warning('delta time offset between new records file and database is not consistent... could mean that something bad has happened')
    keyboard
  end
  if delta_time_offset_est ~= 0
    warning('delta time offset is not zero, is that okay?');
    keyboard
  end
  
  if delta_time_offset_est ~= 0
    %% Determine which points are going to be deleted from that paths (actually all will be deleted and then new points inserted)
    point_paths_search_str = sprintf('FROM %s_point_paths WHERE segment_id = %f AND (gps_time < %f OR gps_time > %f)', system_name, segment_id, records.gps_time(1), records.gps_time(end));
    [status,data] = ops_query(sprintf('SELECT point_path_id,segment_id,season_id,gps_time,heading,roll,pitch,ST_X(point_path),ST_Y(point_path),location_id %s;\n', point_paths_search_str));
    if length(data) > 1 % Make sure we got at least one valid return
      fprintf('Deleting %d points from %s_point_paths\n', size(data,2), system_name);
      hold on;
      plot(cell2mat(data(8,:)),cell2mat(data(9,:)),'rx');
      hold off;
    end
  end
  
  %% Determine what the new point path indices are going to be
  if delta_time_offset_est > 0
    new_idxs = find(records.gps_time >= point_path.gps_time(end));
    decimation_spacing = settings.path_decimation_spacing;
    spacing_idxs = get_equal_alongtrack_spacing_idxs(records.along_track(new_idxs),decimation_spacing);
    if length(spacing_idxs) > 2
      spacing_idxs = spacing_idxs(2:end);
    end
    new_idxs = new_idxs(spacing_idxs);
    
    hold on;
    for idx = 1:length(new_idxs)
      plot(records.lon(new_idxs), records.lat(new_idxs),'kx');
    end
    hold off;
  elseif delta_time_offset_est < 0
    new_idxs = find(records.gps_time < point_path.gps_time(1));
    decimation_spacing = settings.path_decimation_spacing;
    spacing_idxs = get_equal_alongtrack_spacing_idxs(records.along_track(new_idxs),decimation_spacing);
    new_idxs = new_idxs(spacing_idxs);
    
    hold on;
    for idx = 1:length(new_idxs)
      plot(records.lon(new_idxs), records.lat(new_idxs),'kx');
    end
    hold off;
  end
  
  %% Create filename of file to store SQL commands in
  sql_fn = ct_filename_tmp(param, '', 'update_gps_time','sql.txt');
  fprintf('Creating output file %s (%s)\n', sql_fn, datestr(now,'HH:MM:SS'));
  sql_fn_dir = fileparts(sql_fn);
  if ~exist(sql_fn_dir,'dir')
    mkdir(sql_fn_dir);
  end
  [sql_fid msg] = fopen(sql_fn,'w');
  if sql_fid == -1
    error('Failed to open %s: %s', sql_fn, msg);
  end
  
  fprintf(sql_fid,'BEGIN;\n');
  
  %% Update the point_paths table
  %   Delete the old point_paths
  %   Insert the new point paths
  fprintf(sql_fid, 'DELETE FROM %s_point_paths WHERE segment_id = %d;\n', system_name, segment_id);

  decimation_spacing = settings.path_decimation_spacing;
  spacing_idxs = get_equal_alongtrack_spacing_idxs(records.along_track,decimation_spacing);
  for idx = 1:length(spacing_idxs)
    if ~mod(idx-1,1000)
      fprintf('Updating point path %d of %d (%s)\n', idx, length(spacing_idxs), datestr(now,'HH:MM:SS'));
    end
    fprintf(sql_fid, 'INSERT INTO %s_point_paths (segment_id,season_id,gps_time,heading,roll,pitch,point_path,location_id) VALUES (%d,%d,%f,%f,%f,%f,ST_GeomFromText(''POINTZ(%.15f %.15f %.6f)'',4326),%d);\n', ...
      system_name, segment_id, season_id,records.gps_time(spacing_idxs(idx)), ...
      records.heading(spacing_idxs(idx)),records.roll(spacing_idxs(idx)), ...
      records.pitch(spacing_idxs(idx)),records.lon(spacing_idxs(idx)), ...
      records.lat(spacing_idxs(idx)), records.gps_time(spacing_idxs(idx)), location_id);
  end
  
  %% Update the segments table: update start_gps_time, stop_gps_time, line_path
  if delta_time_offset_est ~= 0
    fprintf(sql_fid,'UPDATE %s_segments SET start_gps_time=%f,stop_gps_time=%f WHERE segment_id = %d;\n', ...
      system_name, records.gps_time(1), ...
      records.gps_time(end), segment_id);
  end
  
  fprintf(sql_fid,'UPDATE %s_segments SET line_path = ST_GeomFromText(''LINESTRINGZ(\n', system_name);
  for idx = 1:length(spacing_idxs)
    if ~mod(idx-1,1000)
      fprintf('Updating line path %d of %d (%s)\n', idx, length(spacing_idxs), datestr(now,'HH:MM:SS'));
    end
    if idx == length(spacing_idxs)
      fprintf(sql_fid, '  %f %f %f\n', records.lon(spacing_idxs(idx)), ...
        records.lat(spacing_idxs(idx)), records.gps_time(spacing_idxs(idx)));
    else
      fprintf(sql_fid, '  %f %f %f,\n', records.lon(spacing_idxs(idx)), ...
        records.lat(spacing_idxs(idx)), records.gps_time(spacing_idxs(idx)));
    end
  end
  fprintf(sql_fid,')'',4326) WHERE segment_id = %d;\n', segment_id);
  
  %% Update the frames table: update start_gps_time, stop_gps_time
  if delta_time_offset_est ~= 0
    for frame_idx = 1:length(frames.frame_idxs)
      fprintf(sql_fid,'UPDATE %s_frames SET start_gps_time=%f,stop_gps_time=%f WHERE frame_id = %d;\n', ...
        system_name, frame_start_gps_time(frame_idx), ...
        frame_stop_gps_time(frame_idx), frame_id(frame_idx));
    end
  end
  
  %% Update crossovers table: gps_time_1 or gps_time_2 based on segment_1_id and segment_2_id matching
  if delta_time_offset_est ~= 0
    % Technically we could lose or gain a cross over... could be a problem?
    fprintf(sql_fid,'UPDATE %s_crossovers SET gps_time_1=gps_time_1+%f WHERE segment_1_id= %d;\n', ...
      system_name, delta_time_offset_est, segment_id);
    fprintf(sql_fid,'UPDATE %s_crossovers SET gps_time_2=gps_time_2+%f WHERE segment_2_id= %d;\n', ...
      system_name, delta_time_offset_est, segment_id);
  end
  
  %% Update layer_points table: update gps_time and geometry
  if delta_time_offset_est ~= 0
    fprintf(sql_fid,'UPDATE %s_layer_points SET gps_time=gps_time+%f WHERE segment_id= %d;\n', ...
      system_name, delta_time_offset_est, segment_id);
  end
  
  %% Update each layer point geometry (one point at a time...)
  [status,data] = ops_query(sprintf('SELECT layer_points_id,gps_time FROM %s_layer_points where segment_id = %d', system_name, segment_id));
  lon = interp1(records.gps_time, records.lon, cell2mat(data(2,:)), 'linear','extrap');
  lat = interp1(records.gps_time, records.lat, cell2mat(data(2,:)), 'linear','extrap');
  elev = interp1(records.gps_time, records.elev, cell2mat(data(2,:)), 'linear','extrap');
  for data_idx = 1:size(data,2)
    if ~mod(data_idx-1,1000)
      fprintf('Updating layer point %d of %d (%s)\n', data_idx, size(data,2), datestr(now,'HH:MM:SS'));
    end
    fprintf(sql_fid, 'UPDATE %s_layer_points SET layer_point=ST_GeomFromText(''POINTZ(%.15f %.15f %.6f)'',4326) WHERE layer_points_id=%d;\n', ...
      system_name, lon(data_idx), lat(data_idx), elev(data_idx), data{1,data_idx});
  end
  
  %% Cleanup
  fprintf(sql_fid,'COMMIT;\n');
  fclose(sql_fid);
  
  fprintf('Please review sql file %s before running (%s)\n', sql_fn, datestr(now,'HH:MM:SS'));
end

return;








