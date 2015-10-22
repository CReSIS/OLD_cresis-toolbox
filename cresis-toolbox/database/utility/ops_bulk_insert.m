function ops_bulk_insert(settings)
%
% ops_bulk_insert(settings)
%
% Loads paths and/or layers into the database.
% REQUIRES RECORDS AND FRAMES FILES FOR THE SEGMENTS TO LOAD
%
% Input:
%   settings: structure with fields
%     process_type = string ('path', 'layer' or 'both')
%     layer_filter = function which takes one character array argument
%       and returns true for layers which should be inserted (may be
%       undefined or empty if all layers are to be inserted)
%     param_fn = string (absolute path and filename to param spreadsheet)
%     log_base_path = string (absolute path to logs directory)
%     sys_name = string (rds, snow, accum, kuband)
%     path_decimation_spacing = double (decimation spacing for the gps data (100m is standard for CReSIS))
%     location = string (arctic, antarctic)
%     echo = logical (true, false)
%
% Output:
%   none
%
% Author: Kyle W. Purdon
%
% see also run_cr_insert.m

%% PARAM SETUP FOR PATH/LAYER INSERTION
% =========================================================================

if ~isfield(settings,'layer_filter') || isempty(settings.layer_filter)
  settings.layer_filter = inline('~isempty(regexp(x,''.*''))');
end

%% GET username
if ispc
  user_name = getenv('USERNAME');
else
  [~,user_name] = system('whoami');
  user_name = user_name(1:end-1);
end

insert_path_cmd = any(strcmpi(settings.process_type,{'path','both','all'}));
insert_layer_cmd = any(strcmpi(settings.process_type,{'layer','both','all'}));
insert_atm_cmd = any(strcmpi(settings.process_type,{'atm','all'}));

%% Load the param sheet and perform basic input checks
params = read_param_xls(settings.param_fn);
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
    if ~exist(records_fn,'file') && insert_path_cmd
      error('  %s: missing %s\n', param.day_seg, records_fn);
    elseif ~exist(frames_fn,'file') && insert_path_cmd
      error('  %s: missing %s\n', param.day_seg, frames_fn);
    elseif ~exist(layer_dir,'dir') && insert_layer_cmd
      error('  %s: missing %s\n', param.day_seg, layer_dir);
    else
      fprintf('  %s checked\n', param.day_seg);
    end
  end
end

%% CONFIRMATION DIALOG
ops_sys_cmd;
confirm_params = {sprintf('SEASON NAME: \t %s',season_name),'',sprintf('RADAR NAME: \t %s',radar_name),'',...
  sprintf('SYSTEM CONNECTION: \t %s',server_url),'',sprintf('PROCESS TYPE: \t %s',settings.process_type),'',...
  sprintf('SYSTEM NAME: \t %s',settings.sys_name),'',sprintf('DECIMATION SPACING: \t %d meters',settings.path_decimation_spacing),'',...
  sprintf('LOCATION: \t %s',settings.location),''};

confirm_button = questdlg(confirm_params,'Confirm Settings','YES:LOAD','NO:QUIT','NO:QUIT');

switch confirm_button
  case 'NO:QUIT'
    error('PROCESS STOPPED BY USER.');
end


%% PATH INSERTION
% =========================================================================

% Store any days that fail (dont try and load layers for these days)
failed_segments = {};

if insert_path_cmd
  
  % Set Up Status/Time logging
  log_fn = strcat(settings.log_base_path,strcat('pathinsertlog_',datestr(now,'dd.mm.yyyy'),'.',datestr(now,'HH.MM.SS'),'.txt'));
  diary(log_fn);
  
  for param_idx = 1:length(params)
    param = params(param_idx);
    if param.cmd.generic ~= 1
      continue;
    end
    
    segment = param.day_seg;
    
    fprintf('Loading path for segment %s ... ',segment);
    start = tic;
    stop = [];
    
    % Get cell of records and frames files for the current segment
    records_fn = ct_filename_support(param,'','records');
    frames_fn = ct_filename_support(param,'','frames');
    
    % Load the records/frames file
    load(frames_fn);
    records = load(records_fn);
    if isfield(records,'records')
      % Adding this check to support old record (pre-cr1) file format
      records = records.records;
    end
    
    % Set up needed pre-allocations
    frame_start_gps_time = zeros(1,length(frames.frame_idxs));
    frame_stop_gps_time = zeros(1,length(frames.frame_idxs));
    echogram_names = cell(length(frames.frame_idxs),1);
    
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
    
    % GET THE ECHOGRAM PATH (DATA.CRESIS.KU.EDU)
    for echo_idx = 1:length(echogram_names)
      echogram_names{echo_idx} = strcat('ftp://data.cresis.ku.edu/data/',settings.sys_name,'/',season_name,'/images/',segment,'/',segment,sprintf('_%03d_1echo.jpg',echo_idx));
    end
    
    % Save MATLAB computation time
    stop(end+1) = toc(start);
    
    % Get sampling indexes for decimation
    decimation_spacing = settings.path_decimation_spacing;
    spacing_idxs = get_equal_alongtrack_spacing_idxs(geodetic_to_along_track(records.lat,records.lon,records.elev),decimation_spacing);
    
    % CONSTRUCT STRUCTURE FOR JSON
    cr_json.geometry.coordinates = [records.lon(spacing_idxs)' records.lat(spacing_idxs)' records.gps_time(spacing_idxs)'];
    cr_json.properties.location = settings.location;
    cr_json.properties.season = season_name;
    cr_json.properties.radar = radar_name;
    cr_json.properties.segment = segment;
    cr_json.properties.segment_start_gps_time = records.gps_time(1);
    cr_json.properties.segment_stop_gps_time = records.gps_time(end);
    cr_json.properties.elev = records.elev(spacing_idxs);
    cr_json.properties.roll = records.roll(spacing_idxs);
    cr_json.properties.pitch = records.pitch(spacing_idxs);
    cr_json.properties.heading = records.heading(spacing_idxs);
    cr_json.properties.echo_img_url = echogram_names;
    cr_json.properties.frame_count = length(frames.frame_idxs);
    cr_json.properties.frame_start_gps_time = frame_start_gps_time;
    cr_json.properties.frame_stop_gps_time = frame_stop_gps_time;
    
    % SEND DATA FOR INSERT
    [~,message] = ops_create_path(settings.sys_name,cr_json);
    
    % Save Python/SQL computation time
    stop(end+1) = toc(start);
    
    fprintf('%2.2f s\n',sum(stop));
    fprintf('\t-> Status: %s\n',message);
    
    if ~isempty(failed_segments)
      fprintf('\n\n FAILED SEGMENTS:\n')
      for fail_idx = 1:length(failed_segments)
        fprintf('\t %s\n',failed_segments{fail_idx});
      end
      fprintf('Layers will not be loaded for any segment of these days.\n\n');
    end
    clear cr_json;
  end
end


%% LAYER INSERTION

if insert_layer_cmd
  
  % Set Up Status/Time logging
  log_fn = strcat(settings.log_base_path,strcat('layerinsertlog_',datestr(now,'dd.mm.yyyy'),'.',datestr(now,'HH.MM.SS'),'.txt'));
  fid = fopen(log_fn,'w+');
  
  for param_idx = 1:length(params)
    param = params(param_idx);
    if param.cmd.generic ~= 1
      continue;
    end
    
    segment = param.day_seg;
    % Exclude any failing days
    if any(strcmp(segment,failed_segments))
      continue;
    end
    
    fprintf('Loading layers for segment %s ... \n',segment);
    
    % Load the frames file
    frames_fn = ct_filename_support(param,'','frames');
    load(frames_fn);
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
    
    layer_fn_dir = ct_filename_out(param,settings.layer_post_directory,param.day_seg);
    end_gps = 0;
    stop = zeros(length(param.cmd.frms));
    for frame_idx = 1:length(param.cmd.frms)
      start = tic;
      frm = param.cmd.frms(frame_idx);
      layer_fn_name = sprintf('Data_%s_%03d.mat', param.day_seg, frm);
      layer_fn = fullfile(layer_fn_dir, layer_fn_name);
      
      fprintf('\tLoading frame %s ... ',layer_fn(end-18:end-4));
      
      if ~exist(layer_fn,'file')
        warning('Frame %s layer file does not exist', layer_fn_name);
        continue;
      end
      
      ops_layer_points_param = layerdata_to_ops(layer_fn,settings);
      frame = layer_fn(end-18:end-4);
      
      % Remove the overlap in ops_layer_points_param
      if frame_idx >1
        for layer_idx = 1:length(ops_layer_points_param);
          good_idxs = find(ops_layer_points_param(layer_idx).properties.gps_time>end_gps);
          ops_layer_points_param(layer_idx).properties.gps_time = ops_layer_points_param(layer_idx).properties.gps_time(good_idxs);
          ops_layer_points_param(layer_idx).properties.twtt = ops_layer_points_param(layer_idx).properties.twtt(good_idxs);
          ops_layer_points_param(layer_idx).properties.type = ops_layer_points_param(layer_idx).properties.type(good_idxs);
          ops_layer_points_param(layer_idx).properties.quality = ops_layer_points_param(layer_idx).properties.quality(good_idxs);
          lon = ops_layer_points_param(layer_idx).geometry.coordinates(good_idxs',1);
          lat = ops_layer_points_param(layer_idx).geometry.coordinates(good_idxs',2);
          elev = ops_layer_points_param(layer_idx).geometry.coordinates(good_idxs',3);
          ops_layer_points_param(layer_idx).geometry.coordinates = [lon lat elev];
        end
      end
      
      lyr = load(layer_fn,'GPS_time');
      end_gps = lyr.GPS_time(end);
      clear lyr;
      
      % Check if there are complete empty layers
      bad_idxs = [];
      for layer_idx = 1:length(ops_layer_points_param)
        if isempty(ops_layer_points_param(layer_idx).properties.twtt)
          bad_idxs(end+1) = layer_idx;
        end
      end
      
      ops_layer_points_param(bad_idxs) = [];
            
      if isempty(ops_layer_points_param)
        fprintf('STATUS: ALL layers empty. Nothing to insert.\n');
        message = '';
      else
        
        for layer_idx = 1:length(ops_layer_points_param);
          
          % Fill in needed properties
          ops_layer_points_param(layer_idx).properties.username = user_name;
          ops_layer_points_param(layer_idx).properties.location = settings.location;
          ops_layer_points_param(layer_idx).properties.segment = segment;
          ops_layer_points_param(layer_idx).properties.radar = radar_name;
          ops_layer_points_param(layer_idx).properties.yyyymmdd = layer_fn(end-18:end-11);
          
          djtic= tic;
          [status,message] = ops_create_layer_points(settings.sys_name,ops_layer_points_param(layer_idx));
          djtime = toc(djtic);
        
        end
      end
      stop(frame_idx) = toc(start);
      
      fprintf('%2.2f s\n',stop(frame_idx));
      fprintf('\t\t-> Status: %s\n',message);
      
    end
  end
end

%% ATM INSERTION

if insert_atm_cmd
  
  % Set Up Status/Time logging
  log_fn = strcat(settings.log_base_path,strcat('atminsertlog_',datestr(now,'dd.mm.yyyy'),'.',datestr(now,'HH.MM.SS'),'.txt'));
  fid = fopen(log_fn,'w+');
  
  for param_idx = 1:length(params)
    
    start = tic;
    
    param = params(param_idx);
    if param.cmd.generic ~= 1
      continue;
    end
    
    segment = param.day_seg;
    % Exclude any failing days
    if any(strcmp(segment,failed_segments))
      continue;
    end
    
    fprintf('Loading atm layer for segment %s ... \n',segment);
    stop = [];
    
    % Get a cell of atm files for the current segment
    atm_fns = get_filenames_atm(settings.location,segment);
    if isempty(atm_fns)
      warning('No atm data. Skipping segment %s',segment);
      continue;
    end
    
    % Convert atm data to ops format
    ops_atm_param = atm_to_ops(atm_fns,ct_filename_support(param,'','records'));
    
    if isempty(ops_atm_param.properties.gps_time)
      warning('No ATM data aligned with this segment\n');
    else
      % Fill in needed properties
      ops_atm_param.properties.username = user_name;
      ops_atm_param.properties.location = settings.location;
      ops_atm_param.properties.segment = segment;
      ops_atm_param.properties.radar = radar_name;
      ops_atm_param.properties.yyyymmdd = segment(1:8);
      [~,~] = ops_create_layer_points(settings.sys_name,ops_atm_param);   
    end
  end
end

diary off
end