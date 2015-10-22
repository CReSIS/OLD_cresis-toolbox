% ops_update_files_from_database
%
% Updates data files (Surface and Bottom) and records files (Surface)
% from the database.


% Script run_create_posting
%
% Loads the "post" worksheet from the parameter spreadsheet and then calls
% create_posting with this information.
%
% Authors: Theresa Stumpf, John Paden
%
% See also: create_posting.m

fprintf('\n\n========================================================\n');
fprintf('run create posting\n');
fprintf('========================================================\n');

%% User Settings

params = read_param_xls('/users/paden/scripts/branch/params-cr1/snow_param_2012_Antarctica_DC8.xls',[],'post');

update_data = false;
update_records = true;

out_path = 'CSARP_qlook';

%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================
tic;
global gRadar;

% Input checking
if ~exist('params','var')
  error('Use run_master: A struct array of parameters must be passed in\n');
end
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% =====================================================================
% =====================================================================
% Process each of the days
% =====================================================================
% =====================================================================
for param_idx = 1:length(params)
  param = merge_structs(params(param_idx),param_override);
  if isfield(param.cmd,'generic') && param.cmd.generic
    fprintf('ops_update_files_from_database %s\n', param.day_seg);
    
    if isempty(param.post.ops.layers)
      warning('No layers specified');
      continue;
    end
    
    %% Get the layer data for this segment
    ops_sys = ct_output_dir(param.radar_name);
    ops_param = [];
    ops_param.properties.location = param.post.ops.location;
    ops_param.properties.season = param.season_name;
    ops_param.properties.segment = param.day_seg;
    ops_param.properties.return_geom = 'geog';
    for layer_idx = 1:length(param.post.ops.layers)
      ops_param.properties.lyr_name = param.post.ops.layers{layer_idx};
      [~,ops_layer{layer_idx}] = ops_get_layer_points(ops_sys,ops_param);
      ops_layer{layer_idx} = ops_layer{layer_idx}.properties;
    end
    
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
    
    data_fn_dir = ct_filename_out(param,'',out_path);
    
    if update_data
      for frm = reshape(param.cmd.frms,[1 length(param.cmd.frms)])
        frm_id = sprintf('%s_%03d',param.day_seg,frm);
        data_fn = fullfile(data_fn_dir,sprintf('Data_%s.mat',frm_id));
        fprintf('Updating %s\n', data_fn);
        
        %% Interpolate each layer for this segment onto the master layer's gps time
        
        % Get start/stop time of the frame
        ops_param = [];
        ops_param.properties.search_str = frm_id;
        ops_param.properties.location = param.post.ops.location;
        [~,frame_info] = ops_search_frames(ops_sys,ops_param);
        
        % Get the master layer geographic data
        idxs_in_frm = find(ops_layer{1}.gps_time >= frame_info.properties.start_gps_time ...
          & ops_layer{1}.gps_time <= frame_info.properties.stop_gps_time ...
          & ops_layer{1}.type == 2);
        
        lay = load(data_fn,'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom');
        
        master_along_track = geodetic_to_along_track(lay.Latitude, ...
          lay.Longitude, lay.Elevation);
        % HACK FOR NOW:
        if length(ops_layer) >= 2
          output_map = [2 1 3:length(ops_layer)];
        else
          output_map = 1;
        end
        % Interpolate each layer onto the echogram
        for layer_idx = 1:length(ops_layer)
          % Get just the automated points
          auto_idxs = find(ops_layer{layer_idx}.type == 2);
          [unique_vals unique_idxs] = unique(ops_layer{layer_idx}.gps_time(auto_idxs));
          auto_idxs = auto_idxs(unique_idxs);
          
          slave = ops_layer{layer_idx}.gps_time(auto_idxs);
          tie_idx = find(slave >= lay.GPS_time(1) ...
            & slave <= lay.GPS_time(end),2);
          % Interpolate the points onto the echogram GPS time
          if length(tie_idx) < 2
            % Not enough points to interpolate, so set to all NaN
            lay.layerData{output_map(layer_idx)}.value{2}.data ...
              = NaN*zeros(size(lay.GPS_time));
            lay.layerData{output_map(layer_idx)}.quality ...
              = NaN*zeros(size(lay.GPS_time));
            
          else
            lay.layerData{output_map(layer_idx)}.value{2}.data ...
              = interp1(ops_layer{layer_idx}.gps_time(auto_idxs), ...
              ops_layer{layer_idx}.twtt(auto_idxs), lay.GPS_time, 'linear');
            lay.layerData{output_map(layer_idx)}.quality ...
              = interp1(ops_layer{layer_idx}.gps_time(auto_idxs), ...
              ops_layer{layer_idx}.quality(auto_idxs), lay.GPS_time, 'nearest');
            
            % Remove data gaps (i.e. no data within a certain amount of space)
            slave_along_track = geodetic_to_along_track(ops_layer{layer_idx}.lat(auto_idxs), ...
              ops_layer{layer_idx}.lon(auto_idxs), ...
              ops_layer{layer_idx}.elev(auto_idxs));
            
            % Use the first slave point inside the master domain (tie_idx(1))
            % to align the two dist axes.
            slave_along_track_tie_master_val = interp1(lay.GPS_time, ...
              master_along_track, slave(tie_idx(1)));
            slave_along_track = slave_along_track  + (slave_along_track_tie_master_val - slave_along_track(tie_idx(1)));
            gap_data_idxs = data_gaps_check_mex(master_along_track, slave_along_track, param.post.ops.gaps_dist(1), param.post.ops.gaps_dist(2));
            lay.layerData{output_map(layer_idx)}.value{2}.data(gap_data_idxs) = NaN;
            lay.layerData{output_map(layer_idx)}.quality(gap_data_idxs) = NaN;
          end
          
          lay.layerData{output_map(layer_idx)}.value{1}.data = NaN*zeros(size(lay.GPS_time));
          manual_idxs = find(ops_layer{layer_idx}.type == 1);
          manual_idxs_map = interp1(lay.GPS_time, 1:length(lay.GPS_time), ...
            ops_layer{layer_idx}.gps_time(manual_idxs),'nearest');
          manual_idxs = manual_idxs(~isnan(manual_idxs_map));
          manual_idxs_map = manual_idxs_map(~isnan(manual_idxs_map));
          lay.layerData{output_map(layer_idx)}.value{1}.data(manual_idxs_map) ...
            = ops_layer{layer_idx}.twtt(manual_idxs);
          
          if strcmpi(param.post.ops.layers{layer_idx},'surface')
            lay.Elevation = interp1(slave, ops_layer{layer_idx}.elev(auto_idxs), lay.GPS_time);
            lay.Elevation(isnan(lay.Elevation)) = interp1(slave, ...
              ops_layer{layer_idx}.elev(auto_idxs), ...
              lay.GPS_time(isnan(lay.Elevation)), 'nearest','extrap');
          end
        end
        Surface = lay.layerData{1}.value{2}.data;
        if length(lay.layerData) > 1
          Bottom = lay.layerData{2}.value{2}.data;
          save(data_fn,'-append','Surface','Bottom');
        else
          save(data_fn,'-append','Surface');
        end
        
      end
    end
    
    if update_records
      
      records_fn = ct_filename_support(param,'','records');
      
      %% Interpolate each layer for this segment onto the master layer's gps time
      
      records = load(records_fn,'gps_time','lat','lon','elev','surface');
      lay.Latitude = records.lat;
      lay.Longitude = records.lon;
      lay.Elevation = records.elev;
      lay.GPS_time = records.gps_time;
      lay.Surface = records.surface;
      
      master_along_track = geodetic_to_along_track(lay.Latitude, ...
        lay.Longitude, lay.Elevation);
      % HACK FOR NOW:
      if length(ops_layer) >= 2
        output_map = [2 1 3:length(ops_layer)];
      else
        output_map = 1;
      end
      % Interpolate each layer onto the echogram
      for layer_idx = 1:length(ops_layer)
        % Get just the automated points
        auto_idxs = find(ops_layer{layer_idx}.type == 2);
        [unique_vals unique_idxs] = unique(ops_layer{layer_idx}.gps_time(auto_idxs));
        auto_idxs = auto_idxs(unique_idxs);
        
        slave = ops_layer{layer_idx}.gps_time(auto_idxs);
        tie_idx = find(slave >= lay.GPS_time(1) ...
          & slave <= lay.GPS_time(end),2);
        % Interpolate the points onto the echogram GPS time
        if length(tie_idx) < 2
          % Not enough points to interpolate, so set to all NaN
          lay.layerData{output_map(layer_idx)}.value{2}.data ...
            = NaN*zeros(size(lay.GPS_time));
          lay.layerData{output_map(layer_idx)}.quality ...
            = NaN*zeros(size(lay.GPS_time));
          
        else
          lay.layerData{output_map(layer_idx)}.value{2}.data ...
            = interp1(ops_layer{layer_idx}.gps_time(auto_idxs), ...
            ops_layer{layer_idx}.twtt(auto_idxs), lay.GPS_time, 'linear');
          lay.layerData{output_map(layer_idx)}.quality ...
            = interp1(ops_layer{layer_idx}.gps_time(auto_idxs), ...
            ops_layer{layer_idx}.quality(auto_idxs), lay.GPS_time, 'nearest');
          
          % Remove data gaps (i.e. no data within a certain amount of space)
          slave_along_track = geodetic_to_along_track(ops_layer{layer_idx}.lat(auto_idxs), ...
            ops_layer{layer_idx}.lon(auto_idxs), ...
            ops_layer{layer_idx}.elev(auto_idxs));
          
          % Use the first slave point inside the master domain (tie_idx(1))
          % to align the two dist axes.
          slave_along_track_tie_master_val = interp1(lay.GPS_time, ...
            master_along_track, slave(tie_idx(1)));
          slave_along_track = slave_along_track  + (slave_along_track_tie_master_val - slave_along_track(tie_idx(1)));
          gap_data_idxs = data_gaps_check_mex(master_along_track, slave_along_track, param.post.ops.gaps_dist(1), param.post.ops.gaps_dist(2));
          lay.layerData{output_map(layer_idx)}.value{2}.data(gap_data_idxs) = NaN;
          lay.layerData{output_map(layer_idx)}.quality(gap_data_idxs) = NaN;
        end
        
        lay.layerData{output_map(layer_idx)}.value{1}.data = NaN*zeros(size(lay.GPS_time));
        manual_idxs = find(ops_layer{layer_idx}.type == 1);
        manual_idxs_map = interp1(lay.GPS_time, 1:length(lay.GPS_time), ...
          ops_layer{layer_idx}.gps_time(manual_idxs),'nearest');
        manual_idxs = manual_idxs(~isnan(manual_idxs_map));
        manual_idxs_map = manual_idxs_map(~isnan(manual_idxs_map));
        lay.layerData{output_map(layer_idx)}.value{1}.data(manual_idxs_map) ...
          = ops_layer{layer_idx}.twtt(manual_idxs);
        
        if strcmpi(param.post.ops.layers{layer_idx},'surface')
          lay.Elevation = interp1(slave, ops_layer{layer_idx}.elev(auto_idxs), lay.GPS_time);
          lay.Elevation(isnan(lay.Elevation)) = interp1(slave, ...
            ops_layer{layer_idx}.elev(auto_idxs), ...
            lay.GPS_time(isnan(lay.Elevation)), 'nearest','extrap');
        end
      end
      surface = lay.layerData{1}.value{2}.data;
      save(records_fn,'-append','surface');
    end
  end
end

return;
