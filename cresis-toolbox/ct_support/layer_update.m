% function layer_update(param,param_override)
% layer_update(param,param_override)
%
% This function updates the layer according to the param struct. The
% primary purpose is to update an old layer file to the new file format.
%
% param: parameter spreadsheet structure array
%
% Examples: See run_all_layer_update.m
%
% Author: John Paden
%
% See also: run_all_layer_update, layer_update

param = merge_structs(param,param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Input checks

param.layer_update.in_path = 'layerData';
% param.layer_update.in_path = 'layer';
param.layer_update.out_path = 'layer';

param.layer_update.frames_records_en = true;

%% General Setup

if param.layer_update.frames_records_en
  
  % Load records file
  records_fn = ct_filename_support(param,'','records');
  records = load(records_fn);
  
  % Load frames file
  frames_fn = ct_filename_support(param,'','frames');
  frames = load(frames_fn);
  if ~isfield(frames,'frame_idxs')
    warning('Frames file format issue. Is this an old frames file format? Do not continue until frames_update.m run first.');
    if 0
      keyboard
    else
      % Convert to new format
      frames = frames.frames;
      frames.gps_time = [records.gps_time(frames.frame_idxs), records.gps_time(end)];
    end
  end
  
  trajectory_param = struct('gps_source',records.gps_source, ...
    'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
    'tx_weights', [], 'lever_arm_fh', param.radar.lever_arm_fh);
  records = trajectory_with_leverarm(records,trajectory_param);
else
  keyboard
end

% Load layer organizer file:
% * create if it does not exist
% * check that contents are correct
layer_organizer_fn = fullfile(ct_filename_out(param,param.layer_update.in_path,''),sprintf('layer_%s.mat',param.day_seg));
fprintf('Loading layer organizer: %s\n', layer_organizer_fn);
save_layer_organizer_file = false;
if exist(layer_organizer_fn,'file')
  try
    layer_organizer = load(layer_organizer_fn);
  catch
    save_layer_organizer_file = true;
    layer_organizer = [];
  end
else
  save_layer_organizer_file = true;
  layer_organizer = [];
end
layer_organizer.layer_organizer_fn = layer_organizer_fn;

if ~isfield(layer_organizer,'lyr_name')
  layer_organizer.lyr_name = {};
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer,'lyr_age')
  layer_organizer.lyr_age = nan(size(layer_organizer.lyr_name));
  save_layer_organizer_file = true;
else
  len_diff = length(layer_organizer.lyr_name) - length(layer_organizer.lyr_age);
  if len_diff < 0
    layer_organizer.lyr_age = layer_organizer.lyr_age(1:length(layer_organizer.lyr_name));
    save_layer_organizer_file = true;
  elseif len_diff > 0
    layer_organizer.lyr_age(end+1:length(layer_organizer.lyr_name)) = NaN;
    save_layer_organizer_file = true;
  end
end

if ~isfield(layer_organizer,'lyr_age_source')
  layer_organizer.lyr_age_source = cell(size(layer_organizer.lyr_name));
  save_layer_organizer_file = true;
else
  len_diff = length(layer_organizer.lyr_name) - length(layer_organizer.lyr_age_source);
  if len_diff < 0
    layer_organizer.lyr_age_source = layer_organizer.lyr_age_source(1:length(layer_organizer.lyr_name));
    save_layer_organizer_file = true;
  elseif len_diff > 0
    layer_organizer.lyr_age_source(end+1:length(layer_organizer.lyr_name)) = cell([1 len_diff]);
    save_layer_organizer_file = true;
  end
end

if ~isfield(layer_organizer,'lyr_desc')
  layer_organizer.lyr_desc = cellfun(@char,cell(size(layer_organizer.lyr_name)),'UniformOutput',false);
  save_layer_organizer_file = true;
else
  len_diff = length(layer_organizer.lyr_name) - length(layer_organizer.lyr_desc);
  if len_diff < 0
    layer_organizer.lyr_desc = layer_organizer.lyr_desc(1:length(layer_organizer.lyr_name));
    save_layer_organizer_file = true;
  elseif len_diff > 0
    layer_organizer.lyr_desc(end+1:length(layer_organizer.lyr_name)) = cellfun(@char,cell([1 len_diff]),'UniformOutput',false);
    save_layer_organizer_file = true;
  end
end

if ~isfield(layer_organizer,'lyr_group_name')
  layer_organizer.lyr_group_name = cellfun(@char,cell(size(layer_organizer.lyr_name)),'UniformOutput',false);
  save_layer_organizer_file = true;
else
  len_diff = length(layer_organizer.lyr_name) - length(layer_organizer.lyr_group_name);
  if len_diff < 0
    layer_organizer.lyr_group_name = layer_organizer.lyr_group_name(1:length(layer_organizer.lyr_name));
    save_layer_organizer_file = true;
  elseif len_diff > 0
    layer_organizer.lyr_group_name(end+1:length(layer_organizer.lyr_name)) = cellfun(@char,cell([1 len_diff]),'UniformOutput',false);
    save_layer_organizer_file = true;
  end
end

if ~isfield(layer_organizer,'lyr_id')
  layer_organizer.lyr_id = 1:length(layer_organizer.lyr_name);
  save_layer_organizer_file = true;
else
  len_diff = length(layer_organizer.lyr_name) - length(layer_organizer.lyr_id);
  if len_diff < 0
    layer_organizer.lyr_id = layer_organizer.lyr_id(1:length(layer_organizer.lyr_name));
    save_layer_organizer_file = true;
  elseif len_diff > 0
    layer_organizer.lyr_id(end+1:length(layer_organizer.lyr_name)) = max(layer_organizer.lyr_id)+(1:len_diff);
    save_layer_organizer_file = true;
  end
end

if ~isfield(layer_organizer,'lyr_order')
  layer_organizer.lyr_order = 1:length(layer_organizer.lyr_name);
  save_layer_organizer_file = true;
else
  len_diff = length(layer_organizer.lyr_name) - length(layer_organizer.lyr_order);
  if len_diff < 0
    layer_organizer.lyr_order = layer_organizer.lyr_order(1:length(layer_organizer.lyr_name));
    save_layer_organizer_file = true;
  elseif len_diff > 0
    layer_organizer.lyr_order(end+1:length(layer_organizer.lyr_name)) = max(layer_organizer.lyr_order)+(1:len_diff);
    save_layer_organizer_file = true;
  end
end

if ~isfield(layer_organizer,'file_version') || ~strcmp(layer_organizer.file_version,'1')
  layer_organizer.file_version = '1';
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer,'file_type') || ~strcmp(layer_organizer.file_type,'layer_organizer')
  layer_organizer.file_type = 'layer_organizer';
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer,'param') || isempty(layer_organizer.param)
  layer_organizer.param = [];
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer.param,'day_seg') || ~strcmp(layer_organizer.param.day_seg,param.day_seg)
  layer_organizer.param.day_seg = param.day_seg;
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer.param,'radar_name') || ~strcmp(layer_organizer.param.radar_name,param.radar_name)
  layer_organizer.param.radar_name = param.radar_name;
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer.param,'season_name') || ~strcmp(layer_organizer.param.season_name,param.season_name)
  layer_organizer.param.season_name = param.season_name;
  save_layer_organizer_file = true;
end

if ~isfield(layer_organizer.param,'sw_version')
  save_layer_organizer_file = true;
end
layer_organizer.param.sw_version = param.sw_version;

if ~strcmp(param.layer_update.out_path,param.layer_update.in_path)
  save_layer_organizer_file = true;
end

% layer_organizer did not exist or had bad contents, update the file with
% corrected contents
% if save_layer_organizer_file == true
%   fprintf('Updating layer organizer file.\n');
%   layer_organizer_out_fn = fullfile(ct_filename_out(param,param.layer_update.out_path,''),sprintf('layer_%s.mat',param.day_seg));
%   ct_save(layer_organizer_out_fn,'-struct','layer_organizer','file_type','file_version','gps_source','lyr_name','lyr_age','lyr_age_source', ...
%     'lyr_desc','lyr_group_name','lyr_id','lyr_order','param');
% end

layer_fn_dir = ct_filename_out(param,param.layer_update.in_path,'');
layer_out_fn_dir = ct_filename_out(param,param.layer_update.out_path,'');
fprintf('Loading layer files: %s\n', layer_fn_dir);
num_frm = length(frames.frame_idxs);
for frm = 1:num_frm
  layer_fn=fullfile(layer_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  lay = load(layer_fn);
  new_layer_ids = [];
  save_layer_data_file = 0;
  if isfield(lay,'file_version')
    % New file format
    Nx = length(lay.gps_time);
    new_lay = [];
    new_lay.gps_time = lay.gps_time;
    
    % Check if fields are the wrong length
    if length(lay.elev) ~= Nx
      save_layer_data_file = bitor(2,save_layer_data_file);
    end
    if length(lay.lat) ~= Nx
      save_layer_data_file = bitor(2,save_layer_data_file);
    end
    if length(lay.lon) ~= Nx
      save_layer_data_file = bitor(2,save_layer_data_file);
    end
    
    if param.layer_update.frames_records_en
      % Recreate lat, lon, elev data using records file if available
      new_lay.elev = interp1(records.gps_time,records.elev,lay.gps_time);
      new_lay.lat = interp1(records.gps_time,records.lat,lay.gps_time);
      new_lay.lon = interp1(records.gps_time,records.lon,lay.gps_time);
      
    else
      % Check if fields are too short
      if length(lay.elev) < Nx
        lat.elev(end+1:Nx) = NaN;
      end
      if length(lay.lat) < Nx
        lat.lat(end+1:Nx) = NaN;
      end
      if length(lay.lon) < Nx
        lat.lon(end+1:Nx) = NaN;
      end
      
      % Truncate fields incase they are too long
      new_lay.elev(1:Nx) = lay.elev(1:Nx);
      new_lay.lat(1:Nx) = lay.lat(1:Nx);
      new_lay.lon(1:Nx) = lay.lon(1:Nx);
    end
    
    if ~strcmp(lay.file_type,'layer')
      save_layer_data_file = bitor(1,save_layer_data_file);
    end
    if ~strcmp(lay.file_version,'1')
      save_layer_data_file = bitor(1,save_layer_data_file);
    end
    new_lay.file_type = 'layer';
    new_lay.file_version = '1';
    
    new_lay.id = lay.id(:);
    for lay_idx = 1:length(new_lay.id)
      
      match_idx = find(new_lay.id(lay_idx) == layer_organizer.lyr_id);
      if isempty(match_idx)
        fprintf('layer file with id field that is not present in layer organizer file. Layer %d had an invalid .id field or did not have an id field.\n', lay_idx);
        save_layer_data_file = bitor(save_layer_data_file,1);
        % Add the layer to the layer_organizer
        % Ensure a unique name
        % -------------------------------------------------------------------
        duplicate_idx = 1;
        if lay_idx == 1
          % Enforce standard name/group for old file format
          base_name = 'surface';
          name = base_name;
          group_name = 'standard';
        elseif lay_idx == 2
          % Enforce standard name/group for old file format
          base_name = 'bottom';
          name = base_name;
          group_name = 'standard';
        else
          base_name = 'auto';
          name = sprintf('%s_%03d',base_name,duplicate_idx);
          group_name = '';
        end
        while any(strcmp(name,layer_organizer.lyr_name))
          % This is a duplicate layer name, this loop searches for a unique
          % name
          duplicate_idx = duplicate_idx + 1;
          name = sprintf('%s_%03d',base_name,duplicate_idx);
        end
        layer_organizer.lyr_age(end+1) = NaN;
        layer_organizer.lyr_age_source{end+1} = struct('age',{},'source',{},'type',{});
        layer_organizer.lyr_desc{end+1} = '';
        layer_organizer.lyr_group_name{end+1} = group_name;
        layer_organizer.lyr_id(end+1) = new_lay.id(lay_idx);
        layer_organizer.lyr_name{end+1} = name;
        new_order = max(layer_organizer.lyr_order) + 1;
        if isempty(new_order)
          new_order = 1;
        end
        layer_organizer.lyr_order(end+1) = new_order;
        match_idx = length(layer_organizer.lyr_id);
      end
      
      % Ensure a unique id
      % -------------------------------------------------------------------
      while any(lay.id(lay_idx) == new_layer_ids)
        warning('layer file with duplicate id. Layer %d has a duplicate .id field.', lay_idx);
        save_layer_data_file = bitor(save_layer_data_file,1);
        % Add the layer to the layer_organizer
        base_name = 'auto';
        % Ensure a unique name
        % -------------------------------------------------------------------
        duplicate_idx = 1;
        name = sprintf('%s_%03d',base_name,duplicate_idx);
        while any(strcmp(name,layer_organizer.lyr_name))
          % This is a duplicate layer name, this loop searches for a unique
          % name
          duplicate_idx = duplicate_idx + 1;
          name = sprintf('%s_%03d',base_name,duplicate_idx);
        end
        layer_organizer.lyr_age(end+1) = NaN;
        layer_organizer.lyr_age_source(end+1) = struct('age',{},'source',{},'type',{});
        layer_organizer.lyr_desc{end+1} = '';
        layer_organizer.lyr_group_name{end+1} = '';
        layer_organizer.lyr_id(end+1) = max(layer_organizer.lyr_id) + 1;
        layer_organizer.lyr_name{end+1} = name;
        layer_organizer.lyr_order(end+1) = max(layer_organizer.lyr_order) + 1;
        match_idx = length(layer_organizer.lyr_id);
      end
      new_layer_ids(end+1) = layer_organizer.lyr_id(match_idx);
      
      if ~isfield(lay,'quality') || size(lay.quality,1) < lay_idx
        save_layer_data_file = bitor(2,save_layer_data_file);
        new_lay.quality(lay_idx,1:Nx) = ones(1,Nx);
      elseif size(lay.quality,2) == Nx
        new_lay.quality(lay_idx,1:Nx) = lay.quality(lay_idx,:);
      else
        save_layer_data_file = bitor(2,save_layer_data_file);
        Nx_cur = size(lay.quality,2);
        if Nx_cur < Nx
          new_lay.quality(lay_idx,1:Nx) = [lay.quality(lay_idx,Nx), ones(1,Nx-Nx_cur)];
        else
          new_lay.quality(lay_idx,1:Nx) = lay.quality(lay_idx,1:Nx);
        end
      end
      
      if ~isfield(lay,'type') || size(lay.type,1) < lay_idx
        save_layer_data_file = bitor(2,save_layer_data_file);
        new_lay.type(lay_idx,1:Nx) = 2*ones(1,Nx);
      elseif size(lay.type,2) == Nx
        new_lay.type(lay_idx,1:Nx) = lay.type(lay_idx,:);
      else
        save_layer_data_file = bitor(2,save_layer_data_file);
        Nx_cur = size(lay.type,2);
        if Nx_cur < Nx
          new_lay.type(lay_idx,1:Nx) = [lay.type(lay_idx,Nx), 2*ones(1,Nx-Nx_cur)];
        else
          new_lay.type(lay_idx,1:Nx) = lay.type(lay_idx,1:Nx);
        end
      end
      
      if ~isfield(lay,'twtt') || size(lay.twtt,1) < lay_idx
        save_layer_data_file = bitor(2,save_layer_data_file);
        new_lay.twtt(lay_idx,1:Nx) = nan(1,Nx);
      elseif size(lay.twtt,2) == Nx
        new_lay.twtt(lay_idx,1:Nx) = lay.twtt(lay_idx,:);
      else
        save_layer_data_file = bitor(2,save_layer_data_file);
        Nx_cur = size(lay.twtt,2);
        if Nx_cur < Nx
          new_lay.twtt(lay_idx,1:Nx) = [lay.twtt(lay_idx,Nx), nan(1,Nx-Nx_cur)];
        else
          new_lay.twtt(lay_idx,1:Nx) = lay.twtt(lay_idx,1:Nx);
        end
      end
    end
  else
    % Old file format
    save_layer_data_file = bitor(save_layer_data_file,1);
    
    % Remove data that is not contained within frame boundaries
    frms_mask = false(size(lay.GPS_time));
    if frm < length(frames.frame_idxs)
      frms_mask(lay.GPS_time >= frames.gps_time(frm)...
        & lay.GPS_time < frames.gps_time(frm+1)) = true;
    else
      frms_mask(lay.GPS_time >= frames.gps_time(frm)...
        & lay.GPS_time <= frames.gps_time(frm+1)) = true;
    end
    Nx_old = length(lay.GPS_time);
    
    if param.layer_update.frames_records_en
      % Recreate lat, lon, elev data using records file if available
      lay.GPS_time = lay.GPS_time(frms_mask);
      lay.Elevation = interp1(records.gps_time,records.elev,lay.GPS_time);
      lay.Latitude = interp1(records.gps_time,records.lat,lay.GPS_time);
      lay.Longitude = interp1(records.gps_time,records.lon,lay.GPS_time);
      
    else
      % Just apply the frame mask if records are not available
      try
        lay.Elevation = lay.Elevation(frms_mask);
      catch ME
        keyboard
      end
      try
        lay.Latitude = lay.Latitude(frms_mask);
      catch ME
        keyboard
      end
      try
        lay.Longitude = lay.Longitude(frms_mask);
      catch ME
        keyboard
      end
      lay.GPS_time = lay.GPS_time(frms_mask);
    end
    Nx = length(lay.GPS_time);
    
    for lay_idx = 1:length(lay.layerData)
      if ~isfield(lay.layerData{lay_idx},'value')
        lay.layerData{lay_idx}.value = {};
        save_layer_data_file = bitor(save_layer_data_file,1);
      end
      if isempty(lay.layerData{lay_idx}.value)
        lay.layerData{lay_idx}.value{1}.data = nan(1,Nx);
        save_layer_data_file = bitor(save_layer_data_file,1);
      end
      if length(lay.layerData{lay_idx}.value) < 2
        lay.layerData{lay_idx}.value{2}.data = nan(1,Nx);
        save_layer_data_file = bitor(save_layer_data_file,1);
      end
      if ~isfield(lay.layerData{lay_idx},'quality')
        lay.layerData{lay_idx}.quality = ones(1,Nx);
        save_layer_data_file = bitor(save_layer_data_file,1);
      end
      if isfield(lay.layerData{lay_idx},'name') && ~isfield(lay.layerData{lay_idx},'id')
        % Old file format, switch name to id
        match_idx = strmatch(lay.layerData{lay_idx}.name,layer_organizer.lyr_name,'exact');
        save_layer_data_file = bitor(save_layer_data_file,1);
        
        if isempty(match_idx)
          warning('layerData file with name field that is not present in layer organizer file. Adding layer %d as .name=%s.', lay_idx, lay.layerData{lay_idx}.name);
          % Add the layer to the layer_organizer
          base_name = lay.layerData{lay_idx}.name;
          % Ensure a unique name
          % -------------------------------------------------------------------
          duplicate_idx = 1;
          name = base_name;
          while any(strcmp(name,layer_organizer.lyr_name))
            % This is a duplicate layer name, this loop searches for a unique
            % name
            duplicate_idx = duplicate_idx + 1;
            name = sprintf('%s_%03d',base_name,duplicate_idx);
          end
          
          layer_organizer.lyr_age(end+1) = NaN;
          layer_organizer.lyr_age_source{end+1} = struct('age',{},'source',{},'type',{});
          layer_organizer.lyr_desc{end+1} = '';
          layer_organizer.lyr_group_name{end+1} = '';
          if isempty(layer_organizer.lyr_id)
            layer_organizer.lyr_id(end+1) = 1;
          else
            layer_organizer.lyr_id(end+1) = max(layer_organizer.lyr_id) + 1;
          end
          layer_organizer.lyr_name{end+1} = name;
          if isempty(layer_organizer.lyr_order)
            layer_organizer.lyr_order(end+1) = 1;
          else
            layer_organizer.lyr_order(end+1) = max(layer_organizer.lyr_order) + 1;
          end
          
          match_idx = length(layer_organizer.lyr_id);
        end
        lay.layerData{lay_idx}.id = layer_organizer.lyr_id(match_idx);
        
      else
        if ~isfield(lay.layerData{lay_idx},'id')
          lay.layerData{lay_idx}.id = lay_idx;
        end
        match_idx = find(lay.layerData{lay_idx}.id == layer_organizer.lyr_id);
        if isempty(match_idx)
          fprintf('layerData file with id field that is not present in layer organizer file. Layer %d had an invalid .id field or did not have an id field.\n', lay_idx);
          save_layer_data_file = bitor(save_layer_data_file,1);
          % Add the layer to the layer_organizer
          % Ensure a unique name
          % -------------------------------------------------------------------
          duplicate_idx = 1;
          if lay_idx == 1
            % Enforce standard name/group for old file format
            base_name = 'surface';
            name = base_name;
            group_name = 'standard';
          elseif lay_idx == 2
            % Enforce standard name/group for old file format
            base_name = 'bottom';
            name = base_name;
            group_name = 'standard';
          else
            base_name = 'auto';
            name = sprintf('%s_%03d',base_name,duplicate_idx);
            group_name = '';
          end
          while any(strcmp(name,layer_organizer.lyr_name))
            % This is a duplicate layer name, this loop searches for a unique
            % name
            duplicate_idx = duplicate_idx + 1;
            name = sprintf('%s_%03d',base_name,duplicate_idx);
          end
          layer_organizer.lyr_age(end+1) = NaN;
          layer_organizer.lyr_age_source{end+1} = struct('age',{},'source',{},'type',{});
          layer_organizer.lyr_desc{end+1} = '';
          layer_organizer.lyr_group_name{end+1} = group_name;
          layer_organizer.lyr_id(end+1) = lay.layerData{lay_idx}.id;
          layer_organizer.lyr_name{end+1} = name;
          new_order = max(layer_organizer.lyr_order) + 1;
          if isempty(new_order)
            new_order = 1;
          end
          layer_organizer.lyr_order(end+1) = new_order;
          save_layer_organizer_file = true;
          match_idx = length(layer_organizer.lyr_id);
        end
      end
      
      % Ensure a unique id
      % -------------------------------------------------------------------
      while any(lay.layerData{lay_idx}.id == new_layer_ids)
        warning('layerData file with duplicate id. Layer %d has a duplicate .id field.', lay_idx);
        save_layer_data_file = bitor(save_layer_data_file,1);
        % Add the layer to the layer_organizer
        base_name = 'auto';
        % Ensure a unique name
        % -------------------------------------------------------------------
        duplicate_idx = 1;
        name = sprintf('%s_%03d',base_name,duplicate_idx);
        while any(strcmp(name,layer_organizer.lyr_name))
          % This is a duplicate layer name, this loop searches for a unique
          % name
          duplicate_idx = duplicate_idx + 1;
          name = sprintf('%s_%03d',base_name,duplicate_idx);
        end
        layer_organizer.lyr_age(end+1) = NaN;
        layer_organizer.lyr_age_source(end+1) = struct('age',{},'source',{},'type',{});
        layer_organizer.lyr_desc{end+1} = '';
        layer_organizer.lyr_group_name{end+1} = '';
        layer_organizer.lyr_id(end+1) = max(layer_organizer.lyr_id) + 1;
        layer_organizer.lyr_name{end+1} = name;
        layer_organizer.lyr_order(end+1) = max(layer_organizer.lyr_order) + 1;
        save_layer_organizer_file = true;
        match_idx = length(layer_organizer.lyr_id);
      end
      new_layer_ids(end+1) = layer_organizer.lyr_id(match_idx);
      
      % Ensure all fields are consistent in length
      % -------------------------------------------------------------------
      % Too short:
      if length(lay.layerData{lay_idx}.quality) < Nx_old
        save_layer_data_file = bitor(save_layer_data_file,2);
        lay.layerData{lay_idx}.quality(end+1:Nx_old) = NaN;
      end
      if length(lay.layerData{lay_idx}.value{1}.data) < Nx_old
        save_layer_data_file = bitor(save_layer_data_file,2);
        lay.layerData{lay_idx}.value{1}.data(end+1:Nx_old) = NaN;
      end
      if length(lay.layerData{lay_idx}.value{2}.data) < Nx_old
        save_layer_data_file = bitor(save_layer_data_file,2);
        lay.layerData{lay_idx}.value{2}.data(end+1:Nx_old) = NaN;
      end
      % Too long:
      if length(lay.layerData{lay_idx}.quality) > Nx_old
        save_layer_data_file = bitor(save_layer_data_file,2);
        lay.layerData{lay_idx}.quality = lay.layerData{lay_idx}.quality(1:Nx_old);
      end
      if length(lay.layerData{lay_idx}.value{1}.data) > Nx_old
        save_layer_data_file = bitor(save_layer_data_file,2);
        lay.layerData{lay_idx}.value{1}.data = lay.layerData{lay_idx}.value{1}.data(1:Nx_old);
      end
      if length(lay.layerData{lay_idx}.value{2}.data) > Nx_old
        save_layer_data_file = bitor(save_layer_data_file,2);
        lay.layerData{lay_idx}.value{2}.data = lay.layerData{lay_idx}.value{2}.data(1:Nx_old);
      end
      lay.layerData{lay_idx}.quality = lay.layerData{lay_idx}.quality(frms_mask);
      lay.layerData{lay_idx}.value{1}.data = lay.layerData{lay_idx}.value{1}.data(frms_mask);
      lay.layerData{lay_idx}.value{2}.data = lay.layerData{lay_idx}.value{2}.data(frms_mask);
    end
    
    new_lay.elev = lay.Elevation;
    new_lay.file_type = 'layer';
    new_lay.file_version = '1';
    new_lay.gps_time = lay.GPS_time;
    new_lay.lat = lay.Latitude;
    new_lay.lon = lay.Longitude;
    new_lay.quality = zeros(length(lay.layerData),Nx,'uint8');
    new_lay.twtt = zeros(length(lay.layerData),Nx);
    new_lay.type = zeros(length(lay.layerData),Nx,'uint8');
    for lay_idx = 1:length(lay.layerData)
      new_lay.id(lay_idx,1) = lay.layerData{lay_idx}.id;
      new_lay.quality(lay_idx,1:Nx) = lay.layerData{lay_idx}.quality;
      new_lay.twtt(lay_idx,1:Nx) = lay.layerData{lay_idx}.value{2}.data;
      new_lay.type(lay_idx,1:Nx) = 1 + ~isfinite(lay.layerData{lay_idx}.value{1}.data);
    end
  end
  
  if ~isfield(new_lay,'param') || isempty(new_lay.param)
    new_lay.param = [];
    save_layer_data_file = bitor(save_layer_data_file,8);
  end
  
  if ~isfield(new_lay.param,'day_seg') || ~strcmp(new_lay.param.day_seg,param.day_seg)
    new_lay.param.day_seg = param.day_seg;
    save_layer_data_file = bitor(save_layer_data_file,8);
  end
  
  if ~isfield(new_lay.param,'radar_name') || ~strcmp(new_lay.param.radar_name,param.radar_name)
    new_lay.param.radar_name = param.radar_name;
    save_layer_data_file = bitor(save_layer_data_file,8);
  end
  
  if ~isfield(new_lay.param,'season_name') || ~strcmp(new_lay.param.season_name,param.season_name)
    new_lay.param.season_name = param.season_name;
    save_layer_data_file = bitor(save_layer_data_file,8);
  end
  
  if ~isfield(new_lay.param,'sw_version')
    save_layer_data_file = bitor(save_layer_data_file,8);
  end
  new_lay.param.sw_version = param.sw_version;
  
  if ~isfield(new_lay,'gps_source') || ~strcmp(new_lay.gps_source,records.gps_source)
    new_lay.gps_source = records.gps_source;
    save_layer_data_file = bitor(save_layer_data_file,8);
  end
  
  if param.layer_update.frames_records_en
    if ~isfield(new_lay.param,'radar') || isempty(new_lay.param.radar)
      new_lay.param.radar = [];
      save_layer_data_file = bitor(save_layer_data_file,8);
    end
    
    if ~isfield(new_lay.param.radar,'lever_arm_fh') || isempty(new_lay.param.radar.lever_arm_fh)
      new_lay.param.radar.lever_arm_fh = param.radar.lever_arm_fh;
      save_layer_data_file = bitor(save_layer_data_file,8);
    end
    
    if ~isfield(new_lay.param,'records') || isempty(new_lay.param.records)
      new_lay.param.records = [];
      save_layer_data_file = bitor(save_layer_data_file,8);
    end
    
    if ~isfield(new_lay.param.records,'gps') || isempty(new_lay.param.records.gps)
      new_lay.param.records.gps = [];
      save_layer_data_file = bitor(save_layer_data_file,8);
    end
    
    if ~isfield(new_lay.param.records.gps,'time_offset') ...
        || new_lay.param.records.gps.time_offset ~= param.records.gps.time_offset
      new_lay.param.records.gps.time_offset = param.records.gps.time_offset;
      save_layer_data_file = bitor(save_layer_data_file,8);
    end
  end
  
  new_lay.gps_time = double(new_lay.gps_time);
  new_lay.lat = double(new_lay.lat);
  new_lay.lon = double(new_lay.lon);
  new_lay.elev = double(new_lay.elev);
  new_lay.id = double(new_lay.id);
  new_lay.twtt = double(new_lay.twtt);
  new_lay.type = uint8(new_lay.type);
  new_lay.quality = uint8(new_lay.quality);
  
  if ~strcmpi(layer_fn_dir,layer_out_fn_dir)
    save_layer_data_file = bitor(save_layer_data_file,4);
  end
  
  if save_layer_data_file
    if bitand(save_layer_data_file,1)
      fprintf('Corrected layer data file being saved.\n');
    end
    if bitand(save_layer_data_file,2)
      warning('Some layers did not match GPS_time field in length. Saving corrected fields.');
    end
    layer_out_fn = fullfile(layer_out_fn_dir,sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    fprintf('  Saving: %s\n', layer_out_fn);
    ct_save(layer_out_fn,'-struct','new_lay') % saving to layerData file
  end
  
end

% layer_organizer did not exist, had bad contents, or was missing some
% layers, so update the file with corrected contents
if save_layer_organizer_file == true
  fprintf('Updating layer organizer file.\n');
  layer_organizer_out_fn = fullfile(ct_filename_out(param,param.layer_update.out_path,''),sprintf('layer_%s.mat',param.day_seg));
  ct_save(layer_organizer_out_fn,'-struct','layer_organizer');
end
