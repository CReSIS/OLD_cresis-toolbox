function [layers,layer_params] = opsLoadLayers(param, layer_params)
% [layers,layer_params] = opsLoadLayers(param, layer_params)
%
% This function loads layer data for specified frames from a single segment.
% The main differences compared to opsGetLayerPoints are:
% * The source of the layer data can be records, echogram files, layerData files,
%   ATM or AWI Lidar, OPS.
% * Controlled from param spreadsheet
% * Use opsInsertLayerFromGrid and opsInsertLayerFromPointCloud to compare
%   layers to grids and point clouds
% * Use runOpsCopyLayers to copy layers from one radar to another
%
% param = param spreadsheet structure
% layer_params = N element struct array indicating which layers are to be loaded
%   and which source to use for each layer
%  .name: string (e.g. 'surface', 'Surface', 'bottom', 'atm', etc), default
%    is "surface"
%  .source: string
%    'records': Loads layer data from records file
%    'echogram': Loads layer data from echogram files
%    'layerdata': Loads layer data from layer data files (default)
%    'lidar': Loads (ATM, AWI, or DTU) lidar data
%    'ops': Loads layer data from Open Polar Server
%  .echogram_source = string containing ct_filename_out argument if using
%    'echogram' source
%    (e.g. 'qlook', 'mvdr', 'CSARP_post/standard')
%  .layerdata_source = string containing ct_filename_out argument of where
%    the layerdata is (default is "layer")
%    (e.g. 'layer', 'CSARP_post/layer')
%  .lidar_source = string containing 'atm', 'awi', or 'dtu' if using lidar source
%  .existence_check = boolean, default is true and causes an error to be
%    thrown if the layer does not exist. If false, no data points are
%    returned when the layer does not exist and only a warning is given.
%  .fix_broken_layerdata: logical, default is false and causes any layerdata
%    errors to be fixed while loading by setting the layer for the
%    particular data frame that has the error to the default values (NaN
%    for twtt, quality 1, and type 2). If true, an error is thrown if
%    layerdata errors exist.
%  .debug = default is false
%  .eval: Optional structure for performing operations on the layer
%    .cmd: Command string that will be passed to eval
%    .$(custom): Custom fields to be accessed in cmd as es.$(custom)
%     Variables available are:
%       physical_constants
%       "gps_time" (sec)
%       "along_track" (m)
%       "lat" (deg)
%       "lon" (deg)
%       "elev" (m)
%       "source" (twtt in sec)
%       "es" (this eval structure passed in by the user)
%     The cmd string should generally update "source" variable. For example:
%        '[B,A] = butter(0.1,2); source = filtfilt(B,A,source);' % Filter
%        'source = source + 0.1;' % Apply a twtt shift
%        'source = source*2;' % Surface multiple
%  .frms: This field overrides the param.cmd.frms field, but must be the
%    same for all elements of the layer_params struct array since only the
%    first element will be used.
%
% layers: N element struct array with layer information
%  .gps_time
%  .lat
%  .lon
%  .elev
%  .quality (not available for all sources, set to NaN if not available)
%  .type (not available for all sources, set to NaN if not available)
%  .twtt
%  .point_path_id database key (only filled if source is ops)
%  .frm: numeric vector of frame ids (1 to 999)
%
% Authors: John Paden
%
% See also: runOpsLoadLayers.m

%% Input checks

if ~isfield(param,'records')
  param.records = [];
end

% Check layers defined by regular expression regexp
layer_idx = 1;
while layer_idx <= length(layer_params)
  layer_names = {};
  if isfield(layer_params(layer_idx),'name') && iscell(layer_params(layer_idx).name)
    layer_names = layer_params(layer_idx).name;
  elseif ~isfield(layer_params,'name') || isempty(layer_params(layer_idx).name)
    % Name is not specified, so check regular expression field
    if ~isfield(layer_params,'regexp') || isempty(layer_params(layer_idx).regexp)
      % Default is surface
      warning('No name specified for layer %d, defaulting to use layer "surface.', layer_idx);
      layer_params(layer_idx).name = 'surface';
    else
      error('regexp not supported yet. Just need to set layer_names to the list of layer names that match.');
    end
  end
  if ~isempty(layer_names)
    layer_params = layer_params([1:layer_idx-1 layer_idx*ones(1,length(layer_names)) layer_idx+1:end]);
    for new_idx = 1:length(layer_names)
      layer_params(layer_idx+new_idx-1).name = layer_names{new_idx};
      layer_params(layer_idx+new_idx-1).regexp = '';
    end
    layer_idx = layer_idx+new_idx-1;
    
  end
  layer_idx = layer_idx + 1;
end

% Check layer parameters
ops_en = false;
records_en = false;
echogram_en = false;
lidar_layer_idx = [];
layerdata_sources = {};
for layer_idx = 1:length(layer_params)
  if ~isfield(layer_params,'echogram_source_img') || isempty(layer_params(layer_idx).echogram_source_img)
    % Default is a combined file Data_YYYYMMDD_SS.mat
    layer_params(layer_idx).echogram_source_img = 0;
  end
  if ~isfield(layer_params,'source') || isempty(layer_params(layer_idx).source)
    % Default is layerdata
    warning('No source specified for layer %d, defaulting to use layerdata.', layer_idx);
    layer_params(layer_idx).source = 'layerdata';
  end
  if ~isfield(layer_params,'layerdata_source') || isempty(layer_params(layer_idx).layerdata_source)
    % Default is layerData
    layer_params(layer_idx).layerdata_source = 'layer';
  end
  if ~isfield(layer_params,'read_only') || isempty(layer_params(layer_idx).read_only)
    % Default is layerData
    layer_params(layer_idx).read_only = true;
  end
  if ~isfield(layer_params,'lever_arm_en') || isempty(layer_params(layer_idx).lever_arm_en)
    layer_params(layer_idx).lever_arm_en = false;
  end
  if ~isfield(layer_params(layer_idx),'existence_check') || isempty(layer_params(layer_idx).existence_check)
    layer_params(layer_idx).existence_check = true;
  end
  switch lower(layer_params(layer_idx).source)
    case 'ops'
      ops_en = true;
    case 'layerdata'
      if ~any(strcmp(layer_params(layer_idx).layerdata_source,layerdata_sources))
        layerdata_sources{end+1} = layer_params(layer_idx).layerdata_source;
      end
    case 'records'
      records_en = true;
    case 'echogram'
      echogram_en = true;
    case 'lidar'
      lidar_layer_idx(end+1) = layer_idx;
  end
end

% Nothing to do, return
if isempty(layer_params)
  layers = [];
  return;
end

physical_constants;

%% Get all the frames for this segment
if ~isempty(layerdata_sources) || records_en || echogram_en || ~isempty(lidar_layer_idx)
  % Load frames file
  frames = frames_load(param);

  if records_en || ~isempty(lidar_layer_idx) || ~isfield(frames,'frame_idxs')
    records = records_load(param);
  end
  
  % Determine which frames need to be processed
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
end

if ops_en
  opsAuthenticate(param,false);
  sys = ct_output_dir(param.radar_name);
  ops_param = struct('properties',[]);
  ops_param.properties.season = param.season_name;
  ops_param.properties.segment = param.day_seg;
  [status,ops_seg_data] = opsGetSegmentInfo(sys,ops_param);
end

if ~isempty(lidar_layer_idx)
  %% Load LIDAR surface
  lidar_source = {layer_params.lidar_source};
  lidar_source = lidar_source(find(~cellfun(@isempty,lidar_source)));
  if length(unique(lidar_source)) ~= 1
    error('One and only one lidar source must be specified to opsLoadLayers when one of the sources to load is ''lidar''.');
  end
  if any(strcmpi('atm',{layer_params.lidar_source}))
    lidar_fns = get_filenames_atm(param.post.ops.location,param.day_seg(1:8),param.data_support_path);
    
    lidar = read_lidar_atm(lidar_fns);
    
  elseif any(strcmpi('awi',{layer_params.lidar_source}))
    lidar_fns = get_filenames_lidar(param, 'awi', ...
      records.gps_time([1 end]));
    
    lidar_param = struct('time_reference','utc');
    lidar_param.nc_field = {'TIME','LATITUDE','LONGITUDE','ELEVATION','MJD','FOOT_ROUGH'};
    lidar_param.nc_type = {'v','v','v','v','v','v'};
    lidar_param.types = {'sec','lat_deg','lon_deg','elev_m','mjd_18581117','rms'};
    lidar_param.scale = [1 1 1 1 1 1];
    lidar_param.custom_flag = [0 0 0 0 0 1];
    lidar_param.reshape_en = [1 1 1 1 1 1];
    lidar = read_lidar_netcdf(lidar_fns,lidar_param);
    
  elseif any(strcmpi('awi_L2B',{layer_params.lidar_source}))
    lidar_fns = get_filenames_lidar(param, 'awi_L2B', ...
      records.gps_time([1 end]));
    
    [year month day] = datevec(epoch_to_datenum(records.gps_time(1)));
    lidar_param = struct('time_reference','utc');
    lidar_param.nc_field = {'time','latitude','longitude','l1b_elevation'};
    lidar_param.nc_type = {'v','v','v','v','v','v'};
    lidar_param.types = {'sec','lat_deg','lon_deg','elev_m'};
    lidar_param.scale = [1 1 1 1];
    lidar_param.scale = [1 1 1 1];
    lidar_param.custom_flag = [0 0 0 0];
    lidar_param.reshape_en = [1 1 1 1];
    lidar_param.year = year;
    lidar_param.month = month;
    lidar_param.day = day;
    lidar = read_lidar_netcdf(lidar_fns,lidar_param);
    
  elseif any(strcmpi('dtu',{layer_params.lidar_source}))
    lidar_fns = get_filenames_lidar(param, 'dtu', ...
      records.gps_time([1 end]));
    
    lidar = read_lidar_dtu(lidar_fns,param);
    
  else
    error('Invalid LIDAR source %s', layer_params.lidar_source);
  end
  
  % Remove NAN's from LIDAR Data
  good_lidar_idxs = ~isnan(lidar.gps_time);
  lidar.gps_time = lidar.gps_time(good_lidar_idxs);
  lidar.surface = lidar.surface(good_lidar_idxs);
  lidar.lat = lidar.lat(good_lidar_idxs);
  lidar.lon = lidar.lon(good_lidar_idxs);
  if ~isempty(isnan(lidar.surface))
    good_lidar_idxs = ~isnan(lidar.surface);
    lidar.gps_time = lidar.gps_time(good_lidar_idxs);
    lidar.surface = lidar.surface(good_lidar_idxs);
    lidar.lat = lidar.lat(good_lidar_idxs);
    lidar.lon = lidar.lon(good_lidar_idxs);
  end
  
  % Find reference trajectory
  if ~layer_params(lidar_layer_idx).lever_arm_en
    % Just use the records.elev for the radar phase center elevation. This
    % will likely cause an error because the records elevation field is
    % equal to the GPS data file elevation field which is often the GPS or
    % IMU elevation and not the radar phase center elevation.
    lidar.elev = interp1(records.gps_time,records.elev,lidar.gps_time);
  else
    % Create reference trajectory (rx_path == 0, tx_weights = []). Update
    % the records field with this information.
    trajectory_param = struct('gps_source',records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
      'tx_weights', [], 'lever_arm_fh', param.radar.lever_arm_fh);
    records = trajectory_with_leverarm(records,trajectory_param);
    
    % Find the closest lidar point for each point along the reference
    % trajectory. Set the lidar.gps_time and lidar.elev fields to match
    % this reference point's gps_time and elev so that interpolation later
    % with the gps_time fields will work properly.
    % Convert coordinates to ECEF
    [lidar_x,lidar_y,lidar_z] = geodetic2ecef(lidar.lat/180*pi,lidar.lon/180*pi,zeros(size(lidar.lat)),WGS84.ellipsoid);
    [records_x,records_y,records_z] = geodetic2ecef(records.lat/180*pi,records.lon/180*pi,zeros(size(records.lat)),WGS84.ellipsoid);
    start_idx = 1;
    stop_idx = 1;
    lidar.elev = zeros(size(lidar.lat));
    for lidar_idx = 1:length(lidar.lat)
      while start_idx < length(records.gps_time) && records.gps_time(start_idx) < lidar.gps_time(lidar_idx)-10
        start_idx = start_idx + 1;
      end
      while stop_idx < length(records.gps_time) && records.gps_time(stop_idx) < lidar.gps_time(lidar_idx)+10
        stop_idx = stop_idx + 1;
      end
      if abs(records.gps_time(start_idx) - lidar.gps_time(lidar_idx)) < 10
        recs = start_idx:stop_idx;
        dist = abs(lidar_x(lidar_idx)-records_x(recs)).^2 + abs(lidar_y(lidar_idx)-records_y(recs)).^2 + abs(lidar_z(lidar_idx)-records_z(recs)).^2;
        [min_dist,min_idx] = min(dist);
        lidar.gps_time(lidar_idx) = records.gps_time(recs(min_idx));
        lidar.elev(lidar_idx) = records.elev(recs(min_idx));
      else
        lidar.gps_time(lidar_idx) = NaN;
        lidar.elev(lidar_idx) = NaN;
      end
    end
    % Remove NAN's from LIDAR Data
    good_lidar_idxs = ~isnan(lidar.gps_time);
    lidar.gps_time = lidar.gps_time(good_lidar_idxs);
    lidar.surface = lidar.surface(good_lidar_idxs);
    lidar.lat = lidar.lat(good_lidar_idxs);
    lidar.lon = lidar.lon(good_lidar_idxs);
    lidar.elev = lidar.elev(good_lidar_idxs);
  end
end


%% Initialize Outputs
for layer_idx = 1:length(layer_params)
  layers(layer_idx).gps_time = [];
  layers(layer_idx).twtt = [];
  layers(layer_idx).elev = [];
  layers(layer_idx).lat = [];
  layers(layer_idx).lon = [];
  layers(layer_idx).type = [];
  layers(layer_idx).quality = [];
  layers(layer_idx).point_path_id = [];
end

%% Load each of the frames (lidar, echogram, records file layer sources)
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  
  if 0
    fprintf('  Loading %s frame %03d (%d of %d) (%s)\n', param.day_seg, ...
      frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  end
  
  for layer_idx = 1:length(layer_params)
    layer_param = layer_params(layer_idx);
    
    %% Load LIDAR Data
    if strcmpi(layer_param.source,'lidar')
      if strcmpi(layer_param.name,'surface')
        good_idxs = find(lidar.gps_time >= frames.gps_time(frm) ...
          & lidar.gps_time < frames.gps_time(frm+1));
        Nx = length(good_idxs);
        layers(layer_idx).gps_time(end+(1:Nx)) = lidar.gps_time(good_idxs);
        layers(layer_idx).twtt(end+(1:Nx)) = (lidar.elev(good_idxs)-lidar.surface(good_idxs))/(c/2);
        layers(layer_idx).elev(end+(1:Nx)) = lidar.elev(good_idxs);
        layers(layer_idx).lat(end+(1:Nx)) = lidar.lat(good_idxs);
        layers(layer_idx).lon(end+(1:Nx)) = lidar.lon(good_idxs);
        layers(layer_idx).type(end+(1:Nx)) = 2*ones(1,Nx);
        layers(layer_idx).quality(end+(1:Nx)) = ones(1,Nx);
      else
        error('Unsupported layer %s for lidar source.', layer_param.name);
      end
    end
    
    %% Load Records Data
    if strcmpi(layer_param.source,'records')
      if strcmpi(layer_param.name,'surface')
        good_idxs = frames.frame_idxs(frm):last_idx;
        Nx = length(good_idxs);
        layers(layer_idx).gps_time(end+(1:Nx)) = records.gps_time(good_idxs);
        layers(layer_idx).twtt(end+(1:Nx)) = records.surface(good_idxs);
        layers(layer_idx).elev(end+(1:Nx)) = records.elev(good_idxs);
        layers(layer_idx).lat(end+(1:Nx)) = records.lat(good_idxs);
        layers(layer_idx).lon(end+(1:Nx)) = records.lon(good_idxs);
        layers(layer_idx).type(end+(1:Nx)) = 2*ones(1,Nx);
        layers(layer_idx).quality(end+(1:Nx)) = ones(1,Nx);
      else
        error('Unsupported layer %s for records source.', layer_param.name);
      end
    end
    
    %% Load Echogram Data
    if strcmpi(layer_param.source,'echogram')
      if layer_param.echogram_source_img == 0
        data_fn = fullfile(ct_filename_out(param,layer_param.echogram_source,''), ...
          sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      else
        data_fn = fullfile(ct_filename_out(param,layer_param.echogram_source,''), ...
          sprintf('Data_img_%02d_%s_%03d.mat', layer_param.echogram_source_img, param.day_seg, frm));
      end
      
      if ~exist(data_fn,'file')
        if layer_param.existence_check
          error('Echogram file %s does not exist', data_fn);
        else
          warning('Echogram file %s does not exist', data_fn);
          continue;
        end
      end
      
      % Convert layer name to echogram name
      if strcmpi(layer_params(layer_idx).name,'surface')
        echogram_layer_name = 'Surface';
      elseif strcmpi(layer_params(layer_idx).name,'bottom')
        echogram_layer_name = 'Bottom';
      else
        echogram_layer_name = layer_params(layer_idx).name;
      end
      
      % Load layer data
      warning off;
      mdata = load(data_fn,echogram_layer_name,'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom', ...
        'Elevation_Correction','Truncate_Bins','Time');
      warning on;
      mdata = uncompress_echogram(mdata);
      
      % Remove data that is not contained within frame boundaries
      good_idxs = find(mdata.GPS_time >= frames.gps_time(frm) ...
        & mdata.GPS_time < frames.gps_time(frm+1));
      Nx = length(good_idxs);
      layers(layer_idx).gps_time(end+(1:Nx)) = mdata.GPS_time(good_idxs);
      layers(layer_idx).elev(end+(1:Nx)) = mdata.Elevation(good_idxs);
      layers(layer_idx).lat(end+(1:Nx)) = mdata.Latitude(good_idxs);
      layers(layer_idx).lon(end+(1:Nx)) = mdata.Longitude(good_idxs);
      layers(layer_idx).type(end+(1:Nx)) = 2*ones(1,Nx);
      layers(layer_idx).quality(end+(1:Nx)) = ones(1,Nx);
      
      if isfield(mdata,echogram_layer_name)
        % Layer exists in the file, concatenate it on to twtt
        layers(layer_idx).twtt(end+(1:Nx)) = mdata.(echogram_layer_name)(good_idxs);
      else
        % Layer does not exist in the file
        if layer_param.existence_check
          error('Unsupported layer %s for echogram source.', layer_param.name);
        else
          warning('Unsupported layer %s for echogram source.', layer_param.name);
        end
        
        % Concatenate NaN on to twtt
        layers(layer_idx).twtt(end+(1:Nx)) = nan(1,Nx);
      end
      
    end
  end
end

%% Load layerdata
for layerdata_source_idx = 1:length(layerdata_sources)
  layerdata_source = layerdata_sources{layerdata_source_idx};
  
  tmp_layers = layerdata(param, layerdata_source);
  
  for layer_idx = 1:length(layer_params)
    layer_param = layer_params(layer_idx);
    if ~strcmpi(layer_param.source,'layerdata') || ~strcmpi(layer_param.layerdata_source,layerdata_source)
      continue;
    end
    
    layers(layer_idx).gps_time = tmp_layers.gps_time(param.cmd.frms);
    layers(layer_idx).lat = tmp_layers.lat(param.cmd.frms);
    layers(layer_idx).lon = tmp_layers.lon(param.cmd.frms);
    layers(layer_idx).elev = tmp_layers.elev(param.cmd.frms);
    id = tmp_layers.get_id(layer_param.name);
    if isempty(id)
      if layer_param.existence_check
        error('Layer %s does not exist in %s.', layer_param.name, tmp_layers.layer_organizer_fn());
      end
      layer_organizer = [];
      layer_organizer.lyr_name = {layer_param.name};
      if isfield(layer_param,'group_name')
        layer_organizer.lyr_group_name = {layer_param.group_name};
      end
      tmp_layers.insert_layers(layer_organizer);
    end
    [layers(layer_idx).twtt,layers(layer_idx).quality,layers(layer_idx).type] = tmp_layers.get_layer(param.cmd.frms,layer_param.name);
    if ~layer_param.read_only
      tmp_layers.save();
    end
  end
end

%% Load OPS Data
if ops_en
  for layer_idx = 1:length(layer_params)
    layer_param = layer_params(layer_idx);
    
    if isempty(param.cmd.frms)
      start_gps = ops_seg_data.properties.start_gps_time(1);
      stop_gps = ops_seg_data.properties.stop_gps_time(end);
    else
      start_gps = ops_seg_data.properties.start_gps_time(min(param.cmd.frms));
      stop_gps = ops_seg_data.properties.stop_gps_time(max(param.cmd.frms));
    end
    
    found = true;
    if ~layer_param.existence_check
      % If the layer does not exist, we need to determine this before
      % we call opsGetLayerPoints, otherwise we will get an error in
      % that function.
      [status,data] = opsGetLayers(sys);
      if ~any(strcmpi(data.properties.lyr_name,layer_param.name))
        found = false;
        warning('Layer %s does not exist in OPS.', layer_param.name);
      end
    end
    
    if found
      ops_param = struct('properties',[]);
      ops_param.tmp_path = param.tmp_path;
      ops_param.properties.location = param.post.ops.location;
      ops_param.properties.season = param.season_name;
      ops_param.properties.segment = param.day_seg;
      ops_param.properties.start_gps_time = start_gps;
      ops_param.properties.stop_gps_time = stop_gps;
      ops_param.properties.lyr_name = layer_param.name;
      ops_param.properties.return_geom = 'geog';
      [status,data] = opsGetLayerPoints(sys,ops_param);
      
      [data.properties.gps_time,sort_idxs] = sort(data.properties.gps_time);
      data.properties.twtt = data.properties.twtt(sort_idxs);
      data.properties.elev = data.properties.elev(sort_idxs);
      data.properties.lat = data.properties.lat(sort_idxs);
      data.properties.lon = data.properties.lon(sort_idxs);
      data.properties.type = data.properties.type(sort_idxs);
      data.properties.quality = data.properties.quality(sort_idxs);
      data.properties.point_path_id = data.properties.point_path_id(sort_idxs);
      
      layers(layer_idx).gps_time = data.properties.gps_time;
      layers(layer_idx).twtt = data.properties.twtt;
      layers(layer_idx).elev = data.properties.elev;
      layers(layer_idx).lat = data.properties.lat;
      layers(layer_idx).lon = data.properties.lon;
      layers(layer_idx).type = data.properties.type;
      layers(layer_idx).quality = data.properties.quality;
      layers(layer_idx).point_path_id = data.properties.point_path_id;
    end
  end
end

%% Execute the eval statement on the layer data if supplied
for layer_idx = 1:length(layer_params)
  layer_param = layer_params(layer_idx);
  if isfield(layer_param,'eval') && ~isempty(layer_param.eval)
    s = layers(layer_idx).twtt;
    time = layers(layer_idx).gps_time;
    lat = layers(layer_idx).lat;
    lon = layers(layer_idx).lon;
    elev = layers(layer_idx).elev;
    at = geodetic_to_along_track(lat,lon,elev);
    es = layer_param.eval;
    eval(layer_param.eval.cmd);
    layers(layer_idx).twtt = s;
  end
end

%% Add the frame field
for layer_idx = 1:length(layer_params)
  if strcmp(layer_param.source,'ops')
    layers(layer_idx).frm = zeros(size(layers(layer_idx).gps_time));
    for frm = 1:length(ops_seg_data.properties.start_gps_time)
      if frm == 1
        start_gps_time = 0;
      else
        start_gps_time = ops_seg_data.properties.start_gps_time(frm);
      end
      if frm == length(ops_seg_data.properties.start_gps_time)
        stop_gps_time = inf;
      else
        stop_gps_time = ops_seg_data.properties.start_gps_time(frm+1);
      end
      layers(layer_idx).frm(layers(layer_idx).gps_time >= start_gps_time ...
        & layers(layer_idx).gps_time < stop_gps_time) = frm;
    end
  else
    layers(layer_idx).frm = zeros(size(layers(layer_idx).gps_time));
    for frm = 1:length(frames.frame_idxs)
      if frm == 1
        start_gps_time = 0;
      else
        start_gps_time = frames.gps_time(frm);
      end
      if frm == length(frames.frame_idxs)
        stop_gps_time = inf;
      else
        stop_gps_time = frames.gps_time(frm+1);
      end
      layers(layer_idx).frm(layers(layer_idx).gps_time >= start_gps_time ...
        & layers(layer_idx).gps_time < stop_gps_time) = frm;
    end
  end
end

