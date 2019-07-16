function layers = opsLoadLayers(param, layer_params)
% layers = opsLoadLayers(param, layer_params)
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
%  .name: string (e.g. 'surface', 'Surface', 'bottom', 'atm', etc)
%  .source: string
%    'records': Loads layer data from records file
%    'echogram': Loads layer data from echogram files
%    'layerdata': Loads layer data from layer data files
%    'lidar': Loads (ATM, AWI, or DTU) lidar data
%    'ops': Loads layer data from Open Polar Server
%  .echogram_source = string containing ct_filename_out argument if using
%    'echogram' source
%    (e.g. 'qlook', 'mvdr', 'CSARP_post/standard')
%  .layerdata_source = string containing ct_filename_out argument if using
%    'echogram' source
%    (e.g. 'layerData', 'CSARP_post/layerData')
%  .lidar_source = string containing 'atm', 'awi', or 'dtu' if using lidar source
%  .existence_check = boolean, default is true and causes an error to be
%    thrown if the layer does not exist. If false, no data points are
%    returned when the layer does not exist and only a warning is given.
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

if ~isfield(param,'debug') || isempty(param.debug)
  param.debug = false;
end

if ~isfield(param,'records')
  param.records = [];
end

for layer_idx = 1:length(layer_params)
  if ~isfield(layer_params,'echogram_source_img') || isempty(layer_params(layer_idx).echogram_source_img)
    % Default is a combined file Data_YYYYMMDD_SS.mat
    layer_params(layer_idx).echogram_source_img = 0;
  end
  if ~isfield(layer_params,'name') || isempty(layer_params(layer_idx).name)
    % Default is layerdata
    warning('No name specified for layer %d, using surface.', layer_idx);
    layer_params(layer_idx).name = 'surface';
  end
  if ~isfield(layer_params,'source') || isempty(layer_params(layer_idx).source)
    % Default is layerdata
    warning('No source specified for layer %d, using layerdata.', layer_idx);
    layer_params(layer_idx).source = 'layerdata';
  end
end

physical_constants;

% Load frames file
load(ct_filename_support(param,'','frames'));

if isfield(layer_params,'frms')
  param.cmd.frms = layer_params.frms;
end

%% Determine which frames need to be processed
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

if any(strcmpi('ops',{layer_params.source}))
  %% Get all the frames for this segment
  opsAuthenticate(param,false);
  sys = ct_output_dir(param.radar_name);
  ops_param = struct('properties',[]);
  ops_param.properties.season = param.season_name;
  ops_param.properties.segment = param.day_seg;
  [status,ops_seg_data] = opsGetSegmentInfo(sys,ops_param);
end
if ~all(strcmpi('ops',{layer_params.source}))
  % All other sources use records file for
  % framing gps time info
  records_fn = ct_filename_support(param,'','records');
  records = load(records_fn,'gps_time','surface','elev','lat','lon');
end

if any(strcmpi('lidar',{layer_params.source}))
  %% Load LIDAR surface
  if any(strcmpi('atm',{layer_params.lidar_source}))
    lidar_fns = get_filenames_atm(param.post.ops.location,param.day_seg(1:8),param.data_support_path);
    
    lidar = read_lidar_atm(lidar_fns);
    
  elseif any(strcmpi('awi',{layer_params.lidar_source}))
    lidar_fns = get_filenames_lidar(param, layer_params.lidar_source, ...
      records.gps_time([1 end]));
    
    lidar_param = struct('time_reference','utc');
    lidar_param.nc_field = {'TIME','LATITUDE','LONGITUDE','ELEVATION','MJD','FOOT_ROUGH'};
    lidar_param.nc_type = {'v','v','v','v','v','v'};
    lidar_param.types = {'sec','lat_deg','lon_deg','elev_m','mjd_18581117','rms'};
    lidar_param.scale = [1 1 1 1 1 1];
    lidar_param.custom_flag = [0 0 0 0 0 1];
    lidar_param.reshape_en = [1 1 1 1 1 1];
    lidar = read_lidar_netcdf(lidar_fns,lidar_param);
  
  elseif any(strcmpi('dtu',{layer_params.lidar_source}))
    lidar_fns = get_filenames_lidar(param, layer_params.lidar_source, ...
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
  lidar.elev = interp1(records.gps_time,records.elev,lidar.gps_time);
  
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

%% Load each of the frames
for frm_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frm_idx);
  
  if param.debug
    fprintf('  Loading %s frame %03d (%d of %d) (%s)\n', param.day_seg, ...
      frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  end
  
  for layer_idx = 1:length(layer_params)
    layer_param = layer_params(layer_idx);
    if ~isfield(layer_param,'existence_check') || isempty(layer_param.existence_check)
      layer_param.existence_check = true;
    end
    
    %% Load LIDAR Data
    if strcmpi(layer_param.source,'lidar')
      if strcmpi(layer_param.name,'surface')
        if frm == length(frames.frame_idxs)
          last_idx = length(records.gps_time);
        else
          last_idx = frames.frame_idxs(frm+1);
        end
        good_idxs = lidar.gps_time >= records.gps_time(frames.frame_idxs(frm)) ...
          & lidar.gps_time <= records.gps_time(last_idx);
        layers(layer_idx).gps_time = cat(2,layers(layer_idx).gps_time, ...
          lidar.gps_time(good_idxs));
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt, ...
          (lidar.elev(good_idxs)-lidar.surface(good_idxs))/(c/2));
        layers(layer_idx).elev = cat(2,layers(layer_idx).elev, ...
          lidar.elev(good_idxs));
        layers(layer_idx).lat = cat(2,layers(layer_idx).lat, ...
          lidar.lat(good_idxs));
        layers(layer_idx).lon = cat(2,layers(layer_idx).lon, ...
          lidar.lon(good_idxs));
        layers(layer_idx).type = cat(2,layers(layer_idx).type, ...
          NaN*zeros(1,sum(good_idxs)));
        layers(layer_idx).quality = cat(2,layers(layer_idx).quality, ...
          NaN*zeros(1,sum(good_idxs)));
      else
        error('Unsupported layer %s for lidar source.', layer_param.name);
      end
    end
    
    %% Load Records Data
    if strcmpi(layer_param.source,'records')
      if strcmpi(layer_param.name,'surface')
        if frm == length(frames.frame_idxs)
          last_idx = length(records.gps_time);
          good_idxs = records.gps_time >= records.gps_time(frames.frame_idxs(frm)) ...
            & records.gps_time <= records.gps_time(last_idx);
        else
          last_idx = frames.frame_idxs(frm+1);
          good_idxs = records.gps_time >= records.gps_time(frames.frame_idxs(frm)) ...
            & records.gps_time < records.gps_time(last_idx);
        end
        layers(layer_idx).gps_time = cat(2,layers(layer_idx).gps_time, ...
          records.gps_time(good_idxs));
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt, ...
          records.surface(good_idxs));
        layers(layer_idx).elev = cat(2,layers(layer_idx).elev, ...
          records.elev(good_idxs));
        layers(layer_idx).lat = cat(2,layers(layer_idx).lat, ...
          records.lat(good_idxs));
        layers(layer_idx).lon = cat(2,layers(layer_idx).lon, ...
          records.lon(good_idxs));
        layers(layer_idx).type = cat(2,layers(layer_idx).type, ...
          NaN*zeros(1,sum(good_idxs)));
        layers(layer_idx).quality = cat(2,layers(layer_idx).quality, ...
          NaN*zeros(1,sum(good_idxs)));
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
      
      % Load data
      warning off;
      mdata = load(data_fn,echogram_layer_name,'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom', ...
        'Elevation_Correction','Truncate_Bins','Time');
      warning on;
      mdata = uncompress_echogram(mdata);
      
      % Remove data that is not contained within frame boundaries
      frms_mask = logical(zeros(size(mdata.GPS_time)));
      if frm < length(frames.frame_idxs)
        frms_mask(mdata.GPS_time >= records.gps_time(frames.frame_idxs(frm))...
          & mdata.GPS_time < records.gps_time(frames.frame_idxs(frm+1))) = true;
      else
        frms_mask(mdata.GPS_time >= records.gps_time(frames.frame_idxs(frm))...
          & mdata.GPS_time < records.gps_time(end)) = true;
      end
      
      layers(layer_idx).gps_time = cat(2,layers(layer_idx).gps_time, ...
        mdata.GPS_time(frms_mask));
      layers(layer_idx).elev = cat(2,layers(layer_idx).elev, ...
        mdata.Elevation(frms_mask));
      layers(layer_idx).lat = cat(2,layers(layer_idx).lat, ...
        mdata.Latitude(frms_mask));
      layers(layer_idx).lon = cat(2,layers(layer_idx).lon, ...
        mdata.Longitude(frms_mask));
      layers(layer_idx).type = cat(2,layers(layer_idx).type, ...
        NaN*zeros(size(mdata.GPS_time(frms_mask))));
      layers(layer_idx).quality = cat(2,layers(layer_idx).quality, ...
        NaN*zeros(size(mdata.GPS_time(frms_mask))));
      
      if isfield(mdata,echogram_layer_name)
        % Layer exists in the file, concatenate it on to twtt
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt, ...
          mdata.(echogram_layer_name)(frms_mask));
      else
        % Layer does not exist in the file
        if layer_param.existence_check
          error('Unsupported layer %s for echogram source.', layer_param.name);
        else
          warning('Unsupported layer %s for echogram source.', layer_param.name);
        end
        
        % Concatenate NaN on to twtt
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt,...
          NaN*zeros(size(mdata.GPS_time(frms_mask))));
      end
      
    end
    
    %% Load layerData Data
    if strcmpi(layer_param.source,'layerdata')
      % 1. Open the specific layer data file
      if ~isfield(layer_param,'layerdata_source') || isempty(layer_param.layerdata_source)
        layer_param.layerdata_source = 'layerData';
      end
      layer_fn = fullfile(ct_filename_out(param,layer_param.layerdata_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      if ~exist(layer_fn,'file')
        if isfield(layer_param,'echogram_source') && ~isempty(layer_param.echogram_source)
          % Create the layer data file from the echogram source
          warning('Layer data file does not exist. Since echogram_source specified, creating layer data:\n%s.', layer_fn);
          
          % Load the echogram
          echogram_fn_dir = ct_filename_out(param,layer_param.echogram_source);
          if layer_param.echogram_source_img == 0
            echogram_fn = fullfile(echogram_fn_dir, sprintf('Data_%s_%03d.mat', ...
              param.day_seg, frm));
          else
            echogram_fn = fullfile(echogram_fn_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
              layer_param.echogram_source_img, param.day_seg, frm));
          end
          if ~exist(echogram_fn,'file')
            if layer_param.existence_check
              error('Echogram file %s does not exist', echogram_fn);
            else
              warning('Echogram file %s does not exist', echogram_fn);
              continue;
            end
          end
          % Load and uncompress echogram file
          warning off;
          lay = load(echogram_fn,'GPS_time','Latitude','Longitude','Elevation','Surface','Bottom', ...
            'Elevation_Correction','Truncate_Bins','Time');
          warning on;
          lay = uncompress_echogram(lay);
          
          % In case an old echogram file which has overlap, we remove this
          % overlap region.
          if frm == 1
            start_gps_time = -inf;
          else
            start_gps_time = records.gps_time(frames.frame_idxs(frm));
          end
          if frm == length(frames.frame_idxs)
            stop_gps_time = inf;
          else
            stop_gps_time = records.gps_time(frames.frame_idxs(frm+1));
          end
          good_mask = lay.GPS_time >= start_gps_time & lay.GPS_time < stop_gps_time;
          lay.GPS_time = lay.GPS_time(good_mask);
          lay.Latitude = lay.Latitude(good_mask);
          lay.Longitude = lay.Longitude(good_mask);
          lay.Elevation = lay.Elevation(good_mask);
          lay.Surface = lay.Surface(good_mask);
          if isfield(lay,'Bottom')
            lay.Bottom = lay.Bottom(good_mask);
          else
            lay.Bottom = nan(size(lay.GPS_time));
          end
          Nx = length(lay.GPS_time);

          % Create the layer structure
          for lay_idx = 1:2
            % Manually picked points
            %  inf/nan: no pick
            %  finite: propagation time to target
            lay.layerData{lay_idx}.value{1}.data ...
              = inf*ones(1,Nx);
            % Automatically generated points
            %  inf/nan: no pick
            %  finite: propagation time to target
            if lay_idx == 1 && isfield(lay,'Surface')
              lay.layerData{lay_idx}.value{2}.data = lay.Surface;
            elseif lay_idx == 2 && isfield(lay,'Bottom')
              lay.layerData{lay_idx}.value{2}.data = lay.Bottom;
            else
              lay.layerData{lay_idx}.value{2}.data ...
                = inf*ones(1,Nx);
            end
            % Quality control level
            %  1: good
            %  2: moderate
            %  3: derived
            lay.layerData{lay_idx}.quality ...
              = ones(1,Nx);
          end

          layer_fn_dir = fileparts(layer_fn);
          if ~exist(layer_fn_dir,'dir')
            mkdir(layer_fn_dir);
          end
          save(layer_fn,'-v6','-struct','lay','layerData','Latitude','Longitude','Elevation','GPS_time');

        elseif layer_param.existence_check
          error('Layer file %s does not exist', layer_fn);
        else
          warning('Layer file %s does not exist', layer_fn);
          continue;
        end
      end
      % Load the layerData file
      lay = load(layer_fn);

      % Remove data that is not contained within frame boundaries
      frms_mask = logical(zeros(size(lay.GPS_time)));
      if frm < length(frames.frame_idxs)
        frms_mask(lay.GPS_time >= records.gps_time(frames.frame_idxs(frm))...
          & lay.GPS_time < records.gps_time(frames.frame_idxs(frm+1))) = true;
      else
        frms_mask(lay.GPS_time >= records.gps_time(frames.frame_idxs(frm))...
          & lay.GPS_time < records.gps_time(end)) = true;
      end

      found = false;
      if strcmpi(layer_param.name,'surface')
        lay_idx = 1;
        found = true;
      elseif strcmpi(layer_param.name,'bottom')
        lay_idx = 2;
        found = true;
      else
        for lay_idx = 1:length(lay.layerData)
          if isfield(lay.layerData{lay_idx},'name') ...
              && strcmpi(lay.layerData{lay_idx}.name,layer_param.name)
            found = true;
            break;
          end
        end
        if ~found
          if layer_param.existence_check
            error('Reference layer %s not found\n', layer_param.name);
          else
            warning('Reference layer %s not found\n', layer_param.name);
          end
        end
      end
      layers(layer_idx).gps_time = cat(2,layers(layer_idx).gps_time, ...
        lay.GPS_time(frms_mask));
      if ~found
        % Fill with NaN since layer does not exist
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt, ...
          NaN*zeros(size(lay.GPS_time(frms_mask))));
      else
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt, ...
          lay.layerData{lay_idx}.value{2}.data(frms_mask));
      end
      layers(layer_idx).elev = cat(2,layers(layer_idx).elev, ...
        lay.Elevation(frms_mask));
      layers(layer_idx).lat = cat(2,layers(layer_idx).lat, ...
        lay.Latitude(frms_mask));
      layers(layer_idx).lon = cat(2,layers(layer_idx).lon, ...
        lay.Longitude(frms_mask));
      layer_type = 2*ones(size(lay.GPS_time(frms_mask)));
      layer_type(isfinite(lay.layerData{lay_idx}.value{1}.data(frms_mask))) = 1;
      layers(layer_idx).type = cat(2,layers(layer_idx).type, ...
        layer_type);
      if ~found
        % Fill with 1's since layer does not exist
        layers(layer_idx).quality = cat(2,layers(layer_idx).quality, ...
          ones(size(lay.GPS_time(frms_mask))));
      else
        layers(layer_idx).quality = cat(2,layers(layer_idx).quality, ...
          lay.layerData{lay_idx}.quality(frms_mask));
      end
    end
    
    %% Load OPS Data
    if strcmpi(layer_param.source,'ops')
      start_gps = ops_seg_data.properties.start_gps_time(frm);
      stop_gps = ops_seg_data.properties.stop_gps_time(frm);
      
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
        
        layers(layer_idx).gps_time = cat(2,layers(layer_idx).gps_time, ...
          data.properties.gps_time);
        layers(layer_idx).twtt = cat(2,layers(layer_idx).twtt, ...
          data.properties.twtt);
        layers(layer_idx).elev = cat(2,layers(layer_idx).elev, ...
          data.properties.elev);
        layers(layer_idx).lat = cat(2,layers(layer_idx).lat, ...
          data.properties.lat);
        layers(layer_idx).lon = cat(2,layers(layer_idx).lon, ...
          data.properties.lon);
        layers(layer_idx).type = cat(2,layers(layer_idx).type, ...
          data.properties.type);
        layers(layer_idx).quality = cat(2,layers(layer_idx).quality, ...
          data.properties.quality);
        layers(layer_idx).point_path_id = cat(2,layers(layer_idx).point_path_id, ...
          data.properties.point_path_id);
      end
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
        start_gps_time = records.gps_time(frames.frame_idxs(frm));
      end
      if frm == length(frames.frame_idxs)
        stop_gps_time = inf;
      else
        stop_gps_time = records.gps_time(frames.frame_idxs(frm+1));
      end
      layers(layer_idx).frm(layers(layer_idx).gps_time >= start_gps_time ...
        & layers(layer_idx).gps_time < stop_gps_time) = frm;
    end
  end
end

