function layers = opsCopyLayers(param,copy_param)
% layers = opsCopyLayers(param,copy_param)
%
% Copies contents from one layer into another layer. The copy method can
% be merge, overwrite, or fill gaps. If the destination layer does not
% exist, it will be created. The source and destination for the data can be
% OPS, layerData, echogram, or records. Additionally, the source can be
% lidar (ATM) or a custom GPS-time/two-way-travel-time.
%
% An operation on the source data using the eval field can be added. More
% complicated operations can be done using opsLoadLayers, applying the
% complex operation on that result and then using the custom input.
%
% Input:
%   param: Parameter structure from read_param_xls parameter spreadsheet
%   copy_param: Structure which controls copying process
%     .layer_source: structure specifying the source layer
%       .name: string (e.g. 'surface', 'Surface', 'bottom', 'atm', etc)
%       .source: string (e.g. 'records', 'echogram', 'layerdata', 'lidar',
%         'custom', or 'ops')
%       .echogram_source: used only with "echogram" source, string
%         containing file path argument to ct_filename_out.m
%         (e.g. 'qlook', 'mvdr', 'CSARP_post/standard')
%       .echogram_source_img: used only with "echogram" source, scalar
%         integer specifying the image. Zero is default.
%           0: Data_YYYYMMDD_FF.mat
%           II: Data_img_II_YYYYMMDD_FF.mat
%       .layerdata_source: used only with layerdata source, string
%         containing file path argument to ct_filename_out.m
%         (e.g. 'layerData', 'CSARP_post/layerData')
%       .existence_check: used only with ops source, set to false to allow
%         layers that do not exist (default is true)
%       .gps_time: used only with custom source, Nx by 1 vector,
%         GPS time (ANSI-C standard: seconds since Jan 1, 1970)
%       .twtt: used only with custom source, Nx by 1 vector corresponding
%         to gps_time field containing two way travel time to layer
%       .type: used only with custom source, Nx by 1 vector corresponding
%         to gps_time field containing type (1=manual, 2=auto), 2 is
%         default when not specified or NaN
%       .quality: used only with custom source, Nx by 1 vector corresponding
%         to gps_time field containing quality (1=good, 2=moderate,
%         3=derived/poor), 1 is default when not specified or NaN
%     .layer_dest = structure specifying the destination layer
%       .name: string (e.g. 'surface', 'Surface', 'bottom', 'atm', etc)
%       .source: string (e.g. 'records', 'echogram', 'layerdata', or 'ops')
%       .echogram_source: used only with echogram source, string
%         containing file path argument to ct_filename_out.m
%         (e.g. 'qlook', 'mvdr', 'CSARP_post/standard')
%       .layerdata_source: used only with layerdata source, string
%         containing file path argument to ct_filename_out.m
%         (e.g. 'layerData', 'CSARP_post/layerData')
%       .existence_check: used only with ops source, set to false to allow
%         layers that do not exist (default is true). If false, the layer
%         will be created if it does not exist. If true, an error will be
%         thrown if the layer does not exist.
%       .group: used only with ops source and only needed when the layer
%         does not already exist. Should be a string containing the group
%         name that the layer should be added to. Leave blank to use the
%         standard group.
%       .description: used only with ops source and only needed when the layer
%         does not already exist. Should be a string containing a
%         description of the layer contents.
%     .eval: Optional structure for performing operations on the source
%       before it is written to the destination.
%       .cmd: Command string that will be passed to eval
%       .$(custom): Custom fields
%        Variables available are:
%          physical_constants
%          "gps_time" (sec)
%          "along_track" (m)
%          "lat" (deg)
%          "lon" (deg)
%          "elev" (m)
%          "source" (twtt in sec)
%          "eval_struct" (the eval structure passed in by the user)
%        The cmd string should generally update "source" variable. For example:
%           '[B,A] = butter(0.1,2); source = filtfilt(B,A,source);' % Filter
%           'source = source + 0.1;' % Apply a twtt shift
%           'source = source*2;' % Surface multiple
%           'source = interp1(gps_time+15,source,gps_time);' % Apply a GPS time shift
%     .quality: struct controlling the quality level
%       .mode: string specifying one of these methods:
%         'overwrite': Overwrites the quality level with ".quality"
%         'preserve': Quality level preserved from source
%       .value: scalar containing 1, 2, or 3
%     .copy_method = string specifying one of these methods
%       'fillgaps': only gaps in destination data will be written to
%       'overwrite': none of the previous data are kept
%       'merge': anywhere source data exists, destination data will be
%         overwritten
%     .gaps_fill: struct controlling interpolation across gaps in source
%       .method: string specifying one of these methods
%         'preserve_gaps': runs opsInterpLayersToMasterGPSTime which tries
%           to interpolate where there is data and leave gaps where there
%           is not data (method_args controls this behavior)
%         'interp_finite': runs interp_finite which fills EVERY point
%       .method_args: arguments based on chosen method
%         'preserve_gaps': two element vector used by opsInterpLayersToMasterGPSTime
%         'interp_finite': not used
%
% Authors: John Paden, Abbey Whisler
%
% See also: runOpsCopyLayers.m, opsMergeLayerData, opsCopyLayers

physical_constants;

%% Determine if the copy method and gap filling method are valid
if ~any(strcmpi(copy_param.copy_method,{'fillgaps','overwrite','merge'}))
  error('Invalid copy_method %s', copy_param.copy_method);
end

if ~any(strcmpi(copy_param.gaps_fill.method,{'preserve_gaps','interp_finite'}))
  error('Invalid gap_fill.method %s', copy_param.gaps_fill.method);
end

if strcmpi(copy_param.gaps_fill.method,'preserve_gaps') ...
    && (~isfield(copy_param.gaps_fill,'method_args') || isempty(copy_param.gaps_fill.method_args))
  copy_param.gaps_fill.method_args = [300 60];
end

if ~isfield(copy_param,'quality')
  copy_param.quality.mode = 'preserve';
  copy_param.quality.value = 1;
end

if ~any(strcmpi(copy_param.quality.mode,{'preserve','overwrite'}))
  error('Invalid quality mode %s', copy_param.quality.mode);
end

if ~any(copy_param.quality.value == [1 2 3])
  error('Invalid quality value %d', copy_param.quality.value);
end

if ~isfield(copy_param.layer_dest,'group_name') || isempty(copy_param.layer_dest.group_name)
  % Default is no group name
  copy_param.layer_dest.group_name = '';
end

if ~isfield(copy_param.layer_dest,'layerdata_source') || isempty(copy_param.layer_dest.layerdata_source)
  % Default is the CSARP_layer directory
  copy_param.layer_dest.layerdata_source = 'layer';
end

% Force copy_param.layer_dest.name to be a cell array even if there is just
% one name. The source code below loops over the list of names.
if ischar(copy_param.layer_dest.name)
  copy_param.layer_dest.name = {copy_param.layer_dest.name};
end

if ~isfield(copy_param.layer_source,'layerdata_source') || isempty(copy_param.layer_source.layerdata_source)
  % Default is the CSARP_layer directory
  copy_param.layer_source.layerdata_source = 'layer';
end

if ~isfield(copy_param.layer_source,'echogram_source_img') || isempty(copy_param.layer_source.echogram_source_img)
  % Default is a combined file Data_YYYYMMDD_SS.mat
  copy_param.layer_source.echogram_source_img = 0;
end

%% Load "frames" file
frames = frames_load(param);
param.cmd.frms = frames_param_cmd_frms(param,frames);

%% Load source layer data and existing destination layer data for this segment
% Load all frames (for better edge interpolation)
load_param = param;
load_param.cmd.frms = max(1,min(param.cmd.frms)-1) : min(length(frames.frame_idxs),max(param.cmd.frms)+1);
if strcmpi(copy_param.layer_source.source,'custom')
  layer_source.gps_time = copy_param.layer_source.gps_time;
  layer_source.twtt = copy_param.layer_source.twtt;
  if isfield(copy_param.layer_source,'type') ...
      && ~isempty(copy_param.layer_source.type)
    layer_source.type = copy_param.layer_source.type;
  else
    layer_source.type = 2*ones(size(layer_source.gps_time));
  end
  if isfield(copy_param.layer_source,'quality') ...
      && ~isempty(copy_param.layer_source.quality)
    layer_source.quality = copy_param.layer_source.quality;
  else
    layer_source.quality = ones(size(layer_source.gps_time));
  end
else
  layer_source = opsLoadLayers(load_param,copy_param.layer_source);
end
layer_dest = opsLoadLayers(load_param,copy_param.layer_dest);

%% Load framing information (to determine start/stop gps times of each frame)
if strcmpi(copy_param.layer_dest.source,'ops')
  sys = ct_output_dir(param.radar_name);
  ops_param = struct('properties',[]);
  ops_param.properties.season = param.season_name;
  ops_param.properties.segment = param.day_seg;
  [~,ops_frames] = opsGetSegmentInfo(sys,ops_param);
end

%% Copy/Merge based on destination, copy_method, and gaps_fill.method.
%    Lists below describe the steps for each possible combination of
%    layer_dest.source, copy_method, and gaps_fill.method settings.
%
% If layer_dest.source is 'layerdata':
%
% copy_method = 'fillgaps',gaps_fill.method = 'interp_finite'
%   - Find NaN in destination
%   - Interpolate source onto these points using interpfinite style of interpolation
%   - Only update these points in the destination
% fillgaps,preserve_gaps
%   - Find NaN in destination
%   - Interpolate source onto these points using opsInterpLayersToMasterGPSTime
%   - Only update these points in the destination
% overwrite,interp_finite
%   - Find all points in destination
%   - Interpolate source onto these points using interpfinite style of interpolation
%   - Write these points on to the destination (does not require original layer information)
% overwrite,preserve_gaps
%   - Find all points in destination
%   - Interpolate source onto these points using opsInterpLayersToMasterGPSTime
%   - Write these points on to the destination (does not require original layer information)
% merge,interp_finite
%   - This is identical to overwrite,interp_finite (does not require original layer information)
% merge,preserve_gaps
%   - Find all points in destination
%   - Interpolate source onto these points using opsInterpLayersToMasterGPSTime
%   - Keep only the ~isnan points
%   - Only update these points in the destination
%
% If layer_dest.source is 'ops':
%
% fillgaps,interp_finite
%   - Find flightline points that don't have a pick in the destination
%   - Interpolate source onto these points using interpfinite style of interpolation
%   - Only update these points in the destination (does not require original layer information)
% fillgaps,preserve_gaps
%   - Find flightline points that don't have a pick in the destination
%   - Interpolate source onto these points using opsInterpLayersToMasterGPSTime
%   - Only update these points in the destination (does not require original layer information)
% overwrite,interp_finite
%   - Find all points in destination
%   - Interpolate source onto these points using interpfinite style of interpolation
%   - Write these points on to the destination (does not require original layer information)
% overwrite,preserve_gaps
%   - Find all points in destination
%   - Interpolate source onto these points using opsInterpLayersToMasterGPSTime
%   - Write these points on to the destination (does not require original layer information)
% merge,interp_finite
%   - This is identical to overwrite,interp_finite
% merge,preserve_gaps
%   - Find all points in destination
%   - Interpolate source onto these points using opsInterpLayersToMasterGPSTime
%   - Keep only the ~isnan points
%   - Only update these points in the destination (does not require original layer information)

%% Create all_points: a structure containing every point in the destination
% for this segments with these fields
%  .gps_time, .lat, .lon, .elev
%  .ids: only if ops is the destination (database point IDs)
%  .twtt: the current twtt for each point (NaN means does not exist)
%  .twtt_interp: the new twtt for each point
all_points = [];
if strcmpi(copy_param.layer_dest.source,'ops')
  % OPS source only returns points with picks from opsLoadLayers. So to get the
  % points that do not have picks, we have to load the flight lines.
  ops_param = struct('properties',[]);
  ops_param.properties.location = param.post.ops.location;
  ops_param.properties.season = param.season_name;
  ops_param.properties.start_gps_time = ops_frames.properties.start_gps_time(1);
  ops_param.properties.stop_gps_time = ops_frames.properties.stop_gps_time(end);
  ops_param.properties.nativeGeom = true;
  %ops_param.properties.lyr_id = ops_data.properties.id;
  [~,ops_data] = opsGetPath(sys,ops_param);
  
  for layer_idx = 1:length(layer_dest)
    all_points(layer_idx).gps_time = ops_data.properties.gps_time;
    all_points(layer_idx).lat = ops_data.properties.Y;
    all_points(layer_idx).lon = ops_data.properties.X;
    all_points(layer_idx).elev = ops_data.properties.elev;
    
    % Runs through each data point. At each ops_data point id, if there is
    % no layer_dest point id that matches it, that point is assigned NaN
    all_points(layer_idx).ids = ops_data.properties.id;
    all_points(layer_idx).twtt = zeros(size(ops_data.properties.id));
    
    for point_idx = 1:length(all_points(layer_idx).ids)
      match_idx = find(all_points(layer_idx).ids(point_idx) == layer_dest.point_path_id);
      if isempty(match_idx)
        all_points(layer_idx).twtt(point_idx) = NaN;
        all_points(layer_idx).quality(point_idx) = 1;
        all_points(layer_idx).type(point_idx) = 2;
      else
        all_points(layer_idx).twtt(point_idx) = layer_dest.twtt(match_idx);
        all_points(layer_idx).quality(point_idx) = layer_dest.quality(match_idx);
        all_points(layer_idx).type(point_idx) = layer_dest.type(match_idx);
      end
    end
  end
  
else
  % Non-OPS sources automatically include all points since NaN are used to
  % fill in missing points
  for layer_idx = 1:length(layer_dest)
    all_points(layer_idx).gps_time = layer_dest(layer_idx).gps_time;
    all_points(layer_idx).lat = layer_dest(layer_idx).lat;
    all_points(layer_idx).lon = layer_dest(layer_idx).lon;
    all_points(layer_idx).elev = layer_dest(layer_idx).elev;
    all_points(layer_idx).ids = layer_dest(layer_idx).point_path_id;
    all_points(layer_idx).type = layer_dest(layer_idx).type;
    all_points(layer_idx).twtt = layer_dest(layer_idx).twtt;
    all_points(layer_idx).quality = layer_dest(layer_idx).quality;
  end
end
% Check to make sure there are any points to copy
if any(cellfun(@isempty,{all_points.gps_time}))
  error('Destination layer(s) have no GPS time points. No work to do. This might be an error caused by echogram destination source not finding any echograms.');
end

for layer_idx = 1:length(layer_source)
  % Set invalid types to "auto" or 2
  layer_source(layer_idx).type(~isfinite(layer_source(layer_idx).type)) = 2;
  % Set invalid quality levels to "good" or 1
  layer_source(layer_idx).quality(~isfinite(layer_source(layer_idx).quality)) = 1;
  layer_source(layer_idx).quality(layer_source(layer_idx).quality ~= 1 & layer_source(layer_idx).quality ~= 2 & layer_source(layer_idx).quality ~= 3) = 1;
end

%% Apply evaluation operation to source
if isfield(copy_param,'eval') && ~isempty(copy_param.eval)
  for layer_idx = 1:length(layer_source)
    source = layer_source(layer_idx).twtt;
    gps_time = layer_source(layer_idx).gps_time;
    lat = layer_source(layer_idx).lat;
    lon = layer_source(layer_idx).lon;
    elev = layer_source(layer_idx).elev;
    along_track = geodetic_to_along_track(lat,lon,elev);
    eval_struct = copy_param.eval;
    eval(copy_param.eval.cmd);
    layer_source(layer_idx).twtt = source;
  end
end

%% Interpolate two way travel time (twtt) from layer_source to all_points
if strcmpi(copy_param.gaps_fill.method,'preserve_gaps')
  % interpolation preserves_gaps which preserves gaps in the data
  
  for layer_idx = 1:length(layer_source)
    master = [];
    master.GPS_time = all_points(layer_idx).gps_time;
    master.Latitude = all_points(layer_idx).lat;
    master.Longitude = all_points(layer_idx).lon;
    master.Elevation = all_points(layer_idx).elev;
    ops_layer = [];
    ops_layer{1}.gps_time = layer_source(layer_idx).gps_time;
    ops_layer{1}.type = layer_source(layer_idx).type;
    ops_layer{1}.quality = layer_source(layer_idx).quality;
    ops_layer{1}.twtt = layer_source(layer_idx).twtt;
    if length(master.GPS_time) < 2
      error('Too few points in destination to interpolate onto.');
    end
    lay = opsInterpLayersToMasterGPSTime(master,ops_layer,copy_param.gaps_fill.method_args);
    all_points(layer_idx).twtt_interp = lay.layerData{1}.value{2}.data;
    % All points automated (type 2) by default
    all_points(layer_idx).type_interp = 2*ones(size(all_points(layer_idx).twtt_interp));
    % Overwrite manual points (type 1)
    manual_mask = isfinite(lay.layerData{1}.value{1}.data);
    all_points(layer_idx).type_interp(manual_mask) = 1;
    all_points(layer_idx).twtt_interp(manual_mask) = lay.layerData{1}.value{1}.data(manual_mask);
  end
  
elseif strcmpi(copy_param.gaps_fill.method,'interp_finite')
  % interpolation interp_finite which interpolates through gaps
  
  for layer_idx = 1:length(layer_source)
    % Interpolate source onto destination points using linear interpolation
    if length(layer_source(layer_idx).gps_time) >= 2
      all_points(layer_idx).twtt_interp = interp1(layer_source(layer_idx).gps_time, layer_source(layer_idx).twtt, all_points(layer_idx).gps_time);
    else
      all_points(layer_idx).twtt_interp = NaN(size(all_points(layer_idx).gps_time));
    end
    % Fill in NaN gaps using interp_finite (default is 0 twtt if no good data
    % exists)
    all_points(layer_idx).twtt_interp = interp_finite(all_points(layer_idx).twtt_interp,0);
    
    % All points automated (type 2) by default
    all_points(layer_idx).type_interp = 2*ones(size(all_points(layer_idx).twtt_interp));
    % Overwrite manual points (type 1)
    manual_idxs = find(layer_source(layer_idx).type == 1);
    manual_idxs_map = interp1_robust(all_points(layer_idx).gps_time, 1:length(all_points(layer_idx).gps_time), ...
      layer_source(layer_idx).gps_time(manual_idxs),'nearest');
    manual_idxs = manual_idxs(~isnan(manual_idxs_map));
    manual_idxs_map = manual_idxs_map(~isnan(manual_idxs_map));
    all_points(layer_idx).twtt_interp(manual_idxs_map) ...
      = layer_source(layer_idx).twtt(manual_idxs);
    all_points(layer_idx).type_interp(manual_idxs_map) = 1;
  end
end

update_mask = {};
for layer_idx = 1:length(layer_source)
  
  %% Combine the newly interpolated result with the current values
  if strcmpi(copy_param.copy_method,'fillgaps')
    update_mask{layer_idx} = isnan(all_points(layer_idx).twtt) & ~isnan(all_points(layer_idx).twtt_interp);
  elseif strcmpi(copy_param.copy_method,'merge')
    update_mask{layer_idx} = ~isnan(all_points(layer_idx).twtt_interp);
  elseif strcmpi(copy_param.copy_method,'overwrite')
    update_mask{layer_idx} = logical(ones(size(all_points(layer_idx).twtt)));
  end
  
  %% For nonselected frames, set update_mask{layer_idx} to zero for all points in those frames
  all_frms = 1:length(frames.frame_idxs);
  % If the user selects all frames in the segment on the params
  % spreadsheet, then the update mask will not change.
  % Otherwise we will create a new mask where each point within the start
  % and stop gps times of the selected frames is equal to 1 and all other
  % points are equal to 0.
  %
  % We divide up the segment into frames so that we avoid rounding issues and
  % to every point is included in a frame.
  frms_mask = zeros(size(update_mask{layer_idx}));
  if strcmpi(copy_param.layer_dest.source,'ops')
    for frm = param.cmd.frms
      if frm == 1
        if frm == length(ops_frames.properties.start_gps_time)
          frms_mask(:) = 1;
        else
          frms_mask(all_points(layer_idx).gps_time < ops_frames.properties.start_gps_time(frm+1)) = 1;
        end
      elseif frm < length(ops_frames.properties.start_gps_time)
        frms_mask(all_points(layer_idx).gps_time >= ops_frames.properties.start_gps_time(frm)...
          & all_points(layer_idx).gps_time < ops_frames.properties.start_gps_time(frm+1)) = 1;
      else
        frms_mask(all_points(layer_idx).gps_time >= ops_frames.properties.start_gps_time(frm)) = 1;
      end
    end
    
  else
    for frm = param.cmd.frms
      if frm == 1
        if frm == length(frames.frame_idxs)
          frms_mask(:) = 1;
        else
          frms_mask(all_points(layer_idx).gps_time < frames.gps_time(frm+1)) = 1;
        end
      elseif frm < length(frames.frame_idxs)
        frms_mask(all_points(layer_idx).gps_time >= frames.gps_time(frm)...
          & all_points(layer_idx).gps_time < frames.gps_time(frm+1)) = 1;
      else
        frms_mask(all_points(layer_idx).gps_time >= frames.gps_time(frm)) = 1;
      end
    end
  end
  
  update_mask{layer_idx} = frms_mask & update_mask{layer_idx};
  
  %% Quality level overwrite
  if isempty(layer_source(layer_idx).gps_time)
    all_points(layer_idx).quality_interp = NaN*zeros(size(all_points(layer_idx).gps_time));
  elseif length(layer_source(layer_idx).gps_time) == 1
    all_points(layer_idx).quality_interp = layer_source(layer_idx).quality;
  else
    all_points(layer_idx).quality_interp = interp1(layer_source(layer_idx).gps_time,layer_source(layer_idx).quality,all_points(layer_idx).gps_time,'nearest');
  end
  if strcmpi(copy_param.quality.mode,'overwrite')
    all_points(layer_idx).quality_interp(:) = copy_param.quality.value;
  end
  
  %% Mask the output
  all_points(layer_idx).twtt_final = all_points(layer_idx).twtt;
  all_points(layer_idx).twtt_final(update_mask{layer_idx}) = all_points(layer_idx).twtt_interp(update_mask{layer_idx});
  all_points(layer_idx).quality_final = all_points(layer_idx).quality;
  all_points(layer_idx).quality_final(update_mask{layer_idx}) = all_points(layer_idx).quality_interp(update_mask{layer_idx});
  all_points(layer_idx).type_final = all_points(layer_idx).type;
  all_points(layer_idx).type_final(update_mask{layer_idx}) = all_points(layer_idx).type_interp(update_mask{layer_idx});
  
  if 0
    % Debug code
    figure(1); clf;
    plot(all_points(layer_idx).twtt);
    hold on;
    plot(surface(update_mask{layer_idx}),'rx');
    grid on;
    hold off;
    drawnow;
    pause(1);
  end
end

if strcmpi(copy_param.layer_dest.source,'ops')
  %% Save OPS
  
  % Check to see if layer exists
  % -----------------------------------------------------------------------
  [status,data] = opsGetLayers(sys);
  for layer_idx = 1:length(copy_param.layer_dest.name)
    if ~any(strcmpi(data.properties.lyr_name,copy_param.layer_dest.name))
      % Create the layer if it does not exist
      ops_param = [];
      ops_param.properties.lyr_name = copy_param.layer_dest.name;
      ops_param.properties.lyr_group_name = copy_param.layer_dest.group;
      ops_param.properties.lyr_description = copy_param.layer_dest.description;
      ops_param.properties.public = true;
      
      [status,ops_data] = opsCreateLayer(sys,ops_param);
    end
    
    % Remove all_points that are not in the selected frames
    % Use update_mask to exclude all points that are not getting updated
    % -----------------------------------------------------------------------
    ops_param.properties.point_path_id = all_points(layer_idx).ids(update_mask);
    ops_param.properties.twtt = all_points(layer_idx).twtt_final(update_mask);
    ops_param.properties.type = all_points(layer_idx).type_final(update_mask);
    ops_param.properties.quality = all_points(layer_idx).quality_final(update_mask);
    ops_param.properties.lyr_name = copy_param.layer_dest.name;
    
    % Update these points
    % -----------------------------------------------------------------------
    opsCreateLayerPoints(sys,ops_param);
  end
  
elseif strcmpi(copy_param.layer_dest.source,'layerdata')
  %% Save layerdata
  layers = layerdata(param, copy_param.layer_dest.layerdata_source);
  for layer_idx = 1:length(copy_param.layer_dest.name)
    id = layers.get_id(copy_param.layer_dest.name{layer_idx});
    if isempty(id)
      if copy_param.layer_dest.existence_check
        error('Layer %s not found in layer organizer file %s. Set copy_param.layer_dest.existence_check to false to ignore this error and opsCopyLayers will create layers that do not exist already.\n', copy_param.layer_dest.name{layer_idx}, layers.layer_organizer_fn());
      else
        layer_organizer = [];
        layer_organizer.lyr_group_name = {copy_param.layer_dest.group_name};
        layer_organizer.lyr_name = {copy_param.layer_dest.name{layer_idx}};
        id = layers.insert_layers(layer_organizer);
      end
    end
    layers.update_layer(param.cmd.frms, id, all_points(layer_idx).gps_time, ...
      all_points(layer_idx).twtt_final,all_points(layer_idx).quality_final,all_points(layer_idx).type_final);
  end
  layers.save();
  
end
