function load_layers(obj)
% echowin.load_layers(obj)
%
% Load layer information from database and update layer plot handles

%% Determine GPS time domain for layers to load
dx = obj.eg.gps_time(2)-obj.eg.gps_time(1);
min_gps_time = obj.eg.gps_time(1);
max_gps_time = obj.eg.gps_time(end);
if obj.eg.frms(1) == 1
  % First frame in segment so adjust beginning to make sure all layer points
  % will be displayed.
  min_gps_time = obj.eg.start_gps_time(1)-dx;
end
if obj.eg.frms(end) == length(obj.eg.stop_gps_time)
  % Last frame of segment so adjust end to make sure all layer points will
  % be displayed.
  max_gps_time = obj.eg.stop_gps_time(end)+dx;
end

if strcmpi(obj.eg.layers.source,'OPS')
  %% OPS: Load layer points from database
  fprintf(' Loading layer points from database (%s)\n',datestr(now,'HH:MM:SS'));
  ops_param = struct('properties',[]);
  ops_param.properties.location = obj.eg.cur_sel.location;
  ops_param.properties.season = obj.eg.cur_sel.season_name;
  ops_param.properties.segment_id = obj.eg.cur_sel.segment_id;
  ops_param.properties.start_gps_time = min_gps_time;
  ops_param.properties.stop_gps_time = max_gps_time;
  ops_param.properties.lyr_name = obj.eg.layers.lyr_name;
  % query database
  [status,data] = opsGetLayerPoints(obj.eg.system,ops_param);
  
  %% OPS: Preallocate layer arrays
  obj.eg.layers.x = {};
  obj.eg.layers.y = {};
  obj.eg.layers.qual = {};
  obj.eg.layers.type = {};
  for idx = 1:length(obj.eg.layers.lyr_id)
    obj.eg.layers.x{idx} = double(obj.eg.map_gps_time); % gps-time
    obj.eg.layers.y{idx} = nan(size(obj.eg.map_id)); % twtt
    obj.eg.layers.qual{idx} = nan(size(obj.eg.map_id)); % integer 1-3
    obj.eg.layers.type{idx} = nan(size(obj.eg.map_id)); % this is either 1 (manual) or 2 (auto)
  end
  
  %% OPS: Fill layer arrays with points from database
  for idx = 1:length(data.properties.point_path_id)
    layer_idx = obj.eg.layers.lyr_id == data.properties.lyr_id(idx);
    point_path_idx = obj.eg.map_id == data.properties.point_path_id(idx);
    obj.eg.layers.y{layer_idx}(point_path_idx) = data.properties.twtt(idx);
    obj.eg.layers.qual{layer_idx}(point_path_idx) = data.properties.quality(idx);
    obj.eg.layers.type{layer_idx}(point_path_idx) = data.properties.type(idx);
  end
  
  %% OPS: Update echogram surface if there are enough good points from OPS
  % Find good surface points
  if ~any(obj.eg.layers.lyr_id==1)
    error('standard:surface must be added in the "Layers" preference window before loading echograms.');
  end
  good_mask = ~isnan(obj.eg.layers.y{obj.eg.layers.lyr_id==1});
  if sum(good_mask) > 2
    % There are surface layer points in the database, overwrite the surface
    % with these
    obj.eg.surf_twtt = interp1(obj.eg.map_gps_time(good_mask),obj.eg.layers.y{obj.eg.layers.lyr_id==1}(good_mask),obj.eg.gps_time);
    obj.eg.surf_twtt = interp_finite(obj.eg.surf_twtt,0);
  end
  
elseif strcmpi(obj.eg.layers.source,'layerdata')
  %% LayerData: Load layer points from layerdata
  
  %% LayerData: Preallocate layer arrays
  fprintf(' Loading layer points from layerData (%s)\n',datestr(now));
  for idx = 1:length(obj.eg.layers.lyr_id)
    obj.eg.layers.x{idx} = []; % gps time
    obj.eg.layers.y{idx} = []; % twtt
    obj.eg.layers.qual{idx} = []; % integer 1-3
    obj.eg.layers.type{idx} = []; % this is either 1 (manual) or 2 (auto)
  end
  %% LayerData: Load GPS_time, quality, twtt and type of layers from undo_stack
  lGPS = [];
  for frm = obj.eg.frms(1):obj.eg.frms(end);
    GPS_time = obj.undo_stack.user_data.layer_info(frm).GPS_time;
    lGPS = cat(2,lGPS,GPS_time); % concatenates the layer GPS time
    
    for idx=1:length(obj.eg.layers.lyr_id)
      obj.eg.layers.x{idx} = cat(2,obj.eg.layers.x{idx},GPS_time); % gps time
      qual = obj.undo_stack.user_data.layer_info(frm).layerData{idx}.quality;
      obj.eg.layers.qual{idx} = cat(2,obj.eg.layers.qual{idx},qual); % quality (integer 1-3)
      twtt = obj.undo_stack.user_data.layer_info(frm).layerData{idx}.value{2}.data;
      obj.eg.layers.y{idx} = cat(2,obj.eg.layers.y{idx},twtt); % twtt
      obj.eg.layers.type{idx} = cat(2,obj.eg.layers.type{idx},1 + ~isfinite(obj.undo_stack.user_data.layer_info(frm).layerData{idx}.value{1}.data)); % this is either 1 (manual) or 2 (auto)
    end
  end
  
  %% LayerData: Update echogram surface if there are enough good points 
  % Find good surface points
  if ~any(obj.eg.layers.lyr_id==1)
    error('standard:surface must be added in the "Layers" preference window before loading echograms.');
  end
  
  % logical vector with 1 where the twtt of the surface is a number and 0
  % when NaN.
  good_mask = ~isnan(obj.eg.layers.y{obj.eg.layers.lyr_id==1});
  if sum(good_mask) > 2
    obj.eg.surf_twtt = interp1(lGPS(good_mask),obj.eg.layers.y{obj.eg.layers.lyr_id==1}(good_mask),obj.eg.gps_time);
    obj.eg.surf_twtt = interp_finite(obj.eg.surf_twtt,0);
  end
end
