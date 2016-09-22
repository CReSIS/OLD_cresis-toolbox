function load_layers(obj)
% echowin.load_layers(obj)
%
% Load layer information from database and update layer plot handles

dx = obj.eg.gps_time(2)-obj.eg.gps_time(1);
min_gps_time = obj.eg.gps_time(1);
max_gps_time = obj.eg.gps_time(end);
if obj.eg.frame_idxs(1) == 1
  % First frame in segment so adjust beginning to make sure all layer points
  % will be displayed.
  min_gps_time = obj.eg.start_gps_time(1)-dx;
elseif obj.eg.frame_idxs(end) == length(obj.eg.stop_gps_time)
  % Last frame of segment so adjust end to make sure all layer points will 
  % be displayed.
  max_gps_time = obj.eg.stop_gps_time(end)+dx;
end

fprintf(' Loading layer points from database (%s)\n',datestr(now,'HH:MM:SS'));
ops_param = struct('properties',[]);
ops_param.properties.location = obj.eg.cur_sel.location;
ops_param.properties.season = obj.eg.cur_sel.season_name;
ops_param.properties.segment_id = obj.eg.cur_sel.segment_id;
ops_param.properties.start_gps_time = min_gps_time;
ops_param.properties.stop_gps_time = max_gps_time;
ops_param.properties.lyr_name = obj.eg.layer.layer_names;
% query database
[status,data] = opsGetLayerPoints(obj.eg.system,ops_param);

layer_type = data.properties.lyr_id;
% Preallocate layer arrays 
obj.eg.layer.x = {};
obj.eg.layer.y = {};
obj.eg.layer.qual = {};
obj.eg.layer.type = {};
for idx = 1:length(obj.eg.layer_id)
  obj.eg.layer.x{idx} = double(obj.eg.map_gps_time); % gps-time
  obj.eg.layer.y{idx} = NaN*zeros(size(obj.eg.map_id)); % twtt
  obj.eg.layer.qual{idx} = NaN*zeros(size(obj.eg.map_id)); % integer 1-3
  obj.eg.layer.type{idx} = NaN*zeros(size(obj.eg.map_id)); % this is either 1 (manual) or 2 (auto)
end
% Fill layer arrays with points from database
for idx = 1:length(data.properties.point_path_id)
  layer_idx = obj.eg.layers.lyr_id == data.properties.lyr_id(idx);
  point_path_idx = obj.eg.map_id == data.properties.point_path_id(idx);
  obj.eg.layer.y{layer_idx}(point_path_idx) = data.properties.twtt(idx);
  obj.eg.layer.qual{layer_idx}(point_path_idx) = data.properties.quality(idx);
  obj.eg.layer.type{layer_idx}(point_path_idx) = data.properties.type(idx);
end

%% Update echogram surface if there are enough good points from OPS
% Find good surface points
if ~any(obj.eg.layer_id==1)
  error('standard:surface must be added in the "Layers" preference window before loading echograms.');
end
good_mask = ~isnan(obj.eg.layer.y{obj.eg.layer_id==1});
if sum(good_mask) > 2
  % There are surface layer points in the database, overwrite the surface
  % with these
  obj.eg.surface = interp1(obj.eg.map_gps_time(good_mask),obj.eg.layer.y{obj.eg.layer_id==1}(good_mask),obj.eg.gps_time);
  obj.eg.surface = interp_finite(obj.eg.surface,0);
end

end
