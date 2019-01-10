function load_layers(obj)
frm=2;
% echowin.load_layers(obj)
%
% Load layer information from database and update layer plot handles

%Loading frames file  
  %load(ct_filename_support(param,'','frames'));
  %for frm_idx=1:length(param.frms)
  %frm=param.frms(frm_idx);
  
% if strcmpi(obj.source,'OPS')
% OPS: Determine GPS time domain to load from data base
dx = obj.eg.gps_time(2)-obj.eg.gps_time(1);
min_gps_time = obj.eg.gps_time(1);
max_gps_time = obj.eg.gps_time(end);
if obj.eg.frame_idxs(1) == 1
  % First frame in segment so adjust beginning to make sure all layer points
  % will be displayed.
  min_gps_time = obj.eg.start_gps_time(1)-dx;
end
if obj.eg.frame_idxs(end) == length(obj.eg.stop_gps_time)
  % Last frame of segment so adjust end to make sure all layer points will 
  % be displayed.
  max_gps_time = obj.eg.stop_gps_time(end)+dx;
end

% OPS: Load layer points from database
fprintf(' Loading layer points from database (%s)\n',datestr(now,'HH:MM:SS'));
if 0
ops_param = struct('properties',[]);
ops_param.properties.location = obj.eg.cur_sel.location;
ops_param.properties.season = obj.eg.cur_sel.season_name;
ops_param.properties.segment_id = obj.eg.cur_sel.segment_id;
ops_param.properties.start_gps_time = min_gps_time;
ops_param.properties.stop_gps_time = max_gps_time;
ops_param.properties.lyr_name = obj.eg.layer.layer_names;
% query database
[status,data] = opsGetLayerPoints(obj.eg.system,ops_param);

% OPS: Preallocate layer arrays 
obj.eg.layer.x = {};
obj.eg.layer.y = {};
obj.eg.layer.qual = {};
obj.eg.layer.type = {};
for idx = 1:length(obj.eg.layer_id)
  obj.eg.layer.x{idx} = double(obj.eg.map_gps_time); % gps-time
  obj.eg.layer.y{idx} = nan(size(obj.eg.map_id)); % twtt
  obj.eg.layer.qual{idx} = nan(size(obj.eg.map_id)); % integer 1-3
  obj.eg.layer.type{idx} = nan(size(obj.eg.map_id)); % this is either 1 (manual) or 2 (auto)
end

% OPS: Fill layer arrays with points from database
% start with layer Data gps time, layer Data
% twtt, ops gps time and that would give you the twtt interpolated onto the map or OPS
% time. Check for quality and type. 
for idx = 1:length(data.properties.point_path_id)
  layer_idx = obj.eg.layers.lyr_id == data.properties.lyr_id(idx);
  point_path_idx = obj.eg.map_id == data.properties.point_path_id(idx);
  obj.eg.layer.y{layer_idx}(point_path_idx) = data.properties.twtt(idx);
  obj.eg.layer.qual{layer_idx}(point_path_idx) = data.properties.quality(idx);
  obj.eg.layer.type{layer_idx}(point_path_idx) = data.properties.type(idx);
end

elseif 1
  layer_fn=fullfile(ct_filename_out(obj.eg.cur_sel,'layerData',''),sprintf('Data_%s_%03d.mat','20110502_01',frm));
  %layer_fn=fullfile('X:\ct_data\rds\2011_Greenland_P3\CSARP_standard\20110502_01',frm);
  lay=load(layer_fn);
  %Preallocate layer arrays
  obj.eg.layer.x = {};
  obj.eg.layer.y = {};
  obj.eg.layer.qual = {};
  obj.eg.layer.type = {};
  for idx = 1:length(obj.eg.layer_id)
    %obj.eg.layer.x{idx} = double(obj.eg.map_gps_time); % gps-time
    obj.eg.layer.x{idx}=lay.GPS_time;
    obj.eg.layer.y{idx} = nan(size(obj.eg.map_id)); % twtt
    obj.eg.layer.qual{idx} = nan(size(obj.eg.map_id)); % integer 1-3
    obj.eg.layer.type{idx} = nan(size(obj.eg.map_id)); % this is either 1 (manual) or 2 (auto)
  end
  %Filling the layers
  for idx=1:length(obj.eg.layer_id)
    %for lay_idx=1:lay.layerData;
    % y has twtt of the layer data interpolated on the echogram GPS time
      %y=interp1(lay.GPS_time, lay.layerData{idx}.value{2}.data, obj.eg.layer.x{idx});
      %obj.eg.layer.y{idx}=y;
      obj.eg.layer.y{idx}=lay.layerData{idx}.value{2}.data;
      obj.eg.layer.qual{idx}=lay.layerData{idx}.quality;
      %if not finite fill the type with 2 else 1
      if ~isfinite(lay.layerData{idx}.value{2}.data)
        obj.eg.layer.type{idx}= 2*ones(size(lay.GPS_time));
      else
        obj.eg.layer.type{idx}=ones(size(lay.GPS_time)); 
      end
  end
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
