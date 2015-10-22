function ops_layer_points_param = layerdata_to_ops(layerData_fn,settings)
% ops_layer_points_param = layerdata_to_ops(layerData_fn,settings)
%
% Convert CReSIS layerData to the database insert param structure.
%
% Input:
%   layerData_fn: Absolute path to a CReSIS layerData (.m) file.
%   settings: optional settings structure with the following fields
%     .layer_filter = function which takes one character array argument
%       and returns true for layers which should be inserted (may be
%       undefined or empty if all layers are to be inserted)
%
% NOTE: If the layerData file contains more than 2 layers it must have a
% name field for each layer. Surface and Bottom should be named 'surface'
% and 'bottom' and the rest can be any name.
%
% Output:
%   ops_layer_points_param: structure array with fields
%     geometry.coordinates = double array of format ([lon lat elevation])
%     properties.gps_time = double array
%     properties.twtt = double array
%     properties.type = cell of strings ('auto','manual' or'NULL')
%     properties.quality = integer array (1,2 or 3)
%     properties.lyr_name = string ('suface' or 'bottom' or other name (internal layers))
%
% Author: Kyle W. Purdon

% Parse input arguments and set defaults
if ~exist('settings','var')
  settings = struct();
end
if ~isfield(settings,'layer_filter') || isempty(settings.layer_filter)
  settings.layer_filter = inline('~isempty(regexp(x,''.*''))');
end

% LOAD layerData
lyr = load(layerData_fn,'GPS_time','Latitude','Longitude','Elevation','layerData');

if ~isfield(lyr, 'layerData')
  lyr = load(layerData_fn,'GPS_time','Latitude','Longitude','Elevation','Surface','Truncate_Bins','Elevation_Correction','Time');
  lyr = uncompress_echogram(lyr);
  lyr.layerData{1}.value{1}.data = NaN*zeros(size(lyr.Surface));
  lyr.layerData{1}.value{2}.data = lyr.Surface;
  lyr.layerData{1}.quality = ones(size(lyr.Surface));
end

% CHECK if this is the new layerData format (from OPS)
newLd = false;
if ~isfield(lyr.layerData{1},'value')
  newLd = true;
end

insert_idx = 0;
ops_layer_points_param = [];
for layer_idx = 1:length(lyr.layerData)
  gps_time = [];
  lat = [];
  lon = [];
  elev = [];
  layer = [];
  layer_type = [];
  layer_quality = [];

  % Set default layer names
  if ~isfield(lyr.layerData{layer_idx},'name')
    if layer_idx == 1
      lyr.layerData{layer_idx}.name = 'surface';
    elseif layer_idx == 2
      lyr.layerData{layer_idx}.name = 'bottom';
    end
  end

  % Check to see if this layer should be inserted
  if ~settings.layer_filter(lyr.layerData{layer_idx}.name)
    continue;
  end
  insert_idx = insert_idx + 1;
  
  if ~newLd
    %% Add manual layer points (layer_type == 1)
    layer_orig = lyr.layerData{layer_idx}.value{1}.data;
    good_idxs = find(~isnan(layer_orig) & isfinite(layer_orig));
    gps_time = cat(2,gps_time,lyr.GPS_time(good_idxs));
    lat = cat(2,lat,lyr.Latitude(good_idxs));
    lon = cat(2,lon,lyr.Longitude(good_idxs));
    elev = cat(2,elev,lyr.Elevation(good_idxs));
    layer = cat(2,layer,layer_orig(good_idxs));
    layer_quality = cat(2,layer_quality,lyr.layerData{layer_idx}.quality(good_idxs));
    layer_type = cat(2,layer_type,1*ones(size(layer_orig(good_idxs))));
    
    %% Add automated layer points (layer_type == 2)
    layer_orig = lyr.layerData{layer_idx}.value{2}.data;
    good_idxs = find(~isnan(layer_orig) & isfinite(layer_orig));
    gps_time = cat(2,gps_time,lyr.GPS_time(good_idxs));
    lat = cat(2,lat,lyr.Latitude(good_idxs));
    lon = cat(2,lon,lyr.Longitude(good_idxs));
    elev = cat(2,elev,lyr.Elevation(good_idxs));
    layer = cat(2,layer,layer_orig(good_idxs));
    layer_quality = cat(2,layer_quality,lyr.layerData{layer_idx}.quality(good_idxs));
    layer_type = cat(2,layer_type,2*ones(size(layer_orig(good_idxs))));
  else
    % MANUAL POINTS
    good_idxs = find(lyr.layerData{layer_idx}.type == 1);
    gps_time = cat(2,gps_time,lyr.layerData{layer_idx}.gps_time(good_idxs));
    lat = cat(2,lat,lyr.layerData{layer_idx}.latitude(good_idxs));
    lon = cat(2,lon,lyr.layerData{layer_idx}.longitude(good_idxs));
    elev = cat(2,elev,lyr.layerData{layer_idx}.elevation(good_idxs));
    layer = cat(2,layer,lyr.layerData{layer_idx}.twtt(good_idxs)); %twtt
    layer_quality = cat(2,layer_quality,lyr.layerData{layer_idx}.quality(good_idxs));
    layer_type = cat(2,layer_type,lyr.layerData{layer_idx}.type(good_idxs));
    
    % AUTO POINTS
    good_idxs = find(lyr.layerData{layer_idx}.type == 2);
    gps_time = cat(2,gps_time,lyr.layerData{layer_idx}.gps_time(good_idxs));
    lat = cat(2,lat,lyr.layerData{layer_idx}.latitude(good_idxs));
    lon = cat(2,lon,lyr.layerData{layer_idx}.longitude(good_idxs));
    elev = cat(2,elev,lyr.layerData{layer_idx}.elevation(good_idxs));
    layer = cat(2,layer,lyr.layerData{layer_idx}.twtt(good_idxs)); %twtt
    layer_quality = cat(2,layer_quality,lyr.layerData{layer_idx}.quality(good_idxs));
    layer_type = cat(2,layer_type,lyr.layerData{layer_idx}.type(good_idxs));
  end
  
  ops_layer_points_param(insert_idx).properties.lyr_name = lower(lyr.layerData{layer_idx}.name);
  
  % CORRECT NAN QUALITY
  bad_idxs = find(isnan(layer_quality));
  layer_quality(bad_idxs) = 1;
  
  % STORE DETAILED LAYER INFORMATION
  ops_layer_points_param(insert_idx).properties.gps_time = gps_time;
  ops_layer_points_param(insert_idx).geometry.coordinates = [lon; lat; elev]';
  ops_layer_points_param(insert_idx).properties.twtt = double(layer);
  ops_layer_points_param(insert_idx).properties.type = layer_type;
  ops_layer_points_param(insert_idx).properties.quality = layer_quality;
  
end
end