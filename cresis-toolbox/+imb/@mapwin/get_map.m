function get_map(obj,hObj,event)
% get_map(obj,hObj,event)
%
% This is the callback function which is called when the preference
% window "OK" button is pressed and the prefwin "StateChange" event occurs.

%% Check which settings have changed
if ~strcmpi(obj.cur_map_pref_settings.system, obj.map_pref.settings.system)
  system_changed = true;
else
  system_changed = false;
end

seasons_changed = false;
if length(obj.cur_map_pref_settings.seasons) ~= length(obj.map_pref.settings.seasons)
  seasons_changed = true;
else
  for idx = 1:length(obj.cur_map_pref_settings.seasons)
    if ~any(strcmp(obj.cur_map_pref_settings.seasons{idx},obj.map_pref.settings.seasons))
      seasons_changed = true;
    end
  end
end

if ~strcmpi(obj.cur_map_pref_settings.mapname,obj.map_pref.settings.mapname)
  mapname_changed = true;
else
  mapname_changed = false;
end

if ~strcmpi(obj.cur_map_pref_settings.flightlines,obj.map_pref.settings.flightlines)
  flightlines_changed = true;
else
  flightlines_changed = false;
end

if ~strcmpi(obj.cur_map_pref_settings.mapzone,obj.map_pref.settings.mapzone)
  mapzone_changed = true;
else
  mapzone_changed = false;
end

if ~system_changed && ~seasons_changed && ~mapname_changed && ~mapzone_changed && ~flightlines_changed
  % get_map only needs to update source and layers potentially
  obj.cur_map_pref_settings.sources = obj.map_pref.settings.sources;
  obj.cur_map_pref_settings.layers = obj.map_pref.settings.layers;
  return;
end

%% Copy current preference window settings over
obj.cur_map_pref_settings = obj.map_pref.settings;

%% Update map selection (also called at startup)
% =================================================================
flightlines = obj.cur_map_pref_settings.flightlines;
map_name = obj.cur_map_pref_settings.mapname;
map_zone = obj.cur_map_pref_settings.mapzone;
fprintf('Loading and plotting map %s (%s)\n', map_name, datestr(now,'HH:MM:SS'));

opsCmd;

% CONNECT TO THE WMS SERVER AND GET A LAYER OBJECT
wms = WebMapServer(sprintf('%s%s/wms/',gOps.geoServerUrl,map_zone));
cpbs = wms.getCapabilities();
layer = cpbs.Layer;

% REFINE THE LAYER BASED ON THE ACTIVE SELECTIONS
if strcmpi(flightlines,'Regular Flightlines') && ~isempty(layer.refine('line_paths'))
  layers = [];
  % ADD THE LINE PATH LAYER
  layers = cat(2,layers,layer.refine(sprintf('%s_%s_line_paths',map_zone,obj.cur_map_pref_settings.system)));
  % ADD THE BACKGROUND LAYER
  layers = cat(2,layers,layer.refine(map_name,'matchType','exact'));
  layer = layers.';
elseif strcmpi(flightlines,'Quality Flightlines') && ~isempty(layer.refine('data_quality'))
   layers = [];
  % ADD THE QUALITY LAYER
  layers = cat(2,layers,layer.refine(sprintf('%s_%s_data_quality',map_zone,obj.cur_map_pref_settings.system)));
  % ADD THE BACKGROUND LAYER
  layers = cat(2,layers,layer.refine(map_name,'matchType','exact'));
  layer = layers.';
elseif strcmpi(flightlines,'Coverage Flightlines') && ~isempty(layer.refine('data_coverage'))
  layers = [];
  % ADD THE COVERAGE LAYER
  layers = cat(2,layers,layer.refine(sprintf('%s_%s_data_coverage',map_zone,obj.cur_map_pref_settings.system)));
  % ADD THE BACKGROUND LAYER
  layers = cat(2,layers,layer.refine(map_name,'matchType','exact'));
  layer = layers.';
elseif strcmpi(flightlines,'Crossover Errors') && ~isempty(layer.refine('crossover_errors'))
   layers = [];
  % ADD THE Crossr Errors LAYER
  layers = cat(2,layers,layer.refine(sprintf('%s_%s_crossover_errors',map_zone,obj.cur_map_pref_settings.system)));
  % ADD THE LINE PATH LAYER
  layers = cat(2,layers,layer.refine(sprintf('%s_%s_line_paths',map_zone,obj.cur_map_pref_settings.system)));
  % ADD THE BACKGROUND LAYER
  layers = cat(2,layers,layer.refine(map_name,'matchType','exact'));
  layer = layers.';
elseif strcmpi(flightlines,'Bed Elevation') && ~isempty(layer.refine('data_elevation'))
   layers = [];
  % ADD THE Elevation LAYER
  layers = cat(2,layers,layer.refine(sprintf('%s_%s_data_elevation',map_zone,obj.cur_map_pref_settings.system)));
  % ADD THE BACKGROUND LAYER
  layers = cat(2,layers,layer.refine(map_name,'matchType','exact'));
  layer = layers.';
else
  % JUST ADD THE BACKGROUND LAYER
  layer = layer.refine(map_name,'matchType','exact').';
end

request = WMSMapRequest(layer);
if strcmp(map_zone,'arctic')
  request.CoordRefSysCode = 'EPSG:3413';
  % SET THE START-UP DEFAULT BOUNDING BOX
  bb_x = [-1500000 1500000];
  bb_y = [-4000000 0];
  obj.full_xaxis = bb_x/1e3;
  obj.full_yaxis = bb_y/1e3;
else
  request.CoordRefSysCode = 'EPSG:3031';
  % SET THE START-UP DEFAULT BOUNDING BOX
  bb_x = [-3400000 3400000];
  bb_y = [-3400000 3400000];
  obj.full_xaxis = bb_x/1e3;
  obj.full_yaxis = bb_y/1e3;
end

%% Get the bounding box
% BoundingBox contains xlim and ylim for all valid coordinate systems
% this loop finds the limits for only the relevant EPSG coordinate system
% and also ensures that both the map and the flightlines are fully
% displayed
% for idx1 = 1:length(layer)
%   for idx = 1:length(layer(idx1).Details.BoundingBox)
%     if strcmp(layer(idx1).Details.BoundingBox(idx).CoordRefSysCode,...
%         request.CoordRefSysCode)
%       sz(idx1,1:2) = layer(idx1).Details.BoundingBox(idx).XLim;
%       sz(idx1,3:4) = layer(idx1).Details.BoundingBox(idx).YLim;
%       break;
%     end
%   end
% end
% bb_x = [min(sz(:,1)) max(sz(:,2))];
% bb_y = [min(sz(:,3)) max(sz(:,4))];
% obj.full_xaxis = bb_x/1e3;
% obj.full_yaxis = bb_y/1e3;

%% Set the limits for the new map
if mapzone_changed
  request.XLim = obj.full_xaxis*1e3;
  request.YLim = obj.full_yaxis*1e3;
else
  request.XLim = get(obj.map_panel.h_axes,'XLim')*1e3;
  request.YLim = get(obj.map_panel.h_axes,'YLim')*1e3;
end

%% Fix the aspect ratio of the limits to fit properly in our window
old_u = get(obj.map_panel.h_axes,'units');
set(obj.map_panel.h_axes,'Units','pixels')
PixelBounds = round(get(obj.map_panel.h_axes,'Position'));
set(obj.map_panel.h_axes,'Position',PixelBounds);
set(obj.map_panel.h_axes,'units',old_u);

height = round((PixelBounds(4))*1 - 0);
width = round((PixelBounds(3))*1 - 0);
aspect_ratio = height/width;

if aspect_ratio*diff(request.XLim) > diff(request.YLim)
  growth = aspect_ratio*diff(request.XLim) - diff(request.YLim);
  request.YLim(1) = request.YLim(1) - growth/2;
  request.YLim(2) = request.YLim(2) + growth/2;
elseif aspect_ratio*diff(request.XLim) < diff(request.YLim)
  growth = diff(request.YLim)/aspect_ratio - diff(request.XLim);
  request.XLim(1) = request.XLim(1) - growth/2;
  request.XLim(2) = request.XLim(2) + growth/2;
end
request.ImageHeight =  height;
request.ImageWidth  = width;
request.ImageFormat = 'image/jpeg';

%% Store data about the request in class object for use in other functions
obj.proj = regexp(request.CoordRefSysCode,'\d+','match');
obj.proj = str2double(obj.proj{1});
obj.wms = wms;
obj.cur_request = request;
obj.map.projmat = imb.get_proj_info(map_zone);

%% Make WMS Request

modrequest = strcat(request.RequestURL,'&viewparams=');

% create seasons viewparam
if ~isempty(obj.cur_map_pref_settings.seasons)
  season_names = obj.cur_map_pref_settings.seasons;
  
  % convert season_names to a string for concatenation
  for sidx = 1:size(season_names,2)
    if sidx < size(season_names,2) && size(season_names,2) ~= 1
      season_names{sidx}=['''' season_names{sidx} '''%5C,'];
    else
      season_names{sidx}=['''' season_names{sidx} ''''];
    end
  end
  season_names = cell2mat(season_names);
  obj.seasons_as_string = season_names;
  obj.seasons_modrequest = strcat('season_name:',season_names,';');
else
  obj.seasons_modrequest = '';
  obj.seasons_as_string = '';
end

% create season_group_ids viewparam
if ~isempty(obj.map_pref.profile)
  eval(sprintf('season_group_ids = obj.map_pref.profile.%s_season_group_ids'';',obj.cur_map_pref_settings.system))
  
  if isempty(season_group_ids)
    season_group_ids = {'1'};
  end
  
  % convert season_group_ids to a string for concatenation
  for sidx = 1:size(season_group_ids,2)
    if sidx < size(season_group_ids,2) && size(season_group_ids,2) ~= 1
     season_group_ids{sidx}=['' int2str(season_group_ids{sidx}) '%5C,'];
    else
      season_group_ids{sidx}=['' int2str(season_group_ids{sidx}) ''];
    end
  end
  season_group_ids = cell2mat(season_group_ids);
  obj.season_group_ids_as_string = season_group_ids;
  obj.season_group_ids_modrequest = strcat('season_group_ids:',season_group_ids);

else
  obj.season_group_ids_modrequest = '';
  obj.season_group_ids_as_string = '1';
end

% build modified request
modrequest = strcat(modrequest,obj.seasons_modrequest,obj.season_group_ids_modrequest);

% make and post-process request
A = wms.getMap(modrequest);
R = request.RasterRef;
R = R/1e3;

%% Bring the map into focus
figure(obj.h_fig);
xaxis = R(3,1) + [0 size(A,2)*R(2,1)];
yaxis = R(3,2) + [0 size(A,1)*R(1,2)];
set(obj.map_panel.h_image,'XData', xaxis, ...
  'YData', yaxis, ...
  'CData', A, ...
  'Visible', 'on');
set(obj.map_panel.h_axes, 'Xlim', sort(xaxis([1 end])), ...
  'Ylim', sort(yaxis([1 end])), ...
  'YDir', 'normal', ...
  'Visible', 'on');

zoom on; zoom off;

% Redraw table to ensure everything is the right size
table_draw(obj.table);

fprintf('  Done (%s)\n', datestr(now,'HH:MM:SS'));

return;
