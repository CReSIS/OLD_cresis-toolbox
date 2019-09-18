function okPB_callback(obj,hObj,event)

cur_unit = get(obj.h_fig,'Units');
set(obj.h_fig,'Units','pixels')
prefwin_pos = get(obj.h_fig,'Position');
set(obj.h_fig,'Unit',cur_unit);
obj.default_params.x = prefwin_pos(1);
obj.default_params.y = prefwin_pos(2);
obj.default_params.w = prefwin_pos(3);
obj.default_params.h = prefwin_pos(4);
obj.default_params.sources = get(obj.h_gui.sourceLB,'String');
map_names = get(obj.h_gui.mapsPM,'String');
flightlines = get(obj.h_gui.flightlinesPM,'String');
obj.default_params.map_name = map_names{get(obj.h_gui.mapsPM,'Value')};
obj.default_params.flightlines = flightlines{get(obj.h_gui.flightlinesPM,'Value')};
systems = get(obj.h_gui.systemsLB,'String');
obj.default_params.system = systems{get(obj.h_gui.systemsLB,'Value')};
%
LayerSource = get(obj.h_gui.LayerSourcePM,'String');
obj.default_params.LayerSource = LayerSource{get(obj.h_gui.LayerSourcePM,'Value')};
layerDataSource = get(obj.h_gui.layerDataSourcePM,'String');
obj.default_params.layerDataSource = layerDataSource{get(obj.h_gui.layerDataSourcePM,'Value')};
%
obj.default_params.season_names = obj.h_gui.seasons.get_selected_strings();
obj.default_params.layer_names = obj.h_gui.layers.get_selected_strings();

%% Filter layers
selected_layer_names = obj.h_gui.layers.get_selected_strings().';
selected_layers.lyr_name = {};
selected_layers.lyr_group_name = {};
selected_layers.lyr_id = [];
for idx = 1:length(obj.layers.lyr_name)
  layer_name = sprintf('%s:%s', obj.layers.lyr_group_name{idx}, obj.layers.lyr_name{idx});
  match_idx = strmatch(layer_name,selected_layer_names,'exact');
  if ~isempty(match_idx)
    selected_layers.lyr_name{end+1} = obj.layers.lyr_name{idx};
    selected_layers.lyr_group_name{end+1} = obj.layers.lyr_group_name{idx};
    selected_layers.lyr_id(end+1) = obj.layers.lyr_id(idx);
  end
end

%% Check if anything has been selected
selected_seasons = obj.h_gui.seasons.get_selected_strings().';
if isempty(selected_seasons)
  % if not, don't plot a map
  fprintf('Please select one or more seasons before plotting a map.\n');
  return;
end

%% Filter available maps based on selected system and season
systems = get(obj.h_gui.systemsLB,'string');
system = systems{get(obj.h_gui.systemsLB,'value')};
selected_idxs = zeros(size(selected_seasons));
for idx = 1:length(selected_seasons)
  % For each season selected, find its index in obj.seasons
  found = false;
  for search_idx = 1:length(obj.seasons)
    if strcmp(obj.seasons{search_idx},selected_seasons{idx}) ...
        && strcmp(obj.systems{search_idx},system)
      selected_idxs(idx) = search_idx;
      found = true;
      break;
    end
  end
  if ~found
    warning('This should never happen, selection was not found in the list');
    keyboard
  end
end

%% Pass information to map window
% Give information on the season, location (arctic/antarctic) and system
% (rds accum etc.) to the map window
locs = obj.locations(selected_idxs);
if all(strcmp(locs{1},locs));
  obj.settings.sources = get(obj.h_gui.sourceLB,'String');
  map_names = get(obj.h_gui.mapsPM,'String');
  map_name = map_names{get(obj.h_gui.mapsPM,'Value')};
  
  flightlines = get(obj.h_gui.flightlinesPM,'String');
  flightline = flightlines{get(obj.h_gui.flightlinesPM,'Value')};
   %
  layerDataSources = get(obj.h_gui.layerDataSourcePM,'String');
  layerDataSource = layerDataSources{get(obj.h_gui.layerDataSourcePM,'Value')};
  LayerSources = get(obj.h_gui.LayerSourcePM,'String');
  LayerSource = LayerSources{get(obj.h_gui.LayerSourcePM,'Value')};
  %
  % Parse map name from GUI into 'zone:map'
  [obj.settings.mapzone obj.settings.mapname] = strtok(map_name,':');
  obj.settings.mapname = obj.settings.mapname(2:end);
  obj.settings.system = obj.systems(selected_idxs);
  obj.settings.system = obj.settings.system{1};
  obj.settings.seasons = selected_seasons;
  obj.settings.layers = selected_layers;
  obj.settings.flightlines = flightline;
   %
  obj.settings.layerDataSource = layerDataSource;
  obj.settings.LayerSource = LayerSource;
  %
  % Broadcast notice that StateChange event has occurred (plot map)
  notify(obj,'StateChange');
else
  fprintf('Cannot load multiple zones at once. Select only seasons from a single zone and retry.\n');
end

return
