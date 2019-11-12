function season_update(obj,status,event)

flightlines_value = get(obj.h_gui.flightlinesPM,'Value');
flightlines_list = get(obj.h_gui.flightlinesPM,'String');
flightlines = flightlines_list(flightlines_value);

systems_value = get(obj.h_gui.systemsLB,'Value');
systems = get(obj.h_gui.systemsLB,'String');
if ~isempty(systems)
  system_name = systems{systems_value};
else
  system_name = [];
end

if strcmp(flightlines{1}(1:3),'OPS')
  obj.unique_systems = setdiff(unique(obj.systems),{'layerdata'});
  set(obj.h_gui.systemsLB,'String',obj.unique_systems);
  system_value = find(strcmp(system_name,obj.unique_systems),1);
  if isempty(system_value)
    system_value = 1;
  end
  if isempty(obj.unique_systems)
    set(obj.h_gui.systemsLB,'Value',1);
    system_name = [];
  else
    set(obj.h_gui.systemsLB,'Value',system_value);
    system_name = obj.unique_systems{system_value};
  end
else
  system_name = 'layerdata';
  set(obj.h_gui.systemsLB,'String',{'layerdata'});
  set(obj.h_gui.systemsLB,'Value',1);
end

map_value = get(obj.h_gui.mapsPM,'Value');
map_list = get(obj.h_gui.mapsPM,'String');
map = map_list(map_value);

[map_zone,map_name] = strtok(map,':');
map_name = map_name(2:end);

zone_mask = strcmp(map_zone,obj.locations);
systems_mask = strcmp(system_name,obj.systems);

obj.h_gui.h_seasons.set_list(obj.seasons(zone_mask & systems_mask));

%% Layers
if strcmp(system_name,'layerdata')
  obj.h_gui.h_layers.set_enable(false);
elseif ~isempty(system_name)
  obj.h_gui.h_layers.set_enable(true);
end
obj.layers_callback_refresh();
