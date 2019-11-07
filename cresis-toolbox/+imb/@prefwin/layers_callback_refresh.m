function layers_callback_refresh(obj)

% Get the system that is currently selected
system_name = get(obj.h_gui.systemsLB,'String');
if isempty(system_name)
  system_name = [];
else
  system_name = system_name{get(obj.h_gui.systemsLB,'Value')};
end

if ~isempty(system_name) && ~strcmpi(system_name,'layerdata')
  % Find all the layers associated with this system
  [status,data] = opsGetLayers(system_name);
  obj.layers = data.properties;
  menuString = {};
  for idx = 1:length(data.properties.lyr_name)
    menuString{idx} = sprintf('%s:%s', data.properties.lyr_group_name{idx}, data.properties.lyr_name{idx});
  end
  
  % Get the current list of selected layers
  selectedString = obj.h_gui.h_layers.get_selected_strings();
  
  % Create the new list of layers for the system that is selected
  obj.h_gui.h_layers.set_available(menuString);
  
  % Select all the old layers that were selected (this function
  % automatically ignores layers that do not exist in the new list)
  obj.h_gui.h_layers.set_selected(selectedString,true);

else
  obj.ops.layers = [];
  obj.ops.layers.lyr_name = [];
end
