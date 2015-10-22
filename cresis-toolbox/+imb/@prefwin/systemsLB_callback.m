function systemsLB_callback(obj,status,event)
% systemsLB_callback(obj,status,event)

% Get the system that is currently selected
system_name = get(obj.h_gui.systemsLB,'String');
system_name = system_name{get(obj.h_gui.systemsLB,'Value')};

%% Seasons
% Find all seasons associated with that system
systems_mask = strcmp(system_name,obj.systems);
menuString = obj.seasons(systems_mask);

% Get the current list of selected seasons
selectedString = obj.h_gui.seasons.get_selected_strings();

% Create the new list of seasons for the system that is selected
obj.h_gui.seasons.set_available(menuString);

% Select all the old seasons that were selected (this function
% automatically ignores seasons that do not exist in the new list)
obj.h_gui.seasons.set_selected(selectedString,true);

%% Layers
obj.layers_callback_refresh();

return
