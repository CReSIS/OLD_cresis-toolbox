function update_map_selection(obj,param)
% update_map_selection(obj,param)
%
% Updates the currently selected (obj.map.sel) frame for imb.mapwin class.
% Called from mapwin.button_up.
%
% param.x, param.y = scalars defining search point for new selection in km

%% Query database to find the closest frame to param.x,param.y

% Setup query
ops_param = struct('properties',[]);
ops_param.properties.location = obj.cur_map_pref_settings.map_zone;
ops_param.properties.season = obj.cur_map_pref_settings.seasons;

% Query
ops_param.properties.x = param.x*obj.map.scale;
ops_param.properties.y = param.y*obj.map.scale;
% Get closest frame
[status,data] = obj.get_closest_frame(obj.cur_map_pref_settings.system,ops_param);

% Update map selection plot
set(obj.map_panel.h_cur_sel,{'XData','YData'},{data.properties.X/obj.map.scale,data.properties.Y/obj.map.scale});

% Record current frame selection
obj.map.sel.frame_name = data.properties.frame;
obj.map.sel.season_name = data.properties.season;
obj.map.sel.segment_id = data.properties.segment_id;


% Change map title to the currently selected frame
set(obj.top_panel.flightLabel,'String',obj.map.sel.frame_name);
