function update_map_selection(obj,param)
% update_map_selection(obj,param)
%
% Updates the currently selected (cur_sel) frame for imb.mapwin class.
% Called from mapwin.button_up.
%
% param.x, param.y = scalars defining search point for new selection in km

%% Query database to find the closest frame to param.x,param.y

% Setup query
ops_param = struct('properties',[]);
ops_param.properties.location = obj.cur_map_pref_settings.mapzone;
ops_param.properties.season = obj.cur_map_pref_settings.seasons;

% Query
ops_param.properties.x = param.x*obj.map_scale;
ops_param.properties.y = param.y*obj.map_scale;
if obj.map_source == 1
  % Get closest frame
  [status,data] = obj.google_get_frame_closest(obj.cur_map_pref_settings.system,ops_param);
else
  [status,data] = opsGetFrameClosest(obj.cur_map_pref_settings.system,ops_param);
end

% Update map selection plot
set(obj.map_panel.h_cur_sel,{'XData','YData'},{data.properties.X/obj.map_scale,data.properties.Y/obj.map_scale});

% Record current frame selection
obj.cur_sel.frame_name = data.properties.frame;
obj.cur_sel.season_name = data.properties.season;
obj.cur_sel.segment_id = data.properties.segment_id;


% Change map title to the currently selected frame
set(obj.top_panel.flightLabel,'String',obj.cur_sel.frame_name);
