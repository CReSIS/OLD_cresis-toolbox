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
if obj.isGoogle
  ops_param.properties.x = param.x;
  ops_param.properties.y = param.y;
  
  % Get closest frame
  [status,data] = imb.googleGetFrameClosest(obj.cur_map_pref_settings.system,ops_param);
  
  % Update seg plot
  flightline_plot = get(gca,'Children');
  for graphics_obj_idx = 1:length(flightline_plot)
    if strcmp(flightline_plot(graphics_obj_idx).Tag, 'seg')
      set(flightline_plot(graphics_obj_idx), 'XData', [], 'YData', [])
      set(flightline_plot(graphics_obj_idx), 'XData', data.properties.X, 'YData', data.properties.Y)
    end
  end

else
  ops_param.properties.x = param.x*1e3;
  ops_param.properties.y = param.y*1e3;
  [status,data] = opsGetFrameClosest(obj.cur_map_pref_settings.system,ops_param);

  % Update map selection plot
  set(obj.map_panel.h_cur_sel,{'XData','YData'},{data.properties.X/1e3,data.properties.Y/1e3});
end
% Record current frame selection
obj.cur_sel.frame_name = data.properties.frame;
obj.cur_sel.season_name = data.properties.season;
obj.cur_sel.segment_id = data.properties.segment_id;


% Change map title to the currently selected frame
set(obj.top_panel.flightLabel,'String',obj.cur_sel.frame_name);

return;
