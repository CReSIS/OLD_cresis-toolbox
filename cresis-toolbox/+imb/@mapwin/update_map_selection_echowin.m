function update_map_selection_echowin(obj,src,event)
% update_map_selection_echowin(obj,src,event)
%
% Updates the currently selected (cur_sel) frame for imb.mapwin class.
% Called anytime an echowin "update_map_selection" event occurs.
% This function is similar to update_map_selection.

echowin_idx = find(obj.echowin_list == src);
frames = get(obj.echowin_list(echowin_idx).left_panel.frameLB,'String');
frame_name = frames{get(obj.echowin_list(echowin_idx).left_panel.frameLB,'Value')};

ops_param.properties.search_str = frame_name;
ops_param.properties.season = obj.echowin_list(echowin_idx).eg.cur_sel.season_name;
ops_param.properties.location = obj.cur_map_pref_settings.mapzone;

[status,data] = opsGetFrameSearch(obj.cur_map_pref_settings.system,ops_param);
if status==2 || ~status
  % result not found; warning already printed to console, so just exit
  return;
end

% Record current frame selection
obj.cur_sel.frame_name = data.properties.frame;
obj.cur_sel.season_name = data.properties.season;
obj.cur_sel.segment_id = data.properties.segment_id;

% Update map selection plot
set(obj.map_panel.h_cur_sel,{'XData','YData'},{data.properties.X/1e3,data.properties.Y/1e3});

% Change map title to the currently selected frame
set(obj.top_panel.flightLabel,'String',obj.cur_sel.frame_name);

new_xdata = data.properties.X/1e3;
new_ydata = data.properties.Y/1e3;

% Update map limits if necessary
if get(obj.top_panel.trackCB,'Value')
  [changed,pos] = obj.compute_new_map_limits(new_xdata,new_ydata);
  if changed
    obj.query_redraw_map(pos(1),pos(2),pos(3),pos(4));
  end
end

return;