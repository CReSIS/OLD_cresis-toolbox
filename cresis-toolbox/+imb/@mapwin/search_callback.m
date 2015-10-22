function search_callback(obj,src,event)
% mapwin.search_callback(obj,src,event)

% return focus to figure
try
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  javaFrame = get(obj.h_fig,'JavaFrame');
  javaFrame.getAxisComponent.requestFocus;
catch
  fprintf('JavaFrame figure property not available, click inside echogram window after pressing a listbox button before using key shortcuts\n');
end

if strcmpi(get(obj.map_panel.h_axes,'Visible'),'off')
  % No map selected, so just return
  return;
end

ops_param.properties.search_str = get(obj.top_panel.searchTB,'String');
ops_param.properties.location = obj.cur_map_pref_settings.mapzone;

[status,data] = opsGetFrameSearch(obj.cur_map_pref_settings.system,ops_param);
if status==2
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
[changed,pos] = obj.compute_new_map_limits(new_xdata,new_ydata);
if changed
    obj.query_redraw_map(pos(1),pos(2),pos(3),pos(4));
end

return;