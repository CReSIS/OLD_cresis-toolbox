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
ops_param.properties.season = obj.map_pref.settings.seasons;
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
if obj.isGoogle
  % 1. Load file
  S = load(char(strcat('X:csarp_support\season_layerdata_files\',obj.cur_map_pref_settings.system, '_param_',obj.cur_sel.season_name,'_layerdata.mat')));
  % 2. Convert latlon to world
  [wc_x, wc_y] = imb.latlon_to_world(S.lat, S.lon);
  lat = S.lat;
  lon = S.lon;
  
  % 3. Get world coordinates for the frame
  frm = obj.cur_sel.frame_name;
  frm(regexp(frm,'_')) = [];
  frm = str2double(frm);
  
  frm_idx = ismember(S.frm, frm);
  frm_idx = find(frm_idx);
  data.properties.X = wc_x(frm_idx);
  data.properties.Y = wc_y(frm_idx);
  flightline_plot = get(gca,'Children');
  for graphics_obj_idx = 1:length(flightline_plot)
    if strcmp(flightline_plot(graphics_obj_idx).Tag, 'seg')
      set(flightline_plot(graphics_obj_idx), 'XData', [], 'YData', [])
      set(flightline_plot(graphics_obj_idx), 'XData', data.properties.X, 'YData', data.properties.Y)
    end
  end
  obj.googleObj.c_lat = lat(1);
  obj.googleObj.c_lon = lon(1);
  redraw_google_map(obj, 0, 0, 0, 0);
else
  set(obj.map_panel.h_cur_sel,{'XData','YData'},{data.properties.X/1e3,data.properties.Y/1e3});
  new_xdata = data.properties.X/1e3;
  new_ydata = data.properties.Y/1e3;

  % Update map limits if necessary
  [changed,pos] = obj.compute_new_map_limits(new_xdata,new_ydata);
  if changed
      obj.query_redraw_map(pos(1),pos(2),pos(3),pos(4));
  end
end
% Change map title to the currently selected frame
set(obj.top_panel.flightLabel,'String',obj.cur_sel.frame_name);


return;