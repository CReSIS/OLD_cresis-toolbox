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

% Update map selection plot
if obj.map_source==1
  % Find the first frame that matches the search string
  frm_id = get(obj.top_panel.searchTB,'String');
  frm_id(regexp(frm_id,'_')) = [];
  frm_id = str2num(frm_id);
  
  % Get a logical mask indicating all indices that match the frame
  frm_mask = obj.google_fline_frms == frm_id;
  if all(frm_mask==false)
    % No frames match, so just return
    return;
  end

  % Generate search string
  frm_id = num2str(frm_id);
  day = frm_id(1:8);
  seg = frm_id(9:10);
  frame = frm_id(11:13);
  search_str = strcat(day,'_',seg,'_',frame);
  
  % Get segment id from opsGetFrameSearch
  frame_search_param = struct('properties',[]);
  frame_search_param.properties.search_str = search_str;
  frame_search_param.properties.location = obj.cur_map_pref_settings.mapzone;
  [a,b] = opsGetFrameSearch(obj.cur_map_pref_settings.system,frame_search_param);
  
  % Record current frame selection
  obj.cur_sel.frame_name = b.properties.frame;
  obj.cur_sel.season_name = b.properties.season;
  obj.cur_sel.segment_id = b.properties.segment_id;
  
  set(obj.map_panel.h_cur_sel,{'XData','YData'},{obj.google_fline_x(frm_mask),obj.google_fline_y(frm_mask)});
  new_xdata = obj.google_fline_x(frm_mask);
  new_ydata = obj.google_fline_y(frm_mask);

  % Update map limits if necessary
  [changed,pos] = obj.compute_new_map_limits(new_xdata,new_ydata);
  if changed
      obj.query_redraw_map(pos(1),pos(2),pos(3),pos(4));
  end

else
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