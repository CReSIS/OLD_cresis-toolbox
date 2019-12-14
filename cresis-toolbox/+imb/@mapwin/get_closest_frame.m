function get_closest_frame(obj, param)
% get_closest_frame(obj, param)
%
% Updates the currently selected (obj.map.sel) frame for imb.mapwin class
% by finding the closest frame from where the user clicked. Called from
% mapwin.button_up.
%
% param.x, param.y = scalars defining search point for new selection in km

%% Query database to find the closest frame to param.x,param.y
ops_param.properties.season = obj.cur_map_pref_settings.seasons;

if obj.map.fline_source == 1
  % layerdata flineslines selected
  % -----------------------------------------------------------------------
  
  % Find the frame of the closest point
  [~,idx] = min((obj.layerdata.y-param.y).^2+(obj.layerdata.x-param.x).^2);
  frm_id = obj.layerdata.frms(idx);
  season_idx = obj.layerdata.season_idx(idx);
  season_name = obj.cur_map_pref_settings.seasons{season_idx};
  [sys,season_name_short] = strtok(season_name,'_');
  season_name_short = season_name_short(2:end);
  
  % Get a logical mask indicating all indices that match the frame
  frm_mask = obj.layerdata.frms == frm_id;
  
  % Generate frame string YYYYMMDD_SS_FFF
  frm_id = num2str(frm_id);
  day = frm_id(1:8);
  seg = frm_id(9:10);
  frm = frm_id(11:13);
  frame_name = strcat(day,'_',seg,'_',frm);
  
  if strcmpi(obj.cur_map_pref_settings.layer_source,'layerdata')
    status = [];
    
    % Set data properties
    data = struct('properties',[]);
    data.properties.frame = frame_name;
    data.properties.season = season_name;
    data.properties.segment_id = str2num(frm_id(1:10));
    data.properties.X = obj.layerdata.x(frm_mask);
    data.properties.Y = obj.layerdata.y(frm_mask);
    
  else
    % Get segment id from opsGetFrameSearch
    frame_search_param = struct('properties',[]);
    frame_search_param.properties.search_str = frame_name;
    frame_search_param.properties.location = obj.cur_map_pref_settings.map_zone;
    frame_search_param.properties.season = season_name_short;
    [frm_status,frm_data] = opsGetFrameSearch(sys,frame_search_param);
    
    if frm_status ~= 1
      error_str = sprintf('Frame %s does not exist in OPS for %s:%s.', frame_name, sys, season_name_short);
      uiwait(msgbox(error_str,'Search error','modal'));
      error(error_str);
    else
      % Set data properties
      data = struct('properties',[]);
      data.properties.frame = frame_name;
      data.properties.season = frm_data.properties.season;
      data.properties.segment_id = frm_data.properties.segment_id;
      data.properties.X = obj.layerdata.x(frm_mask);
      data.properties.Y = obj.layerdata.y(frm_mask);
      status = frm_status;
    end
  end
  
else
  % OPS flineslines selected
  % -----------------------------------------------------------------------
  sys = obj.cur_map_pref_settings.system;
  frame_search_param.properties.x = param.x*obj.map.scale;
  frame_search_param.properties.y = param.y*obj.map.scale;
  frame_search_param.properties.location = obj.cur_map_pref_settings.map_zone;
  frame_search_param.properties.season = obj.cur_map_pref_settings.seasons;
  if obj.map.source == 1
    [lat,lon] = google_map.world_to_latlon(frame_search_param.properties.x,256-frame_search_param.properties.y);
    [frame_search_param.properties.x,frame_search_param.properties.y] = projfwd(obj.map.proj,lat,lon);
    [status,data] = opsGetFrameClosest(sys,frame_search_param);
    [lat,lon] = projinv(obj.map.proj,data.properties.X,data.properties.Y);
    [data.properties.X,data.properties.Y] = google_map.latlon_to_world(lat,lon);
    data.properties.Y = 256-data.properties.Y;
  else
    [status,data] = opsGetFrameClosest(sys,frame_search_param);
  end
  data.properties.X = data.properties.X/obj.map.scale;
  data.properties.Y = data.properties.Y/obj.map.scale;
end

% Record current frame selection
obj.map.sel.frame_name = data.properties.frame;
obj.map.sel.season_name = data.properties.season;
obj.map.sel.segment_id = data.properties.segment_id;
obj.map.sel.radar_name = sys;

% Update map selection plot
set(obj.map_panel.h_cur_sel,{'XData','YData'},{data.properties.X,data.properties.Y});

% Change map title to the currently selected frame
set(obj.top_panel.flightLabel,'String',obj.map.sel.frame_name);
