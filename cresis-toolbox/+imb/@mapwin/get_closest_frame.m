function [status, data] = get_closest_frame(obj, sys, param)

if obj.map.fline_source == 1
  
  min_diff = inf;
  min_frame = 0;
  lats = [];
  lons = [];
  frms = [];
  
  % Mouse click coordinates
  wc_x = param.properties.x;
  wc_y = param.properties.y;
  
  % Find the frame of the closest point
  [~,idx] = min((obj.layerdata.y-wc_y).^2+(obj.layerdata.x-wc_x).^2);
  frm_id = obj.layerdata.frms(idx);
  season_idx = obj.layerdata.season_idx(idx);
  season_name = obj.cur_map_pref_settings.seasons{season_idx};
  [sys,season_name] = strtok(season_name,'_');
  season_name = season_name(2:end);
  
  % Get a logical mask indicating all indices that match the frame
  frm_mask = obj.layerdata.frms == frm_id;
  
  % Generate frame string YYYYMMDD_SS_FFF
  frm_id = num2str(frm_id);
  day = frm_id(1:8);
  seg = frm_id(9:10);
  frm = frm_id(11:13);
  frm_str = strcat(day,'_',seg,'_',frm);
  
  if strcmpi(obj.cur_map_pref_settings.layer_source,'layerdata')
    status = [];
    
    % Set data properties
    data = struct('properties',[]);
    data.properties.frame = frm_str;
    data.properties.season = obj.cur_map_pref_settings.seasons{season_idx};
    data.properties.segment_id = str2num(frm_id(1:10));
    data.properties.X = obj.layerdata.x(frm_mask);
    data.properties.Y = obj.layerdata.y(frm_mask);
    
  else
    % Get segment id from opsGetFrameSearch
    frame_search_param = struct('properties',[]);
    frame_search_param.properties.search_str = frm_str;
    frame_search_param.properties.location = param.properties.location;
    frame_search_param.properties.season = season_name;
    [frm_status,frm_data] = opsGetFrameSearch(sys,frame_search_param);
    
    if frm_status ~= 1
      error_str = sprintf('Frame %s does not exist in OPS for %s:%s.', frm_str, sys, season_name);
      uiwait(msgbox(error_str,'Search error','modal'));
      error(error_str);
    else
      % Set data properties
      data = struct('properties',[]);
      data.properties.frame = frm_str;
      data.properties.season = frm_data.properties.season;
      data.properties.segment_id = frm_data.properties.segment_id;
      data.properties.X = obj.layerdata.x(frm_mask);
      data.properties.Y = obj.layerdata.y(frm_mask);
      status = frm_status;
    end
  end
  
else
  % OPS Flightline
  if obj.map.source == 1
    [lat,lon] = google_map.world_to_latlon(param.properties.x,256-param.properties.y);
    [param.properties.x,param.properties.y] = projfwd(obj.map.proj,lat,lon);
    [status,data] = opsGetFrameClosest(obj.cur_map_pref_settings.system,param);
    [lat,lon] = projinv(obj.map.proj,data.properties.X,data.properties.Y);
    [data.properties.X,data.properties.Y] = google_map.latlon_to_world(lat,lon);
    data.properties.Y = 256-data.properties.Y;
  else
    [status,data] = opsGetFrameClosest(obj.cur_map_pref_settings.system,param);
  end
end
