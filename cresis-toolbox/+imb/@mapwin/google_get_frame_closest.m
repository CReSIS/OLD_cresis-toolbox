function [status, data] = google_get_frame_closest(obj, sys, param)
  
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
    data.properties.segment_id = [];
    data.properties.X = obj.layerdata.x(frm_mask);
    data.properties.Y = obj.layerdata.y(frm_mask);
    
  else
    % Get segment id from opsGetFrameSearch
    frame_search_param = struct('properties',[]);
    frame_search_param.properties.frm_str = frm_str;
    frame_search_param.properties.location = param.properties.location;
    [frm_status,frm_data] = opsGetFrameSearch(sys,frame_search_param);
    
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