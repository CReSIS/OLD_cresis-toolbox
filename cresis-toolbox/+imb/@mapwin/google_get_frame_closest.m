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
  [~,idx] = min((obj.google_fline_y-wc_y).^2+(obj.google_fline_x-wc_x).^2);
  frm_id = obj.google_fline_frms(idx);
  
  % Get a logical mask indicating all indices that match the frame
  frm_mask = obj.google_fline_frms == frm_id;

  % Generate search string
  frm_id = num2str(frm_id);
  day = frm_id(1:8);
  seg = frm_id(9:10);
  frame = frm_id(11:13);
  search_str = strcat(day,'_',seg,'_',frame);
  
  % Get segment id from opsGetFrameSearch
  frame_search_param = struct('properties',[]);
  frame_search_param.properties.search_str = search_str;
  frame_search_param.properties.location = param.properties.location;
  [a,b] = opsGetFrameSearch(sys,frame_search_param);
  
  % Set data properties
  data = struct('properties',[]);
  data.properties.frame = search_str;
  data.properties.season = b.properties.season;
  data.properties.segment_id = b.properties.segment_id;
  data.properties.X = obj.google_fline_x(frm_mask);
  data.properties.Y = obj.google_fline_y(frm_mask);
  
  status = [];
end