function [status, data] = googleGetFrameClosest(sys, param)
  
  min_diff = 999999;
  min_frame = 0;
  lats = [];
  lons = [];
  frms = [];
  
  for season_idx = length(param.properties.season)
  
    % Load season layerdata files for selected seasons
    S = load(char(strcat('X:csarp_support\season_layerdata_files\',sys,'_param_',param.properties.season{season_idx},'_layerdata.mat')));  

    % Mouse click coordinates
    wc_x = param.properties.x;
    wc_y = param.properties.y;
    
    % Convert to lat lon
    [lat, lon] = imb.world_to_latlon(wc_x, wc_y);
    
    % Find closest lat lon pair
    [diff index] = min(abs(S.lat-lat)+abs(S.lon-lon));

    % Update closest lat lon and frame
    if(diff < min_diff)
      min_frame = S.frm(index);
      min_diff = diff;
      lats = S.lat;
      lons = S.lon;
      frms = S.frm;
    end
  end
  
  frm = min_frame;
  
  % Get all indices having that frame value
  frm_idx = ismember(frms, frm);
  frm_idx = find(frm_idx);

  % Generate search string
  frm = num2str(frm);
  day = frm(1:8);
  seg = frm(9:10);
  frame = frm(11:13);
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
  [data.properties.X, data.properties.Y] = imb.latlon_to_world(lats(frm_idx), lons(frm_idx));
  
  status = [];
end