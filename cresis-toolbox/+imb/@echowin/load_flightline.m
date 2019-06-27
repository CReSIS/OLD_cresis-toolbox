function load_flightline(obj)
% echowin.load_flightline(obj)
%
% Load flightline from database
%% OPS: Loading flight path from database
if strcmpi(obj.eg.LayerSource,'OPS')
  fprintf(' Loading flight path from database (%s)\n', datestr(now,'HH:MM:SS'));
  ops_param = struct('properties',[]);
  ops_param.properties.location = obj.eg.cur_sel.location;
  ops_param.properties.season = obj.eg.cur_sel.season_name;
  ops_param.properties.start_gps_time = obj.eg.start_gps_time(obj.eg.frame_idxs(1));
  ops_param.properties.stop_gps_time = obj.eg.stop_gps_time(obj.eg.frame_idxs(end));
  [status,data] = opsGetPath(obj.eg.system,ops_param);
  
  obj.eg.map_id = double(data.properties.id);
  obj.eg.map_gps_time = double(data.properties.gps_time);
  obj.eg.map_elev = double(data.properties.elev);
  obj.eg.map_x = double(data.properties.X)/1e3;
  obj.eg.map_y = double(data.properties.Y)/1e3;
  
%% LayerData: Loading flight path from layerData
else
  fprintf(' Loading flight path from layerData(%s)\n', datestr(now,'HH:MM:SS'));
  for idx = 1:length(obj.eg.frame_idxs)
    obj.eg.map_gps_time = double(obj.undo_stack.user_data.layer_info(obj.eg.frame_idxs(idx)).GPS_time);
    obj.eg.map_elev = double(obj.undo_stack.user_data.layer_info(obj.eg.frame_idxs(idx)).Elevation);
    [X,Y] = projfwd(obj.eg.projmat,obj.undo_stack.user_data.layer_info(obj.eg.frame_idxs(idx)).Latitude,obj.undo_stack.user_data.layer_info(obj.eg.frame_idxs(idx)).Longitude);
    obj.eg.map_x = double(X)/1e3;
    obj.eg.map_y = double(Y)/1e3;
  end
  
  % Get the unique point path ids for layerData and save it in the
  % point_path_id field of the undo_stack
  LDpoint_path_id=[];
  for frm = obj.eg.frame_idxs(1):obj.eg.frame_idxs(end)
    k = find(obj.undo_stack.user_data.frame==frm);
    points = (k(1):k(end));
    LDpoint_path_id = cat(2,LDpoint_path_id,points);
  end
  obj.undo_stack.user_data.point_path_id = LDpoint_path_id; %contains the unquie point path ids.
end

return;
