function load_flightline(obj)
% echowin.load_flightline(obj)
%
% Load flightline from database

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

return;
