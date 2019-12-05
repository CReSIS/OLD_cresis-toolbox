function set_cursor_by_map(obj,x,y)
% echowin.set_cursor_by_map(obj,x,y)
%
% Finds the closest point on the echowin's flightline to (x,y) and then
% updates the cursor to that position.

[~,min_idx] = min((obj.eg.map_x-x).^2 + (obj.eg.map_y-y).^2);

obj.cursor.gps_time = interp1(1:length(obj.eg.map_gps_time),obj.eg.map_gps_time,...
  min_idx,'linear','extrap');
obj.cursor.x = interp1(obj.eg.map_gps_time,obj.eg.map_x,obj.cursor.gps_time,'linear','extrap');
obj.cursor.y = interp1(obj.eg.map_gps_time,obj.eg.map_y,obj.cursor.gps_time,'linear','extrap');

% (re)draw cursor in this echogram window
obj.plot_cursors();

str = obj.status_text_cursor();
obj.status_text_set(str,'replace');
