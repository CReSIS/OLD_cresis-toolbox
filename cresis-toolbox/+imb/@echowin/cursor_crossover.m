function cursor_crossover(obj,source,event)
% echowin.cursor_crossover(obj,source,event)

% Get the current cross over
idx = obj.eg.crossovers.gui.get_selected();

%% Set cursor
% get absolute cursor position
obj.cursor.gps_time = obj.eg.crossovers.gps_time(idx);
obj.cursor.x = interp1(obj.eg.map_gps_time,obj.eg.map_x,obj.eg.crossovers.gps_time(idx),'linear','extrap');
obj.cursor.y = interp1(obj.eg.map_gps_time,obj.eg.map_y,obj.eg.crossovers.gps_time(idx),'linear','extrap');

% (re)draw cursor in this echogram window
obj.plot_cursors();

str = obj.status_text_cursor();
obj.status_text_set(str,'replace');

% tell the map to redraw all cursors
notify(obj,'update_cursors');

end
