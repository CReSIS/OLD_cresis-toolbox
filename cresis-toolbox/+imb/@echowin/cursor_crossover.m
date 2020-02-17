function cursor_crossover(obj,source,event)
% echowin.cursor_crossover(obj,source,event)

% Get the current cross over
idx = obj.crossovers.gui.get_selected();

% Convert GPS time to current image x-axis units and call update_cursor
[~,rline] = min(abs(obj.eg.image_gps_time-obj.crossovers.gps_time(idx)));
imb.echowin.update_cursor(obj.eg.image_xaxis(rline));

end
