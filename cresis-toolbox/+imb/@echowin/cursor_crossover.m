function cursor_crossover(obj,source,event)
% echowin.cursor_crossover(obj,source,event)

% Get the current cross over
idx = obj.crossovers.gui.get_selected();

% Convert GPS time to current image x-axis units and call update_cursor
[~,rline] = min(abs(obj.eg.image_gps_time-obj.crossovers.gps_time(idx)));

% CONVERT twtt to y
y = obj.crossovers.twtt(idx);
physical_constants;
vel_air = c/2;
vel_ice = c/2/sqrt(er_ice);
yaxis_choice = get(obj.left_panel.yaxisPM,'Value');
if yaxis_choice == 1 % TWTT
  y = y*1e6;
elseif yaxis_choice == 2 % WGS_84 Elevation
  y = obj.eg.image_elev - min(y,obj.eg.image_surf_twtt(rline))*vel_air - min(0,y-obj.eg.image_surf_twtt(rline))*vel_ice;
elseif yaxis_choice == 3 % Range
  y = min(y,obj.eg.image_surf_twtt(rline))*vel_air + min(0,y-obj.eg.image_surf_twtt(rline))*vel_ice;
elseif yaxis_choice == 4 % Range bin
  y = (y-obj.eg.time(1)) / (obj.eg.time(2)-obj.eg.time(1));
elseif yaxis_choice == 5 % Surface flat
  y = min(0,obj.eg.image_surf_twtt(rline)-y)*vel_air - min(0,y-obj.eg.image_surf_twtt(rline))*vel_ice;
end

obj.update_cursor(obj.eg.image_xaxis(rline),y,true);

end
