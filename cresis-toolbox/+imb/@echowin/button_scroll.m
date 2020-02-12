function button_scroll(obj,src,event)
% button_scroll(obj,src, event)
%
% Called by mapwin when a mouse press is released. Handles zoom in/out,
% frame selection, and marker updating.

[x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);

cur_axis = [get(obj.h_axes,'Xlim') ...
  get(obj.h_axes,'YLim')];
if y < cur_axis(3) || y > cur_axis(4) || x < cur_axis(1) || x > cur_axis(2)
  return;
end

%fprintf('Echo scroll: x = %.3f, y = %.3f, but = %d\n', x, y, but);

zooms = -1 + (event.VerticalScrollCount/2);

cur_axis = [get(obj.h_axes,'Xlim') ...
  get(obj.h_axes,'YLim')];
y_extent = cur_axis(4) - cur_axis(3);
x_extent = cur_axis(2) - cur_axis(1);

% Zoom so that the mouse pointer's position in the echogram does not change
x_percent = (x-cur_axis(1))/x_extent;
y_percent = (y-cur_axis(3))/y_extent;
xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];

% Convert x_min, x_max to GPS time
xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,xlims,'linear','extrap');

% Draw data with new axis, but do not allow new data to be loaded (i.e.
% clip new axis to limits of loaded data
obj.redraw(xlims(1),xlims(2),ylims(1),ylims(2),struct('clipped',true,'ylim_force',obj.shift_pressed));

return
