function button_scroll(obj,src,event)
% mapwin.button_scroll(obj,src, event)
%
% Called by echowin when the mouse scroll button is used. Handles zoom in/out.

if strcmpi(get(obj.map_panel.h_axes,'Visible'),'on')
  [x,y,but] = get_mouse_info(obj.h_fig,obj.map_panel.h_axes);
  cur_axis = axis(obj.map_panel.h_axes);
  if y < cur_axis(3) || y > cur_axis(4) || x < cur_axis(1) || x > cur_axis(2)
    return;
  end
  
  %fprintf('Map scroll: x = %.3f, y = %.3f, but = %d\n', x, y, but);
  
  zooms = -1 + (event.VerticalScrollCount/2);
  
  y_extent = cur_axis(4) - cur_axis(3);
  x_extent = cur_axis(2) - cur_axis(1);
  
  % Zoom so that the mouse pointer's position in the echogram does not change
  x_percent = (x-cur_axis(1))/x_extent;
  y_percent = (y-cur_axis(3))/y_extent;
  xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
  ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];
  
  % Draw map with new axis
  obj.query_redraw_map(xlims(1),xlims(2),ylims(1),ylims(2));
end

return
