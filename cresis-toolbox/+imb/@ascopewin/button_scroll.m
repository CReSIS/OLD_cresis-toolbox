function button_scroll(obj,h_obj,event)
[x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);

% Make sure that click is on the right side panel
mouse_pos = get(obj.h_fig,'CurrentPoint');

% Check to make sure mouse clicked inside of obj.h_axes
%   Since extends the full y-length, just check to the right of minimum x
uipanel_pos = get(obj.left_panel.handle,'Position');
if mouse_pos(1) <= uipanel_pos(1)
  return
end

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

if xlims(1) < obj.xlims(1)
  xlims(1) = obj.xlims(1);
end
if xlims(2) > obj.xlims(2)
  xlims(2) = obj.xlims(2);
end
if ylims(1) < obj.ylims(1)
  ylims(1) = obj.ylims(1);
end
if ylims(2) > obj.ylims(2)
  ylims(2) = obj.ylims(2);
end
xlim(obj.h_axes,xlims);
ylim(obj.h_axes,ylims);

end