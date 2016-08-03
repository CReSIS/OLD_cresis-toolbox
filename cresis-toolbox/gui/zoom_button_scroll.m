function zoom_button_scroll(event,param)
% event = event passed into figure's WindowScrollWheelFcn callback function
% param
%  .xlims = max x-limits
%  .ylims = max y-limits
%  .h_axes = axes handle that we are zooming with
%  .h_fig = figure handle that we are zooming in

[x,y,but] = get_mouse_info(param.h_fig,param.h_axes);

zooms = -1 + (event.VerticalScrollCount/2);

cur_axis = [get(param.h_axes,'Xlim') ...
  get(param.h_axes,'YLim')];
y_extent = cur_axis(4) - cur_axis(3);
x_extent = cur_axis(2) - cur_axis(1);

% Zoom so that the mouse pointer's position in the echogram does not change
x_percent = (x-cur_axis(1))/x_extent;
y_percent = (y-cur_axis(3))/y_extent;
xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];

if ~isempty(param.xlims)
  if xlims(1) < param.xlims(1)
    xlims(1) = param.xlims(1);
  end
  if xlims(2) > param.xlims(2)
    xlims(2) = param.xlims(2);
  end
end
if ~isempty(param.ylims)
  if ylims(1) < param.ylims(1)
    ylims(1) = param.ylims(1);
  end
  if ylims(2) > param.ylims(2)
    ylims(2) = param.ylims(2);
  end
end
xlim(param.h_axes,xlims);
ylim(param.h_axes,ylims);

end