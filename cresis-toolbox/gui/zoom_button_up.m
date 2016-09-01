function zoom_button_up(x,y,but,param)
% zoom_button_up(x,y,but,param)
%
% Function to be used in figure's WindowButtonUpFcn callback function.
%
% but,x,y = values returned from:
%   [x,y,but] = get_mouse_info(h_fig,h_axes);
% param
% .x = x-pos of last click from button_down
% .y = y-pos of last click from button_down
% .xlims = max x-limits
% .ylims = max y-limits
% .h_axes = axes handle that we are zooming with
% .axis_equal: logical that forces equal aspect ratio when true
%
% During figure's WindowButtonDownFcn callback function, one should
% generally run:
%   [param.x,param.y,~] = get_mouse_info(h_fig,h_axes);
%   rbbox;
%
% To set the cursor to zoom cursor:
%  zoom_setup(h_fig);
%  set(h_fig,'pointer','custom');
%
% See also: zoom_arrow.m, zoom_button_up.m zoom_button_scroll.m,
%   zoom_setup.m, zoom_figure_setup.m

if isempty(param.x) || isempty(param.y)
  return;
end

if ~isfield(param,'axis_equal') || isempty(param.axis_equal)
  param.axis_equal = 0;
end

if isempty(param.xlims)
  param.xlims = [-inf inf];
end

if isempty(param.ylims)
  param.ylims = [-inf inf];
end

xlims = sort([x param.x]);
ylims = sort([y param.y]);
xlims(1) = max(xlims(1),param.xlims(1));
ylims(1) = max(ylims(1),param.ylims(1));
xlims(2) = min(xlims(2),param.xlims(2));
ylims(2) = min(ylims(2),param.ylims(2));

if param.axis_equal
  cur_units = get(param.h_axes,'Units');
  set(param.h_axes,'Units','pixels');
  axis_pos = get(param.h_axes,'Position');
  set(param.h_axes,'Units',cur_units);
  aspect_ratio = axis_pos(3)/axis_pos(4);
  x_length = diff(xlims);
  y_length = diff(ylims);
  if x_length > y_length*aspect_ratio
    y_mid = mean(ylims);
    ylims(1) = y_mid - x_length/2/aspect_ratio;
    ylims(2) = y_mid + x_length/2/aspect_ratio;
  elseif y_length*aspect_ratio > x_length
    x_mid = mean(xlims);
    xlims(1) = x_mid - y_length/2*aspect_ratio;
    xlims(2) = x_mid + y_length/2*aspect_ratio;
  end
end

if but == 4
  %% Double click: Zoom reset
  if all(isfinite(param.xlims))
    xlim(param.h_axes,param.xlims);
  end
  if all(isfinite(param.ylims))
    ylim(param.h_axes,param.ylims);
  end
  
elseif but == 1 && x~=param.x && y~=param.y
  %% Left click and drag: Zoom to region
  xlim(param.h_axes,xlims);
  ylim(param.h_axes,ylims);
  
elseif but == 1
  %% Left click: Zoom at point
  zooms = -1.5;
  
  cur_axis = [get(param.h_axes,'Xlim') ...
    get(param.h_axes,'YLim')];
  
  % Make sure the first click was within the axes
  if param.x >= cur_axis(1) && param.x <= cur_axis(2) ...
      && param.y >= cur_axis(3) && param.y <= cur_axis(4)
    
    y_extent = cur_axis(4) - cur_axis(3);
    x_extent = cur_axis(2) - cur_axis(1);
    
    % Zoom so that the mouse pointer's position in the echogram does not change
    x_percent = (x-cur_axis(1))/x_extent;
    y_percent = (y-cur_axis(3))/y_extent;
    xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
    ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];
    
    if xlims(1) < param.xlims(1)
      xlims(1) = param.xlims(1);
    end
    if xlims(2) > param.xlims(2)
      xlims(2) = param.xlims(2);
    end
    if ylims(1) < param.ylims(1)
      ylims(1) = param.ylims(1);
    end
    if ylims(2) > param.ylims(2)
      ylims(2) = param.ylims(2);
    end
    xlim(param.h_axes,sort(xlims));
    ylim(param.h_axes,sort(ylims));
  end
  
elseif but == 3
  %% Right click: Zoom out at point
  zooms = -0.5;
  
  cur_axis = [get(param.h_axes,'Xlim') ...
    get(param.h_axes,'YLim')];
  y_extent = cur_axis(4) - cur_axis(3);
  x_extent = cur_axis(2) - cur_axis(1);
  
  % Zoom so that the mouse pointer's position in the echogram does not change
  x_percent = (x-cur_axis(1))/x_extent;
  y_percent = (y-cur_axis(3))/y_extent;
  xlims = [x - x_extent*2^(zooms+1)*x_percent, x + x_extent*2^(zooms+1)*(1-x_percent)];
  ylims = [y - y_extent*2^(zooms+1)*y_percent, y + y_extent*2^(zooms+1)*(1-y_percent)];
  
  if xlims(1) < param.xlims(1)
    xlims(1) = param.xlims(1);
  end
  if xlims(2) > param.xlims(2)
    xlims(2) = param.xlims(2);
  end
  if ylims(1) < param.ylims(1)
    ylims(1) = param.ylims(1);
  end
  if ylims(2) > param.ylims(2)
    ylims(2) = param.ylims(2);
  end
  xlim(param.h_axes,xlims);
  ylim(param.h_axes,ylims);
end

end