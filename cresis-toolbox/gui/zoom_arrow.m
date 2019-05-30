function zoom_arrow(event,param)
% zoom_arrow(event,param)
%
% event = event passed into figure's WindowKeyPressFcn callback function
% param
%  .xlims = max x-limits
%  .ylims = max y-limits
%  .h_axes = axes handle that we are zooming with
%  .h_fig = figure handle that we are zooming in
%
% See also: zoom_arrow.m, zoom_button_motion.m, zoom_button_scroll.m,
%   zoom_button_up.m, zoom_figure_setup.m, zoom_setup.m

if any(strcmp('ctrl',event.Modifier))
  ctrl_pressed = true;
else
  ctrl_pressed = false;
end

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~isempty(event.Key)
  
  cur_axis = [get(param.h_axes,'Xlim') ...
    get(param.h_axes,'YLim')];
  x_extent = cur_axis(2) - cur_axis(1);
  y_extent = cur_axis(4) - cur_axis(3);
  
  % see event.Modifier for modifiers
  switch event.Key
    
    case 'downarrow' % Down-arrow: Move Echogram Down
      if strcmp(get(param.h_axes,'YDir'),'reverse')
        ylims = cur_axis(3:4) + y_extent*0.25;
      else
        ylims = cur_axis(3:4) - y_extent*0.25;
      end
      xlims = cur_axis(1:2);
      
    case 'uparrow' % Up-arrow: Move Echogram Up
      if strcmp(get(param.h_axes,'YDir'),'reverse')
        ylims = cur_axis(3:4) - y_extent*0.25;
      else
        ylims = cur_axis(3:4) + y_extent*0.25;
      end
      xlims = cur_axis(1:2);
      
    case 'rightarrow' % Right arrow
      if strcmp(get(param.h_axes,'XDir'),'reverse')
        xlims = cur_axis(1:2) - x_extent*0.25;
      else
        xlims = cur_axis(1:2) + x_extent*0.25;
      end
      ylims = cur_axis(3:4);
      
    case 'leftarrow' % Left arrow
      if strcmp(get(param.h_axes,'XDir'),'reverse')
        xlims = cur_axis(1:2) + x_extent*0.25;
      else
        xlims = cur_axis(1:2) - x_extent*0.25;
      end
      ylims = cur_axis(3:4);
      
    otherwise
      return;
  end
  
  if xlims(1) < param.xlims(1)
    xlims = param.xlims(1) + [0 x_extent];
  end
  if xlims(2) > param.xlims(2)
    xlims = param.xlims(2) + [-x_extent 0];
  end
  if ylims(1) < param.ylims(1)
    ylims = param.ylims(1) + [0 y_extent];
  end
  if ylims(2) > param.ylims(2)
    ylims = param.ylims(2) + [-y_extent 0];
  end
  xlim(param.h_axes,xlims);
  ylim(param.h_axes,ylims);
end
