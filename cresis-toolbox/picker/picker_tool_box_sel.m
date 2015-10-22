function [x,y,outside] = picker_tool_box_sel(x_1, x_2, y_1, y_2, axis_lims, x_max)
% [x,y,outside] = picker_tool_box_sel(x_1, x_2, y_1, y_2, axis_lims, x_max)
%
% Support function for picker_pick.m.
% This is used for rubberband/box tools that select a region. The raw
% inputs and the current axis is passed in and these are parsed.
%
% Author: John Paden

% Make sure order is correct and switch around if not
if x_1 > x_2
  tmp = x_1;
  x_1 = x_2;
  x_2 = tmp;
end
if y_1 > y_2
  tmp = y_1;
  y_1 = y_2;
  y_2 = tmp;
end

% Check if both button presses were above the axis. Tools have
% special behavior in this case
if y_1 < axis_lims(3) && y_2 < axis_lims(3)
  outside = true;
else
  outside = false;
end

% Constrain picks to boundaries
if x_1 < max(1,axis_lims(1))
  x_1 = max(1,axis_lims(1));
elseif x_1 > min(x_max, axis_lims(2))
  x_1 = min(x_max, axis_lims(2));
end
if x_2 < max(1,axis_lims(1))
  x_2 = max(1,axis_lims(1));
elseif x_2 > min(x_max, axis_lims(2))
  x_2 = min(x_max, axis_lims(2));
end
if y_1 < axis_lims(3)
  y_1 = axis_lims(3);
elseif y_1 > axis_lims(4)
  y_1 = axis_lims(4);
end
if y_2 < axis_lims(3)
  y_2 = axis_lims(3);
elseif y_2 > axis_lims(4)
  y_2 = axis_lims(4);
end

% Prepare outputs
x = [x_1 x_2];
y = [y_1 y_2];

return;
