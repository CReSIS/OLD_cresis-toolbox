function outside_limits = check_limits(obj,xaxis,yaxis,dir)
% outside_limits = check_limits(obj,xaxis,yaxis,dir)
%
% Support function to make sure the current map view field's limits are
% reasonably far away from the maximum x/y limits.
% dir is a char with a value of either 'u', 'd', 'l', 'r' corresponding to
% the direction of the map window's requested movement.

x_extent = diff(xaxis);
y_extent = diff(yaxis);

switch dir
  case 'u'
    outside_limits = abs(diff([yaxis(end) obj.map.yaxis_default(end)])) <= y_extent*0.05;
  case 'd'
    outside_limits = abs(diff([yaxis(1) obj.map.yaxis_default(1)])) <= y_extent*0.05;
  case 'l'
    outside_limits = abs(diff([xaxis(1) obj.map.xaxis_default(1)])) <= x_extent*0.05;
  case 'r'
    outside_limits = abs(diff([xaxis(end) obj.map.xaxis_default(end)])) <= x_extent*0.05;
end

return;