function [changed,pos] = compute_new_map_limits(obj,new_xdata,new_ydata)
% [changed,pos] = mapwin.compute_new_map_limits(obj,new_xdata,new_ydata)
%
% Computes x limits and y limits based on current limits and new x/y data
% that must be completely visible. Attempts to keep the extent of x and y
% the same. Keeps aspect ratios the same.
%
% changed: logical which indicates that pos is different than current
%          limits
% pos: [xmin xmax ymin ymax]

xlims = sort(get(obj.map_panel.h_axes,'XLim'));
ylims = sort(get(obj.map_panel.h_axes,'YLim'));

% Provide a small buffer around the x-data
min_xdata = min(new_xdata);
max_xdata = max(new_xdata);
new_xextent = max_xdata - min_xdata;
min_xdata = min_xdata - 0.01*max(xlims(2)-xlims(1),new_xextent);
max_xdata = max_xdata + 0.01*max(xlims(2)-xlims(1),new_xextent);

% Provide a small buffer around the y-data
min_ydata = min(new_ydata);
max_ydata = max(new_ydata);
new_yextent = max_ydata - min_ydata;
min_ydata = min_ydata - 0.01*max(ylims(2)-ylims(1),new_yextent);
max_ydata = max_ydata + 0.01*max(ylims(2)-ylims(1),new_yextent);

% Check to see if new data exceeds current boundaries
xmin_exceeded = min_xdata < xlims(1);
xmax_exceeded = max_xdata > xlims(2);
ymin_exceeded = min_ydata < ylims(1);
ymax_exceeded = max_ydata > ylims(2);

if xmin_exceeded || xmax_exceeded || ymin_exceeded || ymax_exceeded
  changed = true;
else
  changed = false;
end

% Provide a slightly larger buffer around the x-data for the new
% recommended position "pos"
min_xdata = min(new_xdata);
max_xdata = max(new_xdata);
new_xextent = max_xdata - min_xdata;
min_xdata = min_xdata - 0.05*max(xlims(2)-xlims(1),new_xextent);
max_xdata = max_xdata + 0.05*max(xlims(2)-xlims(1),new_xextent);

% Provide a slightly larger buffer around the y-data for the new
% recommended position "pos"
min_ydata = min(new_ydata);
max_ydata = max(new_ydata);
new_yextent = max_ydata - min_ydata;
min_ydata = min_ydata - 0.05*max(ylims(2)-ylims(1),new_yextent);
max_ydata = max_ydata + 0.05*max(ylims(2)-ylims(1),new_yextent);

% Shift x-limits and don't change the x-extent unless we have to in order
% to fit all the new data
if xmin_exceeded && xmax_exceeded
  new_x_lims = [min_xdata max_xdata];
elseif xmin_exceeded
  new_x_lims = [min_xdata min_xdata+max(xlims(2)-xlims(1),max_xdata-min_xdata)];
elseif xmax_exceeded
  new_x_lims = [max_xdata-max(xlims(2)-xlims(1),max_xdata-min_xdata) max_xdata];
else
  new_x_lims = xlims;
end

% Shift y-limits and don't change the y-extent unless we have to in order
% to fit all the new data
if ymin_exceeded && ymax_exceeded
  new_y_lims = [min_ydata max_ydata];
elseif ymin_exceeded
  new_y_lims = [min_ydata min_ydata+max(ylims(2)-ylims(1),max_ydata-min_ydata)];
elseif ymax_exceeded
  new_y_lims = [max_ydata-max(ylims(2)-ylims(1),max_ydata-min_ydata) max_ydata];
else
  new_y_lims = ylims;
end

% Now force the aspect ratio to remain the same
new_xextent = new_x_lims(2)-new_x_lims(1);
old_xextent = xlims(2)-xlims(1);
new_yextent = new_y_lims(2)-new_y_lims(1);
old_yextent = ylims(2)-ylims(1);
new_aspect_ratio = new_xextent / new_yextent;
old_aspect_ratio = old_xextent / old_yextent;
if new_aspect_ratio > old_aspect_ratio
  % Grow y limits
  y_center = mean(new_y_lims);
  new_yextent = new_xextent / old_aspect_ratio;
  new_y_lims = y_center + new_yextent*[-0.5 0.5];
else
  % Grow x limits
  x_center = mean(new_x_lims);
  new_xextent = new_yextent * old_aspect_ratio;
  new_x_lims = x_center + new_xextent*[-0.5 0.5];
end

pos = [new_x_lims new_y_lims];

end
