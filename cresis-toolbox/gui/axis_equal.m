function AXIS_LIM = axis_equal(h_axes, varargin)
% AXIS_LIM = axis_equal(h_axes, varargin)
%
%  Takes in an existing axis handle and X and Y vectors
%  Returns appropriate axis limits while adding a buffer zone to the edges
%   and maintaining the original aspect ratio
%
%  An additional 'buffer' convenience argument may be passed in.
%
%  Authors: Victor Berger, John Paden
%
% See also: tomo.create_movie.m, tomo.surfdata_to_DEM, tomo.plot_DEM,
%  zoom_button_up.m

p = inputParser;
addRequired(p, 'h_axes');
addRequired(p, 'x');
addRequired(p, 'y');
addOptional(p, 'buffer', 2, @isnumeric);
addOptional(p, 'xminbuffer', 0, @isnumeric);
addOptional(p, 'xmaxbuffer', 0, @isnumeric);
addOptional(p, 'yminbuffer', 0, @isnumeric);
addOptional(p, 'ymaxbuffer', 0, @isnumeric);
parse(p, h_axes, varargin{:});

cur_units = get(h_axes,'Units');
set(h_axes,'Units','pixels');
axis_pos = get(h_axes,'Position');
set(h_axes,'Units',cur_units);
aspect_ratio = axis_pos(3)/axis_pos(4);
axis(h_axes, 'equal');

if p.Results.xminbuffer ~= 0 || p.Results.xmaxbuffer ~= 0 ...
    || p.Results.yminbuffer ~= 0 || p.Results.ymaxbuffer ~= 0
  buffer = 0;
else
  buffer = p.Results.buffer;
end

if max(p.Results.x)-min(p.Results.x) < max(p.Results.y)-min(p.Results.y)
  % Flight line is longer on y axis than it is on x axis
  % No change to y limits
  % Change x limits to match original aspect ratio
  
  lims.YLim(1) = min(p.Results.y)/1e3;
  lims.YLim(2) = max(p.Results.y)/1e3;
  
  if buffer ~= 0
    lim_ymin = lims.YLim(1)-buffer;
    lim_ymax = lims.YLim(2)+buffer;
    lims.YLim = [lim_ymin, lim_ymax];
  end
  
  lim_xmin = (nanmean(p.Results.x)/1e3)-diff(lims.YLim)/2*aspect_ratio;
  lim_xmax = (nanmean(p.Results.x)/1e3)+diff(lims.YLim)/2*aspect_ratio;
  
  if p.Results.xminbuffer == 0 && p.Results.xmaxbuffer == 0 ...
      && p.Results.yminbuffer == 0 && p.Results.ymaxbuffer == 0
    AXIS_LIM = [lim_xmin, lim_xmax, ...
      lims.YLim(1)-(buffer/2), lims.YLim(2)+(buffer/2)];
  else
    AXIS_LIM = [lim_xmin-p.Results.xminbuffer, lim_xmax+p.Results.xmaxbuffer, ...
      lims.YLim(1)-p.Results.yminbuffer, lims.YLim(2)+p.Results.ymaxbuffer];
  end
  
else
  % Flight line is longer on x axis than it is on y axis
  % No change to x limits
  % Change y limits to match original aspect ratio
  
  lims.XLim(1) = min(p.Results.x)/1e3;
  lims.XLim(2) = max(p.Results.x)/1e3;
  
  if buffer ~= 0
    lim_xmin = lims.XLim(1)-buffer;
    lim_xmax = lims.XLim(2)+buffer;
    lims.XLim = [lim_xmin, lim_xmax];
  end
  
  lim_ymin = (nanmean(p.Results.y)/1e3)-diff(lims.XLim)/2/aspect_ratio;
  lim_ymax = (nanmean(p.Results.y)/1e3)+diff(lims.XLim)/2/aspect_ratio;
  
  if p.Results.xminbuffer == 0 && p.Results.xmaxbuffer == 0 ...
      && p.Results.yminbuffer == 0 && p.Results.ymaxbuffer == 0
    AXIS_LIM = [lims.XLim(1)-(buffer/2), lims.XLim(2)+(buffer/2), ...
      lim_ymin, lim_ymax];
  else
    AXIS_LIM = [lims.XLim(1)-p.Results.xminbuffer, lims.XLim(2)+p.Results.xmaxbuffer, ...
      lim_ymin-p.Results.yminbuffer, lim_ymax+p.Results.ymaxbuffer];
  end
  
end
end

