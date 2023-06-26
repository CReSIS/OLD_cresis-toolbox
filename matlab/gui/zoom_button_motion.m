function zoom_button_motion(event,param)
% zoom_button_motion(event,param)
%
% Function to be used in figure's WindowButtonMotionFcn callback function.
%
% event = event passed into figure's WindowScrollWheelFcn callback function
% param
%  .xlims = max x-limits
%  .ylims = max y-limits
%  .h_axes = axes handle that we are zooming with
%  .h_fig = figure handle that we are zooming in
%  .axes = optional, string which specifies which axes can zoom ('xy' is
%  default, 'x' is x only, 'y' is y only)
%
% See also: zoom_arrow.m, zoom_button_motion.m, zoom_button_scroll.m,
%   zoom_button_up.m, zoom_figure_setup.m, zoom_setup.m

if ~isfield(param,'zoom_mode') || isempty(param.zoom_mode)
  param.zoom_mode = true;
end

[x,y,but] = get_mouse_info(param.h_fig,param.h_axes);
xlims = xlim(param.h_axes);
ylims = ylim(param.h_axes);
if param.zoom_mode && x>=xlims(1) && x<=xlims(2) && y>=ylims(1) && y<=ylims(2)
  if ~strcmpi(get(param.h_fig,'pointer'),'custom')
    set(param.h_fig,'Pointer','custom');
  end
else
  if ~strcmpi(get(param.h_fig,'pointer'),'Arrow')
    set(param.h_fig,'Pointer','Arrow');
  end
end
