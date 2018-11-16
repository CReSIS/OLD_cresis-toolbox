function zoom_figure_setup(h_fig, figure_name)
% zoom_figure_setup(h_fig)
%
% Initialization for zoom figure.
%
% See also: zoom_arrow.m, zoom_button_motion.m, zoom_button_scroll.m,
%   zoom_button_up.m, zoom_figure_setup.m, zoom_setup.m


set(h_fig,'DockControls','off')
set(h_fig,'NumberTitle','off');
if strcmpi(class(h_fig),'double')
  set(h_fig,'Name',sprintf('%d: %s',h_fig, figure_name));
else
  set(h_fig,'Name',sprintf('%d: %s',h_fig.Number, figure_name));
end
set(h_fig,'ToolBar','none');
set(h_fig,'MenuBar','none');

zoom_setup(h_fig);
end