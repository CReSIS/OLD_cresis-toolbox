function zoom_setup(h_fig)
% zoom_setup(h_fig)
%
% Initialization for zoom pointer
%
% See also: zoom_arrow.m, zoom_button_motion.m, zoom_button_scroll.m,
%   zoom_button_up.m, zoom_figure_setup.m, zoom_setup.m

set(h_fig,'Pointer','custom');
zoom_pointer = [NaN   NaN   NaN   NaN     1     1     1     1   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN
  NaN   NaN     1     1   NaN     2   NaN     2     1     1   NaN   NaN   NaN   NaN   NaN   NaN
  NaN     1     2   NaN     2     1     1   NaN     2   NaN     1   NaN   NaN   NaN   NaN   NaN
  NaN     1   NaN     2   NaN     1     1     2   NaN     2     1   NaN   NaN   NaN   NaN   NaN
  1   NaN     2   NaN     2     1     1   NaN     2   NaN     2     1   NaN   NaN   NaN   NaN
  1     2     1     1     1     1     1     1     1     1   NaN     1   NaN   NaN   NaN   NaN
  1   NaN     1     1     1     1     1     1     1     1     2     1   NaN   NaN   NaN   NaN
  1     2   NaN     2   NaN     1     1     2   NaN     2   NaN     1   NaN   NaN   NaN   NaN
  NaN     1     2   NaN     2     1     1   NaN     2   NaN     1   NaN   NaN   NaN   NaN   NaN
  NaN     1   NaN     2   NaN     1     1     2   NaN     2     1     2   NaN   NaN   NaN   NaN
  NaN   NaN     1     1     2   NaN     2   NaN     1     1     1     1     2   NaN   NaN   NaN
  NaN   NaN   NaN   NaN     1     1     1     1   NaN     2     1     1     1     2   NaN   NaN
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1     2   NaN
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1     2
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     1     1
  NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN   NaN     2     1     2];
set(h_fig,'PointerShapeCData',zoom_pointer);
set(h_fig,'PointerShapeHotSpot',[6 6])

end