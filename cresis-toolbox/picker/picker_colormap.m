function picker_colormap(fig_handle)
% picker_colormap(fig_handle)
%
% Sets the specified figure's colormap to the control window's
% colormap popup menu setting
%
% Author: John Paden

global hui;

if ishandle(fig_handle)
  colormap_mode = get(hui.fig.ctrl_panel.colormapPM,'Value');
  if colormap_mode == 1
    set(fig_handle,'Colormap',1-gray(256));
  elseif colormap_mode == 2
    set(fig_handle,'Colormap',jet(256));
  elseif colormap_mode == 3
    set(fig_handle,'Colormap',flipud(copper(256)));
  end
end

return;
