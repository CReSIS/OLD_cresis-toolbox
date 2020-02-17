function yaxisPM_callback(obj,hObj,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

% Replot image
obj.plot_echogram(obj.eg.image_gps_time(1),obj.eg.image_gps_time(end),-inf,inf);
obj.plot_layers();
obj.plot_crossovers();
obj.plot_cursors();
  
%% Set layer and cross over visibility
obj.set_visibility();
