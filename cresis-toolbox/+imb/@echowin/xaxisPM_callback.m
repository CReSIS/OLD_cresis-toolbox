function xaxisPM_callback(obj,hObj,event)
% echowin.xaxisPM_callback(obj,hObj,event)
%
% Callback for x-axis popup menu

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

cur_axis = axis(obj.h_axes);
cur_xaxis_gps = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time, ...
  cur_axis([1,2]),'linear','extrap');

obj.plot_echogram(obj.eg.image_gps_time(1),obj.eg.image_gps_time(end),-inf,inf);
obj.plot_layers();
obj.plot_crossovers();
obj.plot_cursors();

%% Set layer and cross over visibility
obj.set_visibility();

% apply original axis limits
new_xlim = interp1(obj.eg.image_gps_time,obj.eg.image_xaxis,...
  cur_xaxis_gps,'linear','extrap');
xlim(obj.h_axes,new_xlim);
ylim(obj.h_axes,cur_axis([3,4]));
