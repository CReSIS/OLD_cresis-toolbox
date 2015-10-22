function xaxisPM_callback(obj,hObj,event)
% echowin.xaxisPM_callback(obj,hObj,event)
%
% Callback for x-axis popup menu

cur_axis = axis(obj.right_panel.axes.handle);
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
xlim(obj.right_panel.axes.handle,new_xlim);
ylim(obj.right_panel.axes.handle,cur_axis([3,4]));

% next 3 lines are a workaround due to limitations of matlab's gui
% after a menu has been accessed, the menu keeps focus unless the user
% clicks elsewhere
% as a result, the next time the user presses a key shortcut, it will
% toggle this callback because it is still selected
% this is notably annoying when spacebar is pressed to toggle layer
% visibility, because this callback gets called in addition to the echogram
% key press function and the result is massive redraw delay
% if this code generates an error, it can be removed, but the user needs to
% click in the echogram figure after every time a button is pressed
try
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  javaFrame = get(obj.h_fig,'JavaFrame');
  javaFrame.getAxisComponent.requestFocus;
catch
  obj.status_text_set(sprintf('Focus error, click inside echogram window before using key shortcuts'),'replace');
end

return
