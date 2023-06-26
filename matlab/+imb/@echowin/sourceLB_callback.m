function sourceLB_callback(obj,hObj,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

if strcmpi(get(obj.h_fig,'SelectionType'),'open')
  % Double click
  %  - On a double click, the callback is called twice: once with
  %    'normal' SelectionType and then with 'open'.
  %  - On double click, load the selection to the picking window
  
  cur_axis = axis(obj.h_axes);
  
  % Convert x_min, x_max to GPS time
  xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
  
  % Draw data with new axis
  obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',3));
  
else
  % Regular click
  %gCtrl.source.cur_src = get(hui.fig.ctrl_panel.sourceLB,'Value');
end
