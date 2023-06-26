function sourceCM_callback(obj,source,event)
% echowin.sourceCM_callback(obj,source,event)
%
% The source context menu callback function

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

sourceMenus = get(obj.left_panel.sourceCM,'Children');
img = length(sourceMenus)-find(sourceMenus == source)-3;

if img >= 0
  set(sourceMenus,'Checked','off');
  set(sourceMenus(length(sourceMenus)-3-img),'Checked','on');

  cur_axis = axis(obj.h_axes);
  
  % Convert x_min, x_max to GPS time
  xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
  
  % Draw data with new axis
  obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',3));

else
  if strcmpi(get(source,'Label'),'Add')
  elseif strcmpi(get(source,'Label'),'Remove')
  elseif strcmpi(get(source,'Label'),'Refresh')
    obj.update_source_fns_existence();
  end
end
