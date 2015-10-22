function sourceCM_callback(obj,source,event)
% echowin.sourceCM_callback(obj,source,event)
%
% The source context menu callback function

sourceMenus = get(obj.left_panel.sourceCM,'Children');
img = length(sourceMenus)-find(sourceMenus == source)-3;

if img >= 0
  set(sourceMenus,'Checked','off');
  set(sourceMenus(length(sourceMenus)-3-img),'Checked','on');

  cur_axis = axis(obj.right_panel.axes.handle);
  
  % Convert x_min, x_max to GPS time
  xlims = interp1(obj.eg.image_xaxis,obj.eg.image_gps_time,cur_axis(1:2),'linear','extrap');
  
  % Draw data with new axis
  obj.redraw(xlims(1),xlims(2),cur_axis(3),cur_axis(4),struct('clipped',3));
  
  try
    warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
    javaFrame = get(obj.h_fig,'JavaFrame');
    javaFrame.getAxisComponent.requestFocus;
  catch
    obj.status_text_set(sprintf('Focus error, click inside echogram window before using key shortcuts'),'replace');
  end

else
  if strcmpi(get(source,'Label'),'Add')
  elseif strcmpi(get(source,'Label'),'Remove')
  elseif strcmpi(get(source,'Label'),'Refresh')
    obj.update_source_fns_existence();
  end
end

end
