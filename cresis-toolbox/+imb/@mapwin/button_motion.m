function button_motion(obj,src,event)
% mapwin.button_motion(obj,src, event)
%
% Called by mapwin when mouse is moved.  Handles tracking mouse position
% and mouse cursor changing when outside of axes.

if obj.busy_mode
  %% Set the mouse pointer graphic
  set(obj.h_fig,'Pointer','watch');
  
elseif strcmpi(get(obj.map_panel.h_axes,'Visible'),'on')
  %% Check to make sure map is in the map panel
  set(obj.h_fig,'Units','normalized');
  mouse_pos = get(obj.h_fig,'CurrentPoint');
  set(obj.map_panel.handle,'Units','normalized');
  uipanel_pos = get(obj.map_panel.handle,'Position');
  
  set(obj.status_panel.handle,'Units','normalized');
  status_pos = get(obj.status_panel.handle,'Position');
  
  if mouse_pos(2) < uipanel_pos(2) || mouse_pos(2) > uipanel_pos(4)
    % Not in the axes
    set(obj.h_fig,'Pointer','Arrow');
    return;
  end
  
  %% Set the mouse pointer graphic
  if obj.zoom_mode
    set(obj.h_fig,'Pointer','custom');
  else
    set(obj.h_fig,'Pointer','Arrow');
  end
  
  % print current mouse position to status bar
  axis_pos = get(obj.map_panel.h_axes,'CurrentPoint');
  x = axis_pos(1,1);
  y = axis_pos(1,2);
  if(obj.isGoogle)
    [lat, lon] = imb.world_to_latlon(x, y);
    set(obj.status_panel.mouseCoordText,'String',sprintf('%fN %fW; World Coordinates: (%f, %f); ',lat,lon, x, y));
  else
    [lat,lon] = projinv(obj.map.projmat,x*1e3,y*1e3);
    set(obj.status_panel.mouseCoordText,'String',sprintf('%8.3fN %8.3fW; X=%8.3fkm Y=%8.3fkm  ',lat,lon,x,y));
  end
end

end
