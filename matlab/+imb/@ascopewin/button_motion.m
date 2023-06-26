function button_motion(obj,src,event)
% button_motion(obj,src, event)
%
% Called by ascopewin when mouse cursor is moved. Handles tracking mouse position
% and mouse cursor changing when outside of axes

mouse_pos = get(obj.h_fig,'CurrentPoint');
uipanel_pos = get(obj.right_panel.handle,'Position');
% check if inside status bar
status_pos = get(obj.right_panel.status_panel.handle,'Position');
status_h = status_pos(4);
if mouse_pos(1) < uipanel_pos(1) || mouse_pos(2) < status_h
  set(obj.h_fig,'Pointer','Arrow');
  return;
elseif obj.zoom_mode
  set(obj.h_fig,'Pointer','custom');
end
