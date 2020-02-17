function button_down(obj,src,event)
% button_down(obj,src,event)
%
% Support function for echowin class. Called when a mouse button is pressed
% when zoom mode is on and records click position and starts drawing the
% rubber band.

% Make sure that click is on the right side panel
mouse_pos = get(obj.h_fig,'CurrentPoint');
uipanel_pos = get(obj.right_panel.handle,'Position');
if mouse_pos(1) < uipanel_pos(1)
  obj.click_x = [];
  obj.click_y = [];
  return
end

[x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
%fprintf('Echo Button Down: x = %.3f, y = %.3f, but = %d\n', x, y, but);
obj.click_x = x;
obj.click_y = y;

if obj.alt_pressed || (obj.zoom_mode && but == 1) || (~obj.zoom_mode && but == 3)
  rbbox;  % start drawing rubber band
end

return
