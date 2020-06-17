function button_down(obj,h_obj,event)
[x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);

% Make sure that click is on the right side panel
mouse_pos = get(obj.h_fig,'CurrentPoint');

% Check to make sure mouse clicked inside of obj.h_axes
%   Since extends the full y-length, just check to the right of minimum x
uipanel_pos = get(obj.left_panel.handle,'Position');
if mouse_pos(1) <= uipanel_pos(1)
  return
end

obj.zoom_mode_x = x;
obj.zoom_mode_y = y;
rbbox;
