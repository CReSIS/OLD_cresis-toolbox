function update_cursors_echowin(obj,src,event)
% mapwin.update_cursors_echowin(obj,src,event)
%
% Draws cursors in the map window. Can be called via a click in the map
% window (through @mapwin/draw_vector_layers) or via a listener to any echowin
% objects open (@mapwin/<callback_name>) that gets triggered by a shift-left click
% in an echogram window (@echowin/button_up + @echowin/echo_update_cursors)

for idx = 1:length(obj.echowin_list)
  % Update the echowin's cursor in the echogram (this finds the closest
  % point to the supplied x,y and that becomes the next cursor position
  % for this echowin)
  obj.echowin_list(idx).set_cursor_by_map(src.cursor.x,src.cursor.y);
  % Update the echowin's cursor on the map
  set(obj.echowin_maps(idx).h_cursor,'XData',obj.echowin_list(idx).cursor.x, ...
    'YData',obj.echowin_list(idx).cursor.y);
  % Make the cursors be plotted on top
  uistack(obj.echowin_maps(idx).h_cursor,'top');
end

return;
