function display_modePM_callback(obj,hObj,event)

obj.change_display_c;
obj.change_dynamic_range;

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
