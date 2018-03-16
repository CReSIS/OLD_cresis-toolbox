function key_release(obj,src,event)
% key_release(obj,src,event)
% 
% Support function for echowin class.

modifiers = get(event.Source,'CurrentModifier');
obj.shift_pressed = false;
obj.control_pressed = false;
obj.alt_pressed = false;

return
