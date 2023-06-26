function key_release(obj,src,event)
% key_release(obj,src,event)
% 
% Support function for echowin class.

modifiers = get(event.Source,'CurrentModifier');
obj.shift_pressed = ismember('shift',   modifiers);  % true/false
obj.control_pressed  = ismember('control', modifiers);  % true/false
obj.alt_pressed   = ismember('alt',     modifiers);  % true/false