function key_release(obj,src,event)
% mapwin.key_release(obj,src,event)
% 
% Support function for mapwin class.

if ~any(strcmp('control',event.Modifier))
  obj.control_pressed = false;
end

if ~any(strcmp('shift',event.Modifier))
  obj.shift_pressed = false;
end

return