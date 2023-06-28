function flightCM_callback(obj,hObj,event)
% mapwin.flightCM_callback(obj,hObj,event)
%
% Copies to clipboard contents of the flight label string if it is not empty.

str = get(obj.top_panel.flightLabel,'String');
if ~isempty(str)
  clipboard('copy',str);
end

end
