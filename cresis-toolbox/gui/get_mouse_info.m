function [x,y,but] = get_mouse_info(h,h_axis)
% [x,y,but] = get_mouse_info(h,h_axis)
%
% Gets information about the latest mouse/key event.
% Useful for callback functions.
%
% h = figure where click happened (can be empty)
% h_axis = axes where click happened
%
% x,y = position in axes
% but = state of button press (not valid when h is empty or
%   when a mouse button was not pressed)

% Get mouse position
mouse_pos = get(h_axis,'CurrentPoint');
x = mouse_pos(1,1);
y = mouse_pos(1,2);
if ~isempty(h)
  % Get mouse button
  if strcmpi(get(h,'SelectionType'),'normal')
    but = 1;
  elseif strcmpi(get(h,'SelectionType'),'extend')
    but = 2;
  elseif strcmpi(get(h,'SelectionType'),'alt')
    but = 3;
  elseif strcmpi(get(h,'SelectionType'),'open')
    but = 4;
  end
end

return;

