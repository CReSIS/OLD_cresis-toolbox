function button_down(obj,src,event)
% mapwin.button_down(obj,src,event)
%
% Support function for mapwin class. Called when a mouse button is pressed
% when zoom mode is on and records click position and starts drawing the
% rubber band.

[x,y,but] = get_mouse_info(obj.h_fig,obj.map_panel.h_axes);
obj.click_x = x;
obj.click_y = y;

if strcmpi(get(obj.map_panel.h_axes,'Visible'),'on')
  modifier_string = '';
  if obj.control_pressed
    modifier_string = strcat(modifier_string,' ctrl');
  end
  if obj.shift_pressed
    modifier_string = strcat(modifier_string,' shift');
  end
  %fprintf('Map Button Down: x = %.3f, y = %.3f, but = %d %s\n', x, y, but, modifier_string);
  
  % Make sure click is in the axis
  xlims = xlim(obj.map_panel.h_axes);
  ylims = ylim(obj.map_panel.h_axes);
  if obj.control_pressed || x < xlims(1) || x > xlims(2) || y < ylims(1) || y > ylims(2)
    % Click outside of axis or control pressed so nothing to do
    return;
  end
  
  if ~obj.control_pressed && ~obj.shift_pressed % no modifiers
    if but == 1
      % Left click and drag: Zoom
      rbbox;  % start drawing rubber band
    end
  end
end

end
