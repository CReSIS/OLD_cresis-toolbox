function display_modePM_callback(obj,hObj,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

obj.change_display_c;
obj.change_dynamic_range;
