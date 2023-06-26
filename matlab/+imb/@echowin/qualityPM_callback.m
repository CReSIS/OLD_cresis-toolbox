function qualityPM_callback(obj,source,event)
% qualityPM_callback(obj,source,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);
