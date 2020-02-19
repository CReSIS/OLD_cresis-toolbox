function framesPM_callback(obj,hObj,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
uicontrol(obj.right_panel.status_panel.statusText);

obj.default_params.max_frames = get(obj.left_panel.framesPM,'Value');
