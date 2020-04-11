function ascopeLB_callback(obj,source,event)
% ascopeLB_callback(obj,source,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
% uicontrol(obj.right_panel.status_panel.statusText);

val = get(source,'Value');

obj.ascope.selected(:)=false;
obj.ascope.selected(val)=true;

obj.plot_update();
