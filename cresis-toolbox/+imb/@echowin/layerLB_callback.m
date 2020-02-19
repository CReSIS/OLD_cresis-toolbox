function layerLB_callback(obj,source,event)
% layerLB_callback(obj,source,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
% uicontrol(obj.right_panel.status_panel.statusText);

val = get(source,'Value');

obj.eg.layers.selected_layers(:)=false;
obj.eg.layers.selected_layers(val)=true;

% Update plot based on selection
obj.set_visibility();
