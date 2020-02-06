function layerLB_callback(obj,source,event)

val = get(source,'Value');

obj.eg.layers.selected_layers(:)=false;
obj.eg.layers.selected_layers(val)=true;

% Update plot based on selection
obj.set_visibility();

end
