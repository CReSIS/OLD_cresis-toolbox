function layerCM_callback(obj,source,event)

if strcmp(source.Label,'&Visible')
  val = get(obj.left_panel.layerLB,'Value');
  
  obj.eg.layers.visible_layers(val)=true;
  
  % Update plot based on selection
  obj.set_visibility();
elseif strcmp(source.Label,'&Hide')
  val = get(obj.left_panel.layerLB,'Value');
  
  obj.eg.layers.visible_layers(val)=false;
  
  % Update plot based on selection
  obj.set_visibility();
elseif strcmp(source.Label,'&New layer')
  keyboard
  % obj.undo_stack.user_data.layer_info
  % obj.eg.layers
elseif strcmp(source.Label,'&Edit layer')
elseif strcmp(source.Label,'&Up')
elseif strcmp(source.Label,'&Down')
elseif strcmp(source.Label,'&Top')
elseif strcmp(source.Label,'&Bottom')
end

end
