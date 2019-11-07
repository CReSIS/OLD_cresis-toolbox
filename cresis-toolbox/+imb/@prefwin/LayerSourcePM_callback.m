function layerSourcePM_callback(obj,status,event)
% layerSourcePM_callback(obj,status,event)
%
% Author: Anjali Pare

temp = get(obj.h_gui.layerSourcePM,'String');

layer_source = temp{get(obj.h_gui.layerSourcePM,'Value')};

if strcmpi(layer_source,'OPS')
  set(obj.h_gui.layerDataSourcePM,'Enable','off');
  obj.h_gui.h_layers.set_enable(true);
elseif strcmpi(layer_source,'layerdata')
  set(obj.h_gui.layerDataSourcePM,'Enable','on');
  obj.h_gui.h_layers.set_enable(false);
end
