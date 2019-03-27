function LayerSourcePM_callback(obj,status,event)
temp = get(obj.h_gui.LayerSourcePM,'String');
LayerSource = temp{get(obj.h_gui.LayerSourcePM,'Value')};
if strcmpi(LayerSource,'OPS')
  set(obj.h_gui.layerDataSourcePM,'Enable','off');
elseif strcmpi(LayerSource,'layerdata')
  set(obj.h_gui.layerDataSourcePM,'Enable','on');
end
return
