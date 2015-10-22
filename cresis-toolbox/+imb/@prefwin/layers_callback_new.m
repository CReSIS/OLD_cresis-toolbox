function layers_callback_new(obj,status,event)
% imb.prefwin.layers_callback_new(obj,status,event)
%
% Support function for layers_callback new layer dialog box

h_children = get(get(status,'UserData'),'Children');
name = get(h_children(7),'String');
group = get(h_children(5),'String');
description = get(h_children(3),'String');

if ~isempty(name)
  % Get the system that is currently selected
  system_name = get(obj.h_gui.systemsLB,'String');
  system_name = system_name{get(obj.h_gui.systemsLB,'Value')};
  
  ops_param.properties.lyr_name = name;
  ops_param.properties.lyr_group_name = group;
  ops_param.properties.lyr_description = description;
  ops_param.properties.public = true;
  
  [status,data] = opsCreateLayer(system_name,ops_param);
  
  obj.layers_callback_refresh();
end

uiresume(gcbf);

end
