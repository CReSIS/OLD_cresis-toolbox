function layerLB_str(obj)

LB_strings = cell(1,length(obj.eg.layers.lyr_name));
for idx = 1:length(obj.eg.layers.lyr_name)
  if obj.eg.layers.visible_layers(idx)
    LB_strings{idx} = sprintf('(%d) %s:%s',idx,obj.eg.layers.lyr_group_name{idx},obj.eg.layers.lyr_name{idx});
  else
    LB_strings{idx} = sprintf('<HTML><FONT color="red">(%d) %s:%s</FONT></HTML>',idx,obj.eg.layers.lyr_group_name{idx},obj.eg.layers.lyr_name{idx});
  end
end
set(obj.left_panel.layerLB,'String',LB_strings);
