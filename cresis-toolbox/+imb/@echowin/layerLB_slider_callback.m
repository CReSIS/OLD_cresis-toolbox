function layerLB_slider_callback(obj,source,event)

% layerLB_slider_callback(obj,source,event)
%
% Layer listbox callback for the scrollbar. Moves the layer listbox up or
% down and redraws the table containing the listbox elements.
%

pointer = get(source,'value');
slider_max = get(obj.left_panel.layer_panel.slider,'Max');
% max is equal to size(1,obj.left_panel.layer_panel.layers.handles)-4

delta = diff([obj.left_panel.layer_panel.last_pointer pointer]);
if ~isequal(delta,round(delta))
  delta = round(delta);
  set(source,'value',obj.left_panel.layer_panel.last_pointer+delta);
  pointer = obj.left_panel.layer_panel.last_pointer+delta;
end
if delta ~= 0
  % scrolled a distance
  for idx = 1:obj.left_panel.layer_panel.MAX_ROW
    layer_idx = slider_max+1-pointer+idx-1;
    set(obj.left_panel.layer_panel.table.handles{idx,1},'Value',obj.left_panel.layer_panel.visible_layers(layer_idx));
    set(obj.left_panel.layer_panel.table.handles{idx,2},'Value',obj.left_panel.layer_panel.selected_layers(layer_idx));
    if length(obj.left_panel.layer_panel.layer_names{layer_idx}) <= 23
      % string fits in layers LB
      set(obj.left_panel.layer_panel.table.handles{idx,3},'String',obj.left_panel.layer_panel.layer_names{layer_idx});
      set(obj.left_panel.layer_panel.table.handles{idx,3},'TooltipString','');
    else
      % string doesn't fit in layers LB
      set(obj.left_panel.layer_panel.table.handles{idx,3},'String',[obj.left_panel.layer_panel.layer_names{layer_idx}(1:20) '...']);
      set(obj.left_panel.layer_panel.table.handles{idx,3},'TooltipString',obj.left_panel.layer_panel.layer_names{layer_idx});
    end
  end
else
  return;
end

% draw the table
clear row col
table_draw(obj.left_panel.layer_panel.table);

obj.left_panel.layer_panel.last_pointer = pointer;

return;
