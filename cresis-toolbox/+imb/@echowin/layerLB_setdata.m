function layerLB_setdata(obj,data)

% layerLB_setdata(obj,data)
% 
% Sets the initial value of layer names in the layers listbox. Called in
% echowin/draw.
%
% data is a struct with:
%   names: 1D cell vector of layer names
%   surface: surface layer index
%   bottom: bottom layer index
%

% layer name data
for idx = 1:length(data.names)
  obj.left_panel.layer_panel.layer_names{idx,1} = data.names{idx};
end
% set up state variables of radio buttons (selected layer) and checkboxes
% (visible layers)
obj.left_panel.layer_panel.selected_layers = logical(zeros(length(data.names),1));
obj.left_panel.layer_panel.visible_layers = logical(ones(length(data.names),1));

obj.left_panel.layer_panel.surface = data.surface;
obj.left_panel.layer_panel.bottom = data.bottom;

% populate listbox with data
for idx = 1:obj.left_panel.layer_panel.MAX_ROW
  % set visibilities and values
  if idx <= length(obj.left_panel.layer_panel.layer_names)
    set(obj.left_panel.layer_panel.table.handles{idx,3},'Visible','on');
    if length(obj.left_panel.layer_panel.layer_names{idx}) <= 23
      % string fits in layers LB
      set(obj.left_panel.layer_panel.table.handles{idx,3},'String',obj.left_panel.layer_panel.layer_names{idx});
      set(obj.left_panel.layer_panel.table.handles{idx,3},'TooltipString','');
    else
      % string doesn't fit in layers LB
      set(obj.left_panel.layer_panel.table.handles{idx,3},'String',[obj.left_panel.layer_panel.layer_names{idx}(1:20) '...']);
      set(obj.left_panel.layer_panel.table.handles{idx,3},'TooltipString',obj.left_panel.layer_panel.layer_names{idx});
    end
    set(obj.left_panel.layer_panel.table.handles{idx,1},'Visible','on');
    set(obj.left_panel.layer_panel.table.handles{idx,1},'Value',obj.left_panel.layer_panel.visible_layers(idx));
    set(obj.left_panel.layer_panel.table.handles{idx,2},'Visible','on');
    set(obj.left_panel.layer_panel.table.handles{idx,2},'Value',obj.left_panel.layer_panel.selected_layers(idx));
  else
    set(obj.left_panel.layer_panel.table.handles{idx,3},'Visible','off');
    set(obj.left_panel.layer_panel.table.handles{idx,1},'Visible','off');
    set(obj.left_panel.layer_panel.table.handles{idx,2},'Visible','off');
  end
end

% configure the slider
if length(obj.left_panel.layer_panel.layer_names) > obj.left_panel.layer_panel.MAX_ROW
  num_el = length(obj.left_panel.layer_panel.layer_names);
  set(obj.left_panel.layer_panel.slider,'Enable','on');
  set(obj.left_panel.layer_panel.slider,'Min',1);
  set(obj.left_panel.layer_panel.slider,'Max',num_el-(obj.left_panel.layer_panel.MAX_ROW-1));
  set(obj.left_panel.layer_panel.slider,'Value',get(obj.left_panel.layer_panel.slider,'Max'));
  set(obj.left_panel.layer_panel.slider,'SliderStep',[1/(num_el-obj.left_panel.layer_panel.MAX_ROW) 5/(num_el-obj.left_panel.layer_panel.MAX_ROW)]);
  obj.left_panel.layer_panel.last_pointer = get(obj.left_panel.layer_panel.slider,'Max');
end

% redraw
clear row col
table_draw(obj.left_panel.layer_panel.table);