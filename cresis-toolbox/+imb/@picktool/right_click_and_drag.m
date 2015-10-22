function cmds = right_click_and_drag(obj,param)
image_x = param.image_x;
cur_layers = param.cur_layers;
x = param.x;
y = param.y;
point_path_ids = param.point_path_id;
cmds = [];

fprintf('Delete points %f to %f, %f to %f\n', x, y);

for layer_idx = 1:length(cur_layers)
  cur_layer = cur_layers(layer_idx);
  
  point_idxs = find(param.layer.x >= x(1) & param.layer.x <= x(2) ...
    & param.layer.y{cur_layer} >= y(1) & param.layer.y{cur_layer} <= y(2));
  
  if ~isempty(point_idxs)
    cmds(end+1).undo_cmd = 'insert';
    cmds(end).undo_args = {cur_layer, point_idxs, ...
      param.layer.y{cur_layer}(point_idxs), ...
      param.layer.type{cur_layer}(point_idxs), ...
      param.layer.qual{cur_layer}(point_idxs)};
    cmds(end).redo_cmd = 'delete';
    cmds(end).redo_args = {cur_layer, [x y], point_idxs([1 end])};
  end
end

end
