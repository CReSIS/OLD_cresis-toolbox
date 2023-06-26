function cmds = left_click_and_drag(obj,param)
% cmds = left_click_and_drag(obj,param)
%
% Quality tool

image_x = param.image_x;
image_y = param.image_y;
cur_layers = param.cur_layers;
x = param.x;
y = param.y;
cmds = [];

fprintf('Quality points %f to %f, %f to %f\n', x, y);

param.x_bounds = 2;
param.y_bounds = 2;
for layer_idx = 1:length(cur_layers)
  cur_layer = cur_layers(layer_idx);
  
  [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer);
  
  if ~isempty(point_idxs)
    cmds(end+1).undo_cmd = 'insert';
    cmds(end).undo_args = {cur_layer, point_idxs, ...
      param.layer.y{cur_layer}(point_idxs), ...
      param.layer.type{cur_layer}(point_idxs), ...
      param.layer.qual{cur_layer}(point_idxs)};
    cmds(end).redo_cmd = 'insert';
    cmds(end).redo_args = {cur_layer, point_idxs, ...
      param.layer.y{cur_layer}(point_idxs), ...
      param.layer.type{cur_layer}(point_idxs), ...
      param.cur_quality*ones(size(point_idxs))};
  end
end

return

