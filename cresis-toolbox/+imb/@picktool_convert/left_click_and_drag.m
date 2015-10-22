function cmds = left_click_and_drag(obj,param)
% cmds = left_click_and_drag(obj,param)
%
% Convert layer tool

image_x = param.image_x;
image_y = param.image_y;
cur_layers = param.cur_layers;
x = param.x;
y = param.y;
cmds = [];

fprintf('Convert points %f to %f, %f to %f\n', x, y);

tool_param1_str = get(obj.top_panel.target_layers_TE,'String');
% Assumes the tool_param1_str is a matlab expression that can be evaluated
try
  source_layer = round(str2double(tool_param1_str));
  if length(source_layer) ~= 1 || source_layer < 1 || source_layer > length(param.layer.y)
    error('');
  end
catch ME
  error('Layer numbers entered are invalid. Enter a single integer.\n');
end

param.x_bounds = 2;
param.y_bounds = 2;
[manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,source_layer);

for layer_idx = 1:length(cur_layers)
  cur_layer = cur_layers(layer_idx);
  
  if ~isempty(point_idxs)
    cmds(end+1).undo_cmd = 'insert';
    cmds(end).undo_args = {cur_layer, point_idxs, ...
      param.layer.y{cur_layer}(point_idxs), ...
      param.layer.type{cur_layer}(point_idxs), ...
      param.layer.qual{cur_layer}(point_idxs)};
    cmds(end).redo_cmd = 'insert';
    cmds(end).redo_args = {cur_layer, point_idxs, ...
      param.layer.y{source_layer}(point_idxs), ...
      param.layer.type{source_layer}(point_idxs), ...
      param.layer.qual{source_layer}(point_idxs)};
    
  end
end

end
