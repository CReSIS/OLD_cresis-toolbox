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

try
  source_layer = round(str2double(get(obj.panel.sourceTB,'String')));
  if length(source_layer) ~= 1 || source_layer < 1 || source_layer > length(param.layer.y)
    error('');
  end
catch ME
  error('Source layer entered is invalid. Enter a single positive integer.\n');
end

diff_mode_en = get(obj.panel.diff_modeCB,'Value');

if diff_mode_en
  try
    correct_layer = round(str2double(get(obj.panel.correctTB,'String')));
    if length(correct_layer) ~= 1 || correct_layer < 1 || correct_layer > length(param.layer.y)
      error('');
    end
  catch ME
    error('Correct layer entered is invalid. Enter a single positive integer.\n');
  end
end

param.x_bounds = 2;
param.y_bounds = 2;
if diff_mode_en
  [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,correct_layer);
else
  [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,source_layer);
end

for layer_idx = 1:length(cur_layers)
  cur_layer = cur_layers(layer_idx);
  
  if ~isempty(point_idxs)
    if diff_mode_en
      cmds(end+1).undo_cmd = 'insert';
      cmds(end).undo_args = {cur_layer, point_idxs, ...
        param.layer.y{cur_layer}(point_idxs), ...
        param.layer.type{cur_layer}(point_idxs), ...
        param.layer.qual{cur_layer}(point_idxs)};
      cmds(end).redo_cmd = 'insert';
      cmds(end).redo_args = {cur_layer, point_idxs, ...
        param.layer.y{cur_layer}(point_idxs) + param.layer.y{correct_layer}(point_idxs) - param.layer.y{source_layer}(point_idxs), ...
        param.layer.type{cur_layer}(point_idxs), ...
        param.layer.qual{cur_layer}(point_idxs)};
    else
      cmds(end+1).undo_cmd = 'insert';
      cmds(end).undo_args = {cur_layer, point_idxs, ...
        param.layer.y{cur_layer}(point_idxs), ...
        param.layer.type{cur_layer}(point_idxs), ...
        param.layer.qual{cur_layer}(point_idxs)};
      cmds(end).redo_cmd = 'insert';
      s = param.layer.y{source_layer}(point_idxs);
      try
        size_s = size(s);
        eval_cmd = get(obj.panel.source_evalTB,'String');
        evalc(eval_cmd);
        % Ensure size is correct
        if ~isequal(size_s, size(s))
          error('s changed in size which is not allowed. Source eval parameter should be modified to prevent this from happening.');
        end
        % Ensure type is correct.
        if ~isa(s,'double')
          error('s changed in type which is not allowed. Source eval parameter should be modified to prevent this from happening.');
        end
      catch ME;
        warning(ME.getReport);
      end;
      cmds(end).redo_args = {cur_layer, point_idxs, ...
        s, ...
        param.layer.type{source_layer}(point_idxs), ...
        param.layer.qual{source_layer}(point_idxs)};
    end
  end
end

end
