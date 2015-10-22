function save_undo_stack(undo_stack)
% imb.save_undo_stack(undo_stack)
%
% Support function for echowin and mapwin classes.  Does the actual
% saving of commands in the undo stack.
%

% =========================================================================
% IMPORTANT NOTE:
% cmds_list: The structure of each command is explained in
%   imb.echowin.cmds.cmds_execute.m/cmds_execute_insert
%   imb.echowin.cmds.cmds_execute.m/cmds_execute_delete
%
% The command list is two lists deep (i.e. a list of lists)
% 1. cmds_list is a cell vector of commands
% 2. Each cell contains a struct vector of commands (sub_idxs)
% 3. The database save operation does not care about the cell vector/struct vector
%    boundaries and so this save function flattens the lists into a single list.
%    (FYI: Undo/redo operations execute all the commands in a cell and
%    that is why the double list is used.)
% =========================================================================

% Get the command list for this undo stack and delete these commands from
% the list.
cmds_list = undo_stack.get_save_cmds(true);

%% Group commands by layer
cmd_type = char([]);
cmd_layers = [];
cmd_args = {};
cmd_idxs = [];
sub_idxs = [];
% Flatten command list struct to make it easy to parse and apply matrix
% operations to it
for cmd_idx = 1:length(cmds_list)
  for sub_idx = 1:length(cmds_list{cmd_idx})
    if strcmpi(cmds_list{cmd_idx}(sub_idx).redo_cmd,'insert')
      cmd_type(end+1) = 'i';
    else
      cmd_type(end+1) = 'd';
    end
    cmd_layers(end+1) = cmds_list{cmd_idx}(sub_idx).redo_args{1};
    cmd_args{end+1} = cmds_list{cmd_idx}(sub_idx).redo_args;
    cmd_idxs(end+1) = cmd_idx;
    sub_idxs(end+1) = sub_idx;
  end
end

% Get a list of unique database layer IDs that are being modified
unique_layers = unique(cmd_layers);
if isempty(unique_layers)
  unique_layers = [];
end

for cur_layer = unique_layers
  % Commit the commands for each layer, one layer at a time
  fprintf('Saving layer %d\n', cur_layer);
  param.properties.lyr_id = cur_layer;
  
  % Get the commands for just this layer
  cur_layer_cmds = find(cur_layer == cmd_layers);
  cur_cmd_type = cmd_type(cur_layer_cmds);
  cur_cmd_idxs = cmd_idxs(cur_layer_cmds);
  cur_sub_idxs = sub_idxs(cur_layer_cmds);
  
  % Loop through commands. If multiple inserts are done, these are
  % combined into a single insert
  cur_layer_cmd_idx = 1;
  while cur_layer_cmd_idx <= length(cur_layer_cmds)
    cmd_idx = cur_cmd_idxs(cur_layer_cmd_idx);
    sub_idx = cur_sub_idxs(cur_layer_cmd_idx);
    if cur_cmd_type(cur_layer_cmd_idx) == 'i'
      param.properties.point_path_id = cmds_list{cmd_idx}(sub_idx).redo_args{2};
      param.properties.twtt = cmds_list{cmd_idx}(sub_idx).redo_args{3};
      param.properties.type = cmds_list{cmd_idx}(sub_idx).redo_args{4};
      param.properties.quality = cmds_list{cmd_idx}(sub_idx).redo_args{5};
      fprintf('  Insert %d points\n', length(cmds_list{cmd_idx}(sub_idx).redo_args{2}));
      % Append any following insert commands onto this command... (if we
      % hit the end of the command list or a delete command we stop appending)
      cur_layer_cmd_idx = cur_layer_cmd_idx + 1;
      while cur_layer_cmd_idx <= length(cur_layer_cmds) ...
          && cur_cmd_type(cur_layer_cmd_idx) == 'i'
        cmd_idx = cur_cmd_idxs(cur_layer_cmd_idx);
        sub_idx = cur_sub_idxs(cur_layer_cmd_idx);
        param.properties.point_path_id = cat(2,param.properties.point_path_id, ...
          cmds_list{cmd_idx}(sub_idx).redo_args{2});
        param.properties.twtt = cat(2, param.properties.twtt, ...
          cmds_list{cmd_idx}(sub_idx).redo_args{3});
        param.properties.type = cat(2, param.properties.type, ...
          cmds_list{cmd_idx}(sub_idx).redo_args{4});
        param.properties.quality = cat(2, param.properties.quality, ...
          cmds_list{cmd_idx}(sub_idx).redo_args{5});
        cur_layer_cmd_idx = cur_layer_cmd_idx + 1;
        fprintf('  Insert %d points appended\n', length(cmds_list{cmd_idx}(sub_idx).redo_args{2}));
      end
      fprintf(' Insert commit\n');
      % Remove duplicate point ID/layerIDs
      valid_mask = logical(ones(size(param.properties.point_path_id)));
      for idx = length(valid_mask):-1:1
        if ~isempty(find(param.properties.point_path_id(idx+1:end) ...
            == param.properties.point_path_id(idx)))
          valid_mask(idx) = 0;
        end
      end
      param.properties.point_path_id = param.properties.point_path_id(valid_mask);
      param.properties.twtt = param.properties.twtt(valid_mask);
      param.properties.type = param.properties.type(valid_mask);
      param.properties.quality = param.properties.quality(valid_mask);
      opsCreateLayerPoints(undo_stack.unique_id{1},param);
    else
      param.properties.start_point_path_id = cmds_list{cmd_idx}(sub_idx).redo_args{3}(1);
      param.properties.stop_point_path_id = cmds_list{cmd_idx}(sub_idx).redo_args{3}(2);
      param.properties.min_twtt = cmds_list{cmd_idx}(sub_idx).redo_args{2}(3);
      param.properties.max_twtt = cmds_list{cmd_idx}(sub_idx).redo_args{2}(4);
      cur_layer_cmd_idx = cur_layer_cmd_idx + 1;
      fprintf(' Delete commit\n');
      opsDeleteLayerPoints(undo_stack.unique_id{1},param);
    end
  end
end
  
%% Notify all the other echo windows using this stack that a save has been
% done.
undo_stack.save();


end
