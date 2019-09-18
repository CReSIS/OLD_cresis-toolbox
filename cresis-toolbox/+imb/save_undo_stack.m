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
frames = [];
for cur_layer = unique_layers
  % Commit the commands for each layer, one layer at a time
  fprintf('Saving layer %d\n', cur_layer);
  param.properties.lyr_id = cur_layer;
  
  % Get the commands for just this layer
  cur_layer_cmds = find(cur_layer == cmd_layers);
  cur_cmd_type = cmd_type(cur_layer_cmds);
  cur_cmd_idxs = cmd_idxs(cur_layer_cmds);
  cur_sub_idxs = sub_idxs(cur_layer_cmds);
  
  %% OPS: Inserting and deleting points
  if strcmpi(undo_stack.user_data.layerSource,'OPS')
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
      end %end of if insert or delete
    end %while cur_layer_cmds end
    
    %% LayerData
  elseif strcmpi(undo_stack.user_data.layerSource,'layerdata')
    frames_insert=[];
    frames_delete=[];
    cur_layer_cmd_idx = 1;
    while cur_layer_cmd_idx <= length(cur_layer_cmds)
      cmd_idx = cur_cmd_idxs(cur_layer_cmd_idx);
      sub_idx = cur_sub_idxs(cur_layer_cmd_idx);
      %% LayerData: Inserting Points
      % gets information about which points are to be inserted and updates
      % the layerData file accordingly. Performs necessary steps to convert
      % from the unique point_path_id to a point number relative to that position
      % in a single frame
      if cur_cmd_type(cur_layer_cmd_idx) == 'i'
        point_idxs = [];
        point_mask = logical(zeros(size(undo_stack.user_data.frame)));
        mask = logical(zeros(size(cmds_list{cmd_idx}(sub_idx).redo_args{2})));
        for point_path_idx = 1:length(cmds_list{cmd_idx}(sub_idx).redo_args{2});
          point_id = cmds_list{cmd_idx}(sub_idx).redo_args{2}(point_path_idx);
          point = find(point_id == undo_stack.user_data.point_path_id);
          if(~isempty(point))
            point_mask(point_id)=true;
            mask(point_path_idx)=true;
            point_idxs(end+1) = point;
          end
        end
        % updates the undo_stack fields with information of the updated
        % points
        undo_stack.user_data.twtt{cur_layer}(point_idxs) = cmds_list{cmd_idx}(sub_idx).redo_args{3}(mask);
        undo_stack.user_data.type{cur_layer}(point_idxs) = cmds_list{cmd_idx}(sub_idx).redo_args{4}(mask);
        undo_stack.user_data.qual{cur_layer}(point_idxs) = cmds_list{cmd_idx}(sub_idx).redo_args{5}(mask);
        frames_changed = undo_stack.user_data.frame(point_mask);
        frms_changed = unique(frames_changed); % giving all frames changed
        frames_insert = cat(2,frames_insert,frms_changed);
        
        for frm = 1:length(frms_changed)
          found_frm = find(undo_stack.user_data.frame == frms_changed(frm));
          changed_mask = point_mask(found_frm);
          found_frm_idx = undo_stack.user_data.frame_idxs(found_frm);
          found_twtt = undo_stack.user_data.twtt{cur_layer}(changed_mask);
          found_type = undo_stack.user_data.type{cur_layer}(changed_mask);
          found_qual = undo_stack.user_data.qual{cur_layer}(changed_mask);
          changed_frm_idx = found_frm_idx(changed_mask);
          % Performs the update in the layer information stored in the
          % undo_stack
          undo_stack.user_data.layer_info(frms_changed(frm)).layerData{cur_layer}.value{2}.data(changed_frm_idx) = found_twtt; % updating the twtt
          undo_stack.user_data.layer_info(frms_changed(frm)).layerData{cur_layer}.value{1}.data(changed_frm_idx) = found_type; % updating the type
          undo_stack.user_data.layer_info(frms_changed(frm)).layerData{cur_layer}.quality(changed_frm_idx) = found_qual; % updating the quality
          
        end % end for loop
        cur_layer_cmd_idx = cur_layer_cmd_idx + 1;
        
        %% LayerData: Deleting Points
      else
        % gets information about which points are to be deleted and updates
        % the layerData file accordingly. Performs necessary steps to convert
        % from the unique point_path_id to a point number relative to that position
        % in a single frame.
        % gets the range of points to be deleted and stores it in
        % start_point_path_id, stop_point_path_id, min_twtt and max_twtt
        % respectively
        start_point_path_id = cmds_list{cmd_idx}(sub_idx).redo_args{3}(1); 
        stop_point_path_id = cmds_list{cmd_idx}(sub_idx).redo_args{3}(2);
        min_twtt = cmds_list{cmd_idx}(sub_idx).redo_args{2}(3);
        max_twtt = cmds_list{cmd_idx}(sub_idx).redo_args{2}(4);
        point_mask = logical(zeros(size(undo_stack.user_data.frame)));
        point_idxs = find(undo_stack.user_data.point_path_id >= start_point_path_id & undo_stack.user_data.point_path_id <= stop_point_path_id ...
          & undo_stack.user_data.twtt{cur_layer} > min_twtt & undo_stack.user_data.twtt{cur_layer} < max_twtt); % getting the range of points that are deleted
        
        point = undo_stack.user_data.point_path_id(point_idxs);
        point_mask(point)=true;
        
        % updates the undo_stack fields with information of the deleted
        % points
        undo_stack.user_data.twtt{cur_layer}(point_idxs) = NaN;
        undo_stack.user_data.type{cur_layer}(point_idxs) = 1;
        undo_stack.user_data.qual{cur_layer}(point_idxs) = 1;
        frames_changed = undo_stack.user_data.frame(point_mask);
        frms_changed = unique(frames_changed); % giving all frames changed
        frames_delete = cat(2,frames_delete,frms_changed);
        
        for frm = 1:length(frms_changed)
          found_frm = find(undo_stack.user_data.frame == frms_changed(frm));
          changed_mask = point_mask(found_frm);
          found_frm_idx = undo_stack.user_data.frame_idxs(found_frm);
          found_twtt = undo_stack.user_data.twtt{cur_layer}(changed_mask);
          found_qual = undo_stack.user_data.qual{cur_layer}(changed_mask);
          found_type = undo_stack.user_data.type{cur_layer}(changed_mask);
          changed_frm_idx = found_frm_idx(changed_mask);
          % Performs the update in the layer information stored in the
          % undo_stack
          undo_stack.user_data.layer_info(frms_changed(frm)).layerData{cur_layer}.value{2}.data(changed_frm_idx) = found_twtt; % updating the twtt
          undo_stack.user_data.layer_info(frms_changed(frm)).layerData{cur_layer}.quality(changed_frm_idx) = found_qual; % updating the quality
          undo_stack.user_data.layer_info(frms_changed(frm)).layerData{cur_layer}.value{1}.data(changed_frm_idx) = found_type; % updating the type
        end % end for loop
        cur_layer_cmd_idx = cur_layer_cmd_idx + 1;
      end% if-else end
    end%(while loop)
  end% layer_data end
end% end for loop

%% Notify all the other echo windows using this stack that a save has been
% done.
% save(FILENAME,'-append','-struct','LAYER_STRUCT_NAME','FIELD_NAME')
% getting the filename with correct frame number and updating the information in the file.
undo_stack.save();
if strcmpi(undo_stack.user_data.layerSource,'layerdata')
  unique_frms_insert = unique(frames_insert); %getting the frames in which points were inserted
  unique_frms_delete = unique(frames_delete); %getting the frames from which points were deleted
  if ~isempty(unique_frms_insert)
    for idx = 1:length(unique_frms_insert)
      layerData = undo_stack.user_data.layer_info(unique_frms_insert(idx)).layerData;
      filename_insert = undo_stack.user_data.filename{unique_frms_insert(idx)};
      save(filename_insert,'-append','layerData') % saving to layerData file
    end
  end
  if ~isempty(unique_frms_delete)
    for idx = 1:length(unique_frms_delete)
      layerData = undo_stack.user_data.layer_info(unique_frms_delete(idx)).layerData;
      filename_delete = undo_stack.user_data.filename{unique_frms_delete(idx)};
      save(filename_delete,'-append','layerData') % saving to layerData file
    end
  end
end
end
