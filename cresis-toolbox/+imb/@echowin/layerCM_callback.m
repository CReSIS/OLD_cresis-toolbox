function layerCM_callback(obj,source,event)

if source == obj.left_panel.layerCM_visible
  val = get(obj.left_panel.layerLB,'Value');
  
  obj.eg.layers.visible_layers(val)=true;
  
  % Update plot based on selection
  obj.set_visibility();
elseif source == obj.left_panel.layerCM_hide
  val = get(obj.left_panel.layerLB,'Value');
  
  obj.eg.layers.visible_layers(val)=false;
  
  % Update plot based on selection
  obj.set_visibility();
elseif source == obj.left_panel.layerCM_new || source == obj.left_panel.layerCM_copy
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers. The new layer will be inserted
    % before the first of the currently selected layers or at the bottom
    % of the listbox.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2);
    if isempty(val)
      if source == obj.left_panel.layerCM_copy
        return;
      end
      val = length(obj.eg.layers.lyr_id)+1;
    else
      val = val(1);
    end
    prompt = {'Layer Name:','Layer Group Name:','Description:'};
    if source == obj.left_panel.layerCM_copy
      def = {obj.eg.layers.lyr_name{val},char(obj.eg.layers.lyr_group_name{val}),''};
      val = val+1;
    else
      def = {'', '', ''};
    end
    dlg_title = 'New Layer';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if length(answer) == 3 && ~isempty(answer{1})
      name = answer{1};
      group_name = answer{2};
      desc = answer{3};
      
      if any(strcmpi(name,obj.eg.layers.lyr_name))
        fprintf('  Layer %s already exists\n', name);
        return;
      end
      
      fprintf('Add layer %s:%s "%s"\n', group_name, name, desc);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_delete';
      cmds(end).undo_args = {val};
      cmds(end).redo_cmd = 'layer_new';
      cmds(end).redo_args = {val,name,group_name,desc};
      
      % Push the new command(s) to the stack
      obj.undo_stack.push(cmds);
    end
  end
  
elseif source == obj.left_panel.layerCM_edit
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2);
    if isempty(val)
      return;
    else
      val = val(1);
    end
    prompt = {'Layer Name:','Layer Group Name:','Description:'};
    old_name = obj.eg.layers.lyr_name{val};
    old_group_name = char(obj.eg.layers.lyr_group_name{val});
    def = {old_name,old_group_name,''};
    dlg_title = 'Edit Layer';
    num_lines = 1;
    answer = inputdlg(prompt,dlg_title,num_lines,def);
    
    if length(answer) == 3 && ~isempty(answer{1})
      name = answer{1};
      group_name = answer{2};
      desc = answer{3};
      new_val = val;
      
      if any(strcmpi(name,obj.eg.layers.lyr_name([1:val-1 val+1:end])))
        fprintf('  Layer %s already exists\n', name);
        return;
      end
      
      fprintf('Edit layer %s:%s to %s:%s "%s"\n', old_group_name, old_name, group_name, name, desc);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {new_val,old_name,old_group_name,'',val};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {val,name,group_name,desc,new_val};
      
      % Push the new command(s) to the stack
      obj.undo_stack.push(cmds);
    end
  end
  
elseif source == obj.left_panel.layerCM_up
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>3);
    if isempty(val)
      return;
    else
      val = val(1);
    end
    name = obj.eg.layers.lyr_name{val};
    group_name = obj.eg.layers.lyr_group_name{val};
    new_val = val-1;
    
    fprintf('Move layer up %s:%s\n', group_name, name);
    
    cmds = [];
    cmds(end+1).undo_cmd = 'layer_edit';
    cmds(end).undo_args = {new_val,name,group_name,'',val};
    cmds(end).redo_cmd = 'layer_edit';
    cmds(end).redo_args = {val,name,group_name,'',new_val};
    
    % Push the new command(s) to the stack
    obj.undo_stack.push(cmds);
  end
  
elseif source == obj.left_panel.layerCM_down
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2 & val<length(obj.eg.layers.lyr_name));
    if isempty(val)
      return;
    else
      val = val(1);
    end
    name = obj.eg.layers.lyr_name{val};
    group_name = obj.eg.layers.lyr_group_name{val};
    new_val = val+1;
    
    fprintf('Move layer down %s:%s\n', group_name, name);
    
    cmds = [];
    cmds(end+1).undo_cmd = 'layer_edit';
    cmds(end).undo_args = {new_val,name,group_name,'',val};
    cmds(end).redo_cmd = 'layer_edit';
    cmds(end).redo_args = {val,name,group_name,'',new_val};
    
    % Push the new command(s) to the stack
    obj.undo_stack.push(cmds);
  end
  
elseif source == obj.left_panel.layerCM_top
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>3);
    if isempty(val)
      return;
    else
      val = val(1);
    end
    name = obj.eg.layers.lyr_name{val};
    group_name = obj.eg.layers.lyr_group_name{val};
    new_val = 3;
    
    fprintf('Move layer top %s:%s\n', group_name, name);
    
    cmds = [];
    cmds(end+1).undo_cmd = 'layer_edit';
    cmds(end).undo_args = {new_val,name,group_name,'',val};
    cmds(end).redo_cmd = 'layer_edit';
    cmds(end).redo_args = {val,name,group_name,'',new_val};
    
    % Push the new command(s) to the stack
    obj.undo_stack.push(cmds);
  end
  
elseif source == obj.left_panel.layerCM_bottom
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2);
    if isempty(val)
      return;
    else
      val = val(1);
    end
    name = obj.eg.layers.lyr_name{val};
    group_name = obj.eg.layers.lyr_group_name{val};
    new_val = length(obj.eg.layers.lyr_name);
    
    fprintf('Move layer bottom %s:%s\n', group_name, name);
    
    cmds = [];
    cmds(end+1).undo_cmd = 'layer_edit';
    cmds(end).undo_args = {new_val,name,group_name,'',val};
    cmds(end).redo_cmd = 'layer_edit';
    cmds(end).redo_args = {val,name,group_name,'',new_val};
    
    % Push the new command(s) to the stack
    obj.undo_stack.push(cmds);
  end
  
elseif source == obj.left_panel.layerCM_delete
  if strcmpi(obj.eg.layer_source,'layerData')
    % Get the currently selected layers.
    vals = get(obj.left_panel.layerLB,'Value');
    vals = vals(vals>2);
    if isempty(vals)
      return
    end
    
    if length(vals) > 1
      
      prompt = questdlg(sprintf('Are you sure you want to delete the %d selected layers?', ...
        length(vals)), ...
        'Delete Layer','Yes','Cancel','Cancel');
      
      switch prompt
        case 'Yes'
          cmds = [];
          % Sort in descending order or vals will change as layers get
          % deleted.
          vals = sort(vals,'descend');
          for val = vals
            name = obj.eg.layers.lyr_name{val};
            group_name = obj.eg.layers.lyr_group_name{val};
            desc = '';
            fprintf('Delete layer %s:%s "%s"\n', group_name, name, '');
            
            cmds(end+1).undo_cmd = 'layer_new';
            cmds(end).undo_args = {val,name,group_name,desc};
            cmds(end).redo_cmd = 'layer_delete';
            cmds(end).redo_args = {val};
            
            % Push the new command(s) to the stack
          end
          obj.undo_stack.push(cmds);
        case 'Cancel'
      end
      
    else
      name = obj.eg.layers.lyr_name{vals};
      group_name = obj.eg.layers.lyr_group_name{vals};
      desc = '';
      
      prompt = questdlg(sprintf('Are you sure you want to deleted layer %s:%s?', ...
        group_name,name), ...
        'Delete Layer','Yes','Cancel','Cancel');
      
      switch prompt
        case 'Yes'
          fprintf('Delete layer %s:%s "%s"\n', group_name, name, '');
          
          cmds = [];
          cmds(end+1).undo_cmd = 'layer_new';
          cmds(end).undo_args = {vals,name,group_name,desc};
          cmds(end).redo_cmd = 'layer_delete';
          cmds(end).redo_args = {vals};
          
          % Push the new command(s) to the stack
          obj.undo_stack.push(cmds);
        case 'Cancel'
      end
    end
    
  end
  
end
