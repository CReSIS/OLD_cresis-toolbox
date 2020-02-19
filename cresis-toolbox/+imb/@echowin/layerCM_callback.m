function layerCM_callback(obj,source,event)
% layerCM_callback(obj,source,event)

% Ensure focus stays on figure to prevent hotkeys registering with this
% uicontrol.
% uicontrol(obj.right_panel.status_panel.statusText);

if source == obj.left_panel.layerCM_visible || source == obj.left_panel.layerCM_hide
  val = get(obj.left_panel.layerLB,'Value');
  
  obj.eg.layers.visible_layers(val) = source == obj.left_panel.layerCM_visible;
  
  obj.layerLB_str();
  
  % Update plot based on selection
  obj.set_visibility();
elseif source == obj.left_panel.layerCM_new || source == obj.left_panel.layerCM_copy
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers. The new layer will be inserted
    % before the first of the currently selected layers or at the bottom
    % of the listbox.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2);
    % Ensure that there is at least one layer selected if copying
    if source ~= obj.left_panel.layerCM_copy || ~isempty(val)
      if isempty(val)
        val = length(obj.eg.layers.lyr_id)+1;
      else
        val = val(1);
      end
      prompt = {'Layer Age:','Description:','Layer Group Name:','Layer Name:'};
      if source == obj.left_panel.layerCM_copy
        def = {obj.eg.layers.lyr_age(val),obj.eg.layers.lyr_desc{val},char(obj.eg.layers.lyr_group_name{val}),obj.eg.layers.lyr_name{val}};
        val = val+1;
      else
        def = {'', '', '',''};
      end
      dlg_title = 'New Layer';
      num_lines = 1;
      answer = inputdlg(prompt,dlg_title,num_lines,def);
      
      if length(answer) == 4 && ~isempty(answer{4})
        try
          age = eval(answer{1});
        catch
          age = NaN;
        end
        desc = answer{2};
        group_name = answer{3};
        name = answer{4};
        
        if any(strcmpi(name,obj.eg.layers.lyr_name))
          fprintf('  Layer %s already exists\n', name);
        else
          
          fprintf('Add layer %s:%s age %g desc "%s"\n', group_name, name, age, desc);
          
          cmds = [];
          cmds(end+1).undo_cmd = 'layer_delete';
          cmds(end).undo_args = {val};
          cmds(end).redo_cmd = 'layer_new';
          cmds(end).redo_args = {val,age,desc,group_name,name};
          
          % Push the new command(s) to the stack
          obj.undo_stack.push(obj.cmds_convert_units(cmds));
        end
      end
    end
  end
  
elseif source == obj.left_panel.layerCM_edit
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2);
    if ~isempty(val)
      val = val(1);
      prompt = {'Layer Age:','Description:','Layer Group Name:','Layer Name:'};
      old_age = obj.eg.layers.lyr_age(val);
      old_desc = obj.eg.layers.lyr_desc{val};
      old_group_name = char(obj.eg.layers.lyr_group_name{val});
      old_name = obj.eg.layers.lyr_name{val};
      def = {num2str(old_age),old_desc,old_group_name,old_name};
      dlg_title = 'Edit Layer';
      num_lines = 1;
      answer = inputdlg(prompt,dlg_title,num_lines,def);
      
      if length(answer) == 4 && ~isempty(answer{4})
        try
          age = eval(answer{1});
        catch
          age = NaN;
        end
        desc = answer{2};
        group_name = answer{3};
        name = answer{4};
        
        if any(strcmpi(name,obj.eg.layers.lyr_name([1:val-1 val+1:end])))
          fprintf('  Layer %s already exists\n', name);
        else
          
          fprintf('Edit layer %s:%s to %s:%s age: %g desc: "%s"\n', old_group_name, old_name, group_name, name, age, desc);
          
          cmds = [];
          cmds(end+1).undo_cmd = 'layer_edit';
          cmds(end).undo_args = {val,old_age,old_desc,old_group_name,old_name,val};
          cmds(end).redo_cmd = 'layer_edit';
          cmds(end).redo_args = {val,age,desc,group_name,name,val};
          
          % Push the new command(s) to the stack
          obj.undo_stack.push(obj.cmds_convert_units(cmds));
        end
      end
    end
  end
  
elseif source == obj.left_panel.layerCM_up
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>3);
    if ~isempty(val)
      val = val(1);
      age = obj.eg.layers.lyr_age(val);
      desc = obj.eg.layers.lyr_desc{val};
      group_name = obj.eg.layers.lyr_group_name{val};
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      
      new_val = val-1;
      new_age = obj.eg.layers.lyr_age(new_val);
      new_desc = obj.eg.layers.lyr_desc{new_val};
      new_group_name = obj.eg.layers.lyr_group_name{new_val};
      new_name = obj.eg.layers.lyr_name{new_val};
      new_order = obj.eg.layers.lyr_order(new_val);
      
      fprintf('Move layer up %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {val,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {val,age,desc,group_name,name,new_order};
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {new_val,new_age,new_desc,new_group_name,new_name,new_order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {new_val,new_age,new_desc,new_group_name,new_name,order};
      
      % Push the new command(s) to the stack
      obj.undo_stack.push(obj.cmds_convert_units(cmds));
    end
  end
  
elseif source == obj.left_panel.layerCM_down
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2 & val<length(obj.eg.layers.lyr_name));
    if ~isempty(val)
      val = val(1);
      age = obj.eg.layers.lyr_age(val);
      desc = obj.eg.layers.lyr_desc{val};
      group_name = obj.eg.layers.lyr_group_name{val};
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      
      new_val = val+1;
      new_age = obj.eg.layers.lyr_age(new_val);
      new_desc = obj.eg.layers.lyr_desc{new_val};
      new_group_name = obj.eg.layers.lyr_group_name{new_val};
      new_name = obj.eg.layers.lyr_name{new_val};
      new_order = obj.eg.layers.lyr_order(new_val);
      
      fprintf('Move layer down %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {val,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {val,age,desc,group_name,name,new_order};
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {new_val,new_age,new_desc,new_group_name,new_name,new_order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {new_val,new_age,new_desc,new_group_name,new_name,order};
      
      % Push the new command(s) to the stack
      obj.undo_stack.push(obj.cmds_convert_units(cmds));
    end
  end
  
elseif source == obj.left_panel.layerCM_top
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>3);
    if ~isempty(val)
      val = val(1);
      age = obj.eg.layers.lyr_age(val);
      desc = obj.eg.layers.lyr_desc{val};
      group_name = obj.eg.layers.lyr_group_name{val};
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      
      new_val = 3;
      new_age = obj.eg.layers.lyr_age(new_val);
      new_desc = obj.eg.layers.lyr_desc{new_val};
      new_group_name = obj.eg.layers.lyr_group_name{new_val};
      new_name = obj.eg.layers.lyr_name{new_val};
      new_order = obj.eg.layers.lyr_order(new_val);
      
      fprintf('Move layer top %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {val,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {val,age,desc,group_name,name,new_order};
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {new_val,new_age,new_desc,new_group_name,new_name,new_order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {new_val,new_age,new_desc,new_group_name,new_name,order};
      
      % Push the new command(s) to the stack
      obj.undo_stack.push(obj.cmds_convert_units(cmds));
    end
  end
  
elseif source == obj.left_panel.layerCM_bottom
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers.
    val = get(obj.left_panel.layerLB,'Value');
    val = val(val>2);
    if ~isempty(val)
      val = val(1);
      age = obj.eg.layers.lyr_age(val);
      desc = obj.eg.layers.lyr_desc{val};
      group_name = obj.eg.layers.lyr_group_name{val};
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      
      new_val = length(obj.eg.layers.lyr_name);
      new_age = obj.eg.layers.lyr_age(new_val);
      new_desc = obj.eg.layers.lyr_desc{new_val};
      new_group_name = obj.eg.layers.lyr_group_name{new_val};
      new_name = obj.eg.layers.lyr_name{new_val};
      new_order = obj.eg.layers.lyr_order(new_val);
      
      fprintf('Move layer bottom %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {val,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {val,age,desc,group_name,name,new_order};
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {new_val,new_age,new_desc,new_group_name,new_name,new_order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {new_val,new_age,new_desc,new_group_name,new_name,order};
      
      % Push the new command(s) to the stack
      obj.undo_stack.push(obj.cmds_convert_units(cmds));
    end
  end
  
elseif source == obj.left_panel.layerCM_delete
  if strcmpi(obj.eg.layers.source,'layerData')
    % Get the currently selected layers.
    vals = get(obj.left_panel.layerLB,'Value');
    vals = vals(vals>2);
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
            age = obj.eg.layers.lyr_age(val);
            desc = obj.eg.layers.lyr_desc{val};
            group_name = obj.eg.layers.lyr_group_name{val};
            name = obj.eg.layers.lyr_name{val};
            fprintf('Delete layer %s:%s\n', group_name, name);
            
            cmds(end+1).undo_cmd = 'layer_new';
            cmds(end).undo_args = {val,age,desc,group_name,name};
            cmds(end).redo_cmd = 'layer_delete';
            cmds(end).redo_args = {val};
            
            % Push the new command(s) to the stack
          end
          obj.undo_stack.push(obj.cmds_convert_units(cmds));
        case 'Cancel'
      end
      
    elseif length(vals) == 1
      age = obj.eg.layers.lyr_age(vals);
      desc = obj.eg.layers.lyr_desc{vals};
      group_name = obj.eg.layers.lyr_group_name{vals};
      name = obj.eg.layers.lyr_name{vals};
      
      prompt = questdlg(sprintf('Are you sure you want to delete layer %s:%s?', ...
        group_name,name), ...
        'Delete Layer','Yes','Cancel','Cancel');
      
      switch prompt
        case 'Yes'
          fprintf('Delete layer %s:%s\n', group_name, name);
          
          cmds = [];
          cmds(end+1).undo_cmd = 'layer_new';
          cmds(end).undo_args = {vals,age,desc,group_name,name};
          cmds(end).redo_cmd = 'layer_delete';
          cmds(end).redo_args = {vals};
          
          % Push the new command(s) to the stack
          obj.undo_stack.push(obj.cmds_convert_units(cmds));
        case 'Cancel'
      end
    end
    
  end
  
end
