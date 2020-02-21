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
        def = {sprintf('%g',obj.eg.layers.lyr_age(val)),obj.eg.layers.lyr_desc{val},char(obj.eg.layers.lyr_group_name{val}),obj.eg.layers.lyr_name{val}};
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
        if isempty(age)
          age = NaN;
        end
        desc = answer{2};
        group_name = answer{3};
        name = answer{4};
        
        if any(strcmpi(name,obj.eg.layers.lyr_name))
          fprintf('  Layer %s already exists\n', name);
        else
          
          % Get id for new layer
          new_lyr_id = max(obj.eg.layers.lyr_id) + 1;
          id = new_lyr_id;
          
          fprintf('Add layer %s:%s age %g desc "%s"\n', group_name, name, age, desc);
          if val > length(obj.eg.layers.lyr_id)
            order = obj.eg.layers.lyr_order(end) + 1;
          else
            order = obj.eg.layers.lyr_order(val);
          end
          
          cmds = [];
          cmds(end+1).undo_cmd = 'layer_delete';
          cmds(end).undo_args = {id};
          cmds(end).redo_cmd = 'layer_new';
          cmds(end).redo_args = {id,age,desc,group_name,name,order};
          
          for val = val:length(obj.eg.layers.lyr_id)
            age = obj.eg.layers.lyr_age(val);
            desc = obj.eg.layers.lyr_desc{val};
            group_name = obj.eg.layers.lyr_group_name{val};
            id = obj.eg.layers.lyr_id(val);
            name = obj.eg.layers.lyr_name{val};
            order = obj.eg.layers.lyr_order(val);
            if val < length(obj.eg.layers.lyr_order)
              new_order = obj.eg.layers.lyr_order(val+1);
            else
              new_order = obj.eg.layers.lyr_order(end) + 1;
            end
            
            cmds(end+1).undo_cmd = 'layer_edit';
            cmds(end).undo_args = {id,age,desc,group_name,name,order};
            cmds(end).redo_cmd = 'layer_edit';
            cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
          end
          cmds = cmds([2:end 1]);
          
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
      old_id = obj.eg.layers.lyr_id(val);
      old_name = obj.eg.layers.lyr_name{val};
      old_order = obj.eg.layers.lyr_order(val);
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
        if isempty(age)
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
          cmds(end).undo_args = {old_id,old_age,old_desc,old_group_name,old_name,old_order};
          cmds(end).redo_cmd = 'layer_edit';
          cmds(end).redo_args = {old_id,age,desc,group_name,name,old_order};
          
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
      id = obj.eg.layers.lyr_id(val);
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      new_order = obj.eg.layers.lyr_order(val-1);
      
      fprintf('Move layer up %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {id,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
      val = val-1;
      age = obj.eg.layers.lyr_age(val);
      desc = obj.eg.layers.lyr_desc{val};
      group_name = obj.eg.layers.lyr_group_name{val};
      id = obj.eg.layers.lyr_id(val);
      name = obj.eg.layers.lyr_name{val};
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {id,age,desc,group_name,name,new_order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {id,age,desc,group_name,name,order};
      cmds = cmds([2:end 1]);
      
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
      id = obj.eg.layers.lyr_id(val);
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      new_order = obj.eg.layers.lyr_order(val+1);
      
      fprintf('Move layer down %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {id,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
      val = val+1;
      age = obj.eg.layers.lyr_age(val);
      desc = obj.eg.layers.lyr_desc{val};
      id = obj.eg.layers.lyr_id(val);
      group_name = obj.eg.layers.lyr_group_name{val};
      name = obj.eg.layers.lyr_name{val};
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {id,age,desc,group_name,name,new_order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {id,age,desc,group_name,name,order};
      cmds = cmds([2:end 1]);
      
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
      id = obj.eg.layers.lyr_id(val);
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      
      new_val = 3;
      new_order = obj.eg.layers.lyr_order(new_val);
      
      fprintf('Move layer top %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {id,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
      for val = new_val:val-1
        age = obj.eg.layers.lyr_age(val);
        desc = obj.eg.layers.lyr_desc{val};
        group_name = obj.eg.layers.lyr_group_name{val};
        id = obj.eg.layers.lyr_id(val);
        name = obj.eg.layers.lyr_name{val};
        order = obj.eg.layers.lyr_order(val);
        new_order = obj.eg.layers.lyr_order(val+1);
        cmds(end+1).undo_cmd = 'layer_edit';
        cmds(end).undo_args = {id,age,desc,group_name,name,order};
        cmds(end).redo_cmd = 'layer_edit';
        cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
      end
      cmds = cmds([2:end 1]);
      
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
      id = obj.eg.layers.lyr_id(val);
      name = obj.eg.layers.lyr_name{val};
      order = obj.eg.layers.lyr_order(val);
      
      new_val = length(obj.eg.layers.lyr_name);
      new_order = obj.eg.layers.lyr_order(new_val);
      
      fprintf('Move layer bottom %s:%s\n', group_name, name);
      
      cmds = [];
      cmds(end+1).undo_cmd = 'layer_edit';
      cmds(end).undo_args = {id,age,desc,group_name,name,order};
      cmds(end).redo_cmd = 'layer_edit';
      cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
      for val = val+1:new_val
        age = obj.eg.layers.lyr_age(val);
        desc = obj.eg.layers.lyr_desc{val};
        group_name = obj.eg.layers.lyr_group_name{val};
        id = obj.eg.layers.lyr_id(val);
        name = obj.eg.layers.lyr_name{val};
        order = obj.eg.layers.lyr_order(val);
        new_order = obj.eg.layers.lyr_order(val-1);
        cmds(end+1).undo_cmd = 'layer_edit';
        cmds(end).undo_args = {id,age,desc,group_name,name,order};
        cmds(end).redo_cmd = 'layer_edit';
        cmds(end).redo_args = {id,age,desc,group_name,name,new_order};
      end
      cmds = cmds([2:end 1]);
      
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
            id = obj.eg.layers.lyr_id(val);
            name = obj.eg.layers.lyr_name{val};
            order = obj.eg.layers.lyr_order(val);
            fprintf('Delete layer %s:%s\n', group_name, name);
            
            cmds(end+1).undo_cmd = 'layer_new';
            cmds(end).undo_args = {id,age,desc,group_name,name,order};
            cmds(end).redo_cmd = 'layer_delete';
            cmds(end).redo_args = {id};
            
            % Push the new command(s) to the stack
          end
          obj.undo_stack.push(obj.cmds_convert_units(cmds));
        case 'Cancel'
      end
      
    elseif length(vals) == 1
      age = obj.eg.layers.lyr_age(vals);
      desc = obj.eg.layers.lyr_desc{vals};
      group_name = obj.eg.layers.lyr_group_name{vals};
      id = obj.eg.layers.lyr_id(vals);
      name = obj.eg.layers.lyr_name{vals};
      order = obj.eg.layers.lyr_order(vals);
      
      prompt = questdlg(sprintf('Are you sure you want to delete layer %s:%s?', ...
        group_name,name), ...
        'Delete Layer','Yes','Cancel','Cancel');
      
      switch prompt
        case 'Yes'
          fprintf('Delete layer %s:%s\n', group_name, name);
          
          cmds = [];
          cmds(end+1).undo_cmd = 'layer_new';
          cmds(end).undo_args = {id,age,desc,group_name,name,order};
          cmds(end).redo_cmd = 'layer_delete';
          cmds(end).redo_args = {id};
          
          % Push the new command(s) to the stack
          obj.undo_stack.push(obj.cmds_convert_units(cmds));
        case 'Cancel'
      end
    end
    
  end
  
end
