classdef (HandleCompatible = true) slicetool_delete< imb.slicetool
  % Slice_browser tool for deleting
  %
  % Deletes data on control layer for selected portion of slice.
  %
  % Authors: John Paden
  
  properties (SetAccess = protected, GetAccess = public)
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_delete()
      obj.tool_name = 'Delete';
      obj.tool_menu_name = 'Delete';
      obj.tool_shortcut = '';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 'Delete points';
    end
    
    function cmd = apply_PB_callback(obj,sb,slices)
      % sb: slice browser object. Use the following fields to create
      %     commands, cmd, that use sb.data to operate on sb.layer. You 
      %     should not modify any fields of sb.
      %  .layer: struct array containing layer information
      %  .data: 3D image
      %  .slice: current slice in 3D image (third index of .data)
      %  .layer_idx: active layer
      % slices: array of slices to operate on (overrides sb.slice)
      control_idx = sb.layer(sb.layer_idx).control_layer;
      active_idx = sb.layer(sb.layer_idx).active_layer;
      surf_idx = sb.layer(sb.layer_idx).surf_layer;
      mask_idx = sb.layer(sb.layer_idx).mask_layer;
      quality_idx = sb.layer(sb.layer_idx).quality_layer;

      % If no control layer defined, then do nothing
      if isempty(control_idx)
        cmd = [];
        return
      end

      % Make sure slices is a valid contiguous list
      slice_range = 0;
      if ~exist('slices','var') || isempty(slices)
        slice_range = min(slice_range):max(slice_range);
        slices = sb.slice+slice_range;
      end
      slices = intersect(slices,1:size(sb.data,3));

      % Determine selection mask
      cols = find(sb.select_mask);

      % Create new quality values
      new_vals = NaN*zeros(length(cols),1);
      
      % Create cmd for layer change
      cmd = [];
      for idx = 1:length(slices)
        slice = slices(idx);
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.layer = control_idx;
        cmd{end}.redo.layer = control_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.layer(control_idx).y(cols,slice);
        cmd{end}.redo.x = cols;
        cmd{end}.redo.y = new_vals;
        cmd{end}.type = 'standard';
      end
      cmd{end+1}.redo.slice = sb.slice;
      cmd{end}.undo.slice = sb.slice;
      cmd{end}.type = 'slice_dummy';
      
    end
  end
end
