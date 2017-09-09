classdef (HandleCompatible = true) slicetool_quality < imb.slicetool
  % Slice_browser tool for quality
  % Change quality parameter for selected portion of slice.
  %
  % Currently only supports 0 (bad data) and 1 (good data).
  %
  % Authors: Victor Berger and Mohanad Al-Ibadi
  
  properties (SetAccess = protected, GetAccess = public)
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_quality()
      obj.tool_name = 'Quality';
      obj.tool_menu_name = '(Q)uality';
      obj.tool_shortcut = 'q';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = [0 1];
      obj.help_string = 'q: Set quality to bad, Q: Set quality to good';
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

      % If no quality layer defined, then do nothing
      if isempty(quality_idx)
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
      if sb.shift_pressed
        % Set to good
        new_vals = ones(length(cols),1);
      else
        % Set to bad
        new_vals = zeros(length(cols),1);
      end
      
      % Create cmd for layer change
      cmd = [];
      for idx = 1:length(slices)
        slice = slices(idx);
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.layer = quality_idx;
        cmd{end}.redo.layer = quality_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.layer(quality_idx).y(cols,slice);
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
