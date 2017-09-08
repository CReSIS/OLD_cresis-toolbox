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
      obj.create_option_ui();
      obj.tool_name = 'Quality';
      obj.tool_menu_name = '(Q)uality';
      obj.tool_shortcut = 'q';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 'q: Quality';
    end
    
    function cmd = apply_PB_callback(obj,sb,slices)
      cmd = [];
      selection_index = find(sb.select_mask);
      surf_layer_idx = find(strcmp('surf_quality',{sb.layer.name}));
      bottom_layer_idx = find(strcmp('bottom_quality',{sb.layer.name}));
      
      if sb.gui.layerLB.Value == 1 || sb.gui.layerLB.Value == surf_layer_idx
        % Surface selected
        y_old = sb.layer(1).y(selection_index,sb.slice);
        isSurf = true;
        for sel_idx = 1:length(selection_index)
          if sb.qualityShift
            sb.layer(surf_layer_idx).y(selection_index(sel_idx),sb.slice) ...
              = 1;
          else
            sb.layer(surf_layer_idx).y(selection_index(sel_idx),sb.slice) ...
              = 0;
          end
        end
      elseif sb.gui.layerLB.Value == 2 || sb.gui.layerLB.Value == bottom_layer_idx
        % Bottom selected
        y_old = sb.layer(2).y(selection_index,sb.slice);
        isSurf = false;
        for sel_idx = 1:length(selection_index)
          if sb.qualityShift
            sb.layer(bottom_layer_idx).y(selection_index(sel_idx),sb.slice) ...
              = 1;
          else
            sb.layer(bottom_layer_idx).y(selection_index(sel_idx),sb.slice) ...
              = 0;
          end
        end
      else return;
      end
      sb.update_slice();
      
      % Select mask
      sb.gui.select_maskCB = uicontrol('parent',sb.h_fig);
      set(sb.gui.select_maskCB,'style','checkbox')
      set(sb.gui.select_maskCB,'string','Select')
      set(sb.gui.select_maskCB,'value',1)
      set(sb.gui.select_maskCB,'TooltipString','Check to operate only on the selected region.');
      
      % Create cmd for layer change
      %         active_idx = sb.layer(sb.layer_idx).active_layer;
      %active_idx = 8;%sb.gui.layerLB.Value;
      
      if isSurf
        active_idx = surf_layer_idx;
      else
        active_idx = bottom_layer_idx;
      end
      
      if get(sb.gui.select_maskCB,'Value')
        cols = find(sb.select_mask);
      else
        cols = 1:size(sb.data,2);
      end
      cmd{end+1}.undo.slice = sb.slice;
      cmd{end}.redo.slice = sb.slice;
      cmd{end}.undo.layer = active_idx;
      cmd{end}.redo.layer = active_idx;
      cmd{end}.undo.x = cols;
      cmd{end}.undo.y = y_old;
      cmd{end}.redo.x = cols;
      cmd{end}.redo.y = sb.layer(active_idx).y(cols,sb.slice);
      cmd{end}.type = 'standard';
      
      %         cmd{end+1}.redo.slice = sb.slice;
      %         cmd{end}.undo.slice = sb.slice;
      %         cmd{end}.type = 'slice_dummy';
    end
  end
end
