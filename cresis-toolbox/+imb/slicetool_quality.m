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
      if sb.layer(sb.layer_idx).active_layer == 1
        % Surface selected
        for sel_idx = 1:length(selection_index)
          if sb.qualityShift
            sb.layer(surf_layer_idx).y(selection_index(sel_idx),sb.slice) ...
              = 1;
          else
            sb.layer(surf_layer_idx).y(selection_index(sel_idx),sb.slice) ...
              = 0;
          end
        end
      elseif sb.layer(sb.layer_idx).active_layer == 2
        % Bottom selected
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
    end
  end
end
