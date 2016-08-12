classdef (HandleCompatible = true) slicetool_detect < imb.slicetool
  % Slice_browser tool which calls detect.cpp (HMM)
  
  properties
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
  end
  
  methods
    function obj = slicetool_detect()
      obj.create_option_ui();
      obj.tool_name = 'Detect';
      obj.tool_menu_name = '(D)etect';
      obj.tool_shortcut = 'd';
    end
    
    function cmd = apply_PB_callback(obj,sb)
      % Read and set sb.layer
      % Read sb.data
      % Read sb.slice
      % Read sb.layer_idx
      fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, sb.layer_idx, sb.slice);
      
      control_idx = sb.layer(sb.layer_idx).control_layer;
      surf_idx = sb.layer(sb.layer_idx).surf_layer;
      
      slices = sb.slice+[0];
      slices = intersect(slices,1:size(sb.data,3));
      
      % Create ground truth input
      % 1. Each column is one ground truth input
      % 2. Row 1: relative slice/range-line, Row 2: x, Row 3: y
      gt = [];
      for idx = 1:length(slices)
        slice = slices(idx);
        mask = isfinite(sb.layer(control_idx).x(:,slice)) ...
          & isfinite(sb.layer(control_idx).y(:,slice));
        gt = cat(2,gt,[(idx-1)*ones(1,sum(mask)); ...
          sb.layer(control_idx).x(mask,slice).'; ...
          sb.layer(control_idx).y(mask,slice).']);
      end
      gt = gt(2:3,:);
      
      surf_bins = sb.layer(surf_idx).y(:,sb.slice);
      surf_bins(isnan(surf_bins)) = -1;
      
      %bottom_bin = obj.custom_data.bottom(sb.slice);
      bottom_bin = sb.layer(control_idx).y(33,sb.slice);
      bottom_bin(isnan(bottom_bin)) = -1;
      
      labels = tomo.detect(sb.data(:,:,sb.slice), ...
        double(surf_bins), double(bottom_bin), ...
        double(gt), double(sb.layer(sb.layer(sb.layer_idx).mask_layer).y(:,sb.slice)), ...
        double(obj.custom_data.mu), double(obj.custom_data.sigma));
      
      % Create cmd for layer change
      cmd = [];
      cmd{1}.undo.slice = sb.slice;
      cmd{1}.redo.slice = sb.slice;
      cmd{1}.undo.layer = sb.layer_idx;
      cmd{1}.redo.layer = sb.layer_idx;
      cmd{1}.undo.x = 1:size(sb.layer(sb.layer_idx).y,1);
      cmd{1}.undo.y = sb.layer(sb.layer_idx).y(:,sb.slice);
      cmd{1}.redo.x = 1:size(sb.layer(sb.layer_idx).y,1);
      cmd{1}.redo.y = labels;
      cmd{1}.type = 'standard';
    end
    
    function set_custom_data(obj,custom_data)
      obj.custom_data.mu = mean(custom_data.mu);
      obj.custom_data.sigma = mean(custom_data.sigma);
      obj.custom_data.ice_mask = custom_data.ice_mask;
      obj.custom_data.bottom = custom_data.bottom;
    end
    
    function create_option_ui(obj)
    end

  end
  
end


