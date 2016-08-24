classdef (HandleCompatible = true) slicetool_extract < imb.slicetool
  % Slice_browser tool which calls detect.cpp (HMM)
  
  properties
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
  end
  
  methods
    function obj = slicetool_extract()
      obj.create_option_ui();
      obj.tool_name = 'Extract';
      obj.tool_menu_name = '(E)xtract';
      obj.tool_shortcut = 'e';
    end
    
    function cmd = apply_PB_callback(obj,sb)
      % Read and set sb.layer
      % Read sb.data
      % Read sb.slice
      % Read sb.layer_idx
      fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, sb.layer_idx, sb.slice);
      
      control_idx = sb.layer(sb.layer_idx).control_layer;
      surf_idx = sb.layer(sb.layer_idx).surf_layer;
      try
        extract_range = eval(get(obj.gui.extentLE,'String'));
      catch ME
        error('Error in slice range: %s', ME.getReport);
      end
      try
        num_loops = eval(get(obj.gui.numloopsLE,'String'));
      catch ME
        error('Error in number of loops: %s', ME.getReport);
      end
      
      slices = sb.slice+extract_range;
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
      
      surf_bins = sb.layer(surf_idx).y(:,slices);
      surf_bins(isnan(surf_bins)) = -1;
      
      %bottom_bin = obj.custom_data.bottom(slices);
      bottom_bin = sb.layer(control_idx).y(33,slices);
      bottom_bin(isnan(bottom_bin)) = -1;

      correct_surface = tomo.extract(double(sb.data(:,:,slices)), ...
        double(surf_bins), double(bottom_bin), ...
        double(gt), double(obj.custom_data.ice_mask(:,slices)), ...
        double(obj.custom_data.mu), double(obj.custom_data.sigma));
      correct_surface = reshape(correct_surface, [size(sb.data,2) length(slices)]);
      % Create cmd for layer change
      cmd = [];
      for idx = 2:length(slices)-1
        slice = slices(idx);
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.layer = sb.layer_idx;
        cmd{end}.redo.layer = sb.layer_idx;
        cmd{end}.undo.x = 1:size(sb.layer(sb.layer_idx).y,1);
        cmd{end}.undo.y = sb.layer(sb.layer_idx).y(:,slice);
        cmd{end}.redo.x = 1:size(sb.layer(sb.layer_idx).y,1);
        cmd{end}.redo.y = correct_surface(:,idx);
        cmd{end}.type = 'standard';
      end
        cmd{end+1}.redo.slice = sb.slice;
        cmd{end}.undo.slice = sb.slice;
        cmd{end}.type = 'slice_dummy';
      
    end
    
    function set_custom_data(obj,custom_data)
      obj.custom_data.mu = mean(custom_data.mu);
      obj.custom_data.sigma = mean(custom_data.sigma);
      obj.custom_data.ice_mask = custom_data.ice_mask;
      obj.custom_data.bottom = custom_data.bottom;
    end
    
    function create_option_ui(obj)
      obj.h_fig = figure('Visible','off','DockControls','off', ...
        'NumberTitle','off','ToolBar','none','MenuBar','none','Resize','off');
      if strcmpi(class(obj.h_fig),'double')
        set(obj.h_fig,'Name',sprintf('%d: extent tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: extent tool prefs',obj.h_fig.Number));
      end
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      pos = get(obj.h_fig,'Position');
      pos(3) = 200;
      pos(4) = 50;
      set(obj.h_fig,'Position',pos);
      
      % Number of loops
      obj.gui.numloopsTXT = uicontrol('Style','text','string','Loops');
      
      obj.gui.numloopsLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.numloopsLE,'style','edit')
      set(obj.gui.numloopsLE,'string','70')
      set(obj.gui.numloopsLE,'TooltipString','Number of iterations.');
      
      % Extent
      obj.gui.extentTXT = uicontrol('Style','text','string','Slice range');
      
      obj.gui.extentLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.extentLE,'style','edit')
      set(obj.gui.extentLE,'string','-5:5')
      set(obj.gui.extentLE,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:5".');
      
      % GUI container table
      obj.gui.table.ui = obj.h_fig;
      obj.gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.table.height_margin = NaN*zeros(30,30);
      obj.gui.table.false_width = NaN*zeros(30,30);
      obj.gui.table.false_height = NaN*zeros(30,30);
      obj.gui.table.offset = [0 0];
      row = 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.numloopsTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.numloopsLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.extentTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.extentLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.table);      
      
    end

  end
  
end


