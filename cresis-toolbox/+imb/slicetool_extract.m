classdef (HandleCompatible = true) slicetool_extract < imb.slicetool
  % Slice_browser tool which calls detect.cpp (HMM)
  
  properties (SetAccess = protected, GetAccess = public)
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_extract()
      obj.create_option_ui();
      obj.tool_name = 'Extract';
      obj.tool_menu_name = '(E)xtract';
      obj.tool_shortcut = 'e';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 'e: Extract/refine tools which run TRWS solution to HMM inference model to find best layer. Neighboring slices effect cost function to improve solution.';
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
      
      try
        slice_range = eval(get(obj.gui.slice_rangeLE,'String'));
      catch ME
        error('Error in slice range: %s', ME.getReport);
      end
      try
        num_loops = eval(get(obj.gui.numloopsLE,'String'));
      catch ME
        error('Error in number of loops: %s', ME.getReport);
      end
      try
        threshold = eval(get(obj.gui.thresholdLE,'String'));
      catch ME
        error('Error in threshold: %s', ME.getReport);
      end
      try
        smooth = eval(get(obj.gui.smoothLE,'String'));
      catch ME
        error('Error in smooth: %s', ME.getReport);
      end
      if get(obj.gui.select_maskCB,'Value')
        cols = find(sb.select_mask);
      else
        cols = 1:size(sb.data,2);
      end
      left_edge_en = get(obj.gui.leftEdgeCB,'Value');
      right_edge_en = get(obj.gui.rightEdgeCB,'Value');
      
      if ~exist('slices','var') || isempty(slices)
        slice_range = min(slice_range):max(slice_range);
        slices = sb.slice+slice_range;
      end
      slices = intersect(slices,1:size(sb.data,3));
      if ~left_edge_en
        start_slice_idx = 1;
      else
        start_slice_idx = 2;
      end
      if ~right_edge_en
        end_slice_idx = length(slices);
      else
        end_slice_idx = length(slices)-1;
      end
      fprintf('Apply %s to layer %d slices %d - %d\n', obj.tool_name, active_idx, slices(1), slices(end));
      
      gt = [];
      if ~isempty(control_idx)
        % Create ground truth input
        % 1. Each column is one ground truth input
        % 2. Row 1: relative slice/range-line, Row 2: x, Row 3: y
        for idx = 1:length(slices)
          slice = slices(idx);
          mask = isfinite(sb.layer(control_idx).x(:,slice)) ...
            & isfinite(sb.layer(control_idx).y(:,slice));
          gt = cat(2,gt,[(idx-1)*ones(1,sum(mask)); ...
            sb.layer(control_idx).x(mask,slice).'-1; ...
            sb.layer(control_idx).y(mask,slice).'+0.5]);
          bottom_bin = sb.layer(control_idx).y(33,slices);
        end
      else
        bottom_bin = NaN*zeros(1,length(slices));
      end
      
      if isempty(surf_idx)
        %error('extract cannot be run without a surface layer');
        surf_bins = NaN*sb.layer(active_idx).y(:,slices);
      else
        surf_bins = sb.layer(surf_idx).y(:,slices);
      end
      surf_bins(isnan(surf_bins)) = -1;
      
      bottom_bin(isnan(bottom_bin)) = -1;
      
      if isempty(mask_idx)
        mask = ones(size(sb.data,2),length(slices));
      else
        mask = sb.layer(mask_idx).y(:,slices);
      end
      
      begin_slice = max(1, min(slices)-1);
      end_slice = min(size(sb.data,3), max(slices)+1);
      edge = [sb.layer(active_idx).y(:,begin_slice), sb.layer(active_idx).y(:,end_slice)];
      edge(mask(:,1) == 0,1) = surf_bins(mask(:,1) == 0,1);
      edge(mask(:,end) == 0,2) = surf_bins(mask(:,end) == 0,end);
      
      extract_data = sb.data(:,:,slices);
      extract_data(extract_data>threshold) = threshold;
      if ~left_edge_en
        edge(:,1) = -1;
      end
      if ~right_edge_en
        edge(:,end) = -1;
      end
      smooth_slope = [];
      smooth_weight = smooth(1:2);
      smooth_var = smooth(3);
      mu = [-3 -1 3 3 3 4 12 20 12 4 3 3 3 -1 -3];
      sigma = 3*ones(1,15);
      % mu = obj.custom_data.mu;
      % sigma = obj.custom_data.sigma;
      mask = 90*fir_dec(fir_dec(double(shrink(mask,2)),ones(1,5)/3.7).',ones(1,5)/3.7).';
      mask(mask>=90) = inf;
      if 0
        %% DEBUG: For running mex function in debug mode
        save('/tmp/mex_inputs.mat','extract_data','surf_bins','bottom_bin','gt','mask','mu','sigma','smooth_slope','smooth_weight','smooth_var','edge');
        load('/tmp/mex_inputs.mat');
        correct_surface = tomo.extract(double(extract_data), ...
          double(surf_bins), double(bottom_bin), ....
          double(gt), double(mask), ...
          double(mu), double(sigma), smooth_weight, smooth_var, double(smooth_slope));
      end
      
      correct_surface = tomo.refine(double(extract_data), ...
        double(surf_bins), double(bottom_bin), ....
        double(gt), double(mask), ...
        double(mu), double(sigma), smooth_weight, smooth_var, double(smooth_slope),...
        double(edge));
        
     % Create cmd for layer change
      cmd = [];
      for idx = start_slice_idx:end_slice_idx
        slice = slices(idx);
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.layer = active_idx;
        cmd{end}.redo.layer = active_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.layer(active_idx).y(cols,slice);
        cmd{end}.redo.x = cols;
        cmd{end}.redo.y = correct_surface(cols,idx);
        cmd{end}.type = 'standard';
      end
      cmd{end+1}.redo.slice = sb.slice;
      cmd{end}.undo.slice = sb.slice;
      cmd{end}.type = 'slice_dummy';
      
    end
    
    function set_custom_data(obj,custom_data)
      obj.custom_data.mu = mean(custom_data.mu);
      obj.custom_data.sigma = mean(custom_data.sigma);
    end
    
    function create_option_ui(obj)
      obj.h_fig = figure('Visible','off','DockControls','off', ...
        'NumberTitle','off','ToolBar','none','MenuBar','none','Resize','off');
      if strcmpi(class(obj.h_fig),'double')
        set(obj.h_fig,'Name',sprintf('%d: extract tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: extract tool prefs',obj.h_fig.Number));
      end
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      pos = get(obj.h_fig,'Position');
      pos(3) = 200;
      pos(4) = 140;
      set(obj.h_fig,'Position',pos);
      
      % Number of loops
      obj.gui.numloopsTXT = uicontrol('Style','text','string','Loops');
      
      obj.gui.numloopsLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.numloopsLE,'style','edit')
      set(obj.gui.numloopsLE,'string','50')
      set(obj.gui.numloopsLE,'TooltipString','Number of iterations. IGNORED.');
      
      % Slice range
      obj.gui.slice_rangeTXT = uicontrol('Style','text','string','Slice range');
      set(obj.gui.slice_rangeTXT,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:7".');
      
      obj.gui.slice_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.slice_rangeLE,'style','edit')
      set(obj.gui.slice_rangeLE,'string','-5:7')
      set(obj.gui.slice_rangeLE,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:7".');
      
      % Threshold
      obj.gui.thresholdTXT = uicontrol('Style','text','string','Threshold');
      set(obj.gui.thresholdTXT,'TooltipString','Specify an image threshold.');
      
      obj.gui.thresholdLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.thresholdLE,'style','edit')
      set(obj.gui.thresholdLE,'string','13.5')
      set(obj.gui.thresholdLE,'TooltipString','Specify an image threshold.');
      
      % Smooth
      obj.gui.smoothTXT = uicontrol('Style','text','string','Smoothness');
      set(obj.gui.smoothTXT,'TooltipString','Specify layer smoothness ["slice weight" "DOA weight" "edge weight"].');
      
      obj.gui.smoothLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.smoothLE,'style','edit')
      set(obj.gui.smoothLE,'string','[32 16 12]')
      set(obj.gui.smoothLE,'TooltipString','Specify layer smoothness ["slice weight" "DOA weight" "edge weight"].');
      
      % Select mask
      obj.gui.select_maskCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.select_maskCB,'style','checkbox')
      set(obj.gui.select_maskCB,'string','Select')
      set(obj.gui.select_maskCB,'value',1)
      set(obj.gui.select_maskCB,'TooltipString','Check to operate only on the selected region.');
      
      % Left edge
      obj.gui.leftEdgeCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.leftEdgeCB,'style','checkbox')
      set(obj.gui.leftEdgeCB,'string','Left')
      set(obj.gui.leftEdgeCB,'value',1)
      set(obj.gui.leftEdgeCB,'TooltipString','Check to use a left edge boundary condition.');
      
      % Right edge
      obj.gui.rightEdgeCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.rightEdgeCB,'style','checkbox')
      set(obj.gui.rightEdgeCB,'string','Right')
      set(obj.gui.rightEdgeCB,'value',0)
      set(obj.gui.rightEdgeCB,'TooltipString','Check to use a right edge boundary condition.');
      
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
      obj.gui.table.handles{row,col}   = obj.gui.slice_rangeTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.slice_rangeLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.thresholdTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.thresholdLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;

      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.smoothTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.smoothLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.select_maskCB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.leftEdgeCB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 3;
      obj.gui.table.handles{row,col}   = obj.gui.rightEdgeCB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.table);      
      
    end

  end
  
end


