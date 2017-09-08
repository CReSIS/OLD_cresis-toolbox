classdef (HandleCompatible = true) slicetool_detect < imb.slicetool
  % Slice_browser tool which calls detect.cpp (HMM)
  
  properties (SetAccess = protected, GetAccess = public)
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_detect()
      obj.create_option_ui();
      obj.tool_name = 'Detect';
      obj.tool_menu_name = '(D)etect';
      obj.tool_shortcut = 'd';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 'd: Detect tool which runs viterbi solution to HMM inference model to find best layer. Neighboring slices have no influence on solution.';
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
        threshold = eval(get(obj.gui.thresholdLE,'String'));
      catch ME
        error('Error in threshold: %s', ME.getReport);
      end
      try
        egt_weight = eval(get(obj.gui.gt_weightLE,'String'));
      catch ME
        error('Error in ground truth weight: %s', ME.getReport);
      end
      try
        slope = eval(get(obj.gui.slopeLE,'String'));
      catch ME
        error('Error in slope: %s', ME.getReport);
      end
      if get(obj.gui.select_maskCB,'Value')
        cols = find(sb.select_mask);
      else
        cols = 1:size(sb.data,2);
      end
      
      smooth_weight = -1;
      smooth_var = -1;
      
      if ~exist('slices','var') || isempty(slices)
        slices = sb.slice+slice_range;
      end
      [~,slices_idxs] = intersect(slices,1:size(sb.data,3));
      slices = slices(sort(slices_idxs));
      if numel(slices)==1
        fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, active_idx, sb.slice);
      else
        fprintf('Apply %s to layer %d slices %d - %d\n', obj.tool_name, active_idx, slices(1), slices(end));
      end
      if get(obj.gui.previousCB,'Value')
        start_slice_idx = 2;
      else
        start_slice_idx = 1;
      end
      
      cmd = [];
      for slice_idx = start_slice_idx:length(slices)
        slice = slices(slice_idx);
        slice_prev = slices(slice_idx-1);
        % Create ground truth input
        % 1. Each column is one ground truth input
        % 2. Row 1: x, Row 2: y
        if numel(slices)>1
          fprintf('Slice %d\n',slice);
        end
        if get(obj.gui.previousCB,'Value')
          gt = [sb.layer(active_idx).x(cols(1:end),slice_prev).'-1; ...
            sb.layer(active_idx).y(cols(1:end),slice_prev).'+0.5];
        else
          gt = [];
        end
        if ~isempty(control_idx)
          mask = isfinite(sb.layer(control_idx).x(:,slice)) ...
            & isfinite(sb.layer(control_idx).y(:,slice));
          gt = cat(2,gt,[sb.layer(control_idx).x(mask,slice).'-1; ...
            sb.layer(control_idx).y(mask,slice).'+0.5]);
          [~,unique_idxs] = unique(gt(1,:),'last');
          gt = gt(:,unique_idxs);
          [~,sort_idxs] = sort(gt(1,:));
          gt = gt(:,sort_idxs);
          bottom_bin = sb.layer(control_idx).y(ceil(size(sb.data,2)/2)+1,slice);
        else
          bottom_bin = NaN;
        end
        
        if isempty(surf_idx)
          surf_bins = NaN*sb.layer(active_idx).y(:,slice);
        else
          surf_bins = sb.layer(surf_idx).y(:,slice);
        end
        surf_bins(isnan(surf_bins)) = -1;
        
        bottom_bin(isnan(bottom_bin)) = -1;
        
        if isempty(mask_idx)
          mask = ones(size(sb.data,2),1);
        else
          mask = sb.layer(mask_idx).y(:,slice);
        end
        
        detect_data = sb.data(:,:,slice);
        detect_data(detect_data>threshold) = threshold;
        detect_data = fir_dec(detect_data.',hanning(3).'/3,1).';
        % detect_data(182+(-5:5),35).'
%         obj.custom_data.mu = [8.4745    8.3321    9.7678   11.7998   13.1260   13.0728   11.6279   10.1136    9.2768    8.3387    7.3149];
        

        labels = tomo.detect(double(detect_data), ...
          double(surf_bins), double(bottom_bin), ...
          double(gt), double(mask), ...
          double(obj.custom_data.mu), double(obj.custom_data.sigma),-1,double(egt_weight), ...
          double(smooth_weight), double(smooth_var), double(slope));

        % Create cmd for layer change
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.layer = active_idx;
        cmd{end}.redo.layer = active_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.layer(active_idx).y(cols,slice);
        cmd{end}.redo.x = cols;
        cmd{end}.redo.y = labels(cols);
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
        set(obj.h_fig,'Name',sprintf('%d: detect tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: detect tool prefs',obj.h_fig.Number));
      end
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      pos = get(obj.h_fig,'Position');
      pos(3) = 200;
      pos(4) = 140;
      set(obj.h_fig,'Position',pos);
      
      % Slice range
      obj.gui.slice_rangeTXT = uicontrol('Style','text','string','Slice range');
      set(obj.gui.slice_rangeTXT,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-1:10" or "1:-1:-10".');
      
      obj.gui.slice_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.slice_rangeLE,'style','edit')
      set(obj.gui.slice_rangeLE,'string','-1:0')
      set(obj.gui.slice_rangeLE,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-1:10" or "1:-1:-10".');
      
      % Threshold
      obj.gui.thresholdTXT = uicontrol('Style','text','string','Threshold');
      set(obj.gui.thresholdTXT,'TooltipString','Specify an image threshold.');
      
      obj.gui.thresholdLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.thresholdLE,'style','edit')
      set(obj.gui.thresholdLE,'string','13.5')
      set(obj.gui.thresholdLE,'TooltipString','Specify an image threshold.');
       
      % Ground Truth Weight
      obj.gui.gt_weightTXT = uicontrol('Style','text','string','GT weight');
      set(obj.gui.gt_weightTXT,'TooltipString','Specify weighting of ground truth.');
      
      obj.gui.gt_weightLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.gt_weightLE,'style','edit')
      set(obj.gui.gt_weightLE,'string','10')
      set(obj.gui.gt_weightLE,'TooltipString','Specify weighting of ground truth.');
      
      obj.gui.slopeTXT = uicontrol('Style','text','string','Slope');
      obj.gui.slopeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.slopeLE,'style','edit')
      set(obj.gui.slopeLE,'string','zeros(1,63)')
      
      % Select mask
      obj.gui.select_maskCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.select_maskCB,'style','checkbox')
      set(obj.gui.select_maskCB,'string','Select')
      set(obj.gui.select_maskCB,'value',1)
      set(obj.gui.select_maskCB,'TooltipString','Check to operate only on the selected region.');
      
      % Previous ground truth
      obj.gui.previousCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.previousCB,'style','checkbox')
      set(obj.gui.previousCB,'string','Previous')
      set(obj.gui.previousCB,'value',1)
      set(obj.gui.previousCB,'TooltipString','Use previous slice as ground truth.');
      
      % GUI container table
      obj.gui.table.ui = obj.h_fig;
      obj.gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.table.height_margin = NaN*zeros(30,30);
      obj.gui.table.false_width = NaN*zeros(30,30);
      obj.gui.table.false_height = NaN*zeros(30,30);
      obj.gui.table.offset = [0 0];
      row = 1;
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
      obj.gui.table.handles{row,col}   = obj.gui.gt_weightTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.gt_weightLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.slopeTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.slopeLE;
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
      obj.gui.table.handles{row,col}   = obj.gui.previousCB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.table);     
    end

  end
  
end


