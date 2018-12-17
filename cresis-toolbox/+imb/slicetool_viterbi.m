classdef (HandleCompatible = true) slicetool_viterbi < imb.slicetool
  % Slice_browser tool which calls viterbi.cpp (HMM)
  %
  % Compile C++ program with:
  %   mex -largeArrayDims viterbi.cpp
  
  properties (SetAccess = protected, GetAccess = public)
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_viterbi()
      obj.create_option_ui();
      obj.tool_name = 'viterbi';
      obj.tool_menu_name = '(V)iterbi';
      obj.tool_shortcut = 'v';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 'v: Viterbi tool which runs Viterbi solution to HMM inference model to find best surface. Neighboring slices have no influence on solution.';
    end
    
    function cmd = apply_PB_callback(obj,sb,slices)
      % sb: slice browser object. Use the following fields to create
      %     commands, cmd, that use sb.data to operate on sb.sd. You
      %     should not modify any fields of sb.
      %  .sd: surfdata .surf struct array containing surface information
      %  .data: 3D image
      %  .slice: current slice in 3D image (third index of .data)
      %  .surf_idx: active surface
      % slices: array of slices to operate on (overrides sb.slice)
      control_idx = sb.sd.surf(sb.surf_idx).gt;
      active_idx = sb.sd.surf(sb.surf_idx).active;
      surf_idx = sb.sd.surf(sb.surf_idx).top;
      mask_idx = sb.sd.surf(sb.surf_idx).mask;
      
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
      
      if ~exist('slices','var') || isempty(slices)
        slices = sb.slice+slice_range;
      end
      [~,slices_idxs] = intersect(slices,1:size(sb.data,3));
      slices = slices(sort(slices_idxs));
      if numel(slices)==1
        fprintf('Apply %s to surface %d slice %d\n', obj.tool_name, active_idx, sb.slice);
      else
        fprintf('Apply %s to surface %d slices %d - %d\n', obj.tool_name, active_idx, slices(1), slices(end));
      end
      if get(obj.gui.previousCB,'Value')
        start_slice_idx = 2;
      else
        start_slice_idx = 1;
      end
      
      cmd = [];
      for slice_idx = start_slice_idx:length(slices)
        slice = slices(slice_idx);
        % Create ground truth input
        % 1. Each column is one ground truth input
        % 2. Row 1: x, Row 2: y
        if numel(slices)>1
          fprintf('Slice %d\n',slice);
        end
        if get(obj.gui.previousCB,'Value')
          slice_prev = slices(slice_idx-1);
          if slice_idx == 2
            gt = [sb.sd.surf(active_idx).x(:,slice_prev).'-1; ...
              sb.sd.surf(active_idx).y(:,slice_prev).'+0.5];
          else
            gt = [sb.sd.surf(active_idx).x(:,slice_prev).'-1; ...
              labels(:).'+0.5];
          end
        else
          gt = [];
        end
        if ~isempty(control_idx)
          mask = isfinite(sb.sd.surf(control_idx).x(:,slice)) ...
            & isfinite(sb.sd.surf(control_idx).y(:,slice));
          gt = cat(2,gt,[sb.sd.surf(control_idx).x(mask,slice).'-1; ...
            sb.sd.surf(control_idx).y(mask,slice).'+0.5]);
          [~,unique_idxs] = unique(gt(1,:),'last','legacy');
          gt = gt(:,unique_idxs);
          viterbi_weight = ones([1 (size(sb.data,2))]);
          viterbi_weight(1 + gt(1,:)) = 2;
          [~,sort_idxs] = sort(gt(1,:));
          gt = gt(:,sort_idxs);
          bottom_bin = sb.sd.surf(control_idx).y(ceil(size(sb.data,2)/2)+1,slice);
        else
          bottom_bin = NaN;
          viterbi_weight = ones([1 length(gt)]);
        end
        
        
        %%
        if ~isempty(control_idx)
          % Create ground truth input
          % 1. Each column is one ground truth input
          % 2. Row 1: relative slice/range-line, Row 2: x, Row 3: y
          
          slice_range = 3;
          
          m_slices = slice-slice_range : slice+slice_range;
          for m_slice_idx = 1:length(m_slices)
            mask = isfinite(sb.sd.surf(control_idx).x(:,m_slices(m_slice_idx))) ...
              & isfinite(sb.sd.surf(control_idx).y(:,m_slices(m_slice_idx)));
            mask(1:sb.bounds_relative(1)) = 0;
            mask(end-sb.bounds_relative(2)+1:end) = 0;
          end
        end
        if isempty(mask_idx)
          mask = ones(size(sb.data,2),length(m_slices));
        else
          mask = sb.sd.surf(mask_idx).y(:,m_slices);
        end
        
        
        %%
        if isempty(surf_idx)
          surf_bins = NaN*sb.sd.surf(active_idx).y(:,slice);
        else
          surf_bins = sb.sd.surf(surf_idx).y(:,slice);
        end
        surf_bins(isnan(surf_bins)) = -1;
        bottom_bin(isnan(bottom_bin)) = -1;
        
        viterbi_data = sb.data(:,:,slice);
        viterbi_data(viterbi_data>threshold) = threshold;
        viterbi_data = fir_dec(viterbi_data.',hanning(3).'/3,1).';
        
        bounds = [1 (size(sb.data,2))];
        
        mu_size       = 11;
        mu            = sinc(linspace(-1.5, 1.5, mu_size));
        sigma         = sum(mu)/20*ones(1,mu_size);
        smooth_var    = -1;
        smooth_weight = 10;
        repulsion     = 100;
        ice_bin_thr   = 3;
        mc            = -1 * ones(1, size(sb.data,2));
        mc_weight     = 0;
        
        %%%% TO COMPILE
        if 0
          tmp = pwd;
          cd ~/scripts/cresis-toolbox/cresis-toolbox/+tomo/
          mex -largeArrayDims viterbi.cpp
          cd(tmp);
        end
        %%%%
        
        mask           = 90*fir_dec(fir_dec(double(shrink(mask,2)),ones(1,5)/3.7).',ones(1,5)/3.7).';
        mask(mask>=90) = inf;
        mask           = mask(:, slice_range+1);
        
        CF.sensory_distance = 200;
        CF.max_cost = 50;
        CF.lambda = 0.075;
        
        % Call viterbi.cpp
        tic
        labels = tomo.viterbi(double(viterbi_data), double(surf_bins), ...
          double(bottom_bin), double(gt), double(mask), double(mu), ...
          double(sigma), double(egt_weight), double(smooth_weight), ...
          double(smooth_var), double(slope), int64(bounds), ...
          double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
          double(mc), double(mc_weight), ...
          double(CF.sensory_distance), double(CF.max_cost), double(CF.lambda));
        toc
        
        labels(surf_bins(:) > labels(:)) = surf_bins(surf_bins(:) > labels(:));
        
        % Create cmd for surface change
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.surf = active_idx;
        cmd{end}.redo.surf = active_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.sd.surf(active_idx).y(cols,slice);
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
        set(obj.h_fig,'Name',sprintf('%d: viterbi tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: viterbi tool prefs',obj.h_fig.Number));
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


