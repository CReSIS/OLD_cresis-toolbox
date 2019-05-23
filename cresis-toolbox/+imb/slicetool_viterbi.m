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
        fprintf('Apply %s to surface %d slice %d\n', obj.tool_name, ...
          active_idx, sb.slice);
      else
        fprintf('Apply %s to surface %d slices %d - %d\n', ...
          obj.tool_name, active_idx, slices(1), slices(end));
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
        smooth_weight = 1;
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
        
        %% Distance to ice-mask calculation
        mask_dist = round(bwdist(mask == 0));
        mask_dist = mask_dist(:, 4);
        mask      = mask(:, slice_range+1);
        
        %% Distance-to-Ice-Margin model
        clear DIM DIM_costmatrix;
        global gRadar
        DIM = load(fullfile(gRadar.path, '+tomo', 'Layer_tracking_3D_parameters_Matrix.mat'));
        DIM_costmatrix = DIM.Layer_tracking_3D_parameters;
        DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));
        
        %% DoA-to-DoA transition model
        % Obtained from geostatistical analysis of 2014 Greenland P3
        transition_mu = [0.000000, 0.000000, 2.590611, 3.544282, 4.569263, 5.536577, 6.476430, 7.416807, 8.404554, 9.457255, 10.442658, 11.413710, 12.354409, 13.332689, 14.364614, 15.381671, 16.428969, 17.398906, 18.418794, 19.402757, 20.383026, 21.391834, 22.399259, 23.359765, 24.369957, 25.344982, 26.301805, 27.307530, 28.274756, 28.947572, 29.691010, 32.977387, 34.203212, 34.897994, 35.667128, 36.579019, 37.558978, 38.548659, 39.540715, 40.550138, 41.534781, 42.547407, 43.552700, 44.537758, 45.553618, 46.561057, 47.547331, 48.530976, 49.516588, 50.536075, 51.562886, 52.574938, 53.552979, 54.554206, 55.559657, 56.574029, 57.591999, 58.552986, 59.562937, 60.551616, 61.549909, 62.551092, 63.045791, 63.540490];
        transition_sigma = [0.457749, 0.805132, 1.152514, 1.213803, 1.290648, 1.370986, 1.586141, 1.626730, 1.785789, 1.791043, 1.782936, 1.727153, 1.770210, 1.714973, 1.687484, 1.663294, 1.633185, 1.647318, 1.619522, 1.626555, 1.649593, 1.628138, 1.699512, 1.749184, 1.809822, 1.946782, 2.126822, 2.237959, 2.313358, 2.280555, 1.419753, 1.112363, 1.426246, 2.159619, 2.140899, 2.083267, 1.687420, 1.574745, 1.480296, 1.443887, 1.415708, 1.356100, 1.401891, 1.398477, 1.365730, 1.418647, 1.407810, 1.430151, 1.391357, 1.403471, 1.454194, 1.470535, 1.417235, 1.455086, 1.436509, 1.378037, 1.415834, 1.333177, 1.298108, 1.277559, 1.358260, 1.483521, 1.674642, 1.865764];
        
        %% Call viterbi.cpp
        tic
        labels = tomo.viterbi(double(viterbi_data), double(surf_bins), ...
          double(bottom_bin), double(gt), double(mask), double(mu), ...
          double(sigma), double(egt_weight), double(smooth_weight), ...
          double(smooth_var), double(slope), int64(bounds), ...
          double(viterbi_weight), double(repulsion), double(ice_bin_thr), ...
          double(mask_dist), double(DIM_costmatrix), ...
          double(transition_mu), double(transition_sigma));
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