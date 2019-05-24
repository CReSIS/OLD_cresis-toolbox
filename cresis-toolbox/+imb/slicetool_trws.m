classdef (HandleCompatible = true) slicetool_trws < imb.slicetool
  % Slice_browser tool which calls detect.cpp (HMM)
  
  properties (SetAccess = protected, GetAccess = public)
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_trws()
      obj.create_option_ui();
      obj.tool_name = 'trws';
      obj.tool_menu_name = 'TR(W)S';
      obj.tool_shortcut = 'w';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 'w: tracking tool which runs TRWS solution to HMM inference model to find best surface. Neighboring slices effect cost function to improve solution.';
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
        eval_cmd = get(obj.gui.slice_rangeLE,'String');
        slice_range = eval(eval_cmd);
      catch ME
        error('%s\nError in slice range "%s":\n %s\n', repmat('=',[1 80]), ...
          eval_cmd, ME.message);
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
        range = eval(get(obj.gui.rangeLE,'String'));
      catch ME
        error('Error in range: %s', ME.getReport);
      end
      try
        polynomial = eval(get(obj.gui.polynomialLE,'String'));
      catch ME
        error('Error in polynomial: %s', ME.getReport);
      end
      try
        mu_size = eval(get(obj.gui.correlationLE,'String'));
      catch ME
        error('Error in correlation: %s', ME.getReport);
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
      top_edge_en = get(obj.gui.topEdgeCB,'Value');
      bottom_edge_en = get(obj.gui.bottomEdgeCB,'Value');
      
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
      fprintf('Apply %s to surface %d slices %d - %d\n', obj.tool_name, active_idx, slices(1), slices(end));
      
      gt = [];
      if ~isempty(control_idx)
        % Create ground truth input
        % 1. Each column is one ground truth input
        % 2. Row 1: relative slice/range-line, Row 2: x, Row 3: y
        for idx = 1:length(slices)
          slice = slices(idx);
          mask = isfinite(sb.sd.surf(control_idx).x(:,slice)) ...
            & isfinite(sb.sd.surf(control_idx).y(:,slice));
          mask(1:sb.bounds_relative(1)) = 0;
          mask(end-sb.bounds_relative(2)+1:end) = 0;
          gt = cat(2,gt,[(idx-1)*ones(1,sum(mask)); ...
            sb.sd.surf(control_idx).x(mask,slice).'-1; ...
            sb.sd.surf(control_idx).y(mask,slice).'+0.5]);
          bottom_bin = sb.sd.surf(control_idx).y(33,slices);
        end
      else
        bottom_bin = NaN*zeros(1,length(slices));
      end
      
      if isempty(surf_idx)
        surf_bins = NaN*sb.sd.surf(active_idx).y(:,slices);
      else
        surf_bins = sb.sd.surf(surf_idx).y(:,slices);
      end
      surf_bins(isnan(surf_bins)) = -1;
      
      bottom_bin(isnan(bottom_bin)) = -1;
      
      if isempty(mask_idx)
        mask = ones(size(sb.data,2),length(slices));
      else
        mask = sb.sd.surf(mask_idx).y(:,slices);
      end
      
      begin_slice = max(1, min(slices)-1);
      end_slice = min(size(sb.data,3), max(slices)+1);
      edge = [sb.sd.surf(active_idx).y(:,begin_slice), sb.sd.surf(active_idx).y(:,end_slice)];
      edge(mask(:,1) == 0,1) = surf_bins(mask(:,1) == 0,1);
      edge(mask(:,end) == 0,2) = surf_bins(mask(:,end) == 0,end);
      
      trws_data = sb.data(:,:,slices);
      trws_data(trws_data>threshold(2)) = threshold(2);
      for idx = 1:length(slices)
        for col = 1:size(sb.data,2)
          tmp = trws_data(1:min(end,round(surf_bins(col,idx))+70),col,idx);
          tmp(tmp>threshold(1)) = threshold(1);
          trws_data(1:min(end,round(surf_bins(col,idx))+70),col,idx) = tmp;
        end
      end
      if ~left_edge_en
        edge(:,1) = -1;
      end
      if ~right_edge_en
        edge(:,end) = -1;
      end
      if ~isempty(polynomial)
        smooth_slope = polyval(polynomial, linspace(-1,1,size(sb.data,2)-1));
      else
        smooth_slope = [];
      end
      
      smooth_weight = smooth(1:2);
      smooth_var = smooth(3);
      mu = sinc(linspace(-1.5,1.5,mu_size));
      sigma = sum(mu)/20*ones(1,mu_size);
      mask_dist      = round(bwdist(mask == 0));
      bounds = [sb.bounds_relative(1) size(trws_data,2)-sb.bounds_relative(2)-1 -1 -1];  
      mask_dist = round(mask_dist .* 9);
      
      clear DIM DIM_costmatrix;
      global gRadar;
      DIM = load(fullfile(gRadar.path, '+tomo', 'Layer_tracking_3D_parameters_Matrix.mat'));
      DIM_costmatrix = DIM.Layer_tracking_3D_parameters;
      DIM_costmatrix = DIM_costmatrix .* (200 ./ max(DIM_costmatrix(:)));

      
      %% DoA-to-DoA transition model
      % Obtained from geostatistical analysis of 2014 Greenland P3
      transition_mu = [0.000000, 0.000000, 2.590611, 3.544282, 4.569263, 5.536577, 6.476430, 7.416807, 8.404554, 9.457255, 10.442658, 11.413710, 12.354409, 13.332689, 14.364614, 15.381671, 16.428969, 17.398906, 18.418794, 19.402757, 20.383026, 21.391834, 22.399259, 23.359765, 24.369957, 25.344982, 26.301805, 27.307530, 28.274756, 28.947572, 29.691010, 32.977387, 34.203212, 34.897994, 35.667128, 36.579019, 37.558978, 38.548659, 39.540715, 40.550138, 41.534781, 42.547407, 43.552700, 44.537758, 45.553618, 46.561057, 47.547331, 48.530976, 49.516588, 50.536075, 51.562886, 52.574938, 53.552979, 54.554206, 55.559657, 56.574029, 57.591999, 58.552986, 59.562937, 60.551616, 61.549909, 62.551092, 63.045791, 63.540490];
      transition_sigma = [0.457749, 0.805132, 1.152514, 1.213803, 1.290648, 1.370986, 1.586141, 1.626730, 1.785789, 1.791043, 1.782936, 1.727153, 1.770210, 1.714973, 1.687484, 1.663294, 1.633185, 1.647318, 1.619522, 1.626555, 1.649593, 1.628138, 1.699512, 1.749184, 1.809822, 1.946782, 2.126822, 2.237959, 2.313358, 2.280555, 1.419753, 1.112363, 1.426246, 2.159619, 2.140899, 2.083267, 1.687420, 1.574745, 1.480296, 1.443887, 1.415708, 1.356100, 1.401891, 1.398477, 1.365730, 1.418647, 1.407810, 1.430151, 1.391357, 1.403471, 1.454194, 1.470535, 1.417235, 1.455086, 1.436509, 1.378037, 1.415834, 1.333177, 1.298108, 1.277559, 1.358260, 1.483521, 1.674642, 1.865764];
      
      if length(transition_mu) ~= size(sb.data, 2)
        transition_mu = imresize(transition_mu, [1 size(sb.data, 2)]);
      end
      
      if length(transition_sigma) ~= size(sb.data, 2)
        transition_sigma = imresize(transition_sigma, [1 size(sb.data, 2)]);
      end
      
      if top_edge_en && ~isempty(cols) && cols(1) > bounds(1)+1
        % Add ground truth from top edge as long as top edge is bounded
        gt = cat(2,gt,[slices-slices(1); ...
          sb.sd.surf(active_idx).x(cols(1)-1,slices)-1; ...
          sb.sd.surf(active_idx).y(cols(1)-1,slices)+0.5]);
      end
      
      if bottom_edge_en && ~isempty(cols) && cols(end) < bounds(2)
        % Add ground truth from top edge as long as top edge is bounded
        gt = cat(2,gt,[slices-slices(1); ...
          sb.sd.surf(active_idx).x(cols(end)+1,slices)-1; ...
          sb.sd.surf(active_idx).y(cols(end)+1,slices)+0.5]);
      end
      
      rows = [];
      if isfinite(range) && size(gt,2) > 0
        % Restrict search to range of rows around ground truth
        rows = max(1,min(round(gt(3,:)-range))) : min(size(trws_data,1),max(round(gt(3,:)+range)));
        if length(rows) < length(mu)+1
          error('Error: Range restriction leaves too few rows of data. Increase "Row range" option.');
        end
        trws_data = trws_data(rows,:,:);
        surf_bins = surf_bins - rows(1) + 1;
        bottom_bin = bottom_bin - rows(1) + 1;
        gt(3,:) = gt(3,:) - rows(1) + 1;
      end

      tic;
      correct_surface = tomo.trws(double(trws_data), ...
        double(surf_bins), double(bottom_bin), double(gt), double(mask), ...
        double(mu), double(sigma), double(smooth_weight), double(smooth_var), ...
        double(smooth_slope), double(edge), double(num_loops), int64(bounds), ...
        double(mask_dist), double(DIM_costmatrix), ...
        double(transition_mu), double(transition_sigma));
      toc;
      
      fprintf('  %.2f sec per slice\n', toc/size(trws_data,3));
      
      if ~isempty(rows)
        correct_surface = correct_surface + rows(1) - 1;
      end
      
      % Create cmd for surface change
      cmd = [];
      for idx = start_slice_idx:end_slice_idx
        slice = slices(idx);
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.surf = active_idx;
        cmd{end}.redo.surf = active_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.sd.surf(active_idx).y(cols,slice);
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
        set(obj.h_fig,'Name',sprintf('%d: trws tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: trws tool prefs',obj.h_fig.Number));
      end
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      pos = get(obj.h_fig,'Position');
      pos(3) = 220;
      pos(4) = 240;
      set(obj.h_fig,'Position',pos);
      
      % Number of loops
      obj.gui.numloopsTXT = uicontrol('Style','text','string','Loops');
      
      obj.gui.numloopsLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.numloopsLE,'style','edit')
      set(obj.gui.numloopsLE,'string','10')
      set(obj.gui.numloopsLE,'TooltipString','Number of iterations.');
      
      % Slice range
      obj.gui.slice_rangeTXT = uicontrol('Style','text','string','Slice range');
      set(obj.gui.slice_rangeTXT,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:7".');
      
      obj.gui.slice_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.slice_rangeLE,'style','edit')
      set(obj.gui.slice_rangeLE,'string','-5:6')
      set(obj.gui.slice_rangeLE,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:7".');
      
      % Threshold
      obj.gui.thresholdTXT = uicontrol('Style','text','string','Threshold');
      set(obj.gui.thresholdTXT,'TooltipString','Specify image threshold. First number is surface. Second number is whole column.');
      
      obj.gui.thresholdLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.thresholdLE,'style','edit')
      set(obj.gui.thresholdLE,'string','[13.5 inf]')
      set(obj.gui.thresholdLE,'TooltipString','Specify image threshold. First number is surface. Second number is whole column.');
      
      % Row range
      obj.gui.rangeTXT = uicontrol('Style','text','string','Row range');
      set(obj.gui.rangeTXT,'TooltipString','Specify a number of rows beyond the ground truth to search.');
      
      obj.gui.rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.rangeLE,'style','edit')
      set(obj.gui.rangeLE,'string','inf')
      set(obj.gui.rangeLE,'TooltipString','Specify a number of rows beyond the ground truth to search.');
      
      % Polynomial Shape
      obj.gui.polynomialTXT = uicontrol('Style','text','string','Shape poly');
      set(obj.gui.polynomialTXT,'TooltipString','Specify a polynomial coefficients (e.g. -5x^2 is [-5 0 0]).');
      
      obj.gui.polynomialLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.polynomialLE,'style','edit')
      set(obj.gui.polynomialLE,'string','[]')
      set(obj.gui.polynomialLE,'TooltipString','Specify polynomial coefficients (e.g. -5x^2 is [-5 0 0]).');
      
      % Correlation length
      obj.gui.correlationTXT = uicontrol('Style','text','string','Correlation');
      set(obj.gui.correlationTXT,'TooltipString','Specify a correlation length for surface impulse template.');
      
      obj.gui.correlationLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.correlationLE,'style','edit')
      set(obj.gui.correlationLE,'string','11')
      set(obj.gui.correlationLE,'TooltipString','Specify a correlation length for surface impulse template.');
      
      % Smooth
      obj.gui.smoothTXT = uicontrol('Style','text','string','Smoothness');
      set(obj.gui.smoothTXT,'TooltipString','Specify surface smoothness ["slice weight" "column weight" "edges weight"].');
      
      obj.gui.smoothLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.smoothLE,'style','edit')
      set(obj.gui.smoothLE,'string','[22 22 32]')
      set(obj.gui.smoothLE,'TooltipString','Specify surface smoothness ["slice weight" "column weight" "edges weight"].');
      
      % Select mask
      obj.gui.select_maskCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.select_maskCB,'style','checkbox')
      set(obj.gui.select_maskCB,'string','Select')
      set(obj.gui.select_maskCB,'value',1)
      set(obj.gui.select_maskCB,'TooltipString','Check to operate only on the selected region.');
      
      % Left edge
      obj.gui.leftEdgeCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.leftEdgeCB,'style','checkbox')
      set(obj.gui.leftEdgeCB,'string','L')
      set(obj.gui.leftEdgeCB,'value',1)
      set(obj.gui.leftEdgeCB,'TooltipString','Check to force a left edge boundary condition.');
      
      % Right edge
      obj.gui.rightEdgeCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.rightEdgeCB,'style','checkbox')
      set(obj.gui.rightEdgeCB,'string','R')
      set(obj.gui.rightEdgeCB,'value',1)
      set(obj.gui.rightEdgeCB,'TooltipString','Check to force a right edge boundary condition.');
      
      % Top edge
      obj.gui.topEdgeCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.topEdgeCB,'style','checkbox')
      set(obj.gui.topEdgeCB,'string','T')
      set(obj.gui.topEdgeCB,'value',1)
      set(obj.gui.topEdgeCB,'TooltipString','Check to force a top edge boundary condition.');
      
      % Bottom edge
      obj.gui.bottomEdgeCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.bottomEdgeCB,'style','checkbox')
      set(obj.gui.bottomEdgeCB,'string','B')
      set(obj.gui.bottomEdgeCB,'value',1)
      set(obj.gui.bottomEdgeCB,'TooltipString','Check to force a bottom edge boundary condition.');
      
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
      obj.gui.table.handles{row,col}   = obj.gui.rangeTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.rangeLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.polynomialTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.polynomialLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.correlationTXT;
      obj.gui.table.width(row,col)     = 70;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.correlationLE;
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
      col = col+1;
      obj.gui.table.handles{row,col}   = obj.gui.leftEdgeCB;
      obj.gui.table.width(row,col)     = 25;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = col+1;
      obj.gui.table.handles{row,col}   = obj.gui.rightEdgeCB;
      obj.gui.table.width(row,col)     = 25;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = col+1;
      obj.gui.table.handles{row,col}   = obj.gui.topEdgeCB;
      obj.gui.table.width(row,col)     = 25;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = col+1;
      obj.gui.table.handles{row,col}   = obj.gui.bottomEdgeCB;
      obj.gui.table.width(row,col)     = 25;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.table);
      
    end
    
  end
  
end


