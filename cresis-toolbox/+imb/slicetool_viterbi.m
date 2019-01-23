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
          
          slices = slice-slice_range:slice+slice_range;
          for idx = 1:length(slices)
            slice = slices(idx);
            mask = isfinite(sb.sd.surf(control_idx).x(:,slice)) ...
              & isfinite(sb.sd.surf(control_idx).y(:,slice));
            mask(1:sb.bounds_relative(1)) = 0;
            mask(end-sb.bounds_relative(2)+1:end) = 0;
          end
        end
        if isempty(mask_idx)
          mask = ones(size(sb.data,2),length(slices));
        else
          mask = sb.sd.surf(mask_idx).y(:,slices);
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
        % Obtained from geostatistical analysis of 2014 Greenland P3
        prob.DIM_means = [    16.0582   27.3381   34.0492   40.5714...
          46.9463   52.3826   58.4664   63.7750   70.5143   76.3140...
          81.3519   86.9523   92.2285   98.1430   102.8310  107.0558...
          112.4538  116.3923  121.0109  125.1486  128.7611  133.3286...
          136.3500  139.5058  142.3812  146.1565  148.3781  151.1398...
          153.5642  155.2606  157.4824  159.8529  161.0239  163.1799...
          164.2849  166.1013  166.0728  167.3132  168.1448  169.1323...
          169.7559  170.3869  171.4264  171.8926  171.5201  171.8870...
          171.6407  172.3505  171.3533  171.6161];
        
        prob.DIM_vars  = [     20.5619   24.9082   28.0037   32.0840...
          35.6021   38.9544    42.4034   45.3588   48.9714   52.1360...
          55.0826   57.5144    60.3847   63.1485   65.1199   67.3734...
          69.3662   71.2849    72.9471   74.3759   75.5521   76.9737...
          77.9961   79.3596    79.9999   81.0342   81.6340   82.2424...
          82.9658   83.4794    84.1632   84.4168   85.0014   85.3065...
          85.7757   86.1880    86.3563   86.5577   86.7289   87.0748...
          87.1360   87.2473    87.2828   86.6350   86.5453   86.2989...
          86.4736   86.7318   87.1606   87.6966];
        
        costm = ones(1201, 101);
        fd = 18 * fir_dec(1:101, ones(1,5)/3.7.');
        
        for t = 1:1200
          for i = 1:101
            costm(t,i) = 200 * exp(-0.075 .* t) - 200 * exp(-0.075 * 50);
            if t > fd(i)
              costm(t,i) = 200;
            end
          end
        end
        
        prob.DIM_means = [prob.DIM_means prob.DIM_means(end)* ones(1,50)];
        prob.DIM_vars  = [prob.DIM_vars  prob.DIM_vars(end) * ones(1,50)];
        
        DIM_costmatrix = ones(1200, 100);
        for DIM = 1 : 100
          for T = 1 : 1200
            var = DIM * 0.05 * prob.DIM_vars(end);
            cost = 0.1 ./ normpdf(T, prob.DIM_means(DIM), var);
            cost(cost > 800) = 800;
            DIM_costmatrix(T, DIM) = cost;
          end
        end
        
        for k = 1:100
          DIM_costmatrix(1:50, k) = DIM_costmatrix(1:50, k) + 0.06 * k * costm(1:50, k);
        end
        
        DIM_costmatrix(DIM_costmatrix < 0) = 0; 
        DIM_costmatrix(DIM_costmatrix > 800) = 800;
        
        DIM_costmatrix = DIM_costmatrix(3:end, :);
        DIM_costmatrix(end+1, :) = DIM_costmatrix(end,:);
        DIM_costmatrix(end+1, :) = DIM_costmatrix(end,:);
        DIM_costmatrix = DIM_costmatrix ./ 4;
        
        % Visualization of DIM_costmatrix
        if 0
          figure; (imagesc(DIM_costmatrix)); colorbar; hold on;
          xlabel('Distance to nearest ice-margin [m]');
          ylabel('Ice thickness distribution [range-bin]')
          title('Added cost');
        end

        %% DoA-to-DoA transition model
        % Obtained from geostatistical analysis of 2014 Greenland P3
        transition_mu = [0.000000, 0.000000, 2.590611, 3.544282, 4.569263, 5.536577, 6.476430, 7.416807, 8.404554, 9.457255, 10.442658, 11.413710, 12.354409, 13.332689, 14.364614, 15.381671, 16.428969, 17.398906, 18.418794, 19.402757, 20.383026, 21.391834, 22.399259, 23.359765, 24.369957, 25.344982, 26.301805, 27.307530, 28.274756, 28.947572, 29.691010, 32.977387, 34.203212, 34.897994, 35.667128, 36.579019, 37.558978, 38.548659, 39.540715, 40.550138, 41.534781, 42.547407, 43.552700, 44.537758, 45.553618, 46.561057, 47.547331, 48.530976, 49.516588, 50.536075, 51.562886, 52.574938, 53.552979, 54.554206, 55.559657, 56.574029, 57.591999, 58.552986, 59.562937, 60.551616, 61.549909, 62.551092, 63.045791, 63.540490];
        transition_sigma = [0.457749, 0.805132, 1.152514, 1.213803, 1.290648, 1.370986, 1.586141, 1.626730, 1.785789, 1.791043, 1.782936, 1.727153, 1.770210, 1.714973, 1.687484, 1.663294, 1.633185, 1.647318, 1.619522, 1.626555, 1.649593, 1.628138, 1.699512, 1.749184, 1.809822, 1.946782, 2.126822, 2.237959, 2.313358, 2.280555, 1.419753, 1.112363, 1.426246, 2.159619, 2.140899, 2.083267, 1.687420, 1.574745, 1.480296, 1.443887, 1.415708, 1.356100, 1.401891, 1.398477, 1.365730, 1.418647, 1.407810, 1.430151, 1.391357, 1.403471, 1.454194, 1.470535, 1.417235, 1.455086, 1.436509, 1.378037, 1.415834, 1.333177, 1.298108, 1.277559, 1.358260, 1.483521, 1.674642, 1.865764];
        
        % Visualization of mean and variance vectors
        if 0
          figure; (plot(transition_mu)); hold on;
          plot(transition_sigma); xlim([1 64])
          legend('Mean', 'Variance', 'Location', 'northwest');
          xlabel('DoA bins');
        end
        
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