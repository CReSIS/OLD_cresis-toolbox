classdef (HandleCompatible = true) slicetool_threshold < imb.slicetool
  
  properties    
    img
  end
  
  properties (SetAccess = protected, GetAccess = protected)
  end
  
  events
  end
  
  methods
    function obj = slicetool_threshold()
      obj.create_option_ui();
      obj.tool_name = 'Threshold';
      obj.tool_menu_name = '(T)hreshold';
      obj.tool_shortcut = 't';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
      obj.help_string = 't: Finds first threshold (useful for ice surface).';
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
      active_idx = sb.sd.surf(sb.surf_idx).active;
      
      try
        slice_range = eval(get(obj.gui.slice_rangeLE,'String'));
      catch ME
        error('Error in slice range: %s', ME.getReport);
      end
      try
        Noise_Range_bin = eval(get(obj.gui.Noise_Range_binLE,'String'));
      catch ME
        error('Error in Noise_Range_bin: %s', ME.getReport);
      end 
      try
          row_range = eval(get(obj.gui.row_rangeLE,'String'));
      catch ME
          error('Error in row range: %s', ME.getReport);
      end
      try
          threshold = eval(get(obj.gui.thresholdLE,'String'));
      catch ME
          error('Error in threshold: %s', ME.getReport);
      end
      
      if get(obj.gui.select_maskCB,'Value')
        cols = find(sb.select_mask);
      else
        cols = 1:size(sb.data,2);
      end
      
      cmd = [];
      params.Noise_Range_bin = Noise_Range_bin; 
      params.row_range = row_range;
      params.threshold = threshold;
      
      if ~exist('slices','var') || isempty(slices)
        slice_range = min(slice_range):max(slice_range);
        slices = sb.slice+slice_range;
      end
      slices = intersect(slices,1:size(sb.data,3));
      if numel(slices)==1
        fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, active_idx, sb.slice);
      else
        fprintf('Apply %s to layer %d slices %d - %d\n', obj.tool_name, active_idx, slices(1), slices(end));
      end
      
      for slice = slices(:).'     
        img_slice = obj.img(:,:,slice);
        layer = tomo.threshold(img_slice,params);
        
        if 0
          data = sb.data(:,:,slice);
          %data(data>threshold) = threshold;
          data = fir_dec(data.',hanning(11).'/6,1).';
          
          param = [];
          param.noise_rng = [10 -175 -150];
          param.filter_len = 3;
          param.threshold = 2.5;
          param.threshold_rng = 2;
          param.sidelobe = 2;
          param.init.method = 'snake';
          param.init.search_rng = [-10:10];
          param.max_diff = inf;
          param.search_rng = [0:75];
          %param.noise_override = 4;
          
          layer = tracker_threshold(data,param);
        end
        
        % Create cmd for layer change
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.surf = active_idx;
        cmd{end}.redo.surf = active_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.sd.surf(active_idx).y(cols,slice);
        cmd{end}.redo.x = cols;
        cmd{end}.redo.y = layer(cols);
        cmd{end}.type = 'standard';
      end
      
      % Add dummy command at end to make view go to the current slice
      cmd{end+1}.redo.slice = sb.slice;
      cmd{end}.undo.slice = sb.slice;
      cmd{end}.type = 'slice_dummy';
    end
    
    function set_custom_data(obj,custom_data)
    end
    
    function create_option_ui(obj)
      obj.h_fig = figure('Visible','off','DockControls','off', ...
        'NumberTitle','off','ToolBar','none','MenuBar','none','Resize','off');
      if strcmpi(class(obj.h_fig),'double')
        set(obj.h_fig,'Name',sprintf('%d: threshold tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: threshold tool prefs',obj.h_fig.Number));
      end
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      pos = get(obj.h_fig,'Position');
      pos(3) = 200;
      pos(4) = 130; % Height of Threshold tool options window 
      set(obj.h_fig,'Position',pos);
      
      % Slice range
      obj.gui.slice_rangeTXT = uicontrol('Style','text','string','Slice range');
      set(obj.gui.slice_rangeTXT,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:5".');
      
      obj.gui.slice_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.slice_rangeLE,'style','edit')
      set(obj.gui.slice_rangeLE,'string','0')
      set(obj.gui.slice_rangeLE,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:5".');
      
      % Noise_Range_bin
      obj.gui.Noise_Range_binTXT = uicontrol('Style','text','string','Noise Range bin');
      set(obj.gui.Noise_Range_binTXT,'TooltipString','Specify Noise floor Range bin.');
      
      obj.gui.Noise_Range_binLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.Noise_Range_binLE,'style','edit')
      set(obj.gui.Noise_Range_binLE,'string','1')
      set(obj.gui.Noise_Range_binLE,'TooltipString','Specify Noise Range bin.');    
      
      % Row range
      obj.gui.row_rangeTXT = uicontrol('Style','text','string','Row range');
      set(obj.gui.row_rangeTXT,'TooltipString','Enter a vector specifying relative row range to search in. E.g. "-30:30".');
      
      obj.gui.row_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.row_rangeLE,'style','edit')
      set(obj.gui.row_rangeLE,'string','-30:30')
      set(obj.gui.row_rangeLE,'TooltipString','Enter a vector specifying relative row range to search in. E.g. "-30:30".');
      
      % Threshold
      obj.gui.thresholdTXT = uicontrol('Style','text','string','Threshold');
      set(obj.gui.thresholdTXT,'TooltipString','Specify an image threshold.');
      
      obj.gui.thresholdLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.thresholdLE,'style','edit')
      set(obj.gui.thresholdLE,'string','10')
      set(obj.gui.thresholdLE,'TooltipString','Specify an image threshold.');
           
      % Select mask
      obj.gui.select_maskCB = uicontrol('parent',obj.h_fig);
      set(obj.gui.select_maskCB,'style','checkbox')
      set(obj.gui.select_maskCB,'string','Select')
      set(obj.gui.select_maskCB,'value',1)
      set(obj.gui.select_maskCB,'TooltipString','Check to operate only on the selected region.');
      
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
      obj.gui.table.width(row,col)     = 90;
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
      obj.gui.table.handles{row,col}   = obj.gui.Noise_Range_binTXT;
      obj.gui.table.width(row,col)     = 90;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.Noise_Range_binLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;    
      obj.gui.table.handles{row,col}   = obj.gui.row_rangeTXT;
      obj.gui.table.width(row,col)     = 90;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.row_rangeLE;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
           
      row = row + 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.thresholdTXT;
      obj.gui.table.width(row,col)     = 90;
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
      obj.gui.table.handles{row,col}   = obj.gui.select_maskCB;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = 20;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
          
      clear row col
      table_draw(obj.gui.table);
    end
    
  end
  
end