classdef (HandleCompatible = true) slicetool_max < imb.slicetool
  % Slice_browser tool which find the max value in a neighborhood around
  % the current values
  
  properties
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
  end
  
  methods
    function obj = slicetool_max()
      obj.create_option_ui();
      obj.tool_name = 'Max';
      obj.tool_menu_name = '(M)ax';
      obj.tool_shortcut = 'm';
      obj.ctrl_pressed = 0;
      obj.shift_pressed = 0;
    end
    
    function cmd = apply_PB_callback(obj,sb)
      % sb: slice browser object. Use the following fields to create
      %     commands, cmd, that use sb.data to operate on sb.layer. You 
      %     should not modify any fields of sb.
      %  .layer: struct array containing layer information
      %  .data: 3D image
      %  .slice: current slice in 3D image (third index of .data)
      %  .layer_idx: active layer
      fprintf('Apply %s to layer %d slice %d\n', obj.tool_name, sb.layer_idx, sb.slice);
      
      control_idx = sb.layer(sb.layer_idx).control_layer;
      
      try
        row_range = eval(get(obj.gui.row_rangeLE,'String'));
      catch ME
        error('Error in slice range: %s', ME.getReport);
      end
      try
        slice_range = eval(get(obj.gui.slice_rangeLE,'String'));
      catch ME
        error('Error in slice range: %s', ME.getReport);
      end
      
      slices = sb.slice+slice_range;
      slices = intersect(slices,1:size(sb.data,3));

      cmd = [];
      for slice = slices(:).'
        cols = find(sb.select_mask);
        new_y = [];
        for col = cols(:).'
          if ~isempty(control_idx) && ~isnan(sb.layer(control_idx).y(col,slice))
            new_y(end+1) = sb.layer(control_idx).y(col,slice);
          else
            rows = round(sb.layer(sb.layer_idx).y(col,slice)) + row_range;
            rows = intersect(rows,1:size(sb.data,1));
            [max_val,max_row] = max(sb.data(rows,col,slice));
            new_y(end+1) = rows(max_row);
          end
        end
        
        % Create cmd for layer change
        cmd{end+1}.undo.slice = slice;
        cmd{end}.redo.slice = slice;
        cmd{end}.undo.layer = sb.layer_idx;
        cmd{end}.redo.layer = sb.layer_idx;
        cmd{end}.undo.x = cols;
        cmd{end}.undo.y = sb.layer(sb.layer_idx).y(cols,slice);
        cmd{end}.redo.x = cols;
        cmd{end}.redo.y = new_y;
        cmd{end}.type = 'standard';
      end
      
      % Add dummy command at end to make view go to the current slice
      cmd{end+1}.undo.slice = sb.slice;
      cmd{end}.redo.slice = sb.slice;
      cmd{end}.undo.layer = sb.layer_idx;
      cmd{end}.redo.layer = sb.layer_idx;
      cmd{end}.undo.x = [];
      cmd{end}.undo.y = [];
      cmd{end}.redo.x = [];
      cmd{end}.redo.y = [];
      cmd{end}.type = 'standard';
    end
    
    function set_custom_data(obj,custom_data)
    end
    
    function create_option_ui(obj)
      obj.h_fig = figure('Visible','off','DockControls','off', ...
        'NumberTitle','off','ToolBar','none','MenuBar','none','Resize','off');
      if strcmpi(class(obj.h_fig),'double')
        set(obj.h_fig,'Name',sprintf('%d: max tool prefs',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: max tool prefs',obj.h_fig.Number));
      end
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      pos = get(obj.h_fig,'Position');
      pos(3) = 200;
      pos(4) = 50;
      set(obj.h_fig,'Position',pos);
      
      % Row range
      obj.gui.row_rangeTXT = uicontrol('Style','text','string','Row range');
      set(obj.gui.row_rangeTXT,'TooltipString','Enter a vector specifying relative row range to search in. E.g. "-30:30".');
      
      obj.gui.row_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.row_rangeLE,'style','edit')
      set(obj.gui.row_rangeLE,'string','-30:30')
      set(obj.gui.row_rangeLE,'TooltipString','Enter a vector specifying relative row range to search in. E.g. "-30:30".');
      
      % Extent
      obj.gui.slice_rangeTXT = uicontrol('Style','text','string','Slice range');
      set(obj.gui.slice_rangeTXT,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:5".');
      
      obj.gui.slice_rangeLE = uicontrol('parent',obj.h_fig);
      set(obj.gui.slice_rangeLE,'style','edit')
      set(obj.gui.slice_rangeLE,'string','-5:5')
      set(obj.gui.slice_rangeLE,'TooltipString','Enter a vector specifying relative range in slices. E.g. "-5:5".');
      
      % GUI container table
      obj.gui.table.ui = obj.h_fig;
      obj.gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.table.height_margin = NaN*zeros(30,30);
      obj.gui.table.false_width = NaN*zeros(30,30);
      obj.gui.table.false_height = NaN*zeros(30,30);
      obj.gui.table.offset = [0 0];
      row = 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.row_rangeTXT;
      obj.gui.table.width(row,col)     = 70;
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
      
      clear row col
      table_draw(obj.gui.table);      
      
    end
    
  end
  
end


