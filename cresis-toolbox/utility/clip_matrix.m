% Class clip_matrix
%
% Class which creates a GUI for clipping out bad sections of a matrix. It
% can also be used to browse images.
%
% Examples:
% obj = clip_matrix(data); % where data is 2-D real matrix

classdef (HandleCompatible = true) clip_matrix < handle
  properties
    % data = 2D matrix of original data that is in image
    data
    
    % bad_val = value that represents bad data
    bad_val
    
    % Function handle hooks for customizing clip_matrix
    fh_button_up
    fh_key_press

    % tool_list: structure array of tools
    %  .event_key: string containing the matching event.Key,
    %  .shift_pressed: boolean matching shift_pressed,
    %  .ctrl_pressed: boolean matching ctrl_pressed,
    %  .help_str: string containing tool descriptor for F1/help,
    %  .fh_callback: function handle to manipulate data according to tool selection
    tool_list
    % active_tool_idx: index into tool_list of active tool (zero if no
    %   tool is active)
    active_tool_idx
    
    % Undo stack
    undo_stack
    
    % GUI handles
    h_fig
    h_axes
    h_img
    h_imgwin % NOT USED
    sf % selection figure GUI handles
    
    % zoom_mode: boolean
    zoom_mode
    
    % x,y: the location of the mouse at button_down
    x
    y

    % actions: state structure that tells what the current action is
    %   (used for drawing polygon)
    actions
  end
  
  methods
    function obj = clip_matrix(data,param)
      
      % Check inputs
      if ~exist('param','var')
        param = [];
      end
      if ~isfield(param,'h_img')
        param.h_img = [];
      end
      if ~isfield(param,'bad_val')
        param.bad_val = NaN;
      end
      if ~isfield(param,'fh_button_up')
        param.fh_button_up = [];
      end
      if ~isfield(param,'fh_key_press')
        param.fh_key_press = [];
      end
      if ~isfield(param,'tool_list')
        obj.tool_list = [];
        tool_idx = 0;
        
        tool_idx = tool_idx + 1;
        obj.tool_list(tool_idx).event_key = 'd';
        obj.tool_list(tool_idx).shift_pressed = NaN;
        obj.tool_list(tool_idx).ctrl_pressed = NaN;
        obj.tool_list(tool_idx).type = 'delete';
        obj.tool_list(tool_idx).help_str = 'd: delete (set work area to bad values)';
        obj.tool_list(tool_idx).fh_callback = @obj.tool_delete;
        
        tool_idx = tool_idx + 1;
        obj.tool_list(tool_idx).event_key = 'i';
        obj.tool_list(tool_idx).shift_pressed = NaN;
        obj.tool_list(tool_idx).ctrl_pressed = NaN;
        obj.tool_list(tool_idx).type = 'interpolate';
        obj.tool_list(tool_idx).help_str = 'i: interpolate bad values in work area';
        obj.tool_list(tool_idx).fh_callback = @obj.tool_interpolate;
      else
        obj.tool_list = tool_list;
      end
      
      % Set active tool to the first tool (or zero if no tools)
      obj.active_tool_idx = min(1,length(obj.tool_list));

      % Set parameters
      obj.data = data;
      obj.bad_val = param.bad_val;
      obj.fh_button_up = param.fh_button_up;
      obj.fh_key_press = param.fh_key_press;
      obj.actions.getting_polygon = 0;
      
      if ~isempty(param.h_img)
        % Latch onto image, figure, axes if they do exist
        obj.h_img = param.h_img;
        obj.h_fig = get_parent_figure(obj.h_img);
        obj.h_axes = get(obj.h_img,'Parent');
        
        set(obj.h_img,'CData',obj.data);
        set(obj.h_img,'CDataMapping','scaled');
        set(obj.h_axes,'YDir','reverse');
      else
        % Create image, figure, and axes if they do not exist
        obj.h_fig = figure;
        set(obj.h_fig,'DockControls','off');
        set(obj.h_fig,'ToolBar','none');
        set(obj.h_fig,'MenuBar','none');
        obj.h_axes = axes('Parent',obj.h_fig);
        obj.h_img = imagesc(obj.data,'Parent',obj.h_axes);
      end

      % Set figure call back functions
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);

      % Set up undo stack
      obj.undo_stack = imb.undo_stack(struct('id',1));
      addlistener(obj.undo_stack,'synchronize_event',@obj.cmds_synchronize);
      
      % Set up zoom
      zoom_setup(obj.h_fig);
      obj.zoom_mode = true;
      set(obj.h_fig,'pointer','custom');

      % Set limits to size of data
      xlim([1 size(obj.data,2)]);
      ylim([1 size(obj.data,1)]);
      
      %obj.h_imgwin = imagewin(sprintf('%d: Image Params',obj.h_fig), obj.h_img, false);
      obj.create_ui();
    end
    
    function delete(obj)
      % Delete the map figure handle
      try
        delete(obj.h_fig);
      end
      %delete(obj.h_imgwin);
      try
        delete(obj.sf.h_fig);
      end
    end
    
    function close_win(obj,h_obj,event)
      try
        delete(obj);
      end
    end
    
    function button_down(obj,h_obj,event)
      if obj.actions.getting_polygon
        if obj.actions.getting_polygon == 2
          obj.actions.getting_polygon = 0;
        end
      end
      [obj.x,obj.y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      fprintf('Button Down: x = %.3f, y = %.3f, but = %d\n', obj.x, obj.y, but); % DEBUG ONLY
      rbbox;
    end
    
    function button_up(obj,h_obj,event)
      if obj.actions.getting_polygon
        return;
      end
      
      if ~isempty(obj.fh_button_up)
        status = obj.fh_button_up(obj,h_obj,event);
        if status == 0
          return;
        end
      end
      
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      fprintf('Button Up: x = %.3f, y = %.3f, but = %d\n', x, y, but); % DEBUG ONLY
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.x,'y',obj.y, ...
          'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
      elseif obj.active_tool_idx ~= 0
        % Tool mode and there is an active tool
        
        % Determine work area
        xlims = round(sort([x obj.x]));
        ylims = round(sort([y obj.y]));
        xlims(1) = max(xlims(1),1);
        ylims(1) = max(ylims(1),1);
        xlims(2) = min(xlims(2),size(obj.data,2));
        ylims(2) = min(ylims(2),size(obj.data,1));
        polygon_mask = logical(ones(ylims(2)-ylims(1)+1,xlims(2)-xlims(1)+1));

        % Apply active tool
        obj.tool_list(obj.active_tool_idx).fh_callback(xlims,ylims,polygon_mask);
      end
    end
    
    function button_scroll(obj,h_obj,event)
      zoom_button_scroll(event,struct('h_fig',obj.h_fig, ...
        'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
    end
    
    function key_press(obj,src,event)
      
      if any(strcmp('shift',event.Modifier))
        shift_pressed = true;
      else
        shift_pressed = false;
      end
      
      if any(strcmp('control',event.Modifier))
        ctrl_pressed = true;
      else
        ctrl_pressed = false;
      end
      
      if ~isempty(obj.fh_key_press)
        status = obj.fh_key_press(src,event);
        if status == 0
          return;
        end
      end
      
      % Check to make sure that a key was pressed and not
      % just a modifier (e.g. shift, ctrl, alt)
      if ~isempty(event.Key)
        
        % see event.Modifier for modifiers
        switch event.Key
          
          case 'f1'
            fprintf('Key Short Cuts\n');
            fprintf('f: open the selection filter window\n');
            fprintf('p: draw a polygon to select the working area\n');
            fprintf('u: undo the last command\n');
            fprintf('r: redo a command that has been undone\n');
            fprintf('z: Toggle between zoom mode and tool mode\n');
            fprintf('ctrl-z: Reset zoom to view whole image\n');
            fprintf('arrows: pan in the direction of the arrow\n');
            for tool_idx = 1:length(obj.tool_list)
              fprintf(obj.tool_list(tool_idx).help_str);
            end
            
            fprintf('Tool Mode\n');
            fprintf('left-click and drag: selects a rectangular work area\n');
            fprintf('scroll: zoom in/out at point\n');
            
            fprintf('Zoom Mode\n');
            fprintf('left-click and drag: zoom to selection\n');
            fprintf('left-click: zoom in at point\n');
            fprintf('right-click: zoom out at point\n');
            fprintf('scroll: zoom in/out at point\n');
            
          case 'z'
            if ctrl_pressed
              %% zoom reset
              axis tight;
            else
              %% toggle zoom mode
              obj.zoom_mode = ~obj.zoom_mode;
              if obj.zoom_mode
                set(obj.h_fig,'pointer','custom');
              else
                set(obj.h_fig,'pointer','arrow');
              end
            end
            
          case 'p'
            if obj.active_tool_idx ~= 0
              obj.actions.getting_polygon = 1;
              poly_handle = impoly;
              position = wait(poly_handle);
              obj.actions.getting_polygon = 2;
              [polyPts] = getPosition(poly_handle);
              xPoly = polyPts(:,1);
              yPoly = polyPts(:,2);
              delete(poly_handle);

              xlims = round([max(1,min(xPoly)) min(size(obj.data,2),max(xPoly))]);
              ylims = round([max(1,min(yPoly)) min(size(obj.data,1),max(yPoly))]);
              [x_mat y_mat] = meshgrid(xlims(1):xlims(2),ylims(1):ylims(2));

              polygon_mask = inpolygon(x_mat,y_mat,xPoly,yPoly);
            
              % Apply active tool
              obj.tool_list(obj.active_tool_idx).fh_callback(xlims,ylims,polygon_mask);
            end
            
          case 's'
            %% Make the selection filter window visible
            set(obj.sf.h_fig,'Visible','on')
            
          case 'r'
            %% Redo last tool operation
            if length(obj.undo_stack.stack) > obj.undo_stack.pointer
              obj.undo_stack.redo();
              cmd = obj.undo_stack.peak(); cmd = cmd{1};
              if ~isempty(cmd.sf_mask)
                set(obj.sf.h_gui.caxis_maskCB,'Value',cmd.sf_mask);
              end
              if ~isempty(cmd.sf_min_val)
                obj.sf.h_gui.min_slider.set_value(cmd.sf_min_val);
              end
              if ~isempty(cmd.sf_max_val)
                obj.sf.h_gui.max_slider.set_value(cmd.sf_max_val);
              end
            else
              fprintf('No more commands to redo.\n');
            end
            
          case 'u'
            %% Undo last tool operation
            if ~isempty(obj.undo_stack.peak())
              obj.undo_stack.pop();
              cmd = obj.undo_stack.peak();
              if ~isempty(cmd)
                cmd = cmd{1};
                if ~isempty(cmd.sf_mask)
                  set(obj.sf.h_gui.caxis_maskCB,'Value',cmd.sf_mask);
                end
                if ~isempty(cmd.sf_min_val)
                  obj.sf.h_gui.min_slider.set_value(cmd.sf_min_val);
                end
                if ~isempty(cmd.sf_max_val)
                  obj.sf.h_gui.max_slider.set_value(cmd.sf_max_val);
                end
              end
            else
              fprintf('No more commands to undo.\n');
            end
            
          case 'downarrow' % Down-arrow: Move Echogram Down
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'uparrow' % Up-arrow: Move Echogram Up
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'rightarrow' % Right arrow
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'leftarrow' % Left arrow
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          otherwise
            for tool_idx = 1:length(obj.tool_list)
              if strcmpi(obj.tool_list(tool_idx).event_key,event.Key)
                obj.active_tool_idx = tool_idx;
                break;
              end
            end
            
        end
        
      end
      
      if ~isempty(obj.fh_key_press)
        obj.fh_key_press(src,event)
      end
    end
    
    function set_data(obj,CData)
      obj.data = CData;
      set(obj.h_img,'CData',obj.data);
    end

    function tool_delete(obj,xlims,ylims,polygon_mask)
      % Set values in selection mask to bad value
      fprintf('Adding delete command %d %d %d %d\n', xlims, ylims); % DEBUG
      cmd.type = 'delete';
      cmd.xlims = xlims;
      cmd.ylims = ylims;
      cmd.polygon_mask = polygon_mask;
      cmd.old_vals = obj.data(ylims(1):ylims(2),xlims(1):xlims(2));
      cmd.new_vals = obj.data(ylims(1):ylims(2),xlims(1):xlims(2));
      
      cmd.sf_mask = get(obj.sf.h_gui.caxis_maskCB,'Value');
      if cmd.sf_mask
        auto_setting = get(obj.sf.h_gui.caxis_autoCB,'Value');
        if auto_setting
          mean_val = nanmean(cmd.old_vals(polygon_mask));
          std_val = nanstd(cmd.old_vals(polygon_mask));
          num_std = str2double(get(obj.sf.h_gui.caxis_num_stdLE,'String'));
          if num_std > 0 && ~isnan(mean_val) && ~isnan(std_val)
            obj.sf.h_gui.min_slider.set_value(mean_val-num_std*std_val);
            obj.sf.h_gui.max_slider.set_value(mean_val+num_std*std_val);
          end
        end
        cmd.sf_min_val = obj.sf.h_gui.min_slider.get_value();
        cmd.sf_max_val = obj.sf.h_gui.max_slider.get_value();
        
        selection_mask = cmd.polygon_mask ...
          & (cmd.old_vals < cmd.sf_min_val | cmd.old_vals > cmd.sf_max_val);
        
      else
        cmd.sf_min_val = [];
        cmd.sf_max_val = [];
        selection_mask = cmd.polygon_mask;
      end
      
      cmd.new_vals(selection_mask) = obj.bad_val;
      
      obj.undo_stack.push(cmd);
    end
    
    function tool_interpolate(obj,xlims,ylims,polygon_mask)
      % Interpolate over all bad pixels (i.e. pixels with a value equal to
      %   bad value) in polygon mask
      fprintf('Adding interpolate command %d %d %d %d\n', xlims, ylims); % DEBUG
      cmd.type = 'interpolate';
      cmd.xlims = xlims;
      cmd.ylims = ylims;
      cmd.polygon_mask = polygon_mask;
      cmd.old_vals = obj.data(ylims(1):ylims(2),xlims(1):xlims(2));
      cmd.new_vals = obj.data(ylims(1):ylims(2),xlims(1):xlims(2));

      selection_mask = cmd.polygon_mask;
    
      if isnan(obj.bad_val)
        bad_idxs = find(isnan(cmd.new_vals) & selection_mask);
        good_idxs = find(~isnan(cmd.new_vals) & selection_mask);
      else
        bad_idxs = find(cmd.new_vals == obj.bad_val & selection_mask);
        good_idxs = find(cmd.new_vals ~= obj.bad_val & selection_mask);
      end
      x_idxs = repmat(xlims(1):xlims(2),[size(cmd.new_vals,1) 1]);
      y_idxs = repmat((ylims(1):ylims(2))',[1 size(cmd.new_vals,2)]);
      x_vals = x_idxs(good_idxs);
      y_vals = y_idxs(good_idxs);
      z_vals = cmd.new_vals(good_idxs);
      x_out = x_idxs(bad_idxs);
      y_out = y_idxs(bad_idxs);
      try
        z_out = griddata(x_vals,y_vals,z_vals,x_out,y_out);
      catch ME
        fprintf('Grid data command failed, ignoring\n');
        %ME.getReport
        return;
      end
      cmd.new_vals(bad_idxs) = z_out;
      
      cmd.sf_mask = [];
      cmd.sf_min_val = [];
      cmd.sf_max_val = [];

      obj.undo_stack.push(cmd);
    end
    
    function cmds_synchronize(obj,varargin)
      % Applies commands from undo stack
      
      % Get the list of commands from the undo stack that need to be run to
      % synchronize the echowin
      [cmds_list,cmds_direction] = obj.undo_stack.get_synchronize_cmds();
      
      % Execute the commands
      if strcmpi(cmds_direction,'undo')
        for cmd_idx = 1:length(cmds_list)
          obj.data(cmds_list{cmd_idx}.ylims(1):cmds_list{cmd_idx}.ylims(2), ...
            cmds_list{cmd_idx}.xlims(1):cmds_list{cmd_idx}.xlims(2)) = cmds_list{cmd_idx}.old_vals;
        end
      elseif strcmpi(cmds_direction,'redo')
        for cmd_idx = 1:length(cmds_list)
          obj.data(cmds_list{cmd_idx}.ylims(1):cmds_list{cmd_idx}.ylims(2), ...
            cmds_list{cmd_idx}.xlims(1):cmds_list{cmd_idx}.xlims(2)) = cmds_list{cmd_idx}.new_vals;
        end
      end
      set(obj.h_img,'CData',obj.data);
    end
    
    
    % ====================================================================
    % ====================================================================
    % ====================================================================

    function sf_hide_win(obj,varargin)
      set(obj.sf.h_fig,'Visible','off')
    end
    
    function sf_update(obj,varargin)
      % Get the last command off the undostack
      cmd = obj.undo_stack.peak();
      if isempty(cmd)
        return;
      end
      obj.undo_stack.pop();
      
      % Put the command back on the undostack
      for tool_idx = 1:length(obj.tool_list)
        if strcmpi(obj.tool_list(tool_idx).type,cmd{1}.type)
          obj.tool_list(obj.active_tool_idx).fh_callback(cmd{1}.xlims,cmd{1}.ylims,cmd{1}.polygon_mask);
          break;
        end
      end
    end
    
    function create_ui(obj)
      % clip_matrix.create_ui(obj)
      %
      % Create the selection filter figure and GUI
      
      clims = [min(obj.data(:)) max(obj.data(:))];

      %% Create Figure
      obj.sf.h_fig = figure;
      %obj.sf.h_fig = figure('Visible','off');
      
      set(obj.sf.h_fig,'DockControls','off');
      set(obj.sf.h_fig,'NumberTitle','off');
      if isnumeric(obj.h_fig)
        set(obj.sf.h_fig,'Name',sprintf('%d: Selection Filter', obj.h_fig));
      else
        set(obj.sf.h_fig,'Name',sprintf('%d: Selection Filter', obj.h_fig.Number));
      end
      set(obj.sf.h_fig,'ToolBar','none');
      set(obj.sf.h_fig,'MenuBar','none');
      set(obj.sf.h_fig,'CloseRequestFcn',@obj.sf_hide_win);
      set(obj.sf.h_fig,'Units','Points');
      pos = get(obj.sf.h_fig,'Position');
      set(obj.sf.h_fig,'Position',[pos(1:2) 240 60]);
      
      %% Create GUI Objects
       
      obj.sf.h_gui.caxis_maskCB = uicontrol('Parent',obj.sf.h_fig);
      set(obj.sf.h_gui.caxis_maskCB,'Style','CheckBox');
      set(obj.sf.h_gui.caxis_maskCB,'String','Mask');
      set(obj.sf.h_gui.caxis_maskCB,'Value',false);
      set(obj.sf.h_gui.caxis_maskCB,'Callback',@obj.sf_update);
      
      obj.sf.h_gui.caxis_autoCB = uicontrol('Parent',obj.sf.h_fig);
      set(obj.sf.h_gui.caxis_autoCB,'Style','CheckBox');
      set(obj.sf.h_gui.caxis_autoCB,'String','Auto caxis');
      set(obj.sf.h_gui.caxis_autoCB,'Value',false);
      set(obj.sf.h_gui.caxis_autoCB,'Callback',@obj.sf_update);
      
      obj.sf.h_gui.caxis_num_stdLE = uicontrol('Parent',obj.sf.h_fig);
      set(obj.sf.h_gui.caxis_num_stdLE,'Style','Edit');
      set(obj.sf.h_gui.caxis_num_stdLE,'String','1');
      set(obj.sf.h_gui.caxis_num_stdLE,'Callback',@obj.sf_update);

      obj.sf.h_gui.min_slider = slider(obj.sf.h_fig,'Min',clims,clims(1));
      
      obj.sf.h_gui.max_slider = slider(obj.sf.h_fig,'Max',clims,clims(end));
      
      addlistener(obj.sf.h_gui.min_slider,'slider_changed',@obj.sf_update);
      addlistener(obj.sf.h_gui.max_slider,'slider_changed',@obj.sf_update);
      
      %% Create GUI Table
      obj.sf.h_gui.table.ui = obj.sf.h_fig;
      obj.sf.h_gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.sf.h_gui.table.height_margin = NaN*zeros(30,30);
      obj.sf.h_gui.table.false_width = NaN*zeros(30,30);
      obj.sf.h_gui.table.false_height = NaN*zeros(30,30);
      obj.sf.h_gui.table.offset = [0 0];
      
      row = 1;
      col = 1;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.caxis_maskCB;
      obj.sf.h_gui.table.width(row,col)     = inf;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      col = 2;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.caxis_autoCB;
      obj.sf.h_gui.table.width(row,col)     = inf;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      col = 3;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.caxis_num_stdLE;
      obj.sf.h_gui.table.width(row,col)     = inf;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.min_slider.h_text;
      obj.sf.h_gui.table.width(row,col)     = 50;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 5;
      col = 2;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.min_slider.h_slider;
      obj.sf.h_gui.table.width(row,col)     = inf;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      col = 3;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.min_slider.h_LE;
      obj.sf.h_gui.table.width(row,col)     = 50;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      
      row = row + 1;
      col = 1;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.max_slider.h_text;
      obj.sf.h_gui.table.width(row,col)     = 50;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 5;
      col = 2;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.max_slider.h_slider;
      obj.sf.h_gui.table.width(row,col)     = inf;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      col = 3;
      obj.sf.h_gui.table.handles{row,col}   = obj.sf.h_gui.max_slider.h_LE;
      obj.sf.h_gui.table.width(row,col)     = 50;
      obj.sf.h_gui.table.height(row,col)    = 20;
      obj.sf.h_gui.table.width_margin(row,col) = 1;
      obj.sf.h_gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.sf.h_gui.table);
    end

    
    
  end
end




