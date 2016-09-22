% script imb.slice_browser
%
% Class for browsing 3D imagery one 2D slice at a time and for editing
% layers in that imagery.
%
% Contructor: slice_browser(data,h_control_image,param)
% data: 3D imagery
% h_control_image: optional handle to a Matlab "image" which has an x-axis
%   aligned with the third dimension of the data. Clicks in this figure
%   will then choose difference slices out of data based on the third axis.
% param: structure controlling the operation of slice_browser
%  .layer_fn: filename of .mat file containing layer structure array
%
% Layer file should contain:
%  layer: structure array of layer information
%   .x: x-values of the layer
%   .y: y-values of the layer
%   .plot_name_values: name-value pairs to be passed to this layer's plot
%     function
%   .name: name of this layer
%
% Example:
%  See run_slice_browser.m
%
% Author: Elijah Paden, John Paden

classdef slice_browser < handle
  
  properties
    select_mask
    data % N-dimensional matrix with last dimension equal to Nx
    slice % Integer from 1 to Nx
    layer % Layer structures
    layer_fn
    layer_idx % Active layer

    % Slice GUI handles
    h_fig
    h_axes
    h_image
    x, y % Last click position
    
    % Control GUI handles
    h_control_fig
    h_control_axes
    h_control_image
    h_control_plot
    h_control_layer
    h_control_is_child % logical (true means slice_browser created control)
    control_x, control_y % Last click position
    
    % Layer GUI handles
    h_layer_fig
    h_layer_axes
    h_layer_image
    h_layer_plot
    layer_x, layer_y % Last click position
    
    gui
    
    % slice_tool: .list, .timer, .cluster, .cmds
    slice_tool
    
    % Function handle hooks for customizing clip_matrix
    fh_button_up
    fh_key_press
    fh_button_motion
    
    % zoom_mode: boolean, x,y: used by zoom mode
    zoom_mode
    shift_pressed
    ctrl_pressed
    plot_visibility
    undo_stack
    save_callback
  end
  
  events
    SliceChange
  end
  
  methods
    %% constructor/slice_browser:
    function obj = slice_browser(data,h_control_image,param)
      if ~exist('param','var')
        param = [];
      end
      if ~isfield(param,'fh_button_up')
        param.fh_button_up = [];
      end
      if ~isfield(param,'fh_key_press')
        param.fh_key_press = [];
      end
      if ~isfield(param,'fh_button_motion')
        param.fh_button_motion = [];
      end
      undo_param.id = [];
      obj.undo_stack = imb.undo_stack(undo_param);
      obj.data = data;
      obj.slice = 1;
      obj.plot_visibility = true;
      
      obj.slice_tool.list = [];
      
      obj.slice_tool.timer = timer;
      obj.slice_tool.timer.StartDelay = 2;
      obj.slice_tool.timer.Period = 2;
      obj.slice_tool.timer.TimerFcn = @obj.timer_callback;
      obj.slice_tool.timer.ExecutionMode = 'fixedSpacing';
      %start(obj.slice_tool.timer)
      
      % Load layer data
      if isfield(param,'layer_fn') && ~isempty(param.layer_fn)
        tmp = load(param.layer_fn);
        obj.layer_fn = param.layer_fn;
        obj.layer = tmp.surf;
      else
        obj.layer = struct('x',{},'y',{},'name',{},'plot_name_values',{});
      end
      
      if ~isempty(h_control_image)
        obj.h_control_is_child = false;
        obj.h_control_image = h_control_image;
        obj.h_control_axes = get(obj.h_control_image,'Parent');
        obj.h_control_fig = get(obj.h_control_axes,'Parent');
      else
        obj.h_control_is_child = true;
        obj.h_control_fig = figure;
        obj.h_control_axes = axes('Parent',obj.h_control_fig);
        obj.h_control_image = imagesc(10*log10(squeeze(obj.data(:,floor(size(data,2)/2)+1,:))),'Parent',obj.h_control_axes);
        colormap(obj.h_control_axes,parula(256));
        xlabel(obj.h_control_axes,'Along-track range line');
        ylabel(obj.h_control_axes,'Range bin');
      end
      obj.fh_button_up = param.fh_button_up;
      obj.fh_key_press = param.fh_key_press;
      obj.fh_button_motion = param.fh_button_motion;
      
      obj.h_layer_fig = figure;
      pos = get(obj.h_fig,'Position');
      pos(3) = 750;
      pos(4) = 500;
      set(obj.h_fig,'Position',pos);

      obj.h_layer_axes = axes('Parent',obj.h_layer_fig,'YDir','reverse');
      obj.h_layer_image = imagesc(NaN*zeros(size(obj.data,2),size(obj.data,3)),'parent',obj.h_layer_axes);
      colormap(obj.h_layer_axes, parula(256));
      xlabel(obj.h_layer_axes,'Along-track range line');
      ylabel(obj.h_layer_axes,'Cross-track');
      hold(obj.h_layer_axes,'on');
      obj.h_layer_plot = plot(NaN,NaN,'parent',obj.h_layer_axes,'Marker','x','Color','black','LineWidth',2,'MarkerSize',10);
      
      obj.h_fig = figure;

      pos = get(obj.h_fig,'Position');
      pos(3) = 750;
      pos(4) = 500;
      set(obj.h_fig,'Position',pos);
      
      obj.gui.left_panel = uipanel('parent',obj.h_fig);
      obj.gui.right_panel = uipanel('parent',obj.h_fig);
      obj.h_axes = axes('Parent',obj.gui.right_panel,'YDir','reverse');
      hold(obj.h_axes,'on');
      colormap(obj.h_axes, parula(256));
      xlabel(obj.h_axes,'Cross-track');
      ylabel(obj.h_axes,'Range bin');
      
      obj.h_image = imagesc(obj.data(:,:,obj.slice),'parent',obj.h_axes);
      for layer_idx = 1:numel(obj.layer)
        if islogical(obj.layer(layer_idx).y)
          tmp_y = obj.layer(obj.layer(layer_idx).surf_layer).y(:,obj.slice);
          tmp_y(~obj.layer(layer_idx).y(:,obj.slice)) = NaN;
          
          obj.layer(layer_idx).h_plot ...
          = plot(obj.layer(layer_idx).x(:,obj.slice), ...
          tmp_y, 'parent',obj.h_axes,'color','black', ...
          obj.layer(layer_idx).plot_name_values{:});
        else
          obj.layer(layer_idx).h_plot ...
            = plot(obj.layer(layer_idx).x(:,obj.slice), ...
            obj.layer(layer_idx).y(:,obj.slice), ...
            'parent',obj.h_axes,'color','black', ...
            obj.layer(layer_idx).plot_name_values{:});
        end
      end
      
      addlistener(obj.undo_stack,'synchronize_event',@obj.undo_sync);
      
      obj.gui.h_select_plot = plot(NaN,NaN,'m.');
      
      hold(obj.h_control_axes,'on');
      obj.h_control_plot = plot(NaN,NaN,'parent',obj.h_control_axes,'Marker','x','Color','black','LineWidth',2,'MarkerSize',10);
      obj.h_control_layer = plot(NaN,NaN,'parent',obj.h_control_axes,'Marker','.','Color','red');

      % Set up figure callbacks and zoom
      zoom_figure_setup(obj.h_fig,'slice');
      obj.zoom_mode = true;
      set(obj.h_fig,'pointer','custom');
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowButtonMotionFcn',@obj.button_motion);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig,'WindowKeyReleaseFcn',@obj.key_release);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      
      zoom_figure_setup(obj.h_layer_fig,'layer');
      set(obj.h_layer_fig,'pointer','custom');
      set(obj.h_layer_fig,'WindowButtonUpFcn',@obj.layer_button_up);
      set(obj.h_layer_fig,'WindowButtonDownFcn',@obj.layer_button_down);
%       set(obj.h_layer_fig,'WindowButtonMotionFcn',@obj.button_motion);
      set(obj.h_layer_fig,'WindowScrollWheelFcn',@obj.layer_button_scroll);
      set(obj.h_layer_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_layer_fig,'WindowKeyReleaseFcn',@obj.key_release);
      set(obj.h_layer_fig,'CloseRequestFcn',[]);
      
      if obj.h_control_is_child
        zoom_figure_setup(obj.h_control_fig,'echogram');
        set(obj.h_control_fig,'pointer','custom');
        set(obj.h_control_fig,'WindowButtonUpFcn',@obj.control_button_up);
        set(obj.h_control_fig,'WindowButtonDownFcn',@obj.control_button_down);
%         set(obj.h_control_fig,'WindowButtonMotionFcn',@obj.button_motion);
        set(obj.h_control_fig,'WindowScrollWheelFcn',@obj.control_button_scroll);
        set(obj.h_control_fig,'WindowKeyPressFcn',@obj.key_press);
        set(obj.h_control_fig,'WindowKeyReleaseFcn',@obj.key_release);
        set(obj.h_control_fig,'CloseRequestFcn',[]);
      end
      
      obj.select_mask = logical(zeros(size(obj.data,2),1));
      
      obj.gui.table.ui = obj.h_fig;
      obj.gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.table.height_margin = NaN*zeros(30,30);
      obj.gui.table.false_width = NaN*zeros(30,30);
      obj.gui.table.false_height = NaN*zeros(30,30);
      obj.gui.table.offset = [0 0];
      row = 1;
      col = 1;
      obj.gui.table.handles{row,col}   = obj.gui.left_panel;
      obj.gui.table.width(row,col)     = 130;
      obj.gui.table.height(row,col)    = inf;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      row = 1;
      col = 2;
      obj.gui.table.handles{row,col}   = obj.gui.right_panel;
      obj.gui.table.width(row,col)     = inf;
      obj.gui.table.height(row,col)    = inf;
      obj.gui.table.width_margin(row,col) = 1;
      obj.gui.table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.table);
      
      % Set limits to size of data
      xlim(obj.h_axes, [1 size(obj.data(:,:,obj.slice),2)]);
      ylim(obj.h_axes, [1 size(obj.data(:,:,obj.slice),1)]);
      
      obj.gui.nextPB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.nextPB,'style','pushbutton')
      set(obj.gui.nextPB,'string','>')
      set(obj.gui.nextPB,'Callback',@obj.next_button_callback)
      set(obj.gui.nextPB,'TooltipString','Move forward one slice (.)');
      
      obj.gui.prevPB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.prevPB,'style','pushbutton')
      set(obj.gui.prevPB,'string','<')
      set(obj.gui.prevPB,'Callback',@obj.prev_button_callback)
      set(obj.gui.prevPB,'TooltipString','Move backward one slice (,)');
      
      obj.gui.prev10PB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.prev10PB,'style','pushbutton')
      set(obj.gui.prev10PB,'string','<<')
      set(obj.gui.prev10PB,'Callback',@obj.prev10_button_callback)
      set(obj.gui.prev10PB,'TooltipString','Move backward ten slices (<)');
      
      obj.gui.next10PB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.next10PB,'style','pushbutton')
      set(obj.gui.next10PB,'string','>>')
      set(obj.gui.next10PB,'Callback',@obj.next10_button_callback)
      set(obj.gui.next10PB,'TooltipString','Move forward ten slices (>)');
      
      obj.gui.savePB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.savePB,'style','pushbutton')
      set(obj.gui.savePB,'string','(S)ave')
      set(obj.gui.savePB,'Callback',@obj.save_button_callback)
      set(obj.gui.savePB,'TooltipString','(S)ave layers to file');
      
      obj.gui.helpPB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.helpPB,'style','pushbutton')
      set(obj.gui.helpPB,'string','Help (F1)')
      set(obj.gui.helpPB,'Callback',@obj.help_button_callback)
      set(obj.gui.helpPB,'TooltipString','Print help to stdout (F1)');
      
      obj.gui.layerTXT = uicontrol('Style','text','string','Layer');
      obj.gui.layerLB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.layerLB,'style','listbox')
      set(obj.gui.layerLB,'string',{obj.layer.name})
      set(obj.gui.layerLB,'Callback',@obj.layerLB_callback)
      set(obj.gui.layerLB,'TooltipString','Select active layer (#)');
      obj.gui.layerCM = uicontextmenu('Parent',obj.h_fig);
      % Define the context menu items and install their callbacks
      obj.gui.layerCM_visible = uimenu(obj.gui.layerCM, 'Label', 'Toggle Visible', 'Callback', @obj.layerLB_visibility_toggle);
      set(obj.gui.layerLB,'UIContextMenu',obj.gui.layerCM);
      
      obj.gui.applyPB= uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.applyPB,'style','pushbutton')
      set(obj.gui.applyPB,'string','Apply')
      set(obj.gui.applyPB,'Callback',@obj.applyPB_callback)
      set(obj.gui.applyPB,'TooltipString','Apply selected tool');

      obj.gui.optionsPB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.optionsPB,'style','pushbutton');
      set(obj.gui.optionsPB,'string','Options');
      set(obj.gui.optionsPB,'Callback',@obj.optionsPB_callback);
      set(obj.gui.optionsPB,'TooltipString','Open tool options window');
      
      obj.gui.toolPM = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.toolPM,'style','popup');
      set(obj.gui.toolPM,'string',{''});
      set(obj.gui.toolPM,'TooltipString','Select active tool');
      
      %% Create GUI Table
      obj.gui.left_table.ui =  obj.gui.left_panel;
      obj.gui.left_table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
      obj.gui.left_table.height_margin = NaN*zeros(30,30);
      obj.gui.left_table.false_width = NaN*zeros(30,30);
      obj.gui.left_table.false_height = NaN*zeros(30,30);
      obj.gui.left_table.offset = [0 0];
      
      row = 0;
      col = 0;
      
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.prev10PB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.prevPB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.nextPB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.next10PB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = 0;
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.savePB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.helpPB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = 0;
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.layerTXT;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = 0;
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.layerLB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = inf;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = 0;
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.toolPM;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = 0;
      row = row + 1;
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.applyPB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      col = col + 1;
      obj.gui.left_table.handles{row,col}   = obj.gui.optionsPB;
      obj.gui.left_table.width(row,col)     = inf;
      obj.gui.left_table.height(row,col)    = 20;
      obj.gui.left_table.width_margin(row,col) = 1;
      obj.gui.left_table.height_margin(row,col) = 1;
      
      clear row col
      table_draw(obj.gui.left_table);
      
      obj.update_slice();
    end
    
    %% destructor/delete
    function delete(obj)
      try; delete(obj.h_fig); end;
      try; delete(obj.h_layer_fig); end;
      for tool_idx = 1:length(obj.slice_tool.list)
        try; delete(obj.slice_tool.list{tool_idx}); end;
      end
      try; delete(obj.slice_tool.timer); end;

      % Control figure destructor
      try; set(obj.h_control_fig, 'WindowButtonUpFcn', []); end;
      if obj.h_control_is_child
        try; delete(obj.h_control_fig); end;
      end
      try; delete(obj.h_control_plot); end;
      try; delete(obj.h_control_layer); end;
    end
    
    %% close_win
    function close_win(obj,h_obj,event)
      try
        delete(obj);
      end
    end
    
    %% next_button_callback
    function next_button_callback(obj,source,callbackdata)
      obj.change_slice(obj.slice + 1,false);
    end
    
    %% prev_button_callback
    function prev_button_callback(obj,source,callbackdata)
      obj.change_slice(obj.slice - 1,false);
    end
    
    %% help_button_callback
    function help_button_callback(obj,source,callbackdata)
      obj.help_menu()
    end
    
    %% save_button_callback
    function save_button_callback(obj,source,callbackdata)
      obj.save();
    end
    
    %% next10_button_callback
    function next10_button_callback(obj,source,callbackdata)
      obj.change_slice(obj.slice + 10,false);
    end
    
    %% prev10_button_callback
    function prev10_button_callback(obj,source,callbackdata)
      obj.change_slice(obj.slice - 10,false);
    end
    
    %% undo_sync
    function undo_sync(obj,source,callbackdata)
      [cmds_list,cmds_direction] =  obj.undo_stack.get_synchronize_cmds();
      new_slice = [];
      if strcmp(cmds_direction,'redo')
        for cmd_idx = 1:length(cmds_list)
          for subcmd_idx = 1:length(cmds_list{cmd_idx})
            if strcmp(cmds_list{cmd_idx}{subcmd_idx}.type,'standard')
              layer_idx = cmds_list{cmd_idx}{subcmd_idx}.redo.layer;
              obj.layer(layer_idx).y(round(cmds_list{cmd_idx}{subcmd_idx}.redo.x), ...
                cmds_list{cmd_idx}{subcmd_idx}.redo.slice) ...
                = cmds_list{cmd_idx}{subcmd_idx}.redo.y;
              new_slice = cmds_list{cmd_idx}{subcmd_idx}.redo.slice;
            elseif strcmp(cmds_list{cmd_idx}{subcmd_idx}.type,'slice_dummy')
              new_slice = cmds_list{cmd_idx}{subcmd_idx}.redo.slice;
            end
          end
        end
      else
        for cmd_idx = 1:length(cmds_list)
          for subcmd_idx = 1:length(cmds_list{cmd_idx})
            if strcmp(cmds_list{cmd_idx}{subcmd_idx}.type,'standard')
              layer_idx = cmds_list{cmd_idx}{subcmd_idx}.undo.layer;
              obj.layer(layer_idx).y(round(cmds_list{cmd_idx}{subcmd_idx}.undo.x), ...
                cmds_list{cmd_idx}{subcmd_idx}.undo.slice) ...
                = cmds_list{cmd_idx}{subcmd_idx}.undo.y;
              new_slice = cmds_list{cmd_idx}{subcmd_idx}.undo.slice;
            elseif strcmp(cmds_list{cmd_idx}{subcmd_idx}.type,'slice_dummy')
              new_slice = cmds_list{cmd_idx}{subcmd_idx}.redo.slice;
            end
          end
        end
      end
      obj.change_slice(new_slice,true);

    end  
    
    %% button_down
    function button_down(obj,h_obj,event)
      [obj.x,obj.y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      %fprintf('Button Down: x = %.3f, y = %.3f, but = %d\n', obj.x, obj.y, but); % DEBUG ONLY
      rbbox;
    end
    
    function control_button_down(obj,h_obj,event)
      [obj.control_x,obj.control_y,but] = get_mouse_info(obj.h_control_fig,obj.h_control_axes);
      %fprintf('Button Down: x = %.3f, y = %.3f, but = %d\n', obj.x, obj.y, but); % DEBUG ONLY
      rbbox;
    end
    
    function layer_button_down(obj,h_obj,event)
      [obj.layer_x,obj.layer_y,but] = get_mouse_info(obj.h_layer_fig,obj.h_layer_axes);
      %fprintf('Button Down: x = %.3f, y = %.3f, but = %d\n', obj.x, obj.y, but); % DEBUG ONLY
      rbbox;
    end
    
    %% button_up
    function button_up(obj,h_obj,event)
      % Run user defined button up callback
      if ~isempty(obj.fh_button_up)
        status = obj.fh_button_up(obj,h_obj,event);
        if status == 0
          return;
        end
      end
      
      % Get x,y position of user button release
      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      %fprintf('Button Up: x = %.3f, y = %.3f, but = %d\n', x, y, but); % DEBUG ONLY
      
      set(obj.h_fig,'Units','normalized');
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      set(obj.gui.right_panel,'Units','normalized');
      uipanel_pos = get(obj.gui.right_panel,'Position');
      if mouse_pos(1) < uipanel_pos(1)
        return;
      end
      
      layer_idx = get(obj.gui.layerLB,'value');
      
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.x,'y',obj.y, ...
          'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
      else
        if but == 2 || but == 3
          if obj.x == x
            if x > 1 && x < size(obj.data,2)
              obj.select_mask(round(x)) = ~obj.select_mask(round(x));
            else
              obj.select_mask(:) = false;
            end
          else
            obj.shift_pressed
            if ~obj.shift_pressed
              obj.select_mask(:) = false;
              obj.update_slice;
            end
            
            if islogical(obj.layer(layer_idx).y)
              layer_y = obj.layer(obj.layer(layer_idx).surf_layer).y(:,obj.slice);
            else
              layer_y = obj.layer(layer_idx).y(:,obj.slice);
            end

            obj.select_mask = obj.select_mask | (obj.layer(layer_idx).x(:,obj.slice) >= min(x,obj.x) ...
              & obj.layer(layer_idx).x(:,obj.slice) <= max(x,obj.x) ...
              & layer_y >= min(y,obj.y) ...
              & layer_y <= max(y,obj.y));
          end
        else
          xlims = xlim(obj.h_axes);
          ylims = ylim(obj.h_axes);
          if x >= xlims(1) && x <= xlims(2) && y >= ylims(1) && y <= ylims(2)
            layer_idx = get(obj.gui.layerLB,'value');
            cmd = [];
            cmd{1}.undo.slice = obj.slice;
            cmd{1}.redo.slice = obj.slice;
            cmd{1}.undo.layer = layer_idx;
            cmd{1}.redo.layer = layer_idx;
            cmd{1}.undo.x = round(x);
            cmd{1}.undo.y = obj.layer(layer_idx).y(round(x),obj.slice);
            cmd{1}.redo.x = round(x);
            if islogical(obj.layer(layer_idx).y)
              cmd{1}.redo.y = ~obj.layer(layer_idx).y(round(x),obj.slice);
            else
              cmd{1}.redo.y = y;
            end
            cmd{1}.type = 'standard';
            obj.push(cmd);
          end
        end
        obj.update_slice();
      end
    end
    
    function control_button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_control_fig,obj.h_control_axes);
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.control_x,'y',obj.control_y, ...
          'h_axes',obj.h_control_axes,'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,1)]));
      else
        obj.change_slice(round(x),false);
      end
    end
    
    function layer_button_up(obj,h_obj,event)
      [x,y,but] = get_mouse_info(obj.h_layer_fig,obj.h_layer_axes);
      
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.layer_x,'y',obj.layer_y, ...
          'h_axes',obj.h_layer_axes,'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,2)]));
      else
        if obj.layer_x == x
          obj.change_slice(round(x),false);
        else
          ylims = sort([y obj.layer_y]);
          obj.select_mask(:) = false;
          y_idxs = round(ylims(1)):round(ylims(2));
          y_idxs = y_idxs(y_idxs>=1 & y_idxs<=size(obj.data,2));
          obj.select_mask(y_idxs) = true;
          if but ~= 1
            for tool_idx = 1:length(obj.slice_tool.list)
              tool_name_list{tool_idx} = obj.slice_tool.list{tool_idx}.tool_name;
            end
            xlims = sort([x obj.layer_x]);
            slices = round(xlims(1)):round(xlims(2));
            slices = slices(slices>=1 & slices<=size(obj.data,3));
            title(obj.h_layer_axes,sprintf('Slices %d-%d, DOAs %d-%d\n', slices(1), slices(end), y_idxs(1), y_idxs(end)));
            [tool_idx,ok] = listdlg('PromptString','Choose slicetool:',...
              'SelectionMode','single',...
              'ListString',tool_name_list);
            if ok == 1
              obj.layer_idx = get(obj.gui.layerLB,'Value');
              cmd = obj.slice_tool.list{tool_idx}.apply_PB_callback(obj,slices);
              if ~isempty(cmd)
                obj.undo_stack.push(cmd);
              end
            else
              
            end
          end
          obj.update_slice();
        end
      end
    end

    
    %% button_motion
    function button_motion(obj,h_obj,event)
      % Run user defined button up callback
      if ~isempty(obj.fh_button_motion)
        status = obj.fh_button_motion(obj,h_obj,event);
        if status == 0
          return;
        end
      end
      
      set(obj.h_fig,'Units','normalized');
      mouse_pos = get(obj.h_fig,'CurrentPoint');
      set(obj.gui.right_panel,'Units','normalized');
      uipanel_pos = get(obj.gui.right_panel,'Position');
      if mouse_pos(1) < uipanel_pos(1)
        set(obj.h_fig,'Pointer','Arrow');
        return;
      elseif obj.zoom_mode
        set(obj.h_fig,'Pointer','custom');
      end

      [x,y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
      set(obj.h_control_plot,'XData',obj.slice,'YData',y);
      set(obj.h_layer_plot,'XData',obj.slice,'YData',x);
    end
    
    %% button_scroll
    function button_scroll(obj,h_obj,event)
      zoom_button_scroll(event,struct('h_fig',obj.h_fig, ...
        'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
    end
    
    function control_button_scroll(obj,h_obj,event)
      zoom_button_scroll(event,struct('h_fig',obj.h_control_fig, ...
        'h_axes',obj.h_control_axes,'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,1)]));
    end
    
    function layer_button_scroll(obj,h_obj,event)
      zoom_button_scroll(event,struct('h_fig',obj.h_layer_fig, ...
        'h_axes',obj.h_layer_axes,'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,2)]));
    end
    
    %% key_press
    function key_press(obj,src,event)
      
      if any(strcmp('shift',event.Modifier))
        obj.shift_pressed = true;
      else
        obj.shift_pressed = false;
      end
      
      if any(strcmp('control',event.Modifier))
        obj.ctrl_pressed = true;
      else
        obj.ctrl_pressed = false;
      end
      
      if ~isempty(obj.fh_key_press)
        status = obj.fh_key_press(src,event);
        if status == 0
          return;
        end
      end
      
      % Check to see if this is a slicetool shortcut
      for tool_idx = 1:length(obj.slice_tool.list)
        if strcmpi(obj.slice_tool.list{tool_idx}.tool_shortcut, event.Key) ...
            && obj.slice_tool.list{tool_idx}.ctrl_pressed == obj.ctrl_pressed ...
            && obj.slice_tool.list{tool_idx}.shift_pressed == obj.shift_pressed
          obj.layer_idx = get(obj.gui.layerLB,'Value');
          cmd = obj.slice_tool.list{tool_idx}.apply_PB_callback(obj);
          if ~isempty(cmd)
            obj.undo_stack.push(cmd);
          end
          return;
        end
      end
      
      % Check to make sure that a key was pressed and not
      % just a modifier (e.g. shift, ctrl, alt)
      if ~isempty(event.Key)
        
        if length(event.Key) == 1 && event.Key >= '0' && event.Key <= '9'
          set(obj.gui.layerLB,'value',event.Key-48)
          obj.update_slice();
          return;
        end
        
        % see event.Modifier for modifiers
        switch event.Key
          
          case 'f1'
            obj.help_menu()
            
          case 'z'
            if obj.ctrl_pressed
              %% zoom reset
              if src==obj.h_fig
                axis(obj.h_axes,'tight');
              elseif src==obj.h_control_fig
                axis(obj.h_control_axes,'tight');
              elseif src==obj.h_layer_fig
                axis(obj.h_layer_axes,'tight');
              end
              
            else
              % toggle zoom mode
              obj.zoom_mode = ~obj.zoom_mode;
              if obj.zoom_mode
                set(obj.h_fig,'pointer','custom');
                if obj.h_control_is_child
                  set(obj.h_control_fig,'pointer','custom');
                end
                set(obj.h_layer_fig,'pointer','custom');
              else
                set(obj.h_fig,'pointer','arrow');
                if obj.h_control_is_child
                  set(obj.h_control_fig,'pointer','arrow');
                end
                set(obj.h_layer_fig,'pointer','arrow');
              end
            end
            
          case {'downarrow','uparrow','rightarrow','leftarrow'}
            % Arrows: pan axes
            if src==obj.h_fig
              zoom_arrow(event,struct('h_axes',obj.h_axes, ...
                'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            elseif src==obj.h_control_fig
              zoom_arrow(event,struct('h_axes',obj.h_control_axes, ...
                'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,1)]));
            elseif src==obj.h_layer_fig
              zoom_arrow(event,struct('h_axes',obj.h_layer_axes, ...
                'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,2)]));
            end
            
          case 'period'
            if ~obj.shift_pressed
              obj.change_slice(obj.slice + 1,false);
            else
              obj.change_slice(obj.slice + 10,false);
            end
          case 'comma'
            if ~obj.shift_pressed
              obj.change_slice(obj.slice - 1,false);
            else
              obj.change_slice(obj.slice - 10,false);
            end
            
          case 'delete'
            layer_idx = get(obj.gui.layerLB,'Value');
            cmd = [];
            cmd{1}.undo.slice = obj.slice;
            cmd{1}.redo.slice = obj.slice;
            cmd{1}.undo.layer = layer_idx;
            cmd{1}.redo.layer = layer_idx;
            cmd{1}.undo.y = [];
            cmd{1}.undo.x = [];
            cmd{1}.redo.x = [];
            cmd{1}.redo.y = [];
            
            cmd{1}.undo.x = find(obj.select_mask);
            cmd{1}.undo.y = obj.layer(layer_idx).y(obj.select_mask,obj.slice);
            cmd{1}.redo.x = find(obj.select_mask);
            if islogical(obj.layer(layer_idx).y)
              if any(obj.layer(layer_idx).y(obj.select_mask,obj.slice))
                cmd{1}.redo.y = false * ones(size(cmd{1}.redo.x));
              else
                cmd{1}.redo.y = true * ones(size(cmd{1}.redo.x));
              end
            else
              cmd{1}.redo.y = NaN * ones(size(cmd{1}.redo.x));
            end

            cmd{1}.type = 'standard';
            obj.push(cmd);
            
            obj.update_slice();
            
          case 'space'
            if obj.plot_visibility == true;
              obj.plot_visibility = false;
            else
              obj.plot_visibility = true;
            end
            obj.update_slice();
            
          case 'g'
            prompt = {'Enter slice number:'};
            dlg_title = 'Go to slice';
            num_lines = 1;
            def = {sprintf('%d',obj.slice)};
            answer = inputdlg(prompt,dlg_title,num_lines,def);
            try
              obj.change_slice(str2double(answer),false)
            end
            
          case 'r'
            obj.undo_stack.redo();
            
          case 'u'
            obj.undo_stack.pop();
            
          otherwise
            
        end
        
      end
    end
    
    %% key_release
    function key_release(obj,src,event)
      
      if any(strcmp('shift',event.Modifier))
        obj.shift_pressed = true;
      else
        obj.shift_pressed = false;
      end
      
      if any(strcmp('control',event.Modifier))
        obj.ctrl_pressed = true;
      else
        obj.ctrl_pressed = false;
      end
    end
    
    %% Change slice
    function change_slice(obj, new_slice, force_update)
      if new_slice <= 0
        new_slice = 1;
      end
      if new_slice > size(obj.data,3)
        new_slice = size(obj.data,3);
      end
      if new_slice ~= obj.slice
        obj.slice = new_slice;
        obj.update_slice();
        notify(obj,'SliceChange');
        set(obj.h_control_plot,'XData',obj.slice);
        set(obj.h_layer_plot,'XData',obj.slice);
        
        xlims = xlim(obj.h_control_axes);
        ylims = ylim(obj.h_control_axes);
        if xlims(2) < obj.slice
          new_xlims = xlims + (obj.slice - 0.8*diff(xlims) - xlims(1));
        elseif xlims(1) > obj.slice
          new_xlims = xlims - (xlims(1) - (obj.slice - 0.2*diff(xlims)));
        else
          new_xlims = [];
        end
        if ~isempty(new_xlims)
          zoom_button_up(new_xlims(1),ylims(1),1,struct('x',new_xlims(2),'y',ylims(2), ...
            'h_axes',obj.h_control_axes,'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,1)]));
        end
        
        xlims = xlim(obj.h_layer_axes);
        ylims = ylim(obj.h_layer_axes);
        if xlims(2) < obj.slice
          new_xlims = xlims + (obj.slice - 0.8*diff(xlims) - xlims(1));
        elseif xlims(1) > obj.slice
          new_xlims = xlims - (xlims(1) - (obj.slice - 0.2*diff(xlims)));
        else
          new_xlims = [];
        end
        if ~isempty(new_xlims)
          zoom_button_up(new_xlims(1),ylims(1),1,struct('x',new_xlims(2),'y',ylims(2), ...
            'h_axes',obj.h_layer_axes,'xlims',[1 size(obj.data,3)],'ylims',[1 size(obj.data,2)]));
        end
        
      elseif force_update
        obj.update_slice();
      end
    end
    
    %% Update slice
    function update_slice(obj)
      set(obj.h_image,'CData',obj.data(:,:,obj.slice));
      
      title(sprintf('Slice:%d',obj.slice),'parent',obj.h_axes)

      % Update layer plots
      for layer_idx = 1:numel(obj.layer)
        if islogical(obj.layer(layer_idx).y)
          tmp_y = obj.layer(obj.layer(layer_idx).surf_layer).y(:,obj.slice);
          tmp_y(~obj.layer(layer_idx).y(:,obj.slice)) = NaN;
          set(obj.layer(layer_idx).h_plot, ...
            'XData', obj.layer(layer_idx).x(:,obj.slice), ...
            'YData', tmp_y);
        else
          set(obj.layer(layer_idx).h_plot, ...
            'XData', obj.layer(layer_idx).x(:,obj.slice), ...
            'YData', obj.layer(layer_idx).y(:,obj.slice));
        end
      end
      
      if ~isempty(get(obj.gui.layerLB,'String'))
        % Update layer selection related plots
        layer_idx = get(obj.gui.layerLB,'value');
        x_select = obj.layer(layer_idx).x(:,obj.slice);
        if islogical(obj.layer(layer_idx).y)
          y_select = obj.layer(obj.layer(layer_idx).surf_layer).y(:,obj.slice);
        else
          y_select = obj.layer(layer_idx).y(:,obj.slice);
        end
        set(obj.gui.h_select_plot,'XData',x_select(obj.select_mask), ...
          'YData',y_select(obj.select_mask),'Marker','o','LineWidth',2);
        set(obj.h_layer_image,'CData',obj.layer(layer_idx).y);
        
        % Update layer visibility
        if obj.plot_visibility == true
          for layer_idx = 1:numel(obj.layer)
            if obj.layer(layer_idx).visible
              set(obj.layer(layer_idx).h_plot,'visible','on')
            else
              set(obj.layer(layer_idx).h_plot,'visible','off')
            end
          end
          set(obj.gui.h_select_plot,'visible','on')
        else
          for layer_idx = 1:numel(obj.layer)
            set(obj.layer(layer_idx).h_plot,'visible','off')
          end
          set(obj.gui.h_select_plot,'visible','off')
        end
        
        % Update control figure plots
        layer_idx = get(obj.gui.layerLB,'value');
        if ~islogical(obj.layer(layer_idx).y)
          set(obj.h_control_layer,'XData',1:size(obj.data,3),'YData',obj.layer(layer_idx).y(ceil(size(obj.data,2)/2)+1,:));
        else
          surf_idx = obj.layer(layer_idx).surf_layer;
          new_y = obj.layer(surf_idx).y(ceil(size(obj.data,2)/2)+1,:);
          new_y(obj.layer(layer_idx).y(ceil(size(obj.data,2)/2)+1,:)) = NaN;
          set(obj.h_control_layer,'XData',1:size(obj.data,3),'YData',new_y);
        end
      end
    end
    
    %% layerLB_callback Tool
    function layerLB_callback(obj,src,event)
      obj.update_slice();
    end
    
    %% layerLB_visibility_toggle
    function layerLB_visibility_toggle(obj,src,event)
      layer_idx = get(obj.gui.layerLB,'value');
      obj.layer(layer_idx).visible = ~obj.layer(layer_idx).visible;
      obj.update_slice();
    end
    
    %% optionsPB_callback Tool
    function optionsPB_callback(obj,src,event)
      tool_idx = get(obj.gui.toolPM,'Value');
      obj.slice_tool.list{tool_idx}.open_win();
    end
    
    %% timer_callback
    function timer_callback(obj,src,event)
%      fprintf('Timer\n');
    end
    
    %% applyPB_callback Tool
    function applyPB_callback(obj,src,event)
      tool_idx = get(obj.gui.toolPM,'Value');
      obj.layer_idx = get(obj.gui.layerLB,'Value');
      cmd = obj.slice_tool.list{tool_idx}.apply_PB_callback(obj);
      if ~isempty(cmd)
        obj.undo_stack.push(cmd);
      end
    end
    
    %% Insert Tool
    function insert_tool(obj, slice_browser_tool)
      % slice_browser_tool
      obj.slice_tool.list{end+1} = slice_browser_tool;
      obj.slice_tool.list{end}.add_listener(obj);
      obj.add_listener(obj.slice_tool.list{end});
      
      toolPM_str = {};
      for idx = 1:length(obj.slice_tool.list)
        toolPM_str = [toolPM_str obj.slice_tool.list{idx}.tool_menu_name];
      end
      
      set(obj.gui.toolPM,'String',toolPM_str);
    end
    
    %% Help
    function help_menu(obj)
      fprintf('Key Short Cuts\n');
      
      fprintf('\nZoom Mode\n');
      fprintf('left-click and drag: zoom to selection\n');
      fprintf('left-click: zoom in at point\n');
      fprintf('right-click: zoom out at point\n');
      
      fprintf('\nPointer Mode In "slice" window\n');
      fprintf('left-click: set layer point (or toggle logical value)\n');
      fprintf('right-click and drag: select points (shift-key holds selection)\n');
      
      fprintf('\nPointer Mode In "layer" and "echogram" window\n');
      fprintf('left-click: sets current slice\n');
      
      fprintf('\nAll Modes\n');
      fprintf('scroll: zoom in/out at point\n');
      fprintf('delete: deletes selected points (or toggles logical values)\n');

      if ~isempty(obj.slice_tool.list)
        fprintf('\nInstalled tools:\n');
        for tool_idx = 1:length(obj.slice_tool.list)
          if ~isempty(obj.slice_tool.list{tool_idx}.help_string)
            fprintf('%s\n',obj.slice_tool.list{tool_idx}.help_string);
          end
        end
      end
    end
    
    %% getEventData
    function getEventData(obj,src,~)
      cmd = src.cmd;     
      obj.undo_stack.push(cmd);
    end
    
    %% add_listener
    function add_listener(obj,src)
       evnts = src.get_events();
       
       if ~isempty(evnts)
         for i = 1:numel(evnts);
           addlistener(evnts.src,evnts.evnts{i},@obj.getEventData);
         end
       end
       
    end
     
    %% push
    function push(obj,cmd)
      obj.layer_idx = get(obj.gui.layerLB,'Value');
      for tool_idx = 1:length(obj.slice_tool.list)
        cmd = obj.slice_tool.list{tool_idx}.push_request(cmd);
      end
      obj.undo_stack.push(cmd);
    end
    
    %% save
    function save(obj)
      fprintf('Saving surfData (%s)...\n', datestr(now));
      surf = obj.layer;
      save(obj.layer_fn,'-v7.3','surf');
      fprintf('  Done\n');
      for tool_idx = 1:length(obj.slice_tool.list)
        if ~isempty(obj.slice_tool.list{tool_idx}.save_callback) && ...
          isa(obj.slice_tool.list{tool_idx}.save_callback,'function_handle')
          obj.slice_tool.list{tool_idx}.save_callback();
        end
      end
      obj.undo_stack.save();
    end
    
  end
  
end

