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
    
    % GUI handles
    h_control_fig
    h_control_axes
    h_control_image
    h_control_plot
    
    h_fig
    h_axes
    h_image
    
    h_fig_layer
    h_axes_layer
    h_image_layer
    h_layer_plot
    
    gui
    
    slice_tool_list
    slice_tool_timer
    
    % Function handle hooks for customizing clip_matrix
    fh_button_up
    fh_key_press
    fh_button_motion
    
    % zoom_mode: boolean, x,y: used by zoom mode
    zoom_mode
    x
    y
    shift_pressed
    ctrl_pressed
    plot_visibility
    undo_stack
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
      obj.slice_tool_list = [];
      obj.slice_tool_timer = timer('TimerFcn',@obj.timer_callback,'StartDelay',2,'ExecutionMode','fixedSpacing','Period',2);
      obj.slice_tool_timer = timer;
      obj.slice_tool_timer.StartDelay = 2;
      obj.slice_tool_timer.Period = 2;
      obj.slice_tool_timer.TimerFcn = @obj.timer_callback;
      obj.slice_tool_timer.ExecutionMode = 'fixedSpacing';
      start(obj.slice_tool_timer)
      
      % Load layer data
      if isfield(param,'layer_fn') && ~isempty(param.layer_fn)
        tmp = load(param.layer_fn);
        obj.layer_fn = param.layer_fn;
        obj.layer = tmp.layer;
      else
        obj.layer = [];
        obj.layer.x = [];
        obj.layer.y  = [];
        obj.layer.name = '';
        obj.layer.plot_name_values = [];
      end
      
      obj.h_control_image = h_control_image;
      obj.h_control_axes = get(obj.h_control_image,'Parent');
      obj.h_control_fig = get(obj.h_control_axes,'Parent');
      obj.fh_button_up = param.fh_button_up;
      obj.fh_key_press = param.fh_key_press;
      obj.fh_button_motion = param.fh_button_motion;
      
      obj.h_fig_layer = figure;
      obj.h_axes_layer = axes('Parent',obj.h_fig_layer,'YDir','reverse');
      obj.h_image_layer = imagesc(NaN*zeros(size(obj.data,2),size(obj.data,3)),'parent',obj.h_axes_layer);
      colormap(obj.h_axes_layer, parula(256));
      hold(obj.h_axes_layer,'on');
      obj.h_layer_plot = plot(NaN,NaN,'parent',obj.h_axes_layer,'Marker','x','Color','black','LineWidth',2,'MarkerSize',10);
      
      obj.h_fig = figure;
      set(obj.h_fig,'DockControls','off')
      set(obj.h_fig,'NumberTitle','off');
      if strcmpi(class(obj.h_fig),'double')
        set(obj.h_fig,'Name',sprintf('%d: slice',obj.h_fig));
      else
        set(obj.h_fig,'Name',sprintf('%d: slice',obj.h_fig.Number));
      end
      set(obj.h_fig,'ToolBar','none');
      set(obj.h_fig,'MenuBar','none');
      pos = get(obj.h_fig,'Position');
      pos(3) = 750;
      pos(4) = 500;
      set(obj.h_fig,'Position',pos);

      
      obj.gui.left_panel = uipanel('parent',obj.h_fig);
      obj.gui.right_panel = uipanel('parent',obj.h_fig);
      obj.h_axes = axes('Parent',obj.gui.right_panel,'YDir','reverse');
      hold(obj.h_axes,'on');
      colormap(obj.h_axes, parula(256));
      
      obj.h_image = imagesc(obj.data(:,:,obj.slice),'parent',obj.h_axes);
      for layer_idx = 1:numel(obj.layer)
        obj.layer(layer_idx).h_plot ...
          = plot(obj.layer(layer_idx).x(:,obj.slice), ...
          obj.layer(layer_idx).y(:,obj.slice), ...
          'parent',obj.h_axes,'color','black', ...
          obj.layer(layer_idx).plot_name_values{:});
      end
      
      addlistener(obj.undo_stack,'synchronize_event',@obj.undo_sync);
      
      obj.gui.h_select_plot = plot(NaN,NaN,'m.');
      
      set(obj.h_control_fig, 'WindowButtonUpFcn', @obj.control_button_up);
      set(obj.h_fig_layer, 'WindowButtonUpFcn', @obj.control_button_up);
      
      hold(obj.h_control_axes,'on');
      obj.h_control_plot = plot(NaN,NaN,'parent',obj.h_control_axes,'Marker','x','Color','black','LineWidth',2,'MarkerSize',10);
      
      % Set figure call back functions
      set(obj.h_fig,'WindowButtonUpFcn',@obj.button_up);
      set(obj.h_fig,'WindowButtonDownFcn',@obj.button_down);
      set(obj.h_fig,'WindowButtonMotionFcn',@obj.button_motion);
      set(obj.h_fig,'WindowScrollWheelFcn',@obj.button_scroll);
      set(obj.h_fig,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig_layer,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig,'WindowKeyReleaseFcn',@obj.key_release);
      set(obj.h_fig_layer,'WindowKeyPressFcn',@obj.key_press);
      set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
      set(obj.h_fig_layer,'CloseRequestFcn',[]);
      
      % Set up zoom
      zoom_setup(obj.h_fig);
      obj.zoom_mode = true;
      set(obj.h_fig,'pointer','custom');
      
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
      
      obj.gui.layerLB = uicontrol('parent',obj.gui.left_panel);
      set(obj.gui.layerLB,'style','listbox')
      set(obj.gui.layerLB,'string',{obj.layer.name})
      set(obj.gui.layerLB,'Callback',@obj.layerLB_callback)
      set(obj.gui.layerLB,'TooltipString','Select active layer (#)');
      
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
      % set(obj.gui.toolPM,'Callback',@toolPM_callback);
      set(obj.gui.toolPM,'TooltipString','Select active tool');
      
      % obj.gui.variablePM = uicontrol('parent',obj.gui.left_panel);
      % set(obj.gui.variablePM,'style','popup');
      % set(obj.gui.variablePM,'string',{'data'});
      % set(obj.gui.variablePM,'Callback',@variablePM_callback);
      
      obj.gui.layerTXT = uicontrol('Style','text','string','Layer');
      % obj.gui.plotTXT = uicontrol('Style','text','string','Plot');
      % obj.gui.variableTXT = uicontrol('Style','text','string','Variable');
      
      
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
      
      
      %       col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.layerTXT;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      
      %       col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.LB;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      %     col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.plotTXT;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      %
      %       col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.variableTXT;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      
      
      
      %       col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.layerLB;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      %
      %       col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.plotPM;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      %
      %       col = col + 1;
      %       obj.gui.left_table.handles{row,col}   = obj.gui.variablePM;
      %       obj.gui.left_table.width(row,col)     = inf;
      %       obj.gui.left_table.height(row,col)    = 20;
      %       obj.gui.left_table.width_margin(row,col) = 1;
      %       obj.gui.left_table.height_margin(row,col) = 1;
      %
      clear row col
      table_draw(obj.gui.left_table);
      
      obj.update_slice();
    end
    
    %% destructor/delete
    function delete(obj)
      try; set(obj.h_control_fig, 'WindowButtonUpFcn', []); end;
      try; delete(obj.h_fig); end;
      try; delete(obj.h_fig_layer); end;
      for tool_idx = 1:length(obj.slice_tool_list)
        try; delete(obj.slice_tool_list{tool_idx}); end;
      end
      try; delete(obj.slice_tool_timer); end;
    end
    
    %% close_win
    function close_win(obj,h_obj,event)
      try
        delete(obj);
      end
    end
    
    %% next_button_callback
    function next_button_callback(obj,source,callbackdata)
      obj.slice = obj.slice + 1;
      obj.update_slice();
    end
    
    %% prev_button_callback
    function prev_button_callback(obj,source,callbackdata)
      obj.slice = obj.slice -1;
      obj.update_slice();
    end
    
    %% help_button_callback
    function help_button_callback(obj,source,callbackdata)
      obj.help_menu()
    end
    
    %% save_button_callback
    function save_button_callback(obj,source,callbackdata)
      layer = obj.layer;
      save(obj.layer_fn,'layer')
      obj.undo_stack.save();
    end
    
    function next10_button_callback(obj,source,callbackdata)
      obj.slice = obj.slice + 10;
      obj.update_slice();
    end
    
    function prev10_button_callback(obj,source,callbackdata)
      obj.slice = obj.slice -10;
      obj.update_slice();
    end
    
    function undo_sync(obj,source,callbackdata)
      [cmds_list,cmds_direction] =  obj.undo_stack.get_synchronize_cmds();
      if strcmp(cmds_direction,'redo')
        for cmd_idx = 1:length(cmds_list)
          for subcmd_idx = 1:length(cmds_list{cmd_idx})
            layer_idx = cmds_list{cmd_idx}{subcmd_idx}.redo.layer;
            obj.layer(layer_idx).y(round(cmds_list{cmd_idx}{subcmd_idx}.redo.x), ...
              cmds_list{cmd_idx}{subcmd_idx}.redo.slice) ...
              = cmds_list{cmd_idx}{subcmd_idx}.redo.y;
            obj.slice = cmds_list{cmd_idx}{subcmd_idx}.redo.slice;
          end
        end
      else
        for cmd_idx = 1:length(cmds_list)
          for subcmd_idx = 1:length(cmds_list{cmd_idx})
            layer_idx = cmds_list{cmd_idx}{subcmd_idx}.undo.layer;
            obj.layer(layer_idx).y(round(cmds_list{cmd_idx}{subcmd_idx}.undo.x), ...
              cmds_list{cmd_idx}{subcmd_idx}.undo.slice) ...
              = cmds_list{cmd_idx}{subcmd_idx}.undo.y;
            obj.slice = cmds_list{cmd_idx}{subcmd_idx}.undo.slice;
          end
        end
      end 
      obj.update_slice();

    end  
    
    %% control_button_up
    function control_button_up(obj,h_obj,event)
      if h_obj == obj.h_control_fig
        [x,y,but] = get_mouse_info(obj.h_control_fig,obj.h_control_axes);
      else
        [x,y,but] = get_mouse_info(obj.h_fig_layer,obj.h_axes_layer);
      end
      
      obj.slice = round(x);
      obj.update_slice();
      
    end
    
    %% button_down
    function button_down(obj,h_obj,event)
      [obj.x,obj.y,but] = get_mouse_info(obj.h_fig,obj.h_axes);
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
      
      layer_idx = get(obj.gui.layerLB,'value');
      
      if obj.zoom_mode
        zoom_button_up(x,y,but,struct('x',obj.x,'y',obj.y, ...
          'h_axes',obj.h_axes,'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
      else
        if but == 2 || but == 3
          if obj.x == x
            obj.select_mask(round(x)) = true;
          else
            obj.shift_pressed
            if ~obj.shift_pressed
              obj.select_mask(:) = false;
              obj.update_slice;
            end
            
            obj.select_mask = obj.select_mask | (obj.layer(layer_idx).x(:,obj.slice) >= min(x,obj.x) ...
              & obj.layer(layer_idx).x(:,obj.slice) <= max(x,obj.x) ...
              & obj.layer(layer_idx).y(:,obj.slice) >= min(y,obj.y) ...
              & obj.layer(layer_idx).y(:,obj.slice) <= max(y,obj.y));
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
            cmd{1}.redo.y = y;
            obj.undo_stack.push(cmd);
          end
        end
        obj.update_slice();
      end
    end
    
    %% button_motion
    function button_motion(obj,hObj,event)
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
      
      if obj.ctrl_pressed
        for tool_idx = 1:length(obj.slice_tool_list)
          if strcmpi(obj.slice_tool_list{tool_idx}.tool_shortcut, event.Key)
            obj.layer_idx = get(obj.gui.layerLB,'Value');
            obj.layer_idx = obj.layer(obj.layer_idx).active_layer;
            cmd = obj.slice_tool_list{tool_idx}.apply_PB_callback(obj);
            if ~isempty(cmd)
              obj.undo_stack.push(cmd);
            end
            return;
          end
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
              axis(obj.h_axes,'tight');
            else
              % toggle zoom mode
              obj.zoom_mode = ~obj.zoom_mode;
              if obj.zoom_mode
                set(obj.h_fig,'pointer','custom');
              else
                set(obj.h_fig,'pointer','arrow');
              end
            end
            
          case 'downarrow' % Down-arrow: Pan down
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'uparrow' % Up-arrow: Pan up
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'rightarrow' % Right arrow: Pan right
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'leftarrow' % Left arrow: Pan left
            zoom_arrow(event,struct('h_axes',obj.h_axes, ...
              'xlims',[1 size(obj.data,2)],'ylims',[1 size(obj.data,1)]));
            
          case 'period'
            if ~obj.shift_pressed
              obj.slice = obj.slice + 1;
              obj.update_slice();
            else
              obj.slice = obj.slice + 10;
              obj.update_slice();
            end
          case 'comma'
            if ~obj.shift_pressed
              obj.slice = obj.slice - 1;
              obj.update_slice();
            else
              obj.slice = obj.slice - 10;
              obj.update_slice();
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
            for k = 1:64;
              if obj.select_mask(k,1) == 1;
                cmd{1}.undo.y(end+1) = obj.layer(layer_idx).y(k,obj.slice);
                cmd{1}.undo.x(end+1) = k;
                cmd{1}.redo.x(end+1) = k;
                cmd{1}.redo.y(end+1) = NaN;
              end
            end
            obj.undo_stack.push(cmd);
            
            obj.update_slice();
            obj.select_mask = logical(zeros(size(obj.data,2),1));
            
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
              obj.slice = str2double(answer);
              obj.update_slice()
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
    
    %% Update slice
    function update_slice(obj)
      if obj.slice <= 0
        obj.slice = 1;
      end
      if obj.slice > size(obj.data,3)
        obj.slice = size(obj.data,3);
      end
      
      set(obj.h_image,'CData',obj.data(:,:,obj.slice));
      
      title(sprintf('Slice:%d',obj.slice),'parent',obj.h_axes)

      % Update layer plots
      for layer_idx = 1:numel(obj.layer)
        set(obj.layer(layer_idx).h_plot, ...
          'XData', obj.layer(layer_idx).x(:,obj.slice), ...
          'YData', obj.layer(layer_idx).y(:,obj.slice));
      end
      
      % Update layer selection related plots
      layer_idx = get(obj.gui.layerLB,'value');
      x_select = obj.layer(layer_idx).x(:,obj.slice);
      y_select = obj.layer(layer_idx).y(:,obj.slice);
      set(obj.gui.h_select_plot,'XData',x_select(obj.select_mask), ...
        'YData',y_select(obj.select_mask),'Marker','o','LineWidth',2);
      layer_idx = obj.layer(layer_idx).active_layer;
      set(obj.h_image_layer,'CData',obj.layer(layer_idx).y);

      % Update layer visibility
      for layer_idx = 1:numel(obj.layer)
        if obj.plot_visibility == true
          set (obj.layer(layer_idx).h_plot,'visible','on')
        else
          set (obj.layer(layer_idx).h_plot,'visible','off')
        end
      end
      
    end
    
    %% layerLB_callback Tool
    function layerLB_callback(obj,src,event)
      obj.update_slice();
    end
    
    %% optionsPB_callback Tool
    function optionsPB_callback(obj,src,event)
      tool_idx = get(obj.gui.toolPM,'Value');
      obj.slice_tool_list{tool_idx}.open_win();
    end
    
    %% timer_callback
    function timer_callback(obj,src,event)
%      fprintf('Timer\n');
%       obj.slice_tool_timer.start();
    end
    
    %% applyPB_callback Tool
    function applyPB_callback(obj,src,event)
      tool_idx = get(obj.gui.toolPM,'Value');
      obj.layer_idx = get(obj.gui.layerLB,'Value');
      obj.layer_idx = obj.layer(obj.layer_idx).active_layer;
      cmd = obj.slice_tool_list{tool_idx}.apply_PB_callback(obj);
      if ~isempty(cmd)
        obj.undo_stack.push(cmd);
      end
    end
    
    %% Insert Tool
    function insert_tool(obj, slice_browser_tool)
      % slice_browser_tool
      obj.slice_tool_list{end+1} = slice_browser_tool;
      
      toolPM_str = {};
      for idx = 1:length(obj.slice_tool_list)
        toolPM_str = [toolPM_str obj.slice_tool_list{idx}.tool_menu_name];
      end
      
      set(obj.gui.toolPM,'String',toolPM_str);
    end
    
    %% Help
    function help_menu(obj)
      fprintf('Key Short Cuts\n');
      
      fprintf('? Mode\n');
      fprintf('scroll: zoom in/out at point\n');
      
      fprintf('Zoom Mode\n');
      fprintf('left-click and drag: zoom to selection\n');
      fprintf('left-click: zoom in at point\n');
      fprintf('right-click: zoom out at point\n');
      fprintf('scroll: zoom in/out at point\n');
    end
  end
  
end

