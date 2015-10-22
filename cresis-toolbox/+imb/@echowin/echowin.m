classdef (HandleCompatible = true) echowin < handle
% Echogram Window for imb.picker
%
% General operation:
% 1. Creator function is called
%   Initializes obj.eg fields (eg = echogram)
%   Calls "imb.echowin.create_ui" which creates the GUI objects
% 2. Draw function is called
%   Queries file system to get echograms
%   Queries database to get layers
%   Calls convert_data_to_image which converts layer and echogram data
%     to the proper units and displays them in the right_panel
%     imagesc and layer plots are created here. Note that ALL the data
%     is plotted and xlim/ylim/caxis are used to control what is seen.
% 3. If the user loads a new echogram from the map window, then the
%    draw function is called again.
% 4. If the user applies an operation, update_layer updates all the
%    variables and modifies the layer plots
% 5. If the user changes something redraw is called
%    * If the x-axis or y-axis units change or user loads a different
%    frame from within this echowin (either frameLB_callback or by
%    using the left/right arrow keys) then redraw calls convert_data_to_image
%    * If the user causes redraw to be called for some other reason
%    and new data does not need to be loaded, then convert_data_to_image
%    is not called.
% 6. Normal operation
%    keyboard_down --> used to interpret hot keys and ctrl/shift/alt keys
%    button_down --> sets things up for button_up
%    button_up --> Various different modes of operation
%      ctrl: select layer
%      shift: set cursor
%      zoom-mode: just zooms
%      alt: apply region tool (call obj.left_click_and_drag then obj.update_layers)
%      no modifiers: apply point tool (call obj.left_click then obj.update_layers)
%      right click: apply delete tool (call obj.right_click then obj.update_layers)
%  7. set_visibility: this function sets the visibility for all the objects and
%    is called by convert_data_to_image and update_layers. It is also
%    called by a bunch of other places (like keyboard short cuts and
%    special mouse button clicks) when object visibility might have changed).

  properties
    %% GUI handles and UI related objects + Data
    h_fig
    left_panel
    right_panel
    table
    
    eg % eg = Struct with image information, crossover information, and
       % echogram and layer information
    quality_h % Layer quality plot handles (6 * # of layers)
    layer_h % Layer auto/manual plot handles (2 * # of layers)
    cursor % Cursor handle + state information (.h, .gps_time)
    
    %% Keyboard and mouse function states
    alt_pressed
    control_pressed
    shift_pressed
    zoom_mode
    busy_mode
    click_x
    click_y
    switch_layers % struct with fields for key_press switching layers feature
                  % .old_time % stores time since last key press
                  % .accumulated_event_characters % accumulates consecutive key strokes

    %% General properties
    
    % default_params = Default parameters loaded from default parameters file
    default_params
    
    % Tool state information
    tool
    tool_list
    tool_param
    left_click % Current left click tool function handle
    left_click_and_drag % Current left click and drag tool function handle
    right_click_and_drag % Current right click and drag tool function handle
                  
    show_manual_pts
    show_dots_only
    layer_db_open
    
    tool_accessed
    tool_visible
    old_tool_idx
    crossovers_en

    %% Undo stack properties
    undo_stack
    undo_stack_save_listener
    undo_stack_synchronize_listener
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  properties (SetAccess = immutable, GetAccess = public) %constants
  end
  
  events
    close_window % Signalled when a user closes the window
    update_echowin_flightline % Signalled when the x-axis changes (redraw, convert_data_to_image)
    update_cursors % Signalled when a user changes the cursor (button_up)
    update_map_selection % Signalled when the frame selection changes (frameLB_callback)
    open_crossover_event % Signalled when user requests a cross over opened from the cross over window (open_crossover)
  end
  
  methods
    function obj = echowin(h_fig,default_params)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      % If passed a figure, then plot the echogram in that.  Otherwise
      % create a new figure
      if nargin == 0 || isempty(h_fig)
        h_fig = figure;
      else
        figure(h_fig);
      end
      
      % Set up defaults for the param structure
      if ~exist('default_params','var')
        default_params = struct();
      end
      if ~isfield(default_params,'max_frames')
        default_params.max_frames = 2;
      end
      echowin_pos = get(h_fig,'Position');
      if ~isfield(default_params,'x')
        default_params.x = echowin_pos(1);
      end
      if ~isfield(default_params,'y')
        default_params.y = echowin_pos(2);
      end
      if ~isfield(default_params,'w')
        default_params.w = echowin_pos(3);
      end
      if ~isfield(default_params,'h')
        default_params.h = echowin_pos(4);
      end
      
      %%% Post Initialization %%%
      % Any code, including access to object

      %% General setup
      obj.h_fig = h_fig;
      obj.default_params = default_params;

      %% Keyboard and mouse input state control
      obj.alt_pressed = false;
      obj.control_pressed = false;
      obj.shift_pressed = false;
      obj.busy_mode = false;
      obj.zoom_mode = true;
      obj.switch_layers.old_time = [];
      obj.switch_layers.accumulated_event_characters = [];
      
      %% Initialize tool settings
      
      obj.tool.layer_multiple = 1;
      obj.tool.layer_switch = false;
      
      obj.tool_accessed = false;
      obj.old_tool_idx = 1;         % enter tool by default, only updated when tool switched to has h_fig
      obj.tool_visible = false;     % tool not visible by default
      
      %% Echogram imagesc setup
      obj.eg.cur_sel.day_seg = '';
      obj.eg.frame_idxs = [];
      
      obj.eg.data = [];             % echogram data for current frames
      obj.eg.elevation = [];        % elevation of current frames
      obj.eg.gps_time = [];         % gps_time of current frames
      obj.eg.latitude = [];         % latitude of current frames
      obj.eg.longitude = [];        % longitude of current frames
      obj.eg.surface = [];          % surface of current frames
      obj.eg.time = [];             % time of current frames
      
      obj.eg.h_image = [];          % handle to imagesc object
      obj.eg.image_data = [];       % data for 'imagesc' function
      obj.eg.image_xaxis = [];      % xaxis for 'imagesc' function
      obj.eg.x_label = '';          % string for x-label (assigned in convert_data_to_image)
      obj.eg.image_gps_time = [];   % same size as 'image_xaxis'
      obj.eg.image_yaxis = [];      % yaxis for 'imagesc' function
      obj.eg.y_label = '';          % string for x-label (assigned in convert_data_to_image)
      obj.eg.y_order = '';          % string for 'reverse' or 'normal'
      
      obj.eg.old_frame_idx = -1;
      
      % Cross over information
      obj.crossovers_en = false;
      obj.eg.crossovers.gui = []; % Handle to crossovers class object
      % obj.eg.crossovers.? % What ever properties are returned by
                            % opsGetCrossovers
      obj.eg.crossovers.gps_time = []; % gps-time of cross over in source frame
      obj.eg.crossovers.h = []; % 2*N array of cross over plot handles
      obj.eg.crossovers.x_curUnit = []; % x position of cross over in current image units
      obj.eg.crossovers.y_curUnit = []; % y position of cross over in current image units
      
      % Layer setup
      obj.eg.layer_id = [];         % contains layer ID(s) for currently loaded layers
      obj.eg.layer.x = {};          % Nlayers x 1 cell array containing gps_time of all pnts
      obj.eg.layer.y = {};          % Nlayers x 1 cell array containing twtt of all pnts
      obj.eg.layer.x_curUnit = {};  % layer.x converted to current x-axis units
      obj.eg.layer.y_curUnit = {};  % layer.y converted to current x-axis units
      obj.eg.layer.qual = {};       % Nlayers x 1 cell array containing quality of each layer point
      obj.eg.layer.type = {};       % Nlayers x 1 cell array containing type of each layer point
      obj.eg.layerPB_h = [];
      
      obj.layer_h = [];
      obj.quality_h = [];
      obj.show_manual_pts = true;
      obj.show_dots_only = false;
      obj.tool.quality_only = false;
      obj.tool.quality_en = 0;      % quality view on

      % Cursor setup
      obj.cursor.gps_time = [];     % cursor location
      obj.cursor.x = [];            % cursor location
      obj.cursor.y = [];            % cursor location
      obj.cursor.h = [];            % cursor handle

      % Undo stack setup
      obj.undo_stack = [];
      obj.undo_stack_save_listener = [];
      obj.undo_stack_synchronize_listener = [];
      
      create_ui(obj);
    end
    
    function delete(obj)
      % Delete the figure handle
      try
        delete(obj.h_fig);
      end
      
      % Delete the tools
      try
        for idx = 1:length(obj.tool_list)
          try
            delete(obj.tool_list{idx});
          end
        end
      end
      
      try
        delete(obj.eg.layerPB_h);
      end
      
      try
        delete(obj.left_panel.imagewin);
      end
      
      try
        delete(obj.eg.crossovers.gui);
      end
      
      % Remove echowin from undo_stack
      obj.cmds_set_undo_stack([]);
    end
    
    %% Button, key functions
    button_down(obj,src,event);
    button_motion(obj,src,event);
    button_up(obj,src,event);
    key_press(obj,src,event);
    key_release(obj,src,event);
    button_scroll(obj,src,event);
    
    %% Echowin GUI callback functions
    close_win(obj,varargin);
    delete_layerPB_callback(obj,hObj,event);
    display_modePM_callback(obj,hObj,event);
    framesPM_callback(obj,hObj,event);
    frameLB_callback(obj,hObj,event);
    frameCM_callback(obj,hObj,event);
    layerLB_callback(obj,hObj,event);
    new_layerPB_callback(obj,hObj,event);
    paramPB_callback(obj,hObj,event);
    sourceLB_callback(obj,hObj,event);
    toolPM_callback(obj,hObj,event);
    xaxisPM_callback(obj,hObj,event);
    yaxisPM_callback(obj,hObj,event);
    savePB_callback(obj,hObj,event);
    crossoverPB_callback(obj,hObj,event);
    idx = crossovers_closest(obj,gps_time,twtt);
    open_crossover(obj,source,event);
    cursor_crossover(obj,source,event);
    cur_frame = get_crossover(pbj);
    sourceCM_callback(obj,source,event);
    [fn,comment] = imagewin_fn_callback(obj,fn);
    imagewin_save_mat_callback(obj,h_obj,event);
    layerLB_init(obj);
    layerLB_setdata(obj,data);
    layerLB_check_callback(obj,source,event);
    layerLB_radio_callback(obj,source,event);
    layerLB_slider_callback(obj,source,event);
    status_text_copy_callback(obj,source,event);
    status_text_print(obj,str,type);
    status_str = status_text_cursor(obj,param);
    
    %% Load and plot
    create_ui(obj);
    draw(obj,param);  % Load EG data, frames, etc.
    [x_min,x_max,y_min,y_max] = load_echogram(obj,desire_frame_idxs,clipped,x_min,x_max,y_min,y_max);
    load_crossovers(obj,source,event);
    load_flightline(obj);
    load_layers_init(obj);
    load_layers(obj);
    plot_echogram(obj,x_min,x_max,y_min,y_max);
    plot_crossovers(obj);
    plot_layers(obj);
    plot_cursors(obj);
    set_cursor_by_map(obj,x,y);
    update_source_fns_existence(obj);
    update_frame_and_sourceLB(obj,frm);
    
    set_visibility(obj,varargin); % Layer colors
    convert_data_to_image(obj,x_min,x_max,y_min,y_max,param);
    change_dynamic_range(obj);
    change_display_c(obj);
    new_layerPB_OKbutton_callback(obj,hObj,event);
    new_layerPB_close_callback(obj,hObj,event);
    cancel_operation = undo_stack_modified_check(obj);
    quality_menu_callback(obj,source,event);
    toggle_imagewin_visibility(obj,h_obj,event);
    update_layer_plots(obj); % Update layer plots, called from cmds_execute
    
    %% Commands/Undo stack
    cmds_set_undo_stack(obj,undo_stack); % Attaches and detaches undo stack and listener
    cmds = cmds_convert_units(obj,cmds); % Converts tool commands from current units to gps-time and twtt
    cmds_execute(obj,cmds_list,cmds_direction); % Executes a set of commands on this echowin
    cmds_synchronize(obj,varargin); % for tools (callback for the undo_stack synchronize event)
    cmds_save(obj,varargin); % for tools (callback for the undo_stack save event)
    
  end
  
end


