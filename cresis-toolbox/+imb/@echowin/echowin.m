classdef (HandleCompatible = true) echowin < handle
% Echogram Window for imb.picker
%
% General operation:
% 1. Creator function is called
%   Initializes obj.eg fields (eg = echogram)
%   Calls "imb.echowin.create_ui" which creates the GUI objects
% 2. Draw function is called
%   Queries file system to get echograms obj.update_source_fns_existence();
%   Calls load_flightline, load_layers_init, load_layers, load_crossovers
%   Calls plot_echogram, plot_layers, plot_cursors
%   Plot functions convert to the proper units and displays them in the
%   right_panel imagesc and layer plots are created here. Note that ALL the
%   data is plotted and xlim/ylim/caxis are used to control what is seen.
% 3. If the user loads a new echogram from the map window, then the
%    draw function is called again.
% 4. If the user applies an operation, update_layer updates all the
%    variables and modifies the layer plots
% 5. If the user changes something redraw is called
%    * If the x-axis or y-axis units change or user loads a different
%    frame from within this echowin (either frameLB_callback or by
%    using the left/right arrow keys) then redraw calls 
%    * If the user causes redraw to be called for some other reason
%    and new data does not need to be loaded, then plot_echogram
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
%    is called by redraw. It is also called by a bunch of other places
%    (like keyboard short cuts and special mouse button clicks) when object
%    visibility might have changed).

  properties
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
    
%     %% GUI Properties
%     h_fig
%     left_panel % Structure with fields:
%     % handle
%     % toolPM
%     % paramPB
%     % qualityPM
%     % imagewin
%     % imagePB
%     % yaxisPM
%     % xaxisPM
%     % framesPM
%     % crossoverPB
%     % savePB
%     % topTable
%     % layerLB
%     % layerCM
%     % frameLB
%     % frameCM
%     % sourceLB
%     % sourceCM
%     % table
%     right_panel % Structure with fields:
%     % handle
%     % axes
%     % status_panel
%     % echoCM
%     % echoCM_item1
%     % table
%     table
%     cursor % Cursor handle + state information (.h, .gps_time)
%     % gps_time: GPS time of cursor location
%     % x: map x-position of cursor location
%     % y: map y-position of cursor location
%     % h: Cursor plot handle
% 
%     % Layer plotting information
%     h_image % Image handle
%     h_quality % Layer quality plot handles (6 * # of layers)
%     h_layer % Layer auto/manual plot handles (2 * # of layers)
% 
%     %% Keyboard/Mouse Properties
%     alt_pressed % Logical indicating the state of the alt-key
%     control_pressed % Logical indicating the state of the ctrl-key
%     shift_pressed % Logical indicating the state of the shift-key
%     zoom_mode % Logical indicating the zoom mode
%     busy_mode % Logical indicating the busy mode for mouse pointer
%     click_x % Last mouse click x-position
%     click_y % Last mouse click y-position
%     switch_layers % struct with fields for key_press switching layers feature
%                   % .old_time % stores time since last key press
%                   % .accumulated_event_characters % accumulates consecutive key strokes
%     
%     %% Tool properties
%     tool % Struct with tool information
%     % list: List of tools
%     % left_click_fh % Current: left click tool function handle
%     % left_click_and_drag_fh: Current left click and drag tool function handle
%     % right_click_and_drag_fh: Current right click and drag tool function handle
%     % accessed % True if any pick tool parameter window has been opened
%     % visible % True if tool parameter window is visible
%     % old_idx % Index of last tool
% 
%     %% Echogram properties
%     % eg = Struct with image information, crossover information, and
%     % echogram and layer information
%     eg 
%     % sources: cell array of echogram source file paths
%     % system: string containing accum, kaband, kuband, rds, or snow
%     % cur_sel
%     %  frame_name: '20120330_04_087'
%     %  segment_id: 2.0120e+09
%     %  season_name: '2012_Greenland_P3'
%     %  radar_name: 'snow'
%     %  location: 'arctic'
%     %  day_seg: '20120330_04'
%     %
%     % frame_idxs: frames that are actively loaded
%     % frame_names: Nfrm length cell array of frame names
%     % old_frame_idx: last frame to be loaded (-1 default)
%     % start_gps_time: Nfrm length numeric vector of start GPS times for each frame
%     % stop_gps_time: Nfrm length numeric vector of stop GPS times for each frame
%     %
%     % source_fns_existence: logical array that is Nfrm by Nsrc where Nfrm
%     %   is the number of frames in this segment and Nsrc is the number of
%     %   echogram sources in eg.sources
%     % data: original Nt by Nx single data matrix (a copy must be kept in case operations
%     %   are done on the matrix)
%     % gps_time: GPS time of data matrix, represents seconds since Jan 1,
%     %   1970 (ANSI-C standard), 1 by Nx double vector
%     % elevation: Elevation of data matrix in meters, 1 by Nx double vector
%     % latitude: Latitude of data matrix, 1 by Nx double vector
%     % longitude: Longitude of data matrix, 1 by Nx double vector
%     % surface: Surface of data matrix, 1 by Nx double vector
%     % time: Elevation of data matrix, Nt by 1 double vector
%     %
%     % image_data: modified Nt by Nx single data matrix
%     % image_xaxis
%     % image_gps_time
%     % image_yaxis
%     % x_label
%     % y_label
%     % y_order
%     %
%     % crossovers
%     %   en: logical scalar determining whether or not crossovers are loaded
%     % layers
%     %   show_manual_pts % Logical scalar, Show manual points in layer plots
%     %   show_dots_only % Logical scalar, Show dots only in layer plots
%     % layer_source
%     % layer_data_source
%     %
%     % map
%     % proj
%     % map_gps_time
%     % map_elev
%     % map_x
%     % map_y
%     % map_id
% 
%     %% Undo stack properties
%     undo_stack
%     undo_stack_save_listener
%     undo_stack_synchronize_listener
% 
%     % default_params = Default parameters loaded from default parameters file
%     default_params
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  properties (SetAccess = immutable, GetAccess = public) %constants
  end
  
  events
    close_window % Signalled when a user closes the window
    update_echowin_flightline % Signalled when the x-axis changes (draw and redraw)
    update_cursors % Signalled when a user changes the cursor (button_up)
    update_map_selection % Signalled when the frame selection changes (frameLB_callback)
    open_crossover_event % Signalled when user requests a cross over opened from the cross over window (open_crossover)
  end
  
  methods
    function obj = echowin(h_fig,default_params)
      %% Constructor Input Checks
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

      %% Constructor: General setup
      obj.h_fig = h_fig;
      obj.default_params = default_params;

      %% Constructor: Keyboard and mouse input state control
      obj.alt_pressed = false;
      obj.control_pressed = false;
      obj.shift_pressed = false;
      obj.busy_mode = false;
      obj.zoom_mode = true;
      obj.switch_layers.old_time = [];
      obj.switch_layers.accumulated_event_characters = [];
      
      %% Constructor: Initialize tool settings
      
      obj.tool.layer_multiple = 1;
      
      obj.tool_accessed = false; % True if any pick tool parameter window has been opened
      obj.old_tool_idx = 1;      % Last selected tool with a parameter window before current tool
      obj.tool_visible = false;  % True if current pick tool parameter window should be visible
      
      %% Constructor: Echogram and data (eg) setup
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
      obj.eg.x_label = '';          % string for x-label (assigned in plot_echogram)
      obj.eg.image_gps_time = [];   % same size as 'image_xaxis'
      obj.eg.image_yaxis = [];      % yaxis for 'imagesc' function
      obj.eg.y_label = '';          % string for x-label (assigned in plot_echogram)
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
      obj.eg.layers.x = {};          % Nlayers x 1 cell array containing gps_time of all pnts
      obj.eg.layers.y = {};          % Nlayers x 1 cell array containing twtt of all pnts
      obj.eg.layers.x_curUnit = {};  % layer.x converted to current x-axis units
      obj.eg.layers.y_curUnit = {};  % layer.y converted to current x-axis units
      obj.eg.layers.qual = {};       % Nlayers x 1 cell array containing quality of each layer point
      obj.eg.layers.type = {};       % Nlayers x 1 cell array containing type of each layer point
      
      obj.layer_h = [];
      obj.quality_h = [];
      obj.show_manual_pts = true;
      obj.show_dots_only = false;
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
    layerCM_callback(obj,source,event);
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


