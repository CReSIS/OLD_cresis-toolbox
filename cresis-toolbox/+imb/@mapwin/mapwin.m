classdef (HandleCompatible = true) mapwin < handle
% Map Window for imb.picker
  
  properties
    % Sub classes
    map_pref          % Preference window
    echowin_list      % Vector of imb.echowin class windows
    undo_stack_list   % Vector of imb.undo_stack class objects (one per segment that is modified)

    % GUI handles and UI related objects
    h_fig             % Figure handle
    table             % Figure's top level table
    top_panel         % Top panel and handles
    map_panel         % Axes panel and handles
    status_panel      % Bottom panel and status bar handles
    echowin_maps      % Struct vector of echowin map graphics:
                      % flightlines (h_line), labels (h_text), and cursors (h_cursor)

    % Status information for mouse/key board callbacks
    control_pressed
    shift_pressed
    click_x
    click_y
    zoom_mode
    busy_mode

    % Other properties
    default_params % Contains the default parameters for mapwin, echowin, prefwin, etc
    
    % Map properties
    cur_map_pref_settings % Struct containing currently loaded preference window settings (set in get_map)
    
    map
    % map.source % 0 for OPS map, 1 for blank, 2 for Google
    % map.scale % 1e3 (km to meters for OPS/Blank), 1 for Google
    % map.fline_source % 0 for OPS flight lines, 1 for season layerdata
    % map.proj % Matlab Projection Structure
    % map.xaxis_default % Current xaxis default bounds
    % map.yaxis_default % Current yaxis default bounds
    % map.xaxis % Current xaxis
    % map.yaxis % Current yaxis
    % map.sel % map selection (red line): .frame_name, .season_name, .segment_id
    % map.sel.frame_name
    % map.sel.season_name
    % map.sel.segment_id
    % map.CoordRefSysCode % string containing EPSG:3413, EPSG:3031, or EPSG:3857
    
    % OPS
    ops
    % ops.request % wms current request information set in get_map, query_redraw_map
    % ops.season_group_ids_as_string % Current season_group_ids loaded in string format (for WMS requests)
    % ops.seasons_as_string % Current seasons loaded in string format (for WMS requests)
    
    % Google map
    google
    % google.map % Stores imp info about the google map obj
    
    % Season layerdata:
    % (.../csarp_support/layers/layer_MAPZONE_SYSTEM_SEASON.mat)
    layerdata
    % layerdata.x
    % layerdata.y
    % layerdata.frms
    % layerdata.season_idx
    
  end
  
  properties (SetAccess = private, GetAccess = private)
    priv_prop1
  end
  
  events
  end
  
  methods
    %% mapwin(h_fig,param)
    function obj = mapwin(h_fig,param)
      % h_fig: optional scalar integer parameter which specifies a specific
      %   figure handle to use (leave empty or undefined to have a new
      %   figure handle assigned)
      % param: structure controlling general behavior
      %   .picker_param_fn: filename of picker parameters to load
      %      ---
      
      % GUI: properties
      % -------------------------------------------------------------------
      if nargin == 0 || isempty(h_fig)
        h_fig = figure;
      else
        figure(h_fig);
      end
      obj.h_fig = h_fig; % Figure handle that map is drawn in
      % Each echowin has a cursor, line, and text associated with it
      obj.echowin_maps = [];    % Struct vector holding echowin map graphics:
                                %   .h_cursor: cursor plot handle
                                %   .h_line: flight lines plot handle
                                %   .h_text: label handle
      % Mouse and keyboard callback state information
      obj.click_x = []; % Last click x position
      obj.click_y = []; % Last click y position
      obj.control_pressed = false; % Is control key currently pressed?
      obj.shift_pressed = false; % Is shift key currently pressed?
      obj.zoom_mode = true; % Show zoom mouse pointer
      obj.busy_mode = false; % Show busy mouse pointer

      % Echogram window and undo stack lists
      % -------------------------------------------------------------------
      obj.echowin_list = imb.echowin.empty(1,0);
      obj.undo_stack_list = imb.undo_stack.empty(1,0);

      % Map properties
      % -------------------------------------------------------------------
      obj.map = [];
      obj.map.source = []; % 0 for OPS, 1 for Google
      obj.map.scale = []; % 1e3 (km to meters for OPS), 1 for Google
      % Currently selected flight line information
      obj.map.sel.frame_name = ''; % Current frame name
      obj.map.sel.segment_id = []; % Current segment ID
      obj.map.sel.season_name = ''; % Current season name
      obj.map.xaxis = [];
      obj.map.yaxis = [];
      obj.map.xaxis_default = [];
      obj.map.yaxis_default = [];
      obj.map.CoordRefSysCode = [];

      % Open Polar Server (OPS) properties
      % -------------------------------------------------------------------
      obj.ops = [];
      obj.ops.proj = []; % Matlab Projection Structure
      obj.ops.request = []; % wms current request information set in get_map, query_redraw_map
      obj.ops.season_group_ids_as_string = []; % Current season_group_ids loaded in string format (for WMS requests)
      obj.ops.seasons_as_string = []; % Current seasons loaded in string format (for WMS requests)
      
      % Google map properties
      % -------------------------------------------------------------------
      obj.google = [];
      obj.google.map = []; % Stores imp info about the google map obj
      
      % layerdata properties
      % -------------------------------------------------------------------
      obj.layerdata = [];
      obj.layerdata.x = []; % Flightlines in world coordinates
      obj.layerdata.y = []; % Flightlines in world coordinates
      obj.layerdata.frms = []; % Flightline frame ids (as double)
      obj.layerdata.season_idx = []; % Flightline season index vector
      
      % Set default parameters
      % -------------------------------------------------------------------
      if ~exist('param','var')
        param = struct();
      end
      if ~isfield(param,'picker_param_fn')
        param.picker_param_fn = ct_filename_tmp([],'picker.mat');
      end
      obj.set_default_params(param.picker_param_fn);

      %% mapwin: Create user interface and preference window
      % -------------------------------------------------------------------
      try
        create_ui(obj);
        obj.map_pref = imb.prefwin([],obj.default_params.prefwin);
        obj.cur_map_pref_settings = obj.map_pref.settings;
        addlistener(obj.map_pref,'StateChange',@obj.get_map);
      catch ME
        delete(obj);
        rethrow(ME);
      end
    end
    
    %% delete()
    function delete(obj)
      try
        % Delete the map preferences window
        delete(obj.map_pref);
      end
      
      try
        % Delete the map figure handle
        delete(obj.h_fig);
      end
      
      % Delete all the echogram windows
      for idx = 1:length(obj.echowin_list)
        try
          delete(obj.echowin_list(idx));
        end
      end
      
      % Delelte all the undo_stack classes
      for idx = 1:length(obj.undo_stack_list)
        try
          delete(obj.undo_stack_list(idx));
        end
      end

    end
    
    create_ui(obj); % Create GUI, called from constructor
    query_redraw_map(obj,x_min,x_max,y_min,y_max); % Called whenever map perspective changes
    update_vector_layers(obj); % Draws/updates current selection, echowin flightlines, and cursors
    get_map(obj,hObj,event); % Makes original WMS query request
    update_map_selection(obj,param); % Updates the currect flight line selection (from mouse click)
    [status,data] = get_closest_frame(obj, sys, param);
    outside_limits = check_limits(obj,xaxis,yaxis,dir); % Support function for key_press.m pan functions, checks to see if current request is in the map limits
    set_default_params(obj,picker_param_fn); % Set the default parameters loaded from the default preferences file
    save_default_params(obj);
    [changed,pos] = compute_new_map_limits(obj,new_xdata,new_ydata); % Compute new map axis limits based on new data that must be in view
    
    % Callback Functions
    button_up(obj,src,event);
    button_down(obj,src,event);
    button_motion(obj,src,event);
    button_scroll(obj,src,event);
    flightCM_callback(obj,src,event);
    key_press(obj,src,event);
    key_release(obj,src,event);
    loadPB_callback(obj,src,event);
    prefPB_callback(obj,src,event);
    close_win(obj,varargin); % Handles mapwin close event
    search_callback(obj,src,event); % Updates the currect flight line selection (from search)
    open_crossover_echowin(obj,src,event); % Opens a cross over in a new echowindow

    % Echowin-related Callback Functions
    close_echowin(obj,h_obj,event); % Handles echowin close window event
    update_echowin_flightlines(obj,src,event); % Called when echowin x-axis (i.e. flightline position) changes, updates the map
    update_echowin_cursors(obj,src,event); % Called when a echowin cursor position changes, causes all echowin cursors to be updated
    update_map_selection_echowin(obj,src,event); % Updates the currect flight line selection (from echowin frame selection)
  end
  
end


