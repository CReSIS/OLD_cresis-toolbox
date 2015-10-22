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
    
    % Vector graphics
    echowin_maps      % Struct vector of echowin map graphics:
                      % flightlines (h_line), labels (h_text), and cursors (h_cursor)
    cur_sel           % Current map selection struct (red line):
                      % .frame_name, .season_name, .segment_id
    map               % Map image handle and auxilliary data

    cur_map_pref_settings % Struct containing currently loaded preference window settings (set in get_map)
    wms % wms created in get_map
    proj % wms projection set in get_map
    cur_request % wms current request information set in get_map, query_redraw_map
    full_xaxis % axis limits set in get_map, query_redraw_map
    full_yaxis % axis limits set in get_map, query_redraw_map
    
    season_group_ids_as_string % Current season_group_ids loaded in string format (for WMS requests)
    season_group_ids_modrequest
    seasons_as_string % Current seasons loaded in string format (for WMS requests)
    seasons_modrequest
    
    % Other properties
    default_params % Contains the default parameters for mapwin, echowin, prefwin, etc
    
    frame_text
    
    xaxis
    yaxis
    radar_type
    
    curFrame
    param

    % Status information for mouse/key board callbacks
    control_pressed
    shift_pressed
    click_x
    click_y
    zoom_mode
    busy_mode
    
  end
  
  properties (SetAccess = private, GetAccess = private)
    priv_prop1
  end
  
  events
  end
  
  methods
    function obj = mapwin(h_fig,param)
      % h_fig: optional scalar integer parameter which specifies a specific
      %   figure handle to use (leave empty or undefined to have a new
      %   figure handle assigned)
      % param: structure controlling general behavior
      %   .picker_param_fn: filename of picker parameters to load
      %      ---
      
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      if nargin == 0 || isempty(h_fig)
        h_fig = figure;
      else
        figure(h_fig);
      end
      if ~exist('param','var')
        param = struct();
      end
      if ~isfield(param,'picker_param_fn')
        param.picker_param_fn = ct_filename_tmp([],'picker.mat');
      end
      
      %%% Post Initialization %%%
      % Any code, including access to object

      obj.h_fig = h_fig; % Figure handle that map is drawn in
      obj.set_default_params(param.picker_param_fn);
      
      obj.echowin_list = imb.echowin.empty(1,0);
      obj.undo_stack_list = imb.undo_stack.empty(1,0);
      % Each echowin has a cursor, line, and text associated with it
      obj.echowin_maps = [];    % Struct vector holding echowin map graphics:
                                %   .h_cursor: cursor plot handle
                                %   .h_line: flight lines plot handle
                                %   .h_text: label handle

      % Currently selected flight line information
      obj.cur_sel = [];
      obj.cur_sel.frame_name = ''; % Current frame name
      obj.cur_sel.segment_id = []; % Current segment ID
      obj.cur_sel.season_name = ''; % Current season name
      
      obj.click_x = [];
      obj.click_y = [];
      obj.curFrame = [];
      
      % Mouse and keyboard callback state information
      obj.control_pressed = false;
      obj.shift_pressed = false;
      obj.zoom_mode = true;
      obj.busy_mode = false;

      %% Create user interface and preference window
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
    outside_limits = check_limits(obj,xaxis,yaxis,dir); % Support function for key_press.m pan functions, checks to see if current request is in the map limits
    set_default_params(obj,picker_param_fn); % Set the default parameters loaded from the default preferences file
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


