classdef (HandleCompatible = true) ascopewin < handle
  % Ascope Window for imb.picker
  %
  % General operation:
  % 1. Creator function is called
  %   Calls "imb.ascopewin.create_ui" which creates the GUI objects
  % 2. Waits until imb.ascopewin.update_ascope called
  % 3. Waits until imb.ascopewin.memory called
  % 4. Waits until imb.ascopewin.ascopeLB_callback called
  
  properties
    %% GUI Properties
    h_fig % Main ascopewin figure handle
    h_axes % Main axis handle
    h_ascope % A-scope plots
    h_cursor % Cursor plot
    
    left_panel % Structure with fields:
    % handle
    % ascopeLB
    % ascopeCM
    % ascopeCM_visible
    % ascopeCM_hide
    % ascopeCM_memory
    % ascopeCM_copy
    % ascopeCM_up
    % ascopeCM_down
    % ascopeCM_top
    % ascopeCM_bottom
    % ascopeCM_delete
    % xaxisPM
    % table
    right_panel % Structure with fields:
    % handle
    % axes_panel
    % status_panel
    %   handle
    %   statusText
    %   mouseCoordText
    %   table
    % table
    table
    zoom_mode
    zoom_mode_x
    zoom_mode_y
    xlims
    ylims
    cur_xaxis
    
    %% ascope Properties
    ascope % Structure with ascope information
    % ascope.sys = {}; % N_ascope length cell vector of system strings
    % ascope.season_name = {}; % N_ascope length cell vector of season name strings
    % ascope.frm_str = {}; % N_ascope length cell vector of frame strings
    % ascope.gps_time = []; % N_ascope length double vector of gps times
    % ascope.twtt = {}; % N_ascope length cell vector of twtt vectors
    % ascope.data = {}; % N_ascope length cell vector of a-scope waveforms
    % ascope.surf_twtt = []; % N_ascope length double vector of surface twtt
    % ascope.lat = []; % N_ascope length double vector of platform latitudes
    % ascope.lon = []; % N_ascope length double vector of platform longitude
    % ascope.target_twtt = []; % N_ascope length double vector of target elevations
    % ascope.visible = logical([]); % N_ascope length logical vector (plot is visible)
    % ascope.selected = logical([]); % N_ascope length logical vector (plot is selected)
    % ascope.xlims = [inf -inf];
    % ascope.ylims = [inf -inf];
    
    %% default_params Properties
    % default_params = Default parameters loaded from default parameters file
    default_params
    
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  properties (SetAccess = immutable, GetAccess = public) %constants
  end
  
  events
    StateChange % Signalled when ascopewin ascope properties change
  end
  
  methods
    function obj = ascopewin(h_fig,default_params)
      %% Constructor Input Checks
      % If passed a figure, then plot the echogram in that.  Otherwise
      % create a new figure
      if nargin == 0 || isempty(h_fig)
        h_fig = figure('Visible','off');
      else
        figure(h_fig);
      end
      
      % Set up defaults for the param structure
      if ~exist('default_params','var')
        default_params = struct();
      end
      set(obj.h_fig,'Units','pixels')
      ascopewin_pos = get(h_fig,'Position');
      set(obj.h_fig,'Units','normalized')
      if ~isfield(default_params,'x')
        default_params.x = ascopewin_pos(1);
      end
      if ~isfield(default_params,'y')
        default_params.y = ascopewin_pos(2);
      end
      if ~isfield(default_params,'w')
        default_params.w = ascopewin_pos(3);
      end
      if ~isfield(default_params,'h')
        default_params.h = ascopewin_pos(4);
      end
      
      % default_params = Default parameters loaded from default parameters file
      %% Constructor: default_params
      obj.default_params = default_params;
      
      %% Constructor: GUI
      obj.h_fig = h_fig;
      obj.h_axes = []; % Main axis handle
      obj.h_ascope = []; % Ascope plot handles
      obj.h_cursor = []; % Cursor plot handle
      obj.zoom_mode = true;
      obj.cur_xaxis = 1;
    
      obj.ascope = [];
      obj.ascope.echowin = []; % N_ascope length vector of echo window figure numbers
      obj.ascope.sys = {}; % N_ascope length cell vector of system strings
      obj.ascope.season_name = {}; % N_ascope length cell vector of season name strings
      obj.ascope.frm_str = {}; % N_ascope length cell vector of frame strings
      obj.ascope.gps_time = []; % N_ascope length double vector of gps times
      obj.ascope.twtt = {}; % N_ascope length cell vector of twtt vectors
      obj.ascope.data = {}; % N_ascope length cell vector of a-scope waveforms
      obj.ascope.surf_twtt = []; % N_ascope length double vector of surface twtt
      obj.ascope.lat = []; % N_ascope length double vector of platform latitudes
      obj.ascope.lon = []; % N_ascope length double vector of platform longitude
      obj.ascope.target_twtt = []; % N_ascope length double vector of target twtt
      obj.ascope.visible = logical([]); % N_ascope length logical vector
      obj.ascope.selected = logical([]); % N_ascope length logical vector (plot is selected)
      
      create_ui(obj);
    end
    
    function delete(obj)
      % Delete the figure handle
      try
        delete(obj.h_fig);
      end
    end
    
    %% Button, key Methods
    button_down(obj,src,event);
    button_motion(obj,src,event);
    button_up(obj,src,event);
    key_press(obj,src,event);
    button_scroll(obj,src,event);
    
    %% Echowin GUI callback Methods
    close_win(obj,varargin);
    ascopeLB_callback(obj,hObj,event);
    ascopeCM_callback(obj,hObj,event);
    xaxisPM_callback(obj,hObj,event);
    
    %% Load and plot Methods
    create_ui(obj);
    update_ascope(obj,ascope);
    memory(obj,vals);
    plot_update(obj);
    
  end
  
end
