classdef (HandleCompatible = true) slicetool < handle
  % Parent class of slice_browser tools
  
  properties (SetAccess = protected, GetAccess = public)
    tool_name
    tool_menu_name
    tool_shortcut
    ctrl_pressed
    shift_pressed
    help_string
    save_callback
  end
  
  properties (SetAccess = protected, GetAccess = protected)
    custom_data
    
    h_fig
    gui
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
    hide_param
  end
  
  methods
    function obj = slicetool()
      % obj.create_option_ui();
      % obj.tool_name = 'TOOL NAME';
      % obj.tool_menu_name = 'TOOL NAME IN MENU';
      % obj.tool_shortcut = 'KEYBOARD_EVENT';
      % obj.ctrl_pressed = 0;
      % obj.shift_pressed = 0;
      obj.h_fig = [];
    end
    
    function delete(obj)
      try; delete(obj.h_fig); end;
    end
  
    function cmd = apply_PB_callback(obj,sb,slices)
      % sb: slice browser object. Use the following fields to create
      %     commands, cmd, that use sb.data to operate on sb.layer. You 
      %     should not modify any fields of sb.
      %  .layer: struct array containing layer information
      %  .data: 3D image
      %  .slice: current slice in 3D image (third index of .data)
      %  .layer_idx: active layer
      % slices: array of slices to operate on (overrides sb.slice)
    end
    
    function create_option_ui(obj)
    end
    
    function set_custom_data(obj,custom_data)
      obj.custom_data = custom_data;
    end

    function open_win(obj,varargin)
      if ~isempty(obj.h_fig)
        set(obj.h_fig,'Visible','on');
        figure(obj.h_fig);
      end
    end

    function close_win(obj,varargin)
      set(obj.h_fig,'Visible','off');
      notify(obj,'hide_param');
    end
    
    function add_listener(obj,src)
    end
    
    function evnts = get_events(obj)
      evnts = [];
    end
    
    function cmd = push_request(obj,cmd)
    end

  end
  
end


