classdef (HandleCompatible = true) slicetool < handle
  % Parent class of slice_browser tools
  
  properties
    tool_name
    tool_menu_name
    tool_shortcut
    help_string
    custom_data
    save_callback
    
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
    end
    
    function delete(obj)
      try; delete(obj.h_fig); end;
    end
  
    function cmd = apply_PB_callback(obj,sb)
      % sb: slice browser object. Use the following fields to create
      %     commands, cmd, that use sb.data to operate on sb.layer. You 
      %     should not modify any fields of sb.
      %  .layer: struct array containing layer information
      %  .data: 3D image
      %  .slice: current slice in 3D image (third index of .data)
      %  .layer_idx: active layer
    end
    
    function create_option_ui(obj)
    end
    
    function set_custom_data(obj,custom_data)
      obj.custom_data = custom_data;
    end

    function open_win(obj,varargin)
      if ~isempty(obj.h_fig)
        set(obj.h_fig,'Visible','on');
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


