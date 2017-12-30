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
      %     commands, cmd, that use sb.data to operate on sb.sd. You 
      %     should not modify any fields of sb.
      %  .sd: surfdata .surf struct array containing surface information
      %  .data: 3D image
      %  .slice: current slice in 3D image (third index of .data)
      %  .surf_idx: active surface
      % slices: array of slices to operate on (overrides sb.slice)
      
      % control_idx = sb.sd.surf(sb.surf_idx).gt;
      % active_idx = sb.sd.surf(sb.surf_idx).active;
      % surf_idx = sb.sd.surf(sb.surf_idx).top;
      % mask_idx = sb.sd.surf(sb.surf_idx).mask;
      % quality_idx = sb.sd.surf(sb.surf_idx).quality;
      
      % Create cmd for layer change
      cmd = [];
      
      % cmd{end+1}.undo.slice = slice;
      % cmd{end}.redo.slice = slice;
      % cmd{end}.undo.surf = surf_idx;
      % cmd{end}.redo.surf = surf_idx;
      % cmd{end}.undo.x = cols;
      % cmd{end}.undo.y = sb.sd.surf(surf_idx).y(cols,slice);
      % cmd{end}.redo.x = cols;
      % cmd{end}.redo.y = new_vals;
      % cmd{end}.type = 'standard';
      % cmd{end+1}.redo.slice = sb.slice;
      % cmd{end}.undo.slice = sb.slice;
      % cmd{end}.type = 'slice_dummy';
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


