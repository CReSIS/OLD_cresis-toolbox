classdef picktool_browse < imb.picktool & handle
  
  properties
    % h_fig Inherited
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
    ascope_memory % Signalled by right click to write to a-scope memory 
  end

  methods
    function obj = picktool_browse()
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      % The browse tool does not need a ui
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.h_fig = [];
      obj.tool_name = '(b)rowse';
      obj.tool_name_title = 'browse';
      obj.tool_shortcut = 'b';
      obj.help_string = sprintf('Left click: Open an ascope window at the current marker position (place with shift+click).\nAlt + Left click and drag: No function\n\n');
      
      obj.create_ui;
    end
    
    create_ui(obj);
    cmds = left_click(obj,h_image,layers,cur_layers,x,y,param);
    cmds = left_click_and_drag(obj,h_layer,h_image,layers,cur_layers,x,y,param);
    
    cmds = right_click(obj,h_image,layers,cur_layers,x,y,param);
  end
  
end


