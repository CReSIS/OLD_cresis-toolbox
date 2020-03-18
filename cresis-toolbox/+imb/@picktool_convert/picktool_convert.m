classdef picktool_convert < imb.picktool
  
  properties
    table
    panel
  end
  
  properties (SetAccess = immutable, GetAccess = public) %constants
    w       % window's width
    h       % window's height
  end
  
  methods
    function obj = picktool_convert(h_fig,parent)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      if nargin == 0 || isempty(h_fig)
        h_fig = figure('Visible','off');
      else
        figure(h_fig,'Visible','off');
      end
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.parent = parent;
      obj.h_fig = h_fig;
      obj.tool_name = '(c)onvert';
      obj.tool_name_title = 'convert';
      obj.tool_shortcut = 'c';
      obj.help_string = sprintf('Left click: No function\nLeft click and drag: Convert selected layers to match the layer specified in the param window (p)\n\n');
      obj.w = 190; 
      obj.h = 170;
      
      create_ui(obj);
    end
    
    create_ui(obj);
    cmds = left_click(obj,h_image,layers,cur_layers,x,y,param);
    cmds = left_click_and_drag(obj,h_layer,h_image,layers,cur_layers,x,y,param);
    
    upPB_callback(obj,src,event);
    downPB_callback(obj,src,event);
    
  end
  
end


