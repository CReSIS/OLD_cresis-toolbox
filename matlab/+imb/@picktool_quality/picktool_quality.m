classdef picktool_quality < imb.picktool
% Enter points pick tool
  
  properties
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  methods
    function obj = picktool_quality(h_fig)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      if nargin == 0 || isempty(h_fig)
        h_fig = []; %h_fig = figure('Visible','off');
      else
        h_fig = [];%figure(h_fig,'Visible','off');
      end
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.h_fig = h_fig;
      obj.tool_name = '(q)uality';
      obj.tool_name_title = 'quality';
      obj.tool_shortcut = 'q';
      obj.help_string = sprintf('Left click: No function.\nAlt + Left click and drag: Change quality of layer points in selected region. Quality changes to the selection made in the quality pulldown menu in the left panel\n\n');

      
      %obj.create_ui;
    end
    
    create_ui(obj);
    cmds = left_click(obj,h_image,layers,cur_layers,x,y,param);
    cmds = left_click_and_drag(obj,h_layer,h_image,layers,cur_layers,x,y,param);
    
  end
  
end


