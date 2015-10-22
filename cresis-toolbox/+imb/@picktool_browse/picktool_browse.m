classdef picktool_browse < imb.picktool
  
  properties
    % h_fig Inherited
    h_axes
    h_plot
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  methods
    function obj = picktool_browse(h_fig)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      % The browse tool does not need a ui
      if nargin == 0 || isempty(h_fig)
        h_fig = figure('Visible','off');
      else
        h_fig = figure(h_fig,'Visible','off');
      end
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.h_fig = h_fig;
      obj.h_axes = axes;
      obj.h_plot= plot(NaN,NaN);
      obj.tool_name = '(b)rowse';
      obj.tool_name_title = 'browse';
      obj.tool_shortcut = 'b';
      obj.help_string = sprintf('Left click: Open an ascope window at the current marker position (place with shift+click).\nAlt + Left click and drag: No function\n\n');
      
      obj.create_ui;
    end
    
    create_ui(obj);
    cmds = left_click(obj,h_image,layers,cur_layers,x,y,param);
    cmds = left_click_and_drag(obj,h_layer,h_image,layers,cur_layers,x,y,param);
    
  end
  
end


