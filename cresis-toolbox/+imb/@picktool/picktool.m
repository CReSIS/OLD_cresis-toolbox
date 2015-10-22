classdef (HandleCompatible = true) picktool < handle
  % Parent class of tools
  
  properties
    tool_name
    tool_name_title
    tool_shortcut
    help_string;
    
    parent
    h_fig
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
  
  events
    hide_param
  end
  
  methods
    function obj = picktool()
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      if nargin == 0
      else
      end
      
      %%% Post Initialization %%%
      % Any code, including access to object
      
    end
    
    function delete(obj)
      % Delete the map figure handle
      delete(obj.h_fig);
      %delete(obj);
    end

    
    % param struct fields:
    % x, y = click position or region relative to XData,YData
    %   click produces scalars, region (drag) provides 2x1 sorted bounds
    % cur_layers = currently selected layer indexes
    % cur_quality = current quality setting
    % CData = image handle data
    % XData = image handle data
    % YData = image handle data
    % point_ids = 1 by Nx vector of point path ids
    % layer.x = 1 by Nx vector of point x positions corresponding to XData
    % layer_ids = 1 by Nx vector of layer ids
    % layer.y = Nl length cell vector of 1 by Nx arrays (cell element 1 is surface)
    %   NaN = no value, otherwise value corresponding to YData
    % layer.type = Nl length cell vector of 1 by Nx arrays (cell element 1 is surface)
    %   NaN = no value, otherwise value corresponding to YData
    % layer.qual = Nl length cell vector of 1 by Nx arrays (cell element 1 is surface)
    %   NaN = no value, otherwise value corresponding to YData
    %
    % cmds is a struct array with fields to apply a command (redo) and to
    %   undo the effects of the command
    % undo_cmd = string 'insert', 'delete'
    % undo_args = cell vector of arguments
    % redo_cmd = string 'insert', 'delete'
    % redo_args = cell vector of arguments
    %
    % insert args: layer index, vector of point indexes, vector of y-values
    %   to insert, vector of types, vector of qualities
    % delete args: layer index, 4-element vector of [x-range y-range]
    
    function cmds = left_click(obj,param)
    end
    
    function cmds = left_click_and_drag(obj,param)
    end
    
    cmds = right_click_and_drag(obj,param);
    
    create_ui(obj);
    close_win(obj,varargin);

    % Generic function for inserting single points (useful for child
    % classes to implement for the left_click)
    cmds = insert_pnt(obj,param);
    
    % Utility function to find which point indexes were selected (useful
    % for child classes' left_click_and_drag function)
    [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer)
    
  end
  
end


