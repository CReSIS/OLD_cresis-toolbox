classdef crossover < handle
  % crossover: A class for controlling and displaying crossovers.
  
  properties
    h_fig
    h_gui
    img_title
    hide_only
    crossovers
    sort_idxs
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
     
  events
    update_event
    open_crossover_event
    update_cursor
    refresh_crossovers_event
  end

  methods
    function obj = crossover(img_title, hide_only)
      % img_title = string containing figure titlebar title
      % hide_only = logical, when true, closing window just hides the
      %   window
      
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      %%% Post Initialization %%%
      % Any code, including access to object
      obj.hide_only = hide_only;
      obj.img_title = img_title;
      
      obj.create_ui();
    end
    
    function delete(obj)
      try
        delete(obj.h_fig);
      end
    end
    
    close_win(obj,varargin);
    create_ui(obj);
    
    function reset(obj,source,event)
      set(obj.h_gui.minAngleLE,'String','-inf');
      set(obj.h_gui.maxAngleLE,'String','inf');
      set(obj.h_gui.minErrorLE,'String','0');
      set(obj.h_gui.maxErrorLE,'String','inf');
      set(obj.h_gui.nanCB,'Value',1);
      obj.h_gui.crossoverLB.reset(false);
      notify(obj,'update_event');
    end
    
    function set_cursor(obj,source,event)
      notify(obj,'update_cursor');
    end
    
    function refresh_crossovers(obj,source,event)
      notify(obj,'refresh_crossovers_event');
      notify(obj,'update_event');
    end

    % Function for setting or toggling the visibility of the crossover
    % figure window
    function figure_visibility_toggle(obj,visibility_flag)
      if exist('visibility_flag','var')
        cur_visibility = visibility_flag;
      else
        cur_visibility = strcmpi(get(obj.h_fig,'Visible'),'on');
      end
      if cur_visibility
        set(obj.h_fig,'Visible','on');
        figure(obj.h_fig);
      else
        set(obj.h_fig,'Visible','off');
      end
    end
    
    % Function for setting or toggling the visibility of the crossovers
    function visibility_toggle(obj,visibility_flag)
      if exist('visibility_flag','var')
        new_visibility = visibility_flag;
      else
        new_visibility = ~get(obj.h_gui.visibleCB,'Value');
      end
      set(obj.h_gui.visibleCB,'Value',new_visibility);
      notify(obj,'update_event');
    end
    
    function set_crossovers(obj, crossovers)
      orig_idx = obj.h_gui.crossoverLB.cur_value;
      if isempty(orig_idx) || isempty(obj.crossovers) ...
          || orig_idx > length(obj.crossovers.source_point_path_id)
        orig_idx = [];
      else
        orig_cross_point_path_id = obj.crossovers.cross_point_path_id(orig_idx);
        orig_layer_id = obj.crossovers.layer_id(orig_idx);
      end
      
      obj.crossovers = crossovers;

      % Create menu strings
      menuString = {};
      for idx = 1:length(obj.crossovers.frame_name)
        menuString{end+1} = sprintf('%s: %d %6.4g m %6.4g deg ', obj.crossovers.frame_name{idx}, obj.crossovers.layer_id(idx), obj.crossovers.abs_error(idx), obj.crossovers.angle(idx));
      end
      
      mode = get(obj.h_gui.sort_orderCB,'Value');
      if mode == 1
        mode = 'descend';
      else
        mode = 'ascend';
      end
      
      % Sort crossovers based on sortPM
      sort_method = get(obj.h_gui.sortPM,'Value');
      if sort_method == 1
        obj.sort_idxs = 1:length(obj.crossovers.abs_error);
      elseif sort_method == 2
        [~,obj.sort_idxs] = sort(obj.crossovers.abs_error,1,mode);
      elseif sort_method == 3
        [~,obj.sort_idxs] = sort(obj.crossovers.angle,1,mode);
      elseif sort_method == 4
        [~,obj.sort_idxs] = sort(menuString);
        if strcmpi(mode,'descend')
          obj.sort_idxs = obj.sort_idxs(end:-1:1);
        end
      else
        [~,obj.sort_idxs] = sort(obj.crossovers.gps_time,1,mode);
      end
      
      mask = get_angle_error_mask(obj);
      
      new_idx = [];
      if ~isempty(orig_idx)
        for idx = 1:length(obj.crossovers.cross_point_path_id)
          if obj.crossovers.cross_point_path_id(idx) == orig_cross_point_path_id ...
              && obj.crossovers.layer_id(idx) == orig_layer_id
            new_idx = idx;
            break;
          end
        end
      end
     
      obj.h_gui.crossoverLB.set_list_values(menuString,mask,obj.sort_idxs,new_idx);
      
    end
    
    function val = eval_double(obj,str,default_val)
      try
        % Assumes the value entered is a matlab expression that can be evaluated
        val = str2double(str);
        if length(val) ~= 1
          val = default_val;
        end
      catch ME
        warning('Search range parameter is not valid, using default range');
        val = default_val;
      end
    end

    function mask = get_angle_error_mask(obj)
      %% Get search range from tool param window
      minAngle = obj.eval_double(get(obj.h_gui.minAngleLE,'String'),-inf);
      maxAngle = obj.eval_double(get(obj.h_gui.maxAngleLE,'String'),inf);
      minError = obj.eval_double(get(obj.h_gui.minErrorLE,'String'),0);
      maxError = obj.eval_double(get(obj.h_gui.maxErrorLE,'String'),inf);

      mask = (get(obj.h_gui.nanCB,'Value') | ~isnan(obj.crossovers.twtt)) ...
        & obj.crossovers.angle >= minAngle & obj.crossovers.angle <= maxAngle ...
        & ((obj.crossovers.abs_error >= minError & obj.crossovers.abs_error <= maxError) ...
        | isnan(obj.crossovers.abs_error));
    end

    % Called from echowin.set_visibility. That function is called in
    % response to a number of changes (including changes from the crossover
    % window (e.g. min/max angle)
    function [visible,selected] = get_crossover_visibility(obj)
      % First we need to update the listbox mask in case there have been
      % changes to these fields (which is one of the reasons this function
      % might have been called).
      obj.h_gui.crossoverLB.set_mask(obj.get_angle_error_mask());
      
      % Combine masks with visibility checkbox state
      visible = get(obj.h_gui.visibleCB,'Value') ...
        & obj.h_gui.crossoverLB.get_combined_mask();

      % Get the current selection
      selected = logical(zeros(size(visible)));
      if ~isempty(selected)
        selected(obj.h_gui.crossoverLB.cur_value) = true;
      end
    end
    
    function update_callback(obj,varargin)
      if nargin == 3 && varargin{1} == obj.h_gui.crossoverLB ...
          && obj.h_gui.crossoverLB.double_click
        %% Double click
        obj.h_gui.crossoverLB.double_click = false;
        notify(obj,'open_crossover_event')
      else
        %% Single click
        notify(obj,'update_event');
        notify(obj,'update_cursor');
      end
    end
    
    function sortPM_callback(obj,varargin)
      % Create menu strings
      menuString = {};
      for idx = 1:length(obj.crossovers.frame_name)
        menuString{end+1} = sprintf('%s: %d %6.4g m %6.4g deg', obj.crossovers.frame_name{idx}, obj.crossovers.layer_id(idx), obj.crossovers.abs_error(idx), obj.crossovers.angle(idx));
      end
      
      mode = get(obj.h_gui.sort_orderCB,'Value');
      if mode == 1
        mode = 'descend';
      else
        mode = 'ascend';
      end
      
      % Sort crossovers based on sortPM
      sort_method = get(obj.h_gui.sortPM,'Value');
      if sort_method == 1
        obj.sort_idxs = 1:length(obj.crossovers.abs_error);
      elseif sort_method == 2
        [~,obj.sort_idxs] = sort(obj.crossovers.abs_error,1,mode);
      elseif sort_method == 3
        [~,obj.sort_idxs] = sort(obj.crossovers.angle,1,mode);
      elseif sort_method == 4
        [~,obj.sort_idxs] = sort(menuString);
        if strcmpi(mode,'descend')
          obj.sort_idxs = obj.sort_idxs(end:-1:1);
        end
      else
        [~,obj.sort_idxs] = sort(obj.crossovers.gps_time,1,mode);
      end
      
      obj.h_gui.crossoverLB.set_sort(obj.sort_idxs);
    end
    
    % Get selected crossover and return the crossover's information
    % This is called from mapwin after we ask for a crossover frame to be
    % opened.
    function cur_frame = get_crossover(obj)
      cur_idx = obj.h_gui.crossoverLB.cur_value;
      cur_frame.frame_name = obj.crossovers.frame_name{cur_idx};
      cur_frame.segment_id = obj.crossovers.segment_id(cur_idx);
      cur_frame.season_name = obj.crossovers.season_name{cur_idx};
    end
    
    % Sets the selected cross over based on the original cross over data
    % index (i.e. we can sort the crossover data in different ways and this
    % idx that is passed in uses the original cross over data indexing).
    % Also notifies the echowin that a different cross over has been
    % selected.
    % Called from echowin.button_up (ctrl-click for
    %   setting closest crossover).
    function set_selected(obj,idx)
      obj.h_gui.crossoverLB.set_value(idx);
      notify(obj,'update_event');
    end
    
    function idx = get_selected(obj)
      idx = obj.h_gui.crossoverLB.cur_value;
    end
    
    function enabled = crossovers_en(obj)
      % If figure is hidden and crossover visibility is false, then we
      % say that crossover loading is disabled. Otherwise it is enabled.
      enabled = get(obj.h_gui.visibleCB,'Value') | strcmpi(get(obj.h_fig,'Visible'),'on');
    end

    % Copies the currently selected crossover entry to the clipboard
    function copy_crossover(obj,varargin)
      cur_sel = get(obj.h_gui.visibleCB,'Value');
      crossover_list = get(obj.h_gui.crossoverLB.h_list,'String');
      try
        str = crossover_list{cur_sel};
      catch
        str = '';
      end
      if ~isempty(str)
        clipboard('copy',str);
      end
    end
  end
  
end
