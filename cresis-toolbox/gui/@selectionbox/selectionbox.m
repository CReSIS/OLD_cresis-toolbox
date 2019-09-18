classdef selectionbox < handle
% selectionbox: A class for creating two listboxes with items that can be
% moved back and forth between them.  Double click to select/remove items.
%
% parent: parent GUI handle for each of the GUI objects
% name: label
% list_values: cell vector of strings
% init_value: scalar integer index into list_values vector

  properties
    h_text
    h_list_available
    h_list_selected
    h_list_availableCM
    h_list_selectedCM
    fh_available
    fh_selected
  end
  
  properties (SetAccess = private, GetAccess = private)
  end
   
  events
    selection_changed
  end
  
  methods
    function obj = selectionbox(parent,name,list_values,init_value)
      %%% Pre Initialization %%%
      % Any code not using output argument (obj)
      
      %%% Post Initialization %%%
      % Any code, including access to object
      
      obj.h_text = uicontrol('Parent',parent);
      set(obj.h_text,'Style','Text');
      set(obj.h_text,'String',name);
      
      obj.h_list_available = uicontrol('Parent',parent);
      set(obj.h_list_available,'Style','listbox');
      set(obj.h_list_available,'String',list_values);
      set(obj.h_list_available,'Value',init_value);
      set(obj.h_list_available,'Max',1e9);
      set(obj.h_list_available,'Callback',@obj.list_callback);
      
      obj.h_list_selected = uicontrol('Parent',parent);
      set(obj.h_list_selected,'Style','listbox');
      set(obj.h_list_selected,'String',[]);
      set(obj.h_list_selected,'Value',1);
      set(obj.h_list_selected,'Max',1e9);
      set(obj.h_list_selected,'Callback',@obj.list_callback);

      obj.h_list_availableCM = uicontextmenu('Parent',parent);
      % Define the context menu items and install their callbacks
      uimenu(obj.h_list_availableCM, 'Label', 'Add', 'Callback', @obj.add_callback);
      set(obj.h_list_available,'uicontextmenu',obj.h_list_availableCM)
      
      obj.h_list_selectedCM = uicontextmenu('Parent',parent);
      % Define the context menu items and install their callbacks
      uimenu(obj.h_list_selectedCM, 'Label', 'Remove', 'Callback', @obj.remove_callback);
      set(obj.h_list_selected,'uicontextmenu',obj.h_list_selectedCM)

    end
    
    function delete(obj)
      try
        delete(obj.h_list_available);
      end
      try
        delete(obj.h_list_selected);
      end
    end
    
    function add_callback(obj,status,event)
      h_status = get(status);
      if isfield(h_status,'Label') && strcmp(get(status,'Label'),'Add')
        match_idxs = get(obj.h_list_available,'Value');
        available_list = get(obj.h_list_available,'string');
        if isempty(match_idxs) || isempty(available_list); return; end;
        selected_list = get(obj.h_list_selected,'string');
        selected_list = cat(1,selected_list,available_list{match_idxs});
        selected_list = sort(selected_list);
        moved_mask = zeros(size(available_list));
        moved_mask(match_idxs) = 1;
        available_list = available_list(~moved_mask);
        set(obj.h_list_available,'value',1);
        set(obj.h_list_available,'string',available_list);
        set(obj.h_list_selected,'string',selected_list);
        notify(obj,'selection_changed');
      else
        obj.fh_available(status,event);
      end
    end
    
    function remove_callback(obj,status,event)
      h_status = get(status);
      if isfield(h_status,'Label') && strcmp(get(status,'Label'),'Remove')
        match_idxs = get(obj.h_list_selected,'Value');
        available_list = get(obj.h_list_available,'string');
        selected_list = get(obj.h_list_selected,'string');
        if isempty(match_idxs) || isempty(selected_list); return; end;
        available_list = cat(1,available_list,selected_list{match_idxs});
        available_list = sort(available_list);
        moved_mask = zeros(size(selected_list));
        moved_mask(match_idxs) = 1;
        selected_list = selected_list(~moved_mask);
        set(obj.h_list_selected,'value',1);
        set(obj.h_list_available,'string',available_list);
        set(obj.h_list_selected,'string',selected_list);
        notify(obj,'selection_changed');
      else
        obj.fh_available(status,event);
      end
    end
    
    function list_callback(obj, h_obj, event)
      if strcmpi(get(gcf,'SelectionType'),'Open')
        % Double click causes the currently selected item to switch lists
        match_idx = get(h_obj,'Value');
        available_list = get(obj.h_list_available,'string');
        selected_list = get(obj.h_list_selected,'string');
        if h_obj == obj.h_list_available
          selected_list{end+1} = available_list{match_idx};
          selected_list = sort(selected_list);
          available_list = available_list([1:match_idx-1 match_idx+1:end]);
        else
          available_list{end+1} = selected_list{match_idx};
          available_list = sort(available_list);
          selected_list = selected_list([1:match_idx-1 match_idx+1:end]);
        end
        avail_val = get(obj.h_list_available,'value');
        if avail_val > length(available_list)
          set(obj.h_list_available,'value',1);
        end
        select_val = get(obj.h_list_selected,'value');
        if select_val > length(selected_list)
          set(obj.h_list_selected,'value',1);
        end
        set(obj.h_list_available,'string',available_list);
        set(obj.h_list_selected,'string',selected_list);
        notify(obj,'selection_changed');
      end
    end
    
    function vals = get_selected_values(obj)
      vals = get(obj.h_list_selected,'Value');
    end
    
    function vals = get_selected_strings(obj)
      vals = get(obj.h_list_selected,'string');
    end
    
    function set_available(obj, vals)
      % val = cell vector of character strings (resets selectionbox)
      if ~isempty(vals)
        set(obj.h_list_available,'value',1);
      end
      set(obj.h_list_available,'string',sort(vals));
      set(obj.h_list_selected,'string',[]);
      notify(obj,'selection_changed');
    end
    
    function set_selected(obj, vals, selected)
      % val = cell vector of character strings or array of indices to strings
      % selected = logical (true to select, false to deselect)
      %
      % Strings or indices that match will be (de)selected
      
      selection_changed = false;
      for val_idx = 1:length(vals)
        if iscell(vals)
          val = vals{val_idx};
        else
          val = vals(val_idx);
        end
        
        available_list = get(obj.h_list_available,'string');
        selected_list = get(obj.h_list_selected,'string');
        
        if selected == true
          % User wishes to select val
          string_list = available_list;
        else
          % User wishes to deselect val
          string_list = selected_list;
        end
        if ischar(val)
          % Search through string list for a match
          match_idx = find(strcmp(string_list,val),1);
        else
          % Assume val is index into string list
          match_idx = val;
        end
        if ~isempty(match_idx) && match_idx >= 1 && match_idx <= length(string_list)
          % Move selection from one list to the other depending on selected
          if selected == true
            selected_list{end+1} = available_list{match_idx};
            selected_list = sort(selected_list);
            available_list = available_list([1:match_idx-1 match_idx+1:end]);
          else
            available_list{end+1} = selected_list{match_idx};
            available_list = sort(available_list);
            selected_list = selected_list([1:match_idx-1 match_idx+1:end]);
          end
          avail_val = get(obj.h_list_available,'value');
          if avail_val > length(available_list)
            set(obj.h_list_available,'value',1);
          end
          select_val = get(obj.h_list_selected,'value');
          if select_val > length(selected_list)
            set(obj.h_list_selected,'value',1);
          end
          set(obj.h_list_available,'string',available_list);
          set(obj.h_list_selected,'string',selected_list);
          selection_changed = true;
        end
      end
      
      if selection_changed
        notify(obj,'selection_changed');
      end
    end
    
  end
  
end


